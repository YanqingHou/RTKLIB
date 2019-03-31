/* ------------------------------------------------------------------------------
 * rtkpos.c : precise positioning
 *
 *          Copyright (C) 2007-2013 by T.TAKASU, All rights reserved.
 *
 * version : $Revision: 1.1 $ $Date: 2008/07/17 21:48:06 $
 * history : 2007/01/12 1.0  new
 *           2007/03/13 1.1  add slip detection by LLI flag
 *           2007/04/18 1.2  add antenna pcv correction
 *                           change rtkpos argin
 *           2008/07/18 1.3  refactored
 *           2009/01/02 1.4  modify rtk positioning api
 *           2009/03/09 1.5  support glonass, gallileo and qzs
 *           2009/08/27 1.6  fix bug on numerical exception
 *           2009/09/03 1.7  add check of valid satellite number
 *                           add check time sync for moving-base
 *           2009/11/23 1.8  add api rtkopenstat(),rtkclosestat()
 *                           add receiver h/w bias estimation
 *                           add solution status output
 *           2010/04/04 1.9  support ppp-kinematic and ppp-static modes
 *                           support earth tide correction
 *                           changed api:
 *                               rtkpos()
 *           2010/09/07 1.10 add elevation mask to hold ambiguity
 *           2012/02/01 1.11 add extended receiver error model
 *                           add glonass interchannel bias correction
 *                           add slip detectior by L1-L5 gf jump
 *                           output snr of rover receiver in residuals
 *           2013/03/10 1.12 add otl and pole tides corrections
 *           2014/05/26 1.13 support beidou and galileo
 *                           add output of gal-gps and bds-gps time offset
 *           2014/05/28 1.14 fix bug on memory exception with many sys and freq
 *           2014/08/26 1.15 add functino to swap sol-stat file with keywords
 *           2014/10/21 1.16 fix bug on beidou amb-res with pos2-bdsarmode=0
 *           2014/11/08 1.17 fix bug on ar-degradation by unhealthy satellites
 *----------------------------------------------------------------------------- */
#include <stdarg.h>
#include "rtklib.h"

static const char rcsid[] = "$Id:$";

/* constants/macros ---------------------------------------------------------- */

#define SQR(x)      ((x)*(x))
#define SQRT(x)     ((x)<=0.0?0.0:sqrt(x))
#define MIN(x,y)    ((x)<=(y)?(x):(y))
#define ROUND(x)    (int)floor((x)+0.5)

#define VAR_POS     SQR(30.0) /* initial variance of receiver pos (m^2) */
#define VAR_VEL     SQR(10.0) /* initial variance of receiver vel ((m/s)^2) */
#define VAR_ACC     SQR(10.0) /* initial variance of receiver acc ((m/ss)^2) */
#define VAR_HWBIAS  SQR(1.0)  /* initial variance of h/w bias ((m/MHz)^2) */
#define VAR_GRA     SQR(0.001) /* initial variance of gradient (m^2) */
#define INIT_ZWD    0.15     /* initial zwd (m) */

#define PRN_HWBIAS  1E-6     /* process noise of h/w bias (m/MHz/sqrt(s)) */
#define GAP_RESION  120      /* gap to reset ionosphere parameters (epochs) */
#define MAXACC      30.0     /* max accel for doppler slip detection (m/s^2) */

#define VAR_HOLDAMB 0.001    /* constraint to hold ambiguity (cycle^2) */

#define TTOL_MOVEB  (1.0+2*DTTOL)
/* time sync tolerance for moving-baseline (s) */

/* number of parameters (pos,ionos,tropos,hw-bias,phase-bias,real,estimated) */
#define NF(opt)     ((opt)->ionoopt==IONOOPT_IFLC?1:(opt)->nf)
#define NP(opt)     ((opt)->dynamics==0?3:9)
#define NI(opt)     ((opt)->ionoopt!=IONOOPT_EST?0:MAXSAT)
#define NT(opt)     ((opt)->tropopt<TROPOPT_EST?0:((opt)->tropopt<TROPOPT_ESTG?2:6))
#define NL(opt)     ((opt)->glomodear!=2?0:NFREQGLO)
#define NB(opt)     ((opt)->mode<=PMODE_DGPS?0:MAXSAT*NF(opt))
#define NR(opt)     (NP(opt)+NI(opt)+NT(opt)+NL(opt))
#define NX(opt)     (NR(opt)+NB(opt))

/* state variable index */
#define II(s,opt)   (NP(opt)+(s)-1)                 /* ionos (s:satellite no) */
#define IT(r,opt)   (NP(opt)+NI(opt)+NT(opt)/2*(r)) /* tropos (r:0=rov,1:ref) */
#define IL(f,opt)   (NP(opt)+NI(opt)+NT(opt)+(f))   /* receiver h/w bias */
#define IB(s,f,opt) (NR(opt)+MAXSAT*(f)+(s)-1) /* phase bias (s:satno,f:freq) */

#ifdef EXTGSI

extern int resamb_WLNL(rtk_t *rtk, const obsd_t *obs, const int *sat,
                       const int *iu, const int *ir, int ns, const nav_t *nav, const double *azel);
extern int resamb_TCAR(rtk_t *rtk, const obsd_t *obs, const int *sat,
                       const int *iu, const int *ir, int ns, const nav_t *nav, const double *azel);
#else
/*   commented by Yanqing Hou 2016   */
/**/
extern int resamb_WLNL(rtk_t *rtk, const obsd_t *obs, const int *sat,
                       const int *iu, const int *ir, int ns, const nav_t *nav, const double *azel)
{
    return 0;
}

extern int resamb_TCAR(rtk_t *rtk, const obsd_t *obs, const int *sat,
                       const int *iu, const int *ir, int ns, const nav_t *nav, const double *azel)
{
    return 0;
}
#endif

/* global variables ---------------------------------------------------------- */
static int statlevel = 0; /* rtk status output level (0:off) */
static FILE *fp_stat = NULL; /* rtk status file pointer */
static char file_stat[1024] = ""; /* rtk status file original path */
static gtime_t time_stat = {0}; /* rtk status file time */

/* open solution status file ---------------------------------------------------
 * open solution status file and set output level
 * args   : char     *file   I   rtk status file
 *          int      level   I   rtk status level (0: off)
 * return : status (1:ok,0:error)
 * notes  : file can constain time keywords (%Y,%y,%m...) defined in reppath().
 *          The time to replace keywords is based on UTC of CPU time.
 * output : solution status file record format
 *
 *   $POS,week,tow,stat,posx,posy,posz,posxf,posyf,poszf
 *          week/tow : gps week no/time of week (s)
 *          stat     : solution status
 *          posx/posy/posz    : position x/y/z ecef (m) float
 *          posxf/posyf/poszf : position x/y/z ecef (m) fixed
 *
 *   $VELACC,week,tow,stat,vele,veln,velu,acce,accn,accu,velef,velnf,veluf,accef,accnf,accuf
 *          week/tow : gps week no/time of week (s)
 *          stat     : solution status
 *          vele/veln/velu    : velocity e/n/u (m/s) float
 *          acce/accn/accu    : acceleration e/n/u (m/s^2) float
 *          velef/velnf/veluf : velocity e/n/u (m/s) fixed
 *          accef/accnf/accuf : acceleration e/n/u (m/s^2) fixed
 *
 *   $CLK,week,tow,stat,clk1,clk2,clk3,clk4
 *          week/tow : gps week no/time of week (s)
 *          stat     : solution status
 *          clk1     : receiver clock bias GPS (ns)
 *          clk2     : receiver clock bias GLO-GPS (ns)
 *          clk3     : receiver clock bias GAL-GPS (ns)
 *          clk4     : receiver clock bias BDS-GPS (ns)
 *
 *   $ION,week,tow,stat,sat,az,el,ion,ion-fixed
 *          week/tow : gps week no/time of week (s)
 *          stat     : solution status
 *          sat      : satellite id
 *          az/el    : azimuth/elevation angle(deg)
 *          ion      : vertical ionospheric delay L1 (m) float
 *          ion-fixed: vertical ionospheric delay L1 (m) fixed
 *
 *   $TROP,week,tow,stat,rcv,ztd,ztdf
 *          week/tow : gps week no/time of week (s)
 *          stat     : solution status
 *          rcv      : receiver (1:rover,2:base station)
 *          ztd      : zenith total delay (m) float
 *          ztdf     : zenith total delay (m) fixed
 *
 *   $HWBIAS,week,tow,stat,frq,bias,biasf
 *          week/tow : gps week no/time of week (s)
 *          stat     : solution status
 *          frq      : frequency (1:L1,2:L2,...)
 *          bias     : h/w bias coefficient (m/MHz) float
 *          biasf    : h/w bias coefficient (m/MHz) fixed
 *
 *   $SAT,week,tow,sat,frq,az,el,resp,resc,vsat,snr,fix,slip,lock,outc,slipc,rejc
 *          week/tow : gps week no/time of week (s)
 *          sat/frq  : satellite id/frequency (1:L1,2:L2,...)
 *          az/el    : azimuth/elevation angle (deg)
 *          resp     : pseudorange residual (m)
 *          resc     : carrier-phase residual (m)
 *          vsat     : valid data flag (0:invalid,1:valid)
 *          snr      : signal strength (dbHz)
 *          fix      : ambiguity flag  (0:no data,1:float,2:fixed,3:hold,4:ppp)
 *          slip     : cycle-slip flag (bit1:slip,bit2:parity unknown)
 *          lock     : carrier-lock count
 *          outc     : data outage count
 *          slipc    : cycle-slip count
 *          rejc     : data reject (outlier) count
 *
 *----------------------------------------------------------------------------- */
extern int rtkopenstat(const char *file, int level) {
    gtime_t time = utc2gpst(timeget());
    char path[1024];
    
    RTKtrace(3, "rtkopenstat: file=%s level=%d\n", file, level);
    
    if (level <= 0)
        return 0;
    
    reppath(file, path, time, "", "");
    
    if (!(fp_stat = fopen(path, "w"))) {
        RTKtrace(1, "rtkopenstat: file open error path=%s\n", path);
        return 0;
    }
    strcpy(file_stat, file);
    time_stat = time;
    statlevel = level;
    return 1;
}

/* close solution status file --------------------------------------------------
 * close solution status file
 * args   : none
 * return : none
 *----------------------------------------------------------------------------- */
extern void rtkclosestat(void) {
    RTKtrace(3, "rtkclosestat:\n");
    
    if (fp_stat)
        fclose(fp_stat);
    fp_stat = NULL;
    file_stat[0] = '\0';
    statlevel = 0;
}

/* swap solution status file ------------------------------------------------- */
static void swapsolstat(void) {
    gtime_t time = utc2gpst(timeget());
    char path[1024];
    
    if ((int)(time2gpst(time, NULL) / INT_SWAP_STAT) == (int)
        (time2gpst(time_stat, NULL) / INT_SWAP_STAT)) {
        return;
    }
    time_stat = time;
    
    if (!reppath(file_stat, path, time, "", "")) {
        return;
    }
    if (fp_stat)
        fclose(fp_stat);
    
    if (!(fp_stat = fopen(path, "w"))) {
        RTKtrace(2, "swapsolstat: file open error path=%s\n", path);
        return;
    }
    RTKtrace(3, "swapsolstat: path=%s\n", path);
}

/* output solution status ---------------------------------------------------- */
static void outsolstat(rtk_t *rtk) {
    ssat_t *ssat;
    double tow, pos[3], vel[3], acc[3], vela[3] = {0}, acca[3] = {0}, xa[3];
    int i, j, week, est, nfreq, nf = NF(&rtk->opt);
    char id[32];
    
    if (statlevel <= 0 || !fp_stat)
        return;
    
    RTKtrace(3, "outsolstat:\n");
    
    /* swap solution status file */
    swapsolstat();
    
    est = rtk->opt.mode >= PMODE_DGPS;
    nfreq = est ? nf : 1;
    tow = time2gpst(rtk->sol.time, &week);
    
    /* receiver position */
    if (est) {
        for (i = 0; i < 3; i++)
            xa[i] = i < rtk->na ? rtk->xa[i] : 0.0;
        fprintf(fp_stat, "$POS,%d,%.3f,%d,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f\n",
                week, tow, rtk->sol.stat, rtk->x[0], rtk->x[1], rtk->x[2], xa[0],
                xa[1], xa[2]);
    }
    else {
        fprintf(fp_stat, "$POS,%d,%.3f,%d,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f\n",
                week, tow, rtk->sol.stat, rtk->sol.rr[0], rtk->sol.rr[1],
                rtk->sol.rr[2], 0.0, 0.0, 0.0);
    }
    /* receiver velocity and acceleration */
    if (est && rtk->opt.dynamics) {
        ecef2pos(rtk->sol.rr, pos);
        ecef2enu(pos, rtk->x + 3, vel);
        ecef2enu(pos, rtk->x + 6, acc);
        if (rtk->na >= 6)
            ecef2enu(pos, rtk->xa + 3, vela);
        if (rtk->na >= 9)
            ecef2enu(pos, rtk->xa + 6, acca);
        fprintf(fp_stat,
                "$VELACC,%d,%.3f,%d,%.4f,%.4f,%.4f,%.5f,%.5f,%.5f,%.4f,%.4f,%.4f,%.5f,%.5f,%.5f\n",
                week, tow, rtk->sol.stat, vel[0], vel[1], vel[2], acc[0], acc[1],
                acc[2], vela[0], vela[1], vela[2], acca[0], acca[1], acca[2]);
    }
    else {
        ecef2pos(rtk->sol.rr, pos);
        ecef2enu(pos, rtk->sol.rr + 3, vel);
        fprintf(fp_stat,
                "$VELACC,%d,%.3f,%d,%.4f,%.4f,%.4f,%.5f,%.5f,%.5f,%.4f,%.4f,%.4f,%.5f,%.5f,%.5f\n",
                week, tow, rtk->sol.stat, vel[0], vel[1], vel[2], 0.0, 0.0, 0.0,
                0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
    }
    /* receiver clocks */
    fprintf(fp_stat, "$CLK,%d,%.3f,%d,%d,%.3f,%.3f,%.3f,%.3f\n", week, tow,
            rtk->sol.stat, 1, rtk->sol.dtr[0]*1E9, rtk->sol.dtr[1]*1E9,
            rtk->sol.dtr[2]*1E9, rtk->sol.dtr[3]*1E9);
    
    /* ionospheric parameters */
    if (est && rtk->opt.ionoopt == IONOOPT_EST) {
        for (i = 0; i < MAXSAT; i++) {
            ssat = rtk->ssat + i;
            if (!ssat->vs)
                continue;
            satno2id(i + 1, id);
            j = II(i + 1, &rtk->opt);
            xa[0] = j < rtk->na ? rtk->xa[j] : 0.0;
            fprintf(fp_stat, "$ION,%d,%.3f,%d,%s,%.1f,%.1f,%.4f,%.4f\n", week,
                    tow, rtk->sol.stat, id, ssat->azel[0]*R2D, ssat->azel[1]*R2D,
                    rtk->x[j], xa[0]);
        }
    }
    /* tropospheric parameters */
    if (est && (rtk->opt.tropopt == TROPOPT_EST ||
                rtk->opt.tropopt == TROPOPT_ESTG)) {
        for (i = 0; i < 2; i++) {
            j = IT(i, &rtk->opt);
            xa[0] = j < rtk->na ? rtk->xa[j] : 0.0;
            fprintf(fp_stat, "$TROP,%d,%.3f,%d,%d,%.4f,%.4f\n", week, tow,
                    rtk->sol.stat, i + 1, rtk->x[j], xa[0]);
        }
    }
    /* receiver h/w bias */
    if (est && rtk->opt.glomodear == 2) {
        for (i = 0; i < nfreq; i++) {
            j = IL(i, &rtk->opt);
            xa[0] = j < rtk->na ? rtk->xa[j] : 0.0;
            fprintf(fp_stat, "$HWBIAS,%d,%.3f,%d,%d,%.4f,%.4f\n", week, tow,
                    rtk->sol.stat, i + 1, rtk->x[j], xa[0]);
        }
    }
    if (rtk->sol.stat == SOLQ_NONE || statlevel <= 1)
        return;
    
    /* residuals and status */
    for (i = 0; i < MAXSAT; i++) {
        ssat = rtk->ssat + i;
        if (!ssat->vs)
            continue;
        satno2id(i + 1, id);
        for (j = 0; j < nfreq; j++) {
            fprintf(fp_stat,
                    "$SAT,%d,%.3f,%s,%d,%.1f,%.1f,%.4f,%.4f,%d,%.0f,%d,%d,%d,%d,%d,%d\n",
                    week, tow, id, j + 1, ssat->azel[0]*R2D, ssat->azel[1]*R2D,
                    ssat->resp[j], ssat->resc[j], ssat->vsat[j], ssat->snr[j]*0.25,
                    ssat->fix[j], ssat->slip[j]&3, ssat->lock[j], ssat->outc[j],
                    ssat->slipc[j], ssat->rejc[j]);
        }
    }
}

/* save error message -------------------------------------------------------- */
static void errmsg(rtk_t *rtk, const char *format, ...) {
    char buff[256], tstr[32];
    int n;
    va_list ap;
    time2str(rtk->sol.time, tstr, 2);
    n = sprintf(buff, "%s: ", tstr + 11);
    va_start(ap, format);
    n += vsprintf(buff + n, format, ap);
    va_end(ap);
    n = n < MAXERRMSG - rtk->neb ? n : MAXERRMSG - rtk->neb;
    memcpy(rtk->errbuf + rtk->neb, buff, n);
    rtk->neb += n;
    RTKtrace(2, "%s", buff);
}

/* single-differenced observable --------------------------------------------- */
static double sdobs(const obsd_t *obs, int i, int j, int f) {
    double pi = f < NFREQ ? obs[i].L[f] : obs[i].P[f - NFREQ];
    double pj = f < NFREQ ? obs[j].L[f] : obs[j].P[f - NFREQ];
    return pi == 0.0 || pj == 0.0 ? 0.0 : pi - pj;
}

/* single-differenced geometry-free linear combination of phase -------------- */
static double gfobs_L1L2(const obsd_t *obs, int i, int j, const double *lam) {
    /*SDL1*LAM1-SDL2*LAM2*/
    double pi = sdobs(obs, i, j, 0) * lam[0], pj = sdobs(obs, i, j, 1) * lam[1];
    return pi == 0.0 || pj == 0.0 ? 0.0 : pi - pj;
}

static double gfobs_L1L5(const obsd_t *obs, int i, int j, const double *lam) {
    double pi = sdobs(obs, i, j, 0) * lam[0], pj = sdobs(obs, i, j, 2) * lam[2];
    return pi == 0.0 || pj == 0.0 ? 0.0 : pi - pj;
}

/* single-differenced measurement error variance ----------------------------- */
static double varerr(int sat, int sys, double el, double bl, double dt, int f,
                     const prcopt_t *opt) {
    //TODO: varerr
    double a, b, c = opt->err[3] * bl / 1E4, d = CLIGHT * opt->sclkstab * dt,
    fact = 1.0;
    double sinel = sin(el);
    int i = sys == SYS_GLO ? 1 : (sys == SYS_GAL ? 2 : 0), nf = NF(opt);
    
    /* extended error model */
    if (f >= nf && opt->exterr.ena[0]) { /* code */
        a = opt->exterr.cerr[i][(f - nf) * 2];
        b = opt->exterr.cerr[i][1 + (f - nf) * 2];
        if (sys == SYS_SBS) {
            a *= EFACT_SBS;
            b *= EFACT_SBS;
        }
    }
    else if (f < nf && opt->exterr.ena[1]) { /* phase */
        a = opt->exterr.perr[i][f * 2];
        b = opt->exterr.perr[i][1 + f * 2];
        if (sys == SYS_SBS) {
            a *= EFACT_SBS;
            b *= EFACT_SBS;
        }
    }
    else { /* normal error model */
        if (f >= nf)
            fact = opt->eratio[f - nf];
        if (fact <= 0.0)
            fact = opt->eratio[0];
        fact *= sys == SYS_GLO ? EFACT_GLO :
        (sys == SYS_SBS ? EFACT_SBS : EFACT_GPS);
        a = fact * opt->err[1];
        b = fact * opt->err[2];
    }
    return 2.0 * (opt->ionoopt == IONOOPT_IFLC ? 3.0 : 1.0) *
    (a * a + b * b / sinel / sinel + c * c) + d * d;
}


/* single-differenced measurement error variance ----------------------------- */
static double varerrYQ(int sat, int sys, double el, double bl, double dt, int f,
                       const prcopt_t *opt,int mode) {
    /*mode: 0--normal, 1--exponential, 2--sin*/
    double a, b, c = opt->err[3] * bl / 1E4, d = CLIGHT * opt->sclkstab * dt,
    fact = 1.0;
    double sinel = sin(el);
    int i = sys == SYS_GLO ? 1 : (sys == SYS_GAL ? 2 : 0), nf = NF(opt);
    
    double svar=0.0;
    double gpssig0[2]={0.0232,0.0316},gpssiga[2]={0.4817,0.5943},
    bdssig0[2]={0.0587,0.0383},bdssiga[2]={0.6125,0.4358},
    gpssin0[2]={0.0130,0.0062},gpssina[2]={0.0355,0.0490},
    bdssin0[2]={0.0330,0.0177},bdssina[2]={0.0487,0.0357};
    
    /* extended error model */
    if (f >= nf && opt->exterr.ena[0]) { /* code */
        a = opt->exterr.cerr[i][(f - nf) * 2];
        b = opt->exterr.cerr[i][1 + (f - nf) * 2];
        if (sys == SYS_SBS) {
            a *= EFACT_SBS;
            b *= EFACT_SBS;
        }
    }
    else if (f < nf && opt->exterr.ena[1]) { /* phase */
        a = opt->exterr.perr[i][f * 2];
        b = opt->exterr.perr[i][1 + f * 2];
        if (sys == SYS_SBS) {
            a *= EFACT_SBS;
            b *= EFACT_SBS;
        }
    }
    else { /* normal error model */
        if (f >= nf)
            fact = opt->eratio[f - nf];
        if (fact <= 0.0)
            fact = opt->eratio[0];
        fact *= sys == SYS_GLO ? EFACT_GLO :
        (sys == SYS_SBS ? EFACT_SBS : EFACT_GPS);
        a = fact * opt->err[1];
        b = fact * opt->err[2];
    }
    /*added by Yanqing 2016-09-29  */
    svar= (a * a + b * b / sinel / sinel + c * c) + d * d;
    if (f>=nf) {/*code, use my own elevation model */
        switch (mode) {
            case 1:  /* exponential model */
                
                if (sys==SYS_GPS) {
                    svar=gpssig0[f-nf]+gpssiga[f-nf]*exp(-R2D*el/15);
                }else if(sys==SYS_CMP){
                    svar=bdssig0[f-nf]+bdssiga[f-nf]*exp(-R2D*el/15);
                }
                else{
                    ;
                }
                
                
                
                break;
            case 2:   /*sin model*/
                if (sys==SYS_GPS) {
                    svar=gpssin0[f-nf]+gpssina[f-nf]/sinel;
                }else if(sys==SYS_CMP){
                    svar=bdssin0[f-nf]+bdssina[f-nf]/sinel;
                }
                else{
                    ;
                }
                
                break;
            default:
                ;
                
        }
        
    }		 else
        ;/*phase, use the default model   */
    return 2.0 * (opt->ionoopt == IONOOPT_IFLC ? 3.0 : 1.0) * svar;
    /*	(a * a + b * b / sinel / sinel + c * c) + d * d; */
}

/* baseline length ----------------------------------------------------------- */
static double baseline(const double *ru, const double *rb, double *dr) {
    int i;
    for (i = 0; i < 3; i++)
        dr[i] = ru[i] - rb[i];
    return norm(dr, 3);
}

/* initialize state and covariance ------------------------------------------- */
static void initx(rtk_t *rtk, double xi, double var, int i) {
    int j;
    rtk->x[i] = xi;
    for (j = 0; j < rtk->nx; j++) {
        rtk->P[i + j * rtk->nx] = rtk->P[j + i * rtk->nx] = i == j ? var : 0.0;
    }
}

/* select common satellites between rover and reference station -------------- */
static int selsat(const obsd_t *obs, double *azel, int nu, int nr,
                  const prcopt_t *opt, int *sat, int *iu, int *ir) {
    int i, j, k = 0;
    
    RTKtrace(3, "selsat  : nu=%d nr=%d\n", nu, nr);
    
    for (i = 0, j = nu; i < nu && j < nu + nr; i++, j++) {
        if (obs[i].sat < obs[j].sat)
            j--;
        else if (obs[i].sat > obs[j].sat)
            i--;
        else if (azel[1 + j * 2] >= opt->elmin)
        { /* elevation at base station */
            sat[k] = obs[i].sat;
            iu[k] = i;
            ir[k++] = j;
            RTKtrace(4, "(%2d) sat=%3d iu=%2d ir=%2d\n", k - 1, obs[i].sat, i, j);
        }
    }
    return k;
}

/* temporal update of position/velocity/acceleration ------------------------- */
static void udpos(rtk_t *rtk, double tt) {
    double *F, *FP, *xp, pos[3], Q[9] = {0}, Qv[9], var = 0.0;
    int i, j;
    
    RTKtrace(3, "udpos   : tt=%.3f\n", tt);
    
    /* fixed mode */
    if (rtk->opt.mode == PMODE_FIXED) {
        for (i = 0; i < 3; i++)
            initx(rtk, rtk->opt.ru[i], 1E-8, i);
        return;
    }
    /* initialize position for first epoch */
    if (norm(rtk->x, 3) <= 0.0) {
        for (i = 0; i < 3; i++)
            initx(rtk, rtk->sol.rr[i], VAR_POS, i);
        if (rtk->opt.dynamics) {
            for (i = 3; i < 6; i++)
                initx(rtk, rtk->sol.rr[i], VAR_VEL, i);
            for (i = 6; i < 9; i++)
                initx(rtk, 1E-6, VAR_ACC, i);
        }
    }
    /* static mode */
    if (rtk->opt.mode == PMODE_STATIC)
        return;
    
    /* kinmatic mode without dynamics */
    if (!rtk->opt.dynamics) {
        for (i = 0; i < 3; i++)
            initx(rtk, rtk->sol.rr[i], VAR_POS, i);
        return;
    }
    /* check variance of estimated postion */
    for (i = 0; i < 3; i++)
        var += rtk->P[i + i * rtk->nx];
    var /= 3.0;
    
    if (var > VAR_POS) {
        /* reset position with large variance */
        for (i = 0; i < 3; i++)
            initx(rtk, rtk->sol.rr[i], VAR_POS, i);
        for (i = 3; i < 6; i++)
            initx(rtk, rtk->sol.rr[i], VAR_VEL, i);
        for (i = 6; i < 9; i++)
            initx(rtk, 1E-6, VAR_ACC, i);
        RTKtrace(2, "reset rtk position due to large variance: var=%.3f\n", var);
        return;
    }
    /* state transition of position/velocity/acceleration */
    F = eye(rtk->nx);
    FP = mat(rtk->nx, rtk->nx);
    xp = mat(rtk->nx, 1);
    
    for (i = 0; i < 6; i++) {
        F[i + (i + 3) * rtk->nx] = tt;
    }
    /* x=F*x, P=F*P*F+Q */
    matmul("NN", rtk->nx, 1, rtk->nx, 1.0, F, rtk->x, 0.0, xp);
    matcpy(rtk->x, xp, rtk->nx, 1);
    matmul("NN", rtk->nx, rtk->nx, rtk->nx, 1.0, F, rtk->P, 0.0, FP);
    matmul("NT", rtk->nx, rtk->nx, rtk->nx, 1.0, FP, F, 0.0, rtk->P);
    
    /* process noise added to only acceleration */
    Q[0] = Q[4] = SQR(rtk->opt.prn[3]);
    Q[8] = SQR(rtk->opt.prn[4]);
    ecef2pos(rtk->x, pos);
    covecef(pos, Q, Qv);
    for (i = 0; i < 3; i++)
        for (j = 0; j < 3; j++) {
            rtk->P[i + 6 + (j + 6) * rtk->nx] += Qv[i + j * 3];
        }
    free(F);
    free(FP);
    free(xp);
}

/* temporal update of ionospheric parameters --------------------------------- */
static void udion(rtk_t *rtk, double tt, double bl, const int *sat, int ns) {
    double el, fact;
    int i, j;
    
    RTKtrace(3, "udion   : tt=%.1f bl=%.0f ns=%d\n", tt, bl, ns);
    
    for (i = 1; i <= MAXSAT; i++) {
        j = II(i, &rtk->opt);
        if (rtk->x[j] != 0.0 && rtk->ssat[i - 1].outc[0] >
            GAP_RESION && rtk->ssat[i - 1].outc[1] > GAP_RESION)
            rtk->x[j] = 0.0;
    }
    for (i = 0; i < ns; i++) {
        j = II(sat[i], &rtk->opt);
        
        if (rtk->x[j] == 0.0) {
            initx(rtk, 1E-6, SQR(rtk->opt.std[1]*bl / 1E4), j);
            /* initx(rtk,1E-6,SQR(rtk->opt.std[1]*bl/10),j);*/
        }
        else {
            /* elevation dependent factor of process noise */
            el = rtk->ssat[sat[i] - 1].azel[1];
            fact = cos(el);
            rtk->P[j + j * rtk->nx] +=
            SQR(rtk->opt.prn[1] * bl / 1E4 * fact) * tt;
        }
    }
}

/* temporal update of tropospheric parameters -------------------------------- */
static void udtrop(rtk_t *rtk, double tt, double bl) {
    int i, j, k;
    
    RTKtrace(3, "udtrop  : tt=%.1f\n", tt);
    
    for (i = 0; i < 2; i++) {
        j = IT(i, &rtk->opt);
        
        if (rtk->x[j] == 0.0) {
            initx(rtk, INIT_ZWD, SQR(rtk->opt.std[2]), j); /* initial zwd */
            
            if (rtk->opt.tropopt >= TROPOPT_ESTG) {
                for (k = 0; k < 2; k++)
                    initx(rtk, 1E-6, VAR_GRA, ++j);
            }
        }
        else {
            rtk->P[j + j * rtk->nx] += SQR(rtk->opt.prn[2]) * tt;
            
            if (rtk->opt.tropopt >= TROPOPT_ESTG) {
                for (k = 0; k < 2; k++) {
                    rtk->P[++j * (1 + rtk->nx)] +=
                    SQR(rtk->opt.prn[2] * 0.3) * fabs(rtk->tt);
                }
            }
        }
    }
}

/* temporal update of receiver h/w biases ------------------------------------ */
static void udrcvbias(rtk_t *rtk, double tt) {
    int i, j;
    
    RTKtrace(3, "udrcvbias: tt=%.1f\n", tt);
    
    for (i = 0; i < NFREQGLO; i++) {
        j = IL(i, &rtk->opt);
        
        if (rtk->x[j] == 0.0) {
            initx(rtk, 1E-6, VAR_HWBIAS, j);
        }
        /* hold to fixed solution */
        else if (rtk->nfix >= rtk->opt.minfix && rtk->sol.ratio >
                 rtk->opt.thresar[0]) {
            initx(rtk, rtk->xa[j], rtk->Pa[j + j*rtk->na], j);
        }
        else {
            rtk->P[j + j * rtk->nx] += SQR(PRN_HWBIAS) * tt;
        }
    }
}

/* detect cycle slip by LLI -------------------------------------------------- */
static void detslp_ll(rtk_t *rtk, const obsd_t *obs, int i, int rcv) {
    unsigned char slip, LLI1, LLI2, LLI;
    int f, sat = obs[i].sat;
    //rtk->ssat[sat - 1].slip[f] 的1-2位存储移动站的周跳状态，3-4位存储基准站的周跳状态。
    RTKtrace(3, "detslp_ll: i=%d rcv=%d\n", i, rcv);
    
    for (f = 0; f < rtk->opt.nf; f++) {
        
        if (obs[i].L[f] == 0.0)
            continue;
        
        /* restore previous LLI */ //slip的初值为 1111 1100
        LLI1 = (rtk->ssat[sat - 1].slip[f] >> 6) & 3;//1111 1100>>6=0000 0011; 0000 1100 & 0000 0011 = 0000 0000; 主要看最后两位
        LLI2 = (rtk->ssat[sat - 1].slip[f] >> 4) & 3;//1111 1100>>4=0000 1111; 0000 1111 & 0000 0011 = 0000 0011; 主要看倒数第3-4位
        LLI = rcv == 1 ? LLI1 : LLI2;//如果是移动站（rcv==1），则LLI=LLI1；如果是基准站，则LLI=LLI2.
        
        /* detect slip by cycle slip flag *///取出obs[i].LLI[f]中的最后两bit位，与slip[f]进行或运算，赋给slip
        slip = (rtk->ssat[sat - 1].slip[f] | obs[i].LLI[f]) & 3; //obs[i].LLI[f]中的最后两位bit保存的是LLI的状态，和3与运算之后得到其最后两位bit位；再与slip[f]进行或运算，不影响slip[f]中的原来状态（无周跳时初值为1111 1100），将其最后两bit位赋值。
        
        if (obs[i].LLI[f] & 1) {//如果LLI[f]的最后一个bit位为1，则判定有周跳
            errmsg(rtk, "slip detected (sat=%2d rcv=%d LLI%d=%x)\n", sat, rcv,
                   f + 1, obs[i].LLI[f]);
        }
        /* detect slip by parity unknown flag transition */
        if (((LLI & 2) && !(obs[i].LLI[f] & 2)) ||
            (!(LLI & 2) && (obs[i].LLI[f] & 2))) {
            errmsg(rtk, "slip detected (sat=%2d rcv=%d LLI%d=%x->%x)\n", sat,
                   rcv, f + 1, LLI, obs[i].LLI[f]);
            slip |= 1;//奇偶位检验不通过，判定有slip。
        }
        /* save current LLI and slip flag */
        if (rcv == 1)
            rtk->ssat[sat - 1].slip[f] =
            (obs[i].LLI[f] << 6) | (LLI2 << 4) | slip;//slip[f]中八个bit位存储的分别为：1-2位：obs[i].LLI[f]的最后两位；3-4位：LLI2的最后两位（LLI2的前6位都为0）；5-8位：slip的后四位
        else
            rtk->ssat[sat - 1].slip[f] =
            (obs[i].LLI[f] << 4) | (LLI1 << 6) | slip;
    }
}

/* detect cycle slip by L1-L2 geometry free phase jump ----------------------- */
static void detslp_gf_L1L2(rtk_t *rtk, const obsd_t *obs, int i, int j,
                           const nav_t *nav) {
    int sat = obs[i].sat;
    double g0, g1;
    
    RTKtrace(3, "detslp_gf_L1L2: i=%d j=%d\n", i, j);
    
    if (rtk->opt.nf <= 1 || (g1 = gfobs_L1L2(obs, i, j,
                                             nav->lam[sat - 1])) == 0.0)
        return;
    //g1=SDL1*LAM1-SDL2*LAM2;
    g0 = rtk->ssat[sat - 1].gf;
    rtk->ssat[sat - 1].gf = g1;
    //前后历元的gf相减
    if (g0 != 0.0 && fabs(g1 - g0) > rtk->opt.thresslip) {//threshold=0.05
        
        rtk->ssat[sat - 1].slip[0] |= 1;//slip中的最后两bit位代表周跳状态，有周跳位1，无周跳为0
        rtk->ssat[sat - 1].slip[1] |= 1;
        
        errmsg(rtk, "slip detected (sat=%2d GF_L1_L2=%.3f %.3f)\n", sat,
               g0, g1);
    }
}

/* detect cycle slip by L1-L5 geometry free phase jump ----------------------- */
static void detslp_gf_L1L5(rtk_t *rtk, const obsd_t *obs, int i, int j,
                           const nav_t *nav) {
    int sat = obs[i].sat;
    double g0, g1;
    
    RTKtrace(3, "detslp_gf_L1L5: i=%d j=%d\n", i, j);
    
    if (rtk->opt.nf <= 2 || (g1 = gfobs_L1L5(obs, i, j,
                                             nav->lam[sat - 1])) == 0.0)
        return;
    
    g0 = rtk->ssat[sat - 1].gf2;
    rtk->ssat[sat - 1].gf2 = g1;
    
    if (g0 != 0.0 && fabs(g1 - g0) > rtk->opt.thresslip) {
        
        rtk->ssat[sat - 1].slip[0] |= 1;
        rtk->ssat[sat - 1].slip[2] |= 1;
        
        errmsg(rtk, "slip detected (sat=%2d GF_L1_L5=%.3f %.3f)\n", sat,
               g0, g1);
    }
}

/* detect cycle slip by doppler and phase difference ------------------------- */
static void detslp_dop(rtk_t *rtk, const obsd_t *obs, int i, int rcv,
                       const nav_t *nav) {
    /* detection with doppler disabled because of clock-jump issue (v.2.3.0) */
#if 0
    int f, sat = obs[i].sat;
    double tt, dph, dpt, lam, thres;
    
    RTKtrace(3, "detslp_dop: i=%d rcv=%d\n", i, rcv);
    
    for (f = 0; f < rtk->opt.nf; f++) {
        if (obs[i].L[f] == 0.0 || obs[i].D[f] == 0.0 ||
            rtk->ph[rcv - 1][sat - 1][f] == 0.0) {
            continue;
        }
        if (fabs(tt = timediff(obs[i].time, rtk->pt[rcv - 1][sat - 1][f]))
            < DTTOL)
            continue;
        if ((lam = nav->lam[sat - 1][f]) <= 0.0)
            continue;
        
        /* cycle slip threshold (cycle) */
        thres = MAXACC * tt * tt / 2.0 / lam + rtk->opt.err[4] * fabs(tt) * 4.0;
        
        /* phase difference and doppler x time (cycle) */
        dph = obs[i].L[f] - rtk->ph[rcv - 1][sat - 1][f];
        dpt = -obs[i].D[f] * tt;
        
        if (fabs(dph - dpt) <= thres)
            continue;
        
        rtk->slip[sat - 1][f] |= 1;
        
        errmsg(rtk, "slip detected (sat=%2d rcv=%d L%d=%.3f %.3f thres=%.3f)\n",
               sat, rcv, f + 1, dph, dpt, thres);
    }
#endif
}

/* temporal update of phase biases ------------------------------------------- */
static void udbias(rtk_t *rtk, double tt, const obsd_t *obs, const int *sat,
                   const int *iu, const int *ir, int ns, const nav_t *nav) {
    double cp, pr, cp1, cp2, pr1, pr2, *bias, offset, lami, lam1, lam2, C1, C2;
    int i, j, f, slip, reset, nf = NF(&rtk->opt);
    int debugtmp=0;
    //TODO: udbias
    RTKtrace(3, "udbias  : tt=%.1f ns=%d\n", tt, ns);
    
    for (i = 0; i < ns; i++) {
        
        /* detect cycle slip by LLI */
        for (f = 0; f < rtk->opt.nf; f++)
            rtk->ssat[sat[i] - 1].slip[f] &= 0xFC;//1111 1100
        //detslp_ll(rtk, obs, iu[i], 1);
        detslp_ll(rtk, obs, ir[i], 2);
        
        /* detect cycle slip by geometry-free phase jump */
        detslp_gf_L1L2(rtk, obs, iu[i], ir[i], nav);
        detslp_gf_L1L5(rtk, obs, iu[i], ir[i], nav);
        
        /* detect cycle slip by doppler and phase difference *///目前被disable了
        detslp_dop(rtk, obs, iu[i], 1, nav);
        detslp_dop(rtk, obs, ir[i], 2, nav);
    }
    for (f = 0; f < nf; f++) {
        /* reset phase-bias if instantaneous AR or expire obs outage counter */
        for (i = 1; i <= MAXSAT; i++) {
            
            reset = ++rtk->ssat[i - 1].outc[f] > (unsigned int)rtk->opt.maxout;
            debugtmp=IB(i, f, &rtk->opt);
            if (rtk->opt.modear == ARMODE_INST && rtk->x[IB(i, f,
                                                            &rtk->opt)] != 0.0) {
                initx(rtk, 0.0, 0.0, IB(i, f, &rtk->opt));
            }
            else if (reset && rtk->x[IB(i, f, &rtk->opt)] != 0.0) {
                initx(rtk, 0.0, 0.0, IB(i, f, &rtk->opt));
                RTKtrace(3,
                         "udbias : obs outage counter overflow (sat=%3d L%d n=%d)\n",
                         i, f + 1, rtk->ssat[i - 1].outc[f]);
            }
            if (rtk->opt.modear != ARMODE_INST && reset) {
                rtk->ssat[i - 1].lock[f] = -rtk->opt.minlock;
            }
        }
        /* reset phase-bias if detecting cycle slip */
        for (i = 0; i < ns; i++) {
            j = IB(sat[i], f, &rtk->opt);
            rtk->P[j + j * rtk->nx] += rtk->opt.prn[0] * rtk->opt.prn[0] * tt;
            slip = rtk->ssat[sat[i] - 1].slip[f];
            if (rtk->opt.ionoopt == IONOOPT_IFLC)
                slip |= rtk->ssat[sat[i] - 1].slip[1];
            if (rtk->opt.modear == ARMODE_INST || !(slip & 1))
                continue;
            rtk->x[j] = 0.0;
            rtk->ssat[sat[i] - 1].lock[f] = -rtk->opt.minlock;
        }
        bias = zeros(ns, 1);
        
        /* estimate approximate phase-bias by phase - code */
        for (i = j = 0, offset = 0.0; i < ns; i++) {
            
            if (rtk->opt.ionoopt != IONOOPT_IFLC) {
                cp = sdobs(obs, iu[i], ir[i], f); /* cycle */
                pr = sdobs(obs, iu[i], ir[i], f + NFREQ);
                lami = nav->lam[sat[i] - 1][f];
                if (cp == 0.0 || pr == 0.0 || lami <= 0.0)
                    continue;
                
                bias[i] = cp - pr / lami;
            }
            else {
                cp1 = sdobs(obs, iu[i], ir[i], 0);
                cp2 = sdobs(obs, iu[i], ir[i], 1);
                pr1 = sdobs(obs, iu[i], ir[i], NFREQ);
                pr2 = sdobs(obs, iu[i], ir[i], NFREQ + 1);
                lam1 = nav->lam[sat[i] - 1][0];
                lam2 = nav->lam[sat[i] - 1][1];
                if (cp1 == 0.0 || cp2 == 0.0 || pr1 == 0.0 || pr2 == 0.0 ||
                    lam1 <= 0.0 || lam2 <= 0.0)
                    continue;
                
                C1 = SQR(lam2) / (SQR(lam2) - SQR(lam1));
                C2 = -SQR(lam1) / (SQR(lam2) - SQR(lam1));
                bias[i] = (C1 * lam1 * cp1 + C2 * lam2 * cp2) -
                (C1 * pr1 + C2 * pr2);
            }
            if (rtk->x[IB(sat[i], f, &rtk->opt)] != 0.0) {
                offset += bias[i] - rtk->x[IB(sat[i], f, &rtk->opt)];
                j++;
            }
        }
        /* correct phase-bias offset to enssure phase-code coherency */
        if (j > 0) {
            for (i = 1; i <= MAXSAT; i++) {
                if (rtk->x[IB(i, f, &rtk->opt)] != 0.0)
                    rtk->x[IB(i, f, &rtk->opt)] += offset / j;
            }
        }
        //         RTKtrace(1,"bias=");tracemat(1,bias,1,ns,13,4);
        /* set initial states of phase-bias */
        for (i = 0; i < ns; i++) {
            if (bias[i] == 0.0 || rtk->x[IB(sat[i], f, &rtk->opt)] != 0.0){
                continue;
            }
            initx(rtk, bias[i], SQR(rtk->opt.std[0]), IB(sat[i], f, &rtk->opt));
        }
        free(bias);
    }
}

/* temporal update of states -------------------------------------------------- */
static void udstate(rtk_t *rtk, const obsd_t *obs, const int *sat,
                    const int *iu, const int *ir, int ns, const nav_t *nav) {
    double tt = fabs(rtk->tt), bl=0.0, dr[3];
    //TODO: udstate
    RTKtrace(3, "udstate : ns=%d\n", ns);
    
    /* temporal update of position/velocity/acceleration */
    udpos(rtk, tt);
    
    /* temporal update of ionospheric parameters */
    if (rtk->opt.ionoopt >= IONOOPT_EST) {
        bl = baseline(rtk->x, rtk->rb, dr);
        udion(rtk, tt, bl, sat, ns);
    }
    /* temporal update of tropospheric parameters */
    if (rtk->opt.tropopt >= TROPOPT_EST) {
        udtrop(rtk, tt, bl);
    }
    /* temporal update of eceiver h/w bias */
    if (rtk->opt.glomodear == 2 && (rtk->opt.navsys & SYS_GLO)) {
        udrcvbias(rtk, tt);
    }
    /* temporal update of phase-bias */
    if (rtk->opt.mode > PMODE_DGPS) {
        udbias(rtk, tt, obs, sat, iu, ir, ns, nav);
    }
}

/* undifferenced phase/code residual for satellite --------------------------- */
static void zdres_sat(int base, double r, const obsd_t *obs, const nav_t *nav,
                      const double *azel, const double *dant, const prcopt_t *opt, double *y) {
    const double *lam = nav->lam[obs->sat - 1];
    double f1, f2, C1, C2, dant_if;
    int i, nf = NF(opt);
    
    if (opt->ionoopt == IONOOPT_IFLC) { /* iono-free linear combination */
        if (lam[0] == 0.0 || lam[1] == 0.0)
            return;
        
        if (testsnr(base, 0, azel[1], obs->SNR[0] * 0.25, &opt->snrmask) ||
            testsnr(base, 1, azel[1], obs->SNR[1] * 0.25, &opt->snrmask))
            return;
        
        f1 = CLIGHT / lam[0];
        f2 = CLIGHT / lam[1];
        C1 = SQR(f1) / (SQR(f1) - SQR(f2));
        C2 = -SQR(f2) / (SQR(f1) - SQR(f2));
        dant_if = C1 * dant[0] + C2 * dant[1];
        
        if (obs->L[0] != 0.0 && obs->L[1] != 0.0) {
            y[0] = C1 * obs->L[0] * lam[0] + C2 * obs->L[1] * lam[1] - r -
            dant_if;
        }
        if (obs->P[0] != 0.0 && obs->P[1] != 0.0) {
            y[1] = C1 * obs->P[0] + C2 * obs->P[1] - r - dant_if;
        }
    }
    else {
        for (i = 0; i < nf; i++) {
            if (lam[i] == 0.0)
                continue;
            
            /* check snr mask */
            if (testsnr(base, i, azel[1], obs->SNR[i] * 0.25, &opt->snrmask)) {
                continue;
            }
            /* residuals = observable - pseudorange */
            if (obs->L[i] != 0.0)
                y[i] = obs->L[i] * lam[i] - r - dant[i];
            if (obs->P[i] != 0.0)
                y[i + nf] = obs->P[i] - r - dant[i];
        }
    }
}

/* undifferenced phase/code residuals ---------------------------------------- */
static int zdres(int base, const obsd_t *obs, int n, const double *rs,
                 const double *dts, const int *svh, const nav_t *nav, const double *rr,
                 const prcopt_t *opt, int index, double *y, double *e, double *azel) {
    double r, rr_[3], pos[3], dant[NFREQ] = {0}, disp[3];
    double zhd, zazel[] = {0.0, 90.0 * D2R};
    int i, nf = NF(opt);
    
    RTKtrace(3, "zdres   : n=%d\n", n);
    
    for (i = 0; i < n * nf * 2; i++)
        y[i] = 0.0;
    
    if (norm(rr, 3) <= 0.0)
        return 0; /* no receiver position */
    
    for (i = 0; i < 3; i++)
        rr_[i] = rr[i];
    
    /* earth tide correction */
    if (opt->tidecorr) {
        tidedisp(gpst2utc(obs[0].time), rr_, opt->tidecorr, &nav->erp,
                 opt->odisp[base], disp);
        for (i = 0; i < 3; i++)
            rr_[i] += disp[i];
    }
    ecef2pos(rr_, pos);
    
    for (i = 0; i < n; i++) {
        /* compute geometric-range and azimuth/elevation angle */
        if ((r = geodist(rs + i * 6, rr_, e + i * 3)) <= 0.0)
            continue;
        if (satazel(pos, e + i * 3, azel + i * 2) < opt->elmin)
            continue;
        
        /* excluded satellite? */
        if (satexclude(obs[i].sat, svh[i], opt))
            continue;
        
        /* satellite clock-bias */
        r += -CLIGHT * dts[i * 2];
        
        /* troposphere delay model (hydrostatic) */
        zhd = tropmodel(obs[0].time, pos, zazel, 0.0);
        r += tropmapf(obs[i].time, pos, azel + i * 2, NULL) * zhd;
        
        /* receiver antenna phase center correction */
        antmodel(opt->pcvr + index, opt->antdel[index], azel + i*2,
                 opt->posopt[1], dant);
        
        /* undifferenced phase/code residual for satellite */
        zdres_sat(base, r, obs + i, nav, azel + i*2, dant, opt, y + i*nf*2);
    }
    RTKtrace(4, "rr_=%.3f %.3f %.3f\n", rr_[0], rr_[1], rr_[2]);
    RTKtrace(4, "pos=%.9f %.9f %.3f\n", pos[0]*R2D, pos[1]*R2D, pos[2]);
    for (i = 0; i < n; i++) {
        RTKtrace(4, "sat=%2d %13.3f %13.3f %13.3f %13.10f %6.1f %5.1f\n",
                 obs[i].sat, rs[i*6], rs[1 + i*6], rs[2 + i*6], dts[i*2],
                 azel[i*2]*R2D, azel[1 + i*2]*R2D);
    }
    //    RTKtrace(1, "y=\n");
    //    tracemat(1, y, nf*2, n, 15, 9);
    
    return 1;
}

/* test valid observation data ----------------------------------------------- */
static int validobs(int i, int j, int f, int nf, double *y) {
    /* if no phase observable, psudorange is also unusable */
    return y[f + i * nf * 2] != 0.0 && y[f + j * nf * 2] != 0.0 &&
    (f < nf || (y[f - nf + i * nf * 2] != 0.0 && y[f - nf +
                                                   j * nf * 2] != 0.0));
}

/* double-differenced measurement error covariance --------------------------- */
static void ddcov(const int *nb, int n, const double *Ri, const double *Rj,
                  int nv, double *R) {
    int i, j, k = 0, b;
    
    RTKtrace(3, "ddcov   : n=%d\n", n);
    
    for (i = 0; i < nv * nv; i++)
        R[i] = 0.0;
    for (b = 0; b < n; k += nb[b++]) {
        
        for (i = 0; i < nb[b]; i++)
            for (j = 0; j < nb[b]; j++) {
                R[k + i + (k + j) * nv] =
                Ri[k + i] + (i == j ? Rj[k + i] : 0.0);
            }
    }
    RTKtrace(5, "R=\n");
    tracemat(5, R, nv, nv, 8, 6);
}

/* baseline length constraint ------------------------------------------------ */
static int constbl(rtk_t *rtk, const double *x, const double *P, double *v,
                   double *H, double *Ri, double *Rj, int index) {
    const double thres = 0.1; /* threshold for nonliearity (v.2.3.0) */
    double xb[3], b[3], bb, var = 0.0;
    int i;
    
    RTKtrace(3, "constbl : \n");
    
    /* no constraint */
    if (rtk->opt.baseline[0] <= 0.0)
        return 0;
    
    /* time-adjusted baseline vector and length */
    for (i = 0; i < 3; i++) {
        xb[i] = rtk->rb[i] + rtk->rb[i + 3] * rtk->sol.age;
        b[i] = x[i] - xb[i];
    }
    bb = norm(b, 3);
    
    /* approximate variance of solution */
    if (P) {
        for (i = 0; i < 3; i++)
            var += P[i + i * rtk->nx];
        var /= 3.0;
    }
    /* check nonlinearity */
    if (var > thres * thres * bb * bb) {
        RTKtrace(3, "constbl : equation nonlinear (bb=%.3f var=%.3f)\n", bb, var);
        return 0;
    }
    /* constraint to baseline length */
    v[index] = rtk->opt.baseline[0] - bb;
    if (H) {
        for (i = 0; i < 3; i++)
            H[i + index * rtk->nx] = b[i] / bb;
    }
    Ri[index] = 0.0;
    Rj[index] = SQR(rtk->opt.baseline[1]);
    
    RTKtrace(4, "baseline len   v=%13.3f R=%8.6f %8.6f\n", v[index], Ri[index],
             Rj[index]);
    
    return 1;
}

/* precise tropspheric model ------------------------------------------------- */
static double prectrop(gtime_t time, const double *pos, int r,
                       const double *azel, const prcopt_t *opt, const double *x, double *dtdx) {
    double m_w = 0.0, cotz, grad_n, grad_e;
    int i = IT(r, opt);
    
    /* wet mapping function */
    tropmapf(time, pos, azel, &m_w);
    
    if (opt->tropopt >= TROPOPT_ESTG && azel[1] > 0.0) {
        
        /* m_w=m_0+m_0*cot(el)*(Gn*cos(az)+Ge*sin(az)): ref [6] */
        cotz = 1.0 / tan(azel[1]);
        grad_n = m_w * cotz * cos(azel[0]);
        grad_e = m_w * cotz * sin(azel[0]);
        m_w += grad_n * x[i + 1] + grad_e * x[i + 2];
        dtdx[1] = grad_n * x[i];
        dtdx[2] = grad_e * x[i];
    }
    else
        dtdx[1] = dtdx[2] = 0.0;
    dtdx[0] = m_w;
    return m_w*x[i];
}

/* glonass inter-channel bias correction ------------------------------------- */
static double gloicbcorr(int sat1, int sat2, const prcopt_t *opt, double lam1,
                         double lam2, int f) {
    double dfreq;
    
    if (f >= NFREQGLO || f >= opt->nf || !opt->exterr.ena[2])
        return 0.0;
    
    dfreq = (CLIGHT / lam1 - CLIGHT / lam2) / (f == 0 ? DFRQ1_GLO : DFRQ2_GLO);
    
    return opt->exterr.gloicb[f] * 0.01 * dfreq; /* (m) */
}

/* test navi system (m=0:gps/qzs/sbs,1:glo,2:gal,3:bds) ---------------------- */
static int test_sys(int sys, int m) {
    switch (sys) {
        case SYS_GPS:
            return m == 0;
        case SYS_QZS:
            return m == 0;
        case SYS_SBS:
            return m == 0;
        case SYS_GLO:
            return m == 1;
        case SYS_GAL:
            return m == 2;
        case SYS_CMP:
            return m == 3;
    }
    return 0;
}

/* double-differenced phase/code residuals -----------------------------------
 该函数计算双差观测值，以及双差-单差设计矩阵：即
 yDD=H*xSD+\epsilon
 其中，H左边对应双差观测值，右侧对应单差未知数。*/
static int ddres(rtk_t *rtk, const nav_t *nav, double dt, const double *x,
                 const double *P, const int *sat, double *y, double *e, double *azel,
                 const int *iu, const int *ir, int ns, double *v, double *H, double *R,
                 int *vflg) {
    //TODO:ddres
    prcopt_t *opt = &rtk->opt;
    double bl, dr[3], posu[3], posr[3], didxi = 0.0, didxj = 0.0, *im;
    double *tropr, *tropu, *dtdxr, *dtdxu, *Ri, *Rj, s, lami, lamj, fi, fj, df,
    *Hi = NULL;
    int i, j, k, m, f, ff, nv = 0, nb[NFREQ * 4 * 2 + 2] = {0
    }, b = 0, sysi, sysj, nf = NF(opt);
    int mode=0;
    int tmpi=0, tmpj=0;
    double debugtmp1=0.0, debugtmp2=0.0, tmpv=0.0;
    RTKtrace(3, "ddres   : dt=%.1f nx=%d ns=%d\n", dt, rtk->nx, ns);
    
    bl = baseline(x, rtk->rb, dr);
    ecef2pos(x, posu);
    ecef2pos(rtk->rb, posr);
    
    Ri = mat(ns * nf * 2 + 2, 1);
    Rj = mat(ns * nf * 2 + 2, 1);
    im = mat(ns, 1);
    tropu = mat(ns, 1);
    tropr = mat(ns, 1);
    dtdxu = mat(ns, 3);
    dtdxr = mat(ns, 3);
    /*将所有卫星的所有频率伪距和载波残差初始化为0*/
    for (i = 0; i < MAXSAT; i++)
        for (j = 0; j < NFREQ; j++) {
            rtk->ssat[i].resp[j] = rtk->ssat[i].resc[j] = 0.0;
        }
    /* compute factors of ionospheric and tropospheric delay */
    for (i = 0; i < ns; i++) {
        if (opt->ionoopt >= IONOOPT_EST) {/*如果估计电离层误差，或者估计电离层误差及其梯度，则计算投影系数*/
            im[i] = (ionmapf(posu, azel + iu[i] * 2) + ionmapf(posr,
                                                               azel + ir[i] * 2)) / 2.0;
        }
        if (opt->tropopt >= TROPOPT_EST) {/*如果估计对流层误差，或者估计对流层误差及其梯度，则计算投影系数*/
            tropu[i] = prectrop(rtk->sol.time, posu, 0, azel + iu[i] * 2, opt,
                                x, dtdxu + i * 3);
            tropr[i] = prectrop(rtk->sol.time, posr, 1, azel + ir[i] * 2, opt,
                                x, dtdxr + i * 3);
        }
    }
    for (m = 0; m < 4; m++) /* m=0:gps/qzs/sbs,1:glo,2:gal,3:bds */
        
        for (f = opt->mode > PMODE_DGPS ? 0 : nf; f < nf * 2; f++) {
            
            /* search reference satellite with highest elevation */
            for (i = -1, j = 0; j < ns; j++) {
                sysi = rtk->ssat[sat[j] - 1].sys;
                if (!test_sys(sysi, m))
                    continue;
                if (!validobs(iu[j], ir[j], f, nf, y))
                    continue;
                if (i < 0 || azel[1 + iu[j]*2] >= azel[1 + iu[i]*2])
                    i = j;
            }
            if (i < 0)
                continue;
            
            /* make double difference */
            for (j = 0; j < ns; j++) {
                if (i == j)
                    continue;
                sysi = rtk->ssat[sat[i] - 1].sys;
                sysj = rtk->ssat[sat[j] - 1].sys;
                if (!test_sys(sysj, m))
                    continue;
                if (!validobs(iu[j], ir[j], f, nf, y))
                    continue;
                
                ff = f % nf;
                lami = nav->lam[sat[i] - 1][ff];
                lamj = nav->lam[sat[j] - 1][ff];
                if (lami <= 0.0 || lamj <= 0.0)
                    continue;
                if (H)
                    Hi = H + nv * rtk->nx;/*卫星i的LOS向量*/
                
                /* double-differenced residual */
//                tmpi=f+iu[i]*nf*2;
                v[nv] = (y[f + iu[i] * nf * 2] - y[f + ir[i] * nf * 2]) -
                (y[f + iu[j] * nf * 2] - y[f + ir[j] * nf * 2]);
                
                /* partial derivatives by ====ROVER==== position */
                if (H) {
                    for (k = 0; k < 3; k++) {
                        Hi[k] = -e[k + iu[i] * 3] + e[k + iu[j] * 3];
                    }
                }
                /* double-differenced ionospheric delay term */
                if (opt->ionoopt == IONOOPT_EST) {
                    fi = lami / lam_carr[0];
                    fj = lamj / lam_carr[0];
                    didxi = (f < nf ? -1.0 : 1.0) * fi * fi * im[i];
                    didxj = (f < nf ? -1.0 : 1.0) * fj * fj * im[j];
                    v[nv] -= didxi * x[II(sat[i], opt)] -
                    didxj * x[II(sat[j], opt)];
                    if (H) {
                        Hi[II(sat[i], opt)] = didxi;
                        Hi[II(sat[j], opt)] = -didxj;
                    }
                }
                /* double-differenced tropospheric delay term */
                if (opt->tropopt == TROPOPT_EST || opt->tropopt == TROPOPT_ESTG)
                {
                    v[nv] -= (tropu[i] - tropu[j]) - (tropr[i] - tropr[j]);
                    for (k = 0; k < (opt->tropopt < TROPOPT_ESTG ? 1 : 3); k++)
                    {
                        if (!H)
                            continue;
                        Hi[IT(0, opt) + k] = (dtdxu[k + i * 3] - dtdxu[k + j * 3]);
                        Hi[IT(1, opt) + k] = -(dtdxr[k + i * 3] - dtdxr[k + j * 3]);
                    }
                }
                /* double-differenced phase-bias term */
                if (f < nf) {
                    if (opt->ionoopt != IONOOPT_IFLC) {
                        tmpi=   IB(sat[i], f, opt);
                        tmpj=IB(sat[j], f, opt);
                        debugtmp1=   x[IB(sat[i], f, opt)];
                        debugtmp2=   x[IB(sat[j], f, opt)];
                        tmpv=lami * x[IB(sat[i], f, opt)] - lamj * x[IB(sat[j], f, opt)];
                        v[nv] -= lami * x[IB(sat[i], f, opt)] -
                        lamj * x[IB(sat[j], f, opt)];
                        if (H) {
                            Hi[IB(sat[i], f, opt)] = lami;
                            Hi[IB(sat[j], f, opt)] = -lamj;
                        }
                    }
                    else {
                        v[nv] -= x[IB(sat[i], f, opt)] - x[IB(sat[j], f, opt)];
                        if (H) {
                            Hi[IB(sat[i], f, opt)] = 1.0;
                            Hi[IB(sat[j], f, opt)] = -1.0;
                        }
                    }
                }
                /* glonass receiver h/w bias term */
                if (rtk->opt.glomodear == 2 && sysi == SYS_GLO && sysj ==
                    SYS_GLO && ff < NFREQGLO) {
                    df = (CLIGHT / lami - CLIGHT / lamj) / 1E6;
                    /* freq-difference (MHz) */
                    v[nv] -= df * x[IL(ff, opt)];
                    if (H)
                        Hi[IL(ff, opt)] = df;
                }
                /* glonass interchannel bias correction */
                else if (sysi == SYS_GLO && sysj == SYS_GLO) {
                    
                    v[nv] -= gloicbcorr(sat[i], sat[j], &rtk->opt, lami,
                                        lamj, f);
                }
                if (f < nf)
                    rtk->ssat[sat[j] - 1].resc[f] = v[nv];
                else
                    rtk->ssat[sat[j] - 1].resp[f - nf] = v[nv];
                
                /* test innovation */
                if (opt->maxinno > 0.0 && fabs(v[nv]) > opt->maxinno) {
                    if (f < nf) {
                        rtk->ssat[sat[i] - 1].rejc[f]++;
                        rtk->ssat[sat[j] - 1].rejc[f]++;
                    }
                    errmsg(rtk, "outlier rejected (sat=%3d-%3d %s%d v=%.3f)\n",
                           sat[i], sat[j], f < nf ? "L" : "P", f % nf + 1, v[nv]);
                    continue;
                }
                /* single-differenced measurement error variances */
                mode=0;
                //                Ri[nv] = varerrYQ(sat[i], sysi, azel[1 + iu[i] * 2], bl, dt,
                //                    f, opt,mode);
                //                Rj[nv] = varerrYQ(sat[j], sysj, azel[1 + iu[j] * 2], bl, dt,
                //                    f, opt,mode);
                Ri[nv] = varerr(sat[i], sysi, azel[1 + iu[i] * 2], bl, dt,
                                f, opt);
                Rj[nv] = varerr(sat[j], sysj, azel[1 + iu[j] * 2], bl, dt,
                                f, opt);
                /* set valid data flags */
                if (opt->mode > PMODE_DGPS) {
                    if (f < nf)
                        rtk->ssat[sat[i] - 1].vsat[f] =
                        rtk->ssat[sat[j] - 1].vsat[f] = 1;
                }
                else {
                    rtk->ssat[sat[i] - 1].vsat[f - nf] =
                    rtk->ssat[sat[j] - 1].vsat[f - nf] = 1;
                }
                RTKtrace(4, "sat=%3d-%3d %s%d v=%13.3f R=%8.6f %8.6f\n", sat[i],
                         sat[j], f < nf ? "L" : "P", f % nf + 1, v[nv], Ri[nv],
                         Rj[nv]);
                
                vflg[nv++] = (sat[i] << 16) | (sat[j] << 8) |
                ((f < nf ? 0 : 1) << 4) | (f % nf);
                nb[b]++;
            }
            /* restore single-differenced residuals assuming sum equal zero */
            if (f < nf) {
                for (j = 0, s = 0.0; j < MAXSAT; j++)
                    s += rtk->ssat[j].resc[f];
                s /= nb[b] + 1;
                for (j = 0; j < MAXSAT; j++) {
                    if (j == sat[i] - 1 || rtk->ssat[j].resc[f] != 0.0)
                        rtk->ssat[j].resc[f] -= s;
                }
            }
            else {
                for (j = 0, s = 0.0; j < MAXSAT; j++)
                    s += rtk->ssat[j].resp[f - nf];
                s /= nb[b] + 1;
                for (j = 0; j < MAXSAT; j++) {
                    if (j == sat[i] - 1 || rtk->ssat[j].resp[f - nf] != 0.0)
                        rtk->ssat[j].resp[f - nf] -= s;
                }
            }
            b++;
        }
    /* end of system loop */
    
    /* baseline length constraint for moving baseline */
    if (opt->mode == PMODE_MOVEB && constbl(rtk, x, P, v, H, Ri, Rj, nv)) {
        vflg[nv++] = 3 << 4;
        nb[b++]++;
    }
    if (H) {
        RTKtrace(5, "H=\n");
        tracemat(5, H, rtk->nx, nv, 7, 4);
    }
    
    /* double-differenced measurement error covariance */
    ddcov(nb, b, Ri, Rj, nv, R);
    
    free(Ri);
    free(Rj);
    free(im);
    free(tropu);
    free(tropr);
    free(dtdxu);
    free(dtdxr);
    
    return nv;
}

/* time-interpolation of residuals (for post-mission) ------------------------ */
static double intpres(gtime_t time, const obsd_t *obs, int n, const nav_t *nav,
                      rtk_t *rtk, double *y) {
    static obsd_t obsb[MAXOBS];
    static double yb[MAXOBS * NFREQ * 2], rs[MAXOBS * 6], dts[MAXOBS * 2],
    var[MAXOBS];
    static double e[MAXOBS * 3], azel[MAXOBS * 2];
    static int nb = 0, svh[MAXOBS * 2];
    prcopt_t *opt = &rtk->opt;
    double tt = timediff(time, obs[0].time), ttb, *p, *q;
    int i, j, k, nf = NF(opt);
    
    RTKtrace(3, "intpres : n=%d tt=%.1f\n", n, tt);
    
    if (nb == 0 || fabs(tt) < DTTOL) {
        nb = n;
        for (i = 0; i < n; i++)
            obsb[i] = obs[i];
        return tt;
    }
    ttb = timediff(time, obsb[0].time);
    if (fabs(ttb) > opt->maxtdiff*2.0 || ttb == tt)
        return tt;
    
    satposs(time, obsb, nb, nav, opt->sateph, rs, dts, var, svh);
    
    if (!zdres(1, obsb, nb, rs, dts, svh, nav, rtk->rb, opt, 1, yb, e, azel)) {
        return tt;
    }
    for (i = 0; i < n; i++) {
        for (j = 0; j < nb; j++)
            if (obsb[j].sat == obs[i].sat)
                break;
        if (j >= nb)
            continue;
        for (k = 0, p = y + i * nf * 2, q = yb + j * nf * 2; k < nf * 2;
             k++, p++, q++) {
            if (*p == 0.0 || *q == 0.0)
                * p = 0.0;
            else
                *p = (ttb * (*p) - tt * (*q)) / (ttb - tt);
        }
    }
    return fabs(ttb) > fabs(tt) ? ttb : tt;
}
/* L1L2 ambiguity to L1/WL ambiguity, and single to double difference transformation matrix (D') ------------- */
static int yqddmat(rtk_t *rtk, double *D, double *nbs) {
    /*D matrix is the single to double difference matrix (transpose), the same to ddmat function
     WL (nx-2*m) matrix is the L1-WL transformation matrix, m is the number of available sys;
     nbs (4 x 3) is the number of ambiguities per freq per system;
     wlams(2*m) is the wavelength of the L1 and WL ambiguities ;
     The formula is : [a_L1;a_WL]=WL'*D'*rtk->x*/
    
    int i, j, k, m, f, nb = 0, nx = rtk->nx, na = rtk->na, nf = NF(&rtk->opt);
    /* nx, na---number of float states/fixed states */
    
    RTKtrace(3, "D':\n");
    for (i = 0; i < MAXSAT; i++)
        for (j = 0; j < NFREQ; j++) {
            rtk->ssat[i].fix[j] = 0;
        }
    for (i = 0; i < na; i++)
        D[i + i*nx] = 1.0;   /*the first na*na block of G is identity matrix*/
    
    for (m = 0; m < 4; m++) { /* m=0:gps/qzs/sbs,1:glo,2:gal,3:bds */
        
        if (m == 1 && rtk->opt.glomodear == 0)
            continue;
        if (m == 3 && rtk->opt.bdsmodear == 0)
            continue;
        /*Suppose there are two freqs, the order of observations are L1:sat_1, sat_2,...sat_MAXSAT, L2:sat_1,...*/
        for (f = 0, k = na; f < nf; f++, k += MAXSAT) {
            /*the i loop is to identify which satellite of sys 'm' will be included in AR*/
            /*	nbs[m+f*4]=0.0; */
            for (i = k; i < k + MAXSAT; i++) {
                if (rtk->x[i] == 0.0 || !test_sys(rtk->ssat[i - k].sys, m) ||
                    !rtk->ssat[i - k].vsat[f]) {
                    continue;
                }
                if (rtk->ssat[i - k].lock[f] >
                    0 && !(rtk->ssat[i - k].slip[f] & 2) && rtk->ssat[i -
                                                                      k].azel[1] >= rtk->opt.elmaskar) {
                        rtk->ssat[i - k].fix[f] = 2; /* fix */
                        break; /*find the reference satellite? why the first valid sat?*/
                    }
                else
                    rtk->ssat[i - k].fix[f] = 1; /* ambiguity fix flag (1:float,2:fix,3:hold) ?*/
            }/*end of i loop*/
            for (j = k; j < k + MAXSAT; j++) {
                if (i == j || rtk->x[j] == 0.0 ||
                    !test_sys(rtk->ssat[j - k].sys, m) ||
                    !rtk->ssat[j - k].vsat[f]) {
                    continue;
                } /*for sys 'm', the diagonal element is 0*/
                if (rtk->ssat[j - k].lock[f] >
                    0 && !(rtk->ssat[j - k].slip[f] & 2) && rtk->ssat[i -
                                                                      k].vsat[f] && rtk->ssat[j - k].azel[1] >= rtk->opt.elmaskar)
                {
                    D[i + (na + nb) * nx] = 1.0;
                    D[j + (na + nb) * nx] = -1.0;
                    nb++;
                    nbs[m+f*4]=nbs[m+f*4]+1.0;  /*nbs[m][f]++;*/
                    rtk->ssat[j - k].fix[f] = 2; /* fix */
                }
                else
                    rtk->ssat[j - k].fix[f] = 1;
            }/*end of j loop*/
        }/*end of f loop*/
    }/*end of m loop*/
    RTKtrace(5, "D'=\n");
    tracemat(5, D, nx, na + nb, 2, 0);
    return nb;
}
/* single to double-difference transformation matrix (D') -------------------- */
static int ddmat(rtk_t *rtk, double *D) {
    int i, j, k, m, f, nb = 0, nx = rtk->nx, na = rtk->na, nf = NF(&rtk->opt);
    /*In fact, this D matrix is the G matrix in P165 in RTKLIB manual,
     G=diag([I_{6x6},D,D,D]),
     here, na is the length of the baseline parameter vector.
     The initial value of input D is zero.
     and in fact, it is G', not G */
    
    /* nx, na---number of float states/fixed states */
    RTKtrace(3, "ddmat   :\n");
    
    for (i = 0; i < MAXSAT; i++)
        for (j = 0; j < NFREQ; j++) {
            rtk->ssat[i].fix[j] = 0;
        }
    for (i = 0; i < na; i++)
        D[i + i*nx] = 1.0;   /*the first na*na block of G is identity matrix*/
    
    for (m = 0; m < 4; m++) { /* m=0:gps/qzs/sbs,1:glo,2:gal,3:bds */
        
        if (m == 1 && rtk->opt.glomodear == 0)
            continue;
        if (m == 3 && rtk->opt.bdsmodear == 0)
            continue;
        /*Suppose there are two freqs, the order of observations are L1:sat_1, sat_2,...sat_MAXSAT, L2:sat_1,...*/
        for (f = 0, k = na; f < nf; f++, k += MAXSAT) {
            /*the i loop is to identify which satellite of sys 'm' will be included in AR*/
            for (i = k; i < k + MAXSAT; i++) {
                if (rtk->x[i] == 0.0 || !test_sys(rtk->ssat[i - k].sys, m) ||
                    !rtk->ssat[i - k].vsat[f]) {
                    continue;
                }
                if (rtk->ssat[i - k].lock[f] >
                    0 && !(rtk->ssat[i - k].slip[f] & 2) && rtk->ssat[i -
                                                                      k].azel[1] >= rtk->opt.elmaskar) {
                        rtk->ssat[i - k].fix[f] = 2; /* fix */
                        break; /*find the reference satellite? why the first valid sat?*/
                    }
                else
                    rtk->ssat[i - k].fix[f] = 1; /* ambiguity fix flag (1:float,2:fix,3:hold) ?*/
            }/*end of i loop*/
            for (j = k; j < k + MAXSAT; j++) {
                if (i == j || rtk->x[j] == 0.0 ||
                    !test_sys(rtk->ssat[j - k].sys, m) ||
                    !rtk->ssat[j - k].vsat[f]) {
                    continue;
                } /*for sys 'm', the diagonal element is 0*/
                if (rtk->ssat[j - k].lock[f] >
                    0 && !(rtk->ssat[j - k].slip[f] & 2) && rtk->ssat[i -
                                                                      k].vsat[f] && rtk->ssat[j - k].azel[1] >= rtk->opt.elmaskar)
                {
                    D[i + (na + nb) * nx] = 1.0;
                    D[j + (na + nb) * nx] = -1.0;
                    nb++;
                    rtk->ssat[j - k].fix[f] = 2; /* fix */
                }
                else
                    rtk->ssat[j - k].fix[f] = 1;
            }/*end of j loop*/
        }/*end of f loop*/
    }/*end of m loop*/
    RTKtrace(5, "D'=\n");
    tracemat(5, D, nx, na + nb, 2, 0);
    return nb;
}

/* restore single-differenced ambiguity -------------------------------------- */
static void restamb(rtk_t *rtk, const double *bias, int nb, double *xa) {
    int i, n, m, f, index[MAXSAT], nv = 0, nf = NF(&rtk->opt);
    
    RTKtrace(3, "restamb :\n");
    
    for (i = 0; i < rtk->nx; i++)
        xa[i] = rtk->x[i];
    for (i = 0; i < rtk->na; i++)
        xa[i] = rtk->xa[i];
    
    for (m = 0; m < 4; m++)
        for (f = 0; f < nf; f++) {
            
            for (n = i = 0; i < MAXSAT; i++) {
                if (!test_sys(rtk->ssat[i].sys, m) || rtk->ssat[i].fix[f] != 2)
                {
                    continue;
                }
                index[n++] = IB(i + 1, f, &rtk->opt);
            }
            if (n < 2)
                continue;
            
            xa[index[0]] = rtk->x[index[0]];
            
            for (i = 1; i < n; i++) {
                xa[index[i]] = xa[index[0]] - bias[nv++];
            }
        }
}

/* hold integer ambiguity ---------------------------------------------------- */
static void holdamb(rtk_t *rtk, const double *xa) {
    double *v, *H, *R;
    int i, n, m, f, info, index[MAXSAT], nb = rtk->nx - rtk->na, nv = 0, nf =
    NF(&rtk->opt);
    
    RTKtrace(3, "holdamb :\n");
    
    v = mat(nb, 1);
    H = zeros(nb, rtk->nx);
    
    for (m = 0; m < 4; m++)
        for (f = 0; f < nf; f++) {
            
            for (n = i = 0; i < MAXSAT; i++) {
                if (!test_sys(rtk->ssat[i].sys, m) ||
                    rtk->ssat[i].fix[f] != 2 || rtk->ssat[i].azel[1] <
                    rtk->opt.elmaskhold) {
                    continue;
                }
                index[n++] = IB(i + 1, f, &rtk->opt);
                rtk->ssat[i].fix[f] = 3; /* hold */
            }
            /* constraint to fixed ambiguity */
            for (i = 1; i < n; i++) {
                v[nv] = (xa[index[0]] - xa[index[i]]) -
                (rtk->x[index[0]] - rtk->x[index[i]]);
                
                H[index[0] + nv * rtk->nx] = 1.0;
                H[index[i] + nv * rtk->nx] = -1.0;
                nv++;
            }
        }
    if (nv > 0) {
        R = zeros(nv, nv);
        for (i = 0; i < nv; i++)
            R[i + i * nv] = VAR_HOLDAMB;
        
        /* update states with constraints */
        if ((info = RTKfilter(rtk->x, rtk->P, H, v, R, rtk->nx, nv))) {
            errmsg(rtk, "RTKfilter error (info=%d)\n", info);
        }
        free(R);
    }
    free(v);
    free(H);
}

/* save float solutions ---------------------------------------*/
static int savefloatsolsdualfreq(rtk_t *rtk, double *bias, double *xa)
{
    int i,j,ny,nb,info,nx=rtk->nx,na=rtk->na;
    int nbsl=3,prcm=16,prcn=8;
    double *D,*DP,*y,*Qy,*b,*db,*Qb,*Qab,*Qbsl,*QQ,*bsl,*nbs;
    
    gtime_t time,ts={0};
    double gpst;
    int week,timeu;
    const char *sep=" ";
    char s[64],s1[25];
    solopt_t  solopt;
    solopt.posf=1; solopt.times=0; solopt.timef=1;
    solopt.timeu=3; solopt.degf=0; solopt.outhead=1;
    solopt.outopt=1; solopt.datum=0; solopt.height=0;
    solopt.geoid=0; solopt.solstatic=0; solopt.sstat=0;
    solopt.RTKtrace=1; solopt.nmeaintv[0]=0; solopt.nmeaintv[1]=0;
    /*solopt->sep=" ";
     solopt.prog="Yanqing Hou";     */
    
    
    
    RTKtrace(3,"savefloatsols : nx=%d\n",nx);
    
    timeu=solopt.timeu<0?0:(solopt.timeu>20?20:solopt.timeu);
    time=rtk->sol.time;
    if (solopt.times>=TIMES_UTC) time=gpst2utc(time);
    if (solopt.times==TIMES_JST) time=timeadd(time,9*3600.0);
    
    if (solopt.timef) time2str(time,s,timeu);
    else {
        gpst=time2gpst(time,&week);
        if (86400*7-gpst<0.5/pow(10.0,timeu)) {
            week++;
            gpst=0.0;
        }
        sprintf(s,"%4d%s%*.*f",week,sep,6+(timeu<=0?0:timeu+1),timeu,gpst);
    }
    for (i = 0; i < 25; i++) {
        s1[i]=s[i];
    }
    s1[23]='\r';s1[24]='\0';
    fprintf(stdout,"%s",s1);
    s1[23]='\n';
    /* single to double-difference transformation matrix (D') */
    D=zeros(nx,nx);
    nbs=zeros(4,2);
    if ((nb=yqddmat(rtk,D,nbs))<=0) {
        errmsg(rtk,"no valid double-difference\n");
        free(D);
        return 0;
    }
    ny=na+nb; y=mat(ny,1); Qy=mat(ny,ny); DP=mat(ny,nx);
    b=mat(nb,2); db=mat(nb,1); Qb=mat(nb,nb); Qab=mat(nbsl,nb); QQ=mat(na,nb);
    Qbsl=mat(nbsl,nbsl);  bsl=mat(1,nbsl);
    /* transform single to double-differenced phase-bias (y=D'*x, Qy=D'*P*D) */
    matmul("TN",ny, 1,nx,1.0,D ,rtk->x,0.0,y );
    matmul("TN",ny,nx,nx,1.0,D ,rtk->P,0.0,DP);
    matmul("NN",ny,ny,nx,1.0,DP,D     ,0.0,Qy);
    
    for (i=0;i<nbsl;i++) {
        bsl[i]=rtk->x[i];
        for (j=0;j<nbsl;j++) Qbsl[i+j*nbsl]=rtk->P[i+j*nx];
    }
    /* phase-bias covariance (Qb) and real-parameters to bias covariance (Qab) */
    for (i=0;i<nb;i++) for (j=0;j<nb;j++) Qb [i+j*nb]       =Qy[na+i+(na+j)*ny];
    for (i=0;i<nbsl;i++) for (j=0;j<nb;j++) Qab[i+j*nbsl]   =Qy[   i+(na+j)*ny];
    /* for (i=0;i<nbsl;i++) for (j=0;j<nbsl;j++) Qbsl[i+j*nbsl]=Qy[   i+j*ny]; */
    RTKtrace(1,s1);
    RTKtrace(1,"aflt\n"); tracemat(1,y+na,1,nb,prcm,prcn); /*float ambiguities*/
    RTKtrace(1,"bflt\n"); tracemat(1,bsl,1,nbsl,prcm,prcn);    /* float baseline solution 	*/
    RTKtrace(1,"Qah\n");  tracemat(1,Qb,nb,nb,prcm,prcn);  /* Qahat_ahat */
    RTKtrace(1,"Qbh\n");  tracemat(1,Qbsl,nbsl,nbsl,prcm,prcn);  /* Qahat_ahat */
    RTKtrace(1,"Qba\n");  tracemat(1,Qab,nbsl,nb,prcm,prcn);  /* Qbhat_ahat */
    RTKtrace(1,"nambs\n");tracemat(1,nbs,4,2,prcm,prcn);
    free(D); free(y); free(Qy); free(DP);
    free(b); free(db); free(Qb); free(Qab); free(QQ); free(Qbsl);
    free(nbs);
    return nb; /* number of ambiguities */
    /*					rtk->xa[i]=rtk->x[i];
     for (j=0;j<na;j++) rtk->Pa[i+j*na]=rtk->P[i+j*nx];      */
}


/* save float solutions ---------------------------------------*/
static int savefloatsols(rtk_t *rtk, double *bias, double *xa)
{
    int i,j,ny,nb,info,nx=rtk->nx,na=rtk->na;
    int nbsl=3,prcm=16,prcn=8;
    double *D,*DP,*y,*Qy,*b,*db,*Qb,*Qab,*Qbsl,*QQ,*bsl;
    
    gtime_t time,ts={0};
    double gpst;
    int week,timeu;
    const char *sep=" ";
    char s[64],s1[25];
    solopt_t  solopt;
    solopt.posf=1; solopt.times=0; solopt.timef=1;
    solopt.timeu=3; solopt.degf=0; solopt.outhead=1;
    solopt.outopt=1; solopt.datum=0; solopt.height=0;
    solopt.geoid=0; solopt.solstatic=0; solopt.sstat=0;
    solopt.RTKtrace=1; solopt.nmeaintv[0]=0; solopt.nmeaintv[1]=0;
    /*solopt->sep=" ";
     solopt.prog="Yanqing Hou";     */
    
    
    
    RTKtrace(3,"savefloatsols : nx=%d\n",nx);
    
    timeu=solopt.timeu<0?0:(solopt.timeu>20?20:solopt.timeu);
    time=rtk->sol.time;
    if (solopt.times>=TIMES_UTC) time=gpst2utc(time);
    if (solopt.times==TIMES_JST) time=timeadd(time,9*3600.0);
    
    if (solopt.timef) time2str(time,s,timeu);
    else {
        gpst=time2gpst(time,&week);
        if (86400*7-gpst<0.5/pow(10.0,timeu)) {
            week++;
            gpst=0.0;
        }
        sprintf(s,"%4d%s%*.*f",week,sep,6+(timeu<=0?0:timeu+1),timeu,gpst);
    }
    for (i = 0; i < 25; i++) {
        s1[i]=s[i];
    }
    s1[23]='\n'; s1[24]='\0';
    /* single to double-difference transformation matrix (D') */
    D=zeros(nx,nx);
    if ((nb=ddmat(rtk,D))<=0) {
        errmsg(rtk,"no valid double-difference\n");
        free(D);
        return 0;
    }
    ny=na+nb; y=mat(ny,1); Qy=mat(ny,ny); DP=mat(ny,nx);
    b=mat(nb,2); db=mat(nb,1); Qb=mat(nb,nb); Qab=mat(nbsl,nb); QQ=mat(na,nb);
    Qbsl=mat(nbsl,nbsl);  bsl=mat(1,nbsl);
    /* transform single to double-differenced phase-bias (y=D'*x, Qy=D'*P*D) */
    matmul("TN",ny, 1,nx,1.0,D ,rtk->x,0.0,y );
    matmul("TN",ny,nx,nx,1.0,D ,rtk->P,0.0,DP);
    matmul("NN",ny,ny,nx,1.0,DP,D     ,0.0,Qy);
    
    for (i=0;i<nbsl;i++) {
        bsl[i]=rtk->x[i];
        for (j=0;j<nbsl;j++) Qbsl[i+j*nbsl]=rtk->P[i+j*nx];
    }
    /* phase-bias covariance (Qb) and real-parameters to bias covariance (Qab) */
    for (i=0;i<nb;i++) for (j=0;j<nb;j++) Qb [i+j*nb]       =Qy[na+i+(na+j)*ny];
    for (i=0;i<nbsl;i++) for (j=0;j<nb;j++) Qab[i+j*nbsl]   =Qy[   i+(na+j)*ny];
    /* for (i=0;i<nbsl;i++) for (j=0;j<nbsl;j++) Qbsl[i+j*nbsl]=Qy[   i+j*ny]; */
    RTKtrace(1,s1);
    RTKtrace(1,"aflt\n"); tracemat(1,y+na,1,nb,prcm,prcn); /*float ambiguities*/
    RTKtrace(1,"bflt\n"); tracemat(1,bsl,1,nbsl,prcm,prcn);    /* float baseline solution 	*/
    RTKtrace(1,"Qah\n");  tracemat(1,Qb,nb,nb,prcm,prcn);  /* Qahat_ahat */
    RTKtrace(1,"Qbh\n");  tracemat(1,Qbsl,nbsl,nbsl,prcm,prcn);  /* Qahat_ahat */
    RTKtrace(1,"Qba\n");  tracemat(1,Qab,nbsl,nb,prcm,prcn);  /* Qbhat_ahat */
    
    free(D); free(y); free(Qy); free(DP);
    free(b); free(db); free(Qb); free(Qab); free(QQ); free(Qbsl);
    
    return nb; /* number of ambiguities */
    /*					rtk->xa[i]=rtk->x[i];
     for (j=0;j<na;j++) rtk->Pa[i+j*na]=rtk->P[i+j*nx];      */
}
/* save fixed ambiguities -------------(rtk, bias, nb, xa)--------------------------*/
static int saveAmbsGeometry(rtk_t *rtk)
{
    int i, cnt=0;
    
    gtime_t time,ts={0};
    double gpst;
    int week,timeu;
    const char *sep=" ";
    char s[500], s2[20]={" "}, s3[40]={" "}, s4[100]={" "};
    solopt_t  solopt;
    solopt.posf=1; solopt.times=0; solopt.timef=1;
    solopt.timeu=3; solopt.degf=0; solopt.outhead=1;
    solopt.outopt=1; solopt.datum=0; solopt.height=0;
    solopt.geoid=0; solopt.solstatic=0; solopt.sstat=0;
    solopt.RTKtrace=1; solopt.nmeaintv[0]=0; solopt.nmeaintv[1]=0;
    
    /*initialize time*/
    timeu=solopt.timeu<0?0:(solopt.timeu>20?20:solopt.timeu);
    time=rtk->sol.time;
    if (solopt.times>=TIMES_UTC) time=gpst2utc(time);
    if (solopt.times==TIMES_JST) time=timeadd(time,9*3600.0);
    
    if (solopt.timef) time2str(time,s,timeu);
    else {
        gpst=time2gpst(time,&week);
        if (86400*7-gpst<0.5/pow(10.0,timeu)) {
            week++;
            gpst=0.0;
        }
        sprintf(s,"%4d%s%*.*f",week,sep,6+(timeu<=0?0:timeu+1),timeu,gpst);
    }
    if (rtk->sol.time.time==1424507940) {
        cnt++;
        cnt--;
    }
    /*output interested parameters*/
    for (i = 0; i < MAXSAT; i++) {
        if(rtk->ssat[i].fix[0]==2)  /*fix==2 means it will be fixed by AR. 0-L1, 1-L2*/
        {
            cnt++;
            /*sats*/
            sprintf(s2,"%2d ",i+1);  /*sat prn*/
            strcat(s3,s2);
            /*eles*/
            sprintf(s2,"%4.3f ", rtk->ssat[i].azel[1]);   /*ele of sat*/
            strcat(s4,s2);
        }
        
    }
    sprintf(s2,"%2d ",cnt);
    /*strcat(s,s2,s3,s4);s2--# of sats, s3--sats prns, s4--sats eles*/
    strcat(s,s2);
    strcat(s,s3);
    strcat(s,s4);
    
    
    s2[0]='\n'; s2[1]='\0';
    strcat(s,s2);
    RTKtrace(1, s);
    return cnt; /* number of sats */
}


/* save fixed ambiguities -------------(rtk, bias, nb, xa)--------------------------*/
static int saveMuAndAmbs(rtk_t *rtk, double *bias, int nb, double mu)
{
    int i;
    
    gtime_t time,ts={0};
    double gpst;
    int week,timeu;
    const char *sep=" ";
    char s[500], s2[20]={" "};
    solopt_t  solopt;
    solopt.posf=1; solopt.times=0; solopt.timef=1;
    solopt.timeu=3; solopt.degf=0; solopt.outhead=1;
    solopt.outopt=1; solopt.datum=0; solopt.height=0;
    solopt.geoid=0; solopt.solstatic=0; solopt.sstat=0;
    solopt.RTKtrace=1; solopt.nmeaintv[0]=0; solopt.nmeaintv[1]=0;
    /*solopt->sep=" ";
     solopt.prog="Yanqing Hou";     */
    
    timeu=solopt.timeu<0?0:(solopt.timeu>20?20:solopt.timeu);
    time=rtk->sol.time;
    if (solopt.times>=TIMES_UTC) time=gpst2utc(time);
    if (solopt.times==TIMES_JST) time=timeadd(time,9*3600.0);
    
    if (solopt.timef) time2str(time,s,timeu);
    else {
        gpst=time2gpst(time,&week);
        if (86400*7-gpst<0.5/pow(10.0,timeu)) {
            week++;
            gpst=0.0;
        }
        sprintf(s,"%4d%s%*.*f",week,sep,6+(timeu<=0?0:timeu+1),timeu,gpst);
    }
    sprintf(s2,"%12.5f %4d ",mu,nb);
    strcat(s,s2);
    for (i = 0; i < nb; i++) {
        sprintf(s2,"%8d ", (int)bias[i]);
        strcat(s,s2);
    }
    s2[0]='\n'; s2[1]='\0';
    strcat(s,s2);
    RTKtrace(1, s);
    return nb; /* number of ambiguities */
}
/* save mu --------------------------------------- */
static int savemu(rtk_t *rtk, double mu) {
    /* ,j,ny,nb,info,nx=rtk->nx,na=rtk->na; */
    int i;
    gtime_t time, ts = {0};
    double gpst;
    int week, timeu;
    const char *sep = " ";
    char s[64], s1[41],s2[16];
    solopt_t solopt;
    solopt.posf = 1;
    solopt.times = 0;
    solopt.timef = 1;
    solopt.timeu = 3;
    solopt.degf = 0;
    solopt.outhead = 1;
    solopt.outopt = 1;
    solopt.datum = 0;
    solopt.height = 0;
    solopt.geoid = 0;
    solopt.solstatic = 0;
    solopt.sstat = 0;
    solopt.RTKtrace = 1;
    solopt.nmeaintv[0] = 0;
    solopt.nmeaintv[1] = 0;
    /* solopt->sep=" ";
     solopt.prog="Yanqing Hou"; */
    
    /* RTKtrace(3,"savemu : nx=%d\n",nx); */
    
    timeu = solopt.timeu < 0 ? 0 : (solopt.timeu > 20 ? 20 : solopt.timeu);
    time = rtk->sol.time;
    if (solopt.times >= TIMES_UTC)
        time = gpst2utc(time);
    if (solopt.times == TIMES_JST)
        time = timeadd(time, 9 * 3600.0);
    
    if (solopt.timef)
        time2str(time, s, timeu);
    else {
        gpst = time2gpst(time, &week);
        if (86400 * 7 - gpst < 0.5 / pow(10.0, timeu)) {
            week++;
            gpst = 0.0;
        }
        sprintf(s, "%4d%s%*.*f", week, sep,
                6 + (timeu <= 0 ? 0 : timeu + 1), timeu, gpst);
    }
    sprintf(s2,"%12.5f",mu);
    for (i = 0; i < 25; i++) {
        s1[i] = s[i];
    }
    s1[23]=s1[24]=' ';
    for (i=0; i < 12; i++) {
        s1[i+25]=s2[i];
    }
    s1[36] = '\n';
    s1[37] = '\0';
    RTKtrace(1, s1);
    /* mu=%12.5f*/
    return 1;
}
/*resolve integer ambiguity by WLNL, implemented by Yanqing Hou 2016*/
/*
 static int resamb_WLNL(rtk_t *rtk, const obsd_t *obs, int *sat, int *iu, int *ir, int ns, const nav_t *nav, const double *azel){
 int i, j, k, ny, nb, nbs[12]={0}, nnbs[4]={0}, sysn=0,rowWL, info, nx = rtk->nx, na = rtk->na, n, m, f, index[MAXSAT], nv = 0, nf = NF(&rtk->opt);
 prcopt_t *opt = &rtk->opt;
 double *D, *WL, *DP, *WLD,*y, *Qy, *b, *db, *Qb, *Qab, *QQ,*wlfreqs, *wlams, s[2];
 double  tmpfreqs[12]={0},tmpwlams[12]={0},freqs[12]={FREQ1,FREQ2,FREQ5,FREQ6,FREQ7,FREQ8,FREQ1_GLO,FREQ2_GLO,FREQ3_GLO,FREQ1_CMP,FREQ2_CMP,FREQ3_CMP};
 RTKtrace(3, "resolve amb by WLNL:\n");
 
 if (rtk->opt.mode <= PMODE_DGPS || rtk->opt.modear == ARMODE_OFF ||
 rtk->opt.thresar[0] < 1.0) {
 return 0;
 }
 single to double-difference transformation matrix (D')
 D = zeros(nx, nx);
 if ((nb = yqddmat(rtk, D, nbs)) <= 0) {
 errmsg(rtk, "no valid double-difference\n");
 free(D);
 return 0;
 }
 find available dual frequency systems
 for(m=0;m<4;m++){
 if (m == 1 && rtk->opt.glomodear == 0)
 continue;
 if (m == 3 && rtk->opt.bdsmodear == 0)
 continue;
 if (nbs[m]>0&&nbs[m+4]==nbs[m]){
 nnbs[sysn]=nbs[m];
 tmpfreqs[sysn*3]= freqs[m*3]; Freq L1 or B1
 tmpfreqs[sysn*3+1]= freqs[m*3+1]; Freq L1 or B1
 tmpfreqs[sysn*3+2]= freqs[m*3]+freqs[m*3+1];
 tmplams[sysn*3]=CLIGHT/tmpfreqs[sysn*3];
 tmplams[sysn*3+1]=CLIGHT/tmpfreqs[sysn*3+1];
 tmplams[sysn*3+2]=CLIGHT/tmpfreqs[sysn*3+2];
 sysn++;
 }
 else{
 errmsg(rtk, "no valid widelane ambs\n");
 return 0;
 }
 
 }
 make WL matrix, not WL'
 nb=0;
 for (i = 0; i < sysn; i++) {
 nb+=2*nnbs[i]-2;
 }
 ny=nb+na;
 wlfreqs=zeros(sysn*2,sysn*2);
 wlams=zeros(sysn*2,sysn*2);
 for(i=0;i<sysn;i++){
 wlfreqs[i*2]=tmpfreqs[i*3];
 wlfreqs[i*2+1]=tmpfreqs[i*3+2];
 wlams[i*2]=tmplams[i*3];
 wlams[i*2+1]=tmplams[i*3+2];
 }
 
 WL=zeros(ny,ny);
 for(i=0;i<na;i++){
 WL[i+i*ny]=1; I_na
 }
 WL=[I_na, 0,          0,    0;
 0,    I_dg, 0,    0,    0;
 0,    I_dg, -Idg, 0,    0;
 0,    0,    0,    I_db, 0;
 0,    0,    0,    I_db, -I_db];
 since the ambiguities a are already in cycles, to make wl combination, we should
 multiply the wavelength to the corresponding ambiguities.
 for (i = 0,k=na; i < sysn; k+=2*nnbs[i]-2, i++) {
 for (j=k; j < k+nnbs[i]-1; j++) {
 WL[j+j*ny]=CLIGHT/tmpfreqs[i*3]; I_dg*wavelength_L1
 }
 for (j = k+nnbs[i]-1; j < k+2*nnbs[i]-2; j++) {
 WL[j+j*ny]=-CLIGHT/tmpfreqs[i*3+1]; -I_dg*wavelength_L2
 }
 for (j = k+nnbs[i]-1; j < k+2*nnbs[i]-2; j++) {
 WL[j+ny*(j-nnbs[i]+1)]=CLIGHT/tmpfreqs[i*3]; I_dg*wavelength_L1
 }
 
 }
 Display WL matrix
 RTKtrace(5, "WL*lambda=\n");
 tracemat(5, WL, ny, ny, 2, 0);
 make L1/WL ambiguities
 ny = na + nb;
 y = mat(ny, 1);
 Qy = mat(ny, ny);
 DP = mat(ny, nx);
 b = mat(nb, 2);
 db = mat(nb, 1);
 Qb = mat(nb, nb);
 Qab = mat(na, nb);
 QQ = mat(na, nb);
 WLD=mat(ny,nx);
 
 matmul("NT", ny, nx, ny, 1.0, WL,D,0.0,WLD);
 transform single to double-differenced phase-bias (y=D'*x, Qy=D'*P*D)
 matmul("NN", ny,  1, ny,  1.0, WLD, rtk->x, 0.0, y);
 matmul("NN", ny, nx, ny, 1.0, WLD, rtk->P, 0.0, DP);
 matmul("NT", ny, ny, nx, 1.0, DP,  WLD,    0.0, Qy);
 
 phase-bias covariance (Qb) and real-parameters to bias covariance (Qab)
 for (i = 0; i < nb; i++)
 for (j = 0; j < nb; j++)
 Qb[i + j * nb] = Qy[na + i + (na + j) * ny];
 for (i = 0; i < na; i++)
 for (j = 0; j < nb; j++)
 Qab[i + j * na] = Qy[i + (na + j) * ny];
 
 RTKtrace(4, "N(0)=");
 tracemat(4, y + na, 1, nb, 10, 3);
 if(!info=(wlnlar(nnbs,wlfreqs,wlams,sysn,2,y+na,Qb,b,s,mu,Pf))){
 
 }
 else{
 errmsg(rtk,"WL/NL AR failed!\n");
 nb=0;
 }
 return nb;
 }
 */
/* resolve integer ambiguity by LAMBDA --------------------------------------- */
static int resamb_LAMBDA(rtk_t *rtk, double *bias, double *xa) {
    prcopt_t *opt = &rtk->opt;
    int i, j, ny, nb, info, nx = rtk->nx, na = rtk->na;
    double *D, *DP, *y, *Qy, *b, *db, *Qb, *Qab, *QQ, s[2];
    
    double Pf = 0.01, mu[2]={0,0}, muselect = 0;
    /* int ffrt=0; */
    /*
     if (opt->thresar[0]==100 || opt->thresar[0]==1000) {
     Pf=1/opt->thresar[0];
     ffrt=1;
     } */
    if (opt->thresar[0] == 100) {
        Pf = 0.01;
    }
    else if (opt->thresar[0] == 200) {
        Pf = 0.01;
    }
    else if (opt->thresar[0] == 1000) {
        Pf = 0.001;
    }
    else if (opt->thresar[0] == 2000) {
        Pf = 0.001;
    }
    else {
        Pf = 0.01;
    }
    RTKtrace(3, "resamb_LAMBDA : nx=%d\n", nx);
    
    rtk->sol.ratio = 0.0;
    
    if (rtk->opt.mode <= PMODE_DGPS || rtk->opt.modear == ARMODE_OFF ||
        rtk->opt.thresar[0] < 1.0) {
        return 0;
    }
    /* single to double-difference transformation matrix (D') */
    D = zeros(nx, nx);
    if ((nb = ddmat(rtk, D)) <= 0) {
        errmsg(rtk, "no valid double-difference\n");
        free(D);
        return 0;
    }
    
    // RTKtrace(4, "DDmat=");
    //tracemat(4, D, nb, nb, 10, 6);
    ny = na + nb;
    y = mat(ny, 1);
    Qy = mat(ny, ny);
    DP = mat(ny, nx);
    b = mat(nb, 2);
    db = mat(nb, 1);
    Qb = mat(nb, nb);
    Qab = mat(na, nb);
    QQ = mat(na, nb);
    
    /* transform single to double-differenced phase-bias (y=D'*x, Qy=D'*P*D) */
    matmul("TN", ny, 1, nx, 1.0, D, rtk->x, 0.0, y);
    matmul("TN", ny, nx, nx, 1.0, D, rtk->P, 0.0, DP);
    matmul("NN", ny, ny, nx, 1.0, DP, D, 0.0, Qy);
    
    /* phase-bias covariance (Qb) and real-parameters to bias covariance (Qab) */
    for (i = 0; i < nb; i++)
        for (j = 0; j < nb; j++)
            Qb[i + j * nb] = Qy[na + i + (na + j) * ny];
    for (i = 0; i < na; i++)
        for (j = 0; j < nb; j++)
            Qab[i + j * na] = Qy[i + (na + j) * ny];
    
    RTKtrace(4, "N(0)=");
    tracemat(4, y + na, 1, nb, 10, 6);
    
    
    RTKtrace(4, "Qb=");
    tracemat(4, Qb, nb, nb, 10, 6);
    
    
    RTKtrace(4, "RTKQba=");
    tracemat(4, Qab, na, nb, 10, 6);
    /* if(rtk->sol.time.time>=1424487600){
     Pf=Pf;
     }; */
    /* lambda/mlambda integer least-square estimation */
    if (!(info = lambda(nb, 2, y + na, Qb, b, s, mu, Pf))) {
        
        RTKtrace(4, "N(1)=");
        tracemat(4, b, 1, nb, 10, 3);
        RTKtrace(4, "N(2)=");
        tracemat(4, b + nb, 1, nb, 10, 3);
        
        rtk->sol.ratio = s[0] > 0 ? (float)(s[1] / s[0]) : 0.0f;
        if (rtk->sol.ratio > 999.9)
            rtk->sol.ratio = 999.9f;
        
        if (opt->thresar[0] == 100) {
            muselect = mu[0];
        }
        else if (opt->thresar[0] == 200) {
            muselect = mu[1];
        }
        else if (opt->thresar[0] == 1000) {
            muselect = mu[0];
        }
        else if (opt->thresar[0] == 2000) {
            muselect = mu[1];
        }
        else {
            muselect = opt->thresar[0];
        }
        /* if (!ffrt) *mu  = opt->thresar[0]; */
        /* validation by popular ratio-test */
        /*savemu(rtk, muselect);  */
        /*saveMuAndAmbs(rtk, b, nb, muselect); */
        /*saveAmbsGeometry(rtk);             */
        if (s[0] <= 0.0 || s[1] / s[0] >= muselect) {
            
            /* transform float to fixed solution (xa=xa-Qab*Qb\(b0-b)) */
            for (i = 0; i < na; i++) {
                rtk->xa[i] = rtk->x[i];
                for (j = 0; j < na; j++)
                    rtk->Pa[i + j * na] = rtk->P[i + j * nx];
            }
            for (i = 0; i < nb; i++) {
                bias[i] = b[i];
                y[na + i] -= b[i];
            }
            if (!matinv(Qb, nb)) {
                matmul("NN", nb, 1, nb, 1.0, Qb, y + na, 0.0, db);
                matmul("NN", na, 1, nb, -1.0, Qab, db, 1.0, rtk->xa);
                //                RTKtrace(1,"db=[");         tracemat(1, db, 1, nb, 10, 3); RTKtrace(1,"];");
                //                RTKtrace(1,"RTKx=[");       tracemat(1, rtk->x, 1, 3, 10, 3); RTKtrace(1,"];");
                //                RTKtrace(1,"RTKxa=[");       tracemat(1, rtk->xa, 1, 3, 10, 3); RTKtrace(1,"];");
                
                
                /* covariance of fixed solution (Qa=Qa-Qab*Qb^-1*Qab') */
                matmul("NN", na, nb, nb, 1.0, Qab, Qb, 0.0, QQ);
                matmul("NT", na, na, nb, -1.0, QQ, Qab, 1.0, rtk->Pa);
                
                RTKtrace(3,
                         "resamb : validation ok (nb=%d ratio=%.2f s=%.2f/%.2f)\n",
                         nb, s[0] == 0.0 ? 0.0 : s[1] / s[0], s[0], s[1]);
                
                /* restore single-differenced ambiguity */
                restamb(rtk, bias, nb, xa);
            }
            else
                nb = 0;
        }
        else { /* validation failed */
            errmsg(rtk,
                   "ambiguity validation failed (nb=%d ratio=%.2f s=%.2f/%.2f)\n",
                   nb, s[1] / s[0], s[0], s[1]);
            nb = 0;
        }
    }
    else {
        errmsg(rtk, "lambda error (info=%d)\n", info);
    }
    free(D);
    free(y);
    free(Qy);
    free(DP);
    free(b);
    free(db);
    free(Qb);
    free(Qab);
    free(QQ);
    
    return nb; /* number of ambiguities */
}
static int MatTrans(double *A, int iRow, int iCol, double *B)
{
    int i,j;
    for (j=0;j<iCol;j++)
        for (i=0;i<iRow;i++)
            B[i*iCol+j]=A[j*iRow+i];
    return 1;
}
/* resolve integer ambiguity by LAMBDA --------------------------------------- */
static int resamb_TSRC(rtk_t *rtk, double *bias, double *xa) {
    prcopt_t *opt = &rtk->opt;
    int i, j, ny, nb, info=0, nx = rtk->nx, na = rtk->na,fixno=0;
    double *D, *DP, *y, *Qy, *b, *db, *Qb, *Qab,  *Qaa, s[2],*bchk,*bhat,*Qbchk;
    double Pf = 0.01, mu[2]={0,0}, muselect = 0;
    
    static double dbgQbb[9]={0.35158000014094726,-0.4392067940856300,-0.2175948406996578,
        -0.4392067940856448,1.11971547442439650,0.39141582475067221,
        -0.2175948406996612,0.39141582475066816,0.32311417336082854};
    static double dbgQab[30]={-0.0180544573765851,0.20362046018947033,0.10062408907781301,0.48929812047434440,0.13293040290647182,0.04761389599833459,0.25706816665298643,0.05435282955319777,0.17447623094196174,0.46912377687894813,
        -0.2095810224883106,-0.2308445548425488,-0.5968976404111836,-0.7223094006231740,-0.2639596705780548,-0.2323442210009805,-0.6867070332189844,-0.0161651967144759,-0.2850574615639984,-0.7501016430141588,
        -0.0328416500566006,-0.1273354482914481,-0.1590891768191323,-0.3334601154634608,-0.1660759357849755,0.10987098392484669,-0.3419721029687737,-0.0066219406177061,0.02603616847327139,-0.2234747369992826};
    static double dbgQaa[100]={0.10412815159868744,-0.0161435979940245,0.19048313298841235,0.02386500076645337,0.03524084466093265,0.07681778188124482,0.13646661924792250,-0.0193144791980361,0.02653417031521104,0.04950016161276592,
        -0.0161435979940278,0.12818977562333855,0.04284571863995889,0.28313905303706144,0.07942401906890789,0.01780219992285117,0.14193464333326139,0.03704406260499806,0.09668204500465164,0.26512419738747611,
        0.19048313298840996,0.04284571863996590,0.42845908126976484,0.23564842358237045,0.12215490547003446,0.16232127089037018,0.37806000067728779,-0.0228958009699509,0.11104212562877466,0.27484292374504737,
        0.02386500076644493,0.28313905303706366,0.23564842358235438,0.73055598936198130,0.21094874525970053,0.09601602100836534,0.43771157952791800,0.06890355445241899,0.25324080362446477,0.68427634567205198,
        0.03524084466093042,0.07942401906890952,0.12215490547003093,0.21094874525970314,0.10086390509406487,-0.0163167197536024,0.20510699989928607,0.00946366059719929,0.02415841971252602,0.17033434268281553,
        0.07681778188124344,0.01780219992285360,0.16232127089036905,0.09601602100837058,-0.0163167197536014,0.25576359459930648,0.03531597983201048,0.01313755546099259,0.20706590509412109,0.22599319314931218,
        0.13646661924791775,0.14193464333326714,0.37806000067728185,0.43771157952792922,0.20510699989928760,0.03531597983200933,0.51725255724929708,0.00014511536218223,0.07823734193075304,0.38622896156081699,
        -0.0193144791980369,0.03704406260499740,-0.0228958009699531,0.06890355445241679,0.00946366059719845,0.01313755546099173,0.00014511536217993,0.02284646665421571,0.03961161488721408,0.07173905279570741,
        0.02653417031520775,0.09668204500465281,0.11104212562876861,0.25324080362446555,0.02415841971252518,0.20706590509411921,0.07823734193074911,0.03961161488721502,0.23741678385692594,0.34786638267089498,
        0.04950016161275754,0.26512419738747905,0.27484292374503194,0.68427634567205387,0.17033434268281339,0.22599319314930730,0.38622896156080688,0.07173905279570978,0.34786638267089498,0.75956394039646336 };
    static double dbga[10]={-7.0986337563673736,
        -33.973405774071608,
        24.798613286826853,
        -13.039620112002261,
        98.911336664379334,
        34.887123546749670,
        39.835221010152317,
        2.0153822327987143,
        35.952859072996411,
        30.845380396408935};
    double dbgQba[30]={0};
    MatTrans(dbgQab,10,3,dbgQba);
    /* int ffrt=0; */
    /*
     if (opt->thresar[0]==100 || opt->thresar[0]==1000) {
     Pf=1/opt->thresar[0];
     ffrt=1;
     } */
    if (opt->thresar[0] == 100) {
        Pf = 0.01;
    }
    else if (opt->thresar[0] == 200) {
        Pf = 0.01;
    }
    else if (opt->thresar[0] == 1000) {
        Pf = 0.001;
    }
    else if (opt->thresar[0] == 2000) {
        Pf = 0.001;
    }
    else {
        Pf = 0.01;
    }
    RTKtrace(3, "resamb_TSRC : nx=%d\n", nx);
    
    rtk->sol.ratio = 0.0;
    
    if (rtk->opt.mode <= PMODE_DGPS || rtk->opt.modear == ARMODE_OFF ||
        rtk->opt.thresar[0] < 1.0) {
        return 0;
    }
    /* single to double-difference transformation matrix (D') */
    D = zeros(nx, nx);
    if ((nb = ddmat(rtk, D)) <= 0) {
        errmsg(rtk, "no valid double-difference\n");
        free(D);
        return 0;
    }
    ny = na + nb;
    y = mat(ny, 1);
    Qy = mat(ny, ny);
    DP = mat(ny, nx);
    b = mat(nb, 2);
    db = mat(nb, 1);
    Qb = mat(nb, nb);
    Qab = mat(3, nb);
    //QQ  = mat(na, nb);
    Qaa = mat(3,3);
    
    bchk=mat(1,3);
    // bhat=mat(1,3);
    Qbchk=mat(3,3);
    /* transform single to double-differenced phase-bias (y=D'*x, Qy=D'*P*D) */
    matmul("TN", ny, 1, nx, 1.0, D, rtk->x, 0.0, y);
    matmul("TN", ny, nx, nx, 1.0, D, rtk->P, 0.0, DP);
    matmul("NN", ny, ny, nx, 1.0, DP, D, 0.0, Qy);
    
    /* phase-bias covariance (Qb) and real-parameters to bias covariance (Qab) */
    for (i = 0; i < nb; i++)
        for (j = 0; j < nb; j++)
            Qb[i + j * nb] = Qy[na + i + (na + j) * ny];
    for (i = 0; i < 3; i++)
        for (j = 0; j < nb; j++)
            Qab[i + j * 3] = Qy[i + (na + j) * ny];
    for (i = 0; i < 3; i++) {
        for (j = 0; j < 3; j++)
            Qaa[i + j * 3] = rtk->P[i + j * nx];
    }
    for (i = 0; i < 3; i++)
        bchk[i] = rtk->x[i];
    RTKtrace(4, "N(0)=");
    tracemat(4, y + na, 1, nb, 10, 3);
    /* if(rtk->sol.time.time>=1424487600){
     Pf=Pf;
     }; */
    /* lambda/mlambda integer least-square estimation */
    /*for debug*/
    RTKtrace(1,"ahat=");tracemat(1,y+na,1,nb,10,5);
    RTKtrace(1,"Qahat=");tracemat(1,Qb,nb,nb,10,5);
    RTKtrace(1,"Qbhat=");tracemat(1,Qaa,3,3,10,5);
    RTKtrace(1,"Qba=");tracemat(1,Qab,3,nb,10,5);
    if ((fixno = TSRC(10, 2, dbga, dbgQbb,dbgQba,dbgQaa, b, s,bchk,Qbchk,Pf))) {
        //	if ((fixno = TSRC(nb, 2, y + na, Qaa,Qab,Qb, b, s,bchk,Qbchk,Pf))) {
        //extern int TSRC(int na, int m, const double *a, const double *Qbb, const double *Qba, const double *Qaa,
        //		double *F, double *s, double *bchk, double *Qbchk, double Pf)
        RTKtrace(4, "N(1)=");
        tracemat(4, b, 1, nb, 10, 3);
        //RTKtrace(4, "N(2)=");
        //tracemat(4, b + nb, 1, nb, 10, 3);
        
        /* if (!ffrt) *mu  = opt->thresar[0]; */
        /* validation by popular ratio-test */
        /*savemu(rtk, muselect);  */
        /*saveMuAndAmbs(rtk, b, nb, muselect); */
        /*saveAmbsGeometry(rtk);             */
        if (fixno>0) {
            
            /* transform float to fixed solution (xa=xa-Qab*Qb\(b0-b)) */
            for (i = 0; i < 3; i++) {
                rtk->xa[i] = rtk->x[i];
                for (j = 0; j < 3; j++)
                    rtk->Pa[i + j * 3] = rtk->P[i + j * nx];
            }
            for (i = 0; i < nb; i++) {
                bias[i] = b[i];
                y[na + i] -= b[i];
            }
            
            for(i=0;i<3;i++)
                rtk->xa[i]=bchk[i];
            matcpy(rtk->Pa,Qbchk,3,3);
            restamb(rtk, bias, nb, xa);
        }
        else { /* validation failed */
            errmsg(rtk,
                   "ambiguity validation failed (nb=%d ratio=%.2f s=%.2f/%.2f)\n",
                   nb, s[1] / s[0], s[0], s[1]);
            nb = 0;
        }
    }
    else {
        errmsg(rtk, "lambda error (info=%d)\n", info);
    }
    free(D);
    free(y);
    free(Qy);
    free(DP);
    free(b);
    free(db);
    free(Qb);
    free(Qab);
    //free(QQ);
    free(bchk);
    //    free(bhat);
    free(Qbchk);
    
    return nb; /* number of ambiguities */
}



/* validation of solution ---------------------------------------------------- */
static int valpos(rtk_t *rtk, const double *v, const double *R, const int *vflg,
                  int nv, double thres) {
#if 0
    prcopt_t *opt = &rtk->opt;
    double vv = 0.0;
#endif
    double fact = thres * thres;
    int i, stat = 1, sat1, sat2, type, freq;
    char *stype;
    
    RTKtrace(3, "valpos  : nv=%d thres=%.1f\n", nv, thres);
    
    /* post-fit residual test */
    for (i = 0; i < nv; i++) {
        if (v[i] * v[i] <= fact * R[i + i * nv])
            continue;
        sat1 = (vflg[i] >> 16) & 0xFF;
        sat2 = (vflg[i] >> 8) & 0xFF;
        type = (vflg[i] >> 4) & 0xF;
        freq = vflg[i] & 0xF;
        stype = type == 0 ? "L" : (type == 1 ? "L" : "C");
        errmsg(rtk, "large residual (sat=%2d-%2d %s%d v=%6.3f sig=%.3f)\n",
               sat1, sat2, stype, freq + 1, v[i], SQRT(R[i + i*nv]));
    }
#if 0 /* omitted v.2.4.0 */
    if (stat && nv > NP(opt)) {
        
        /* chi-square validation */
        for (i = 0; i < nv; i++)
            vv += v[i] * v[i] / R[i + i * nv];
        
        if (vv > chisqr[nv - NP(opt) - 1]) {
            errmsg(rtk,
                   "residuals validation failed (nv=%d np=%d vv=%.2f cs=%.2f)\n",
                   nv, NP(opt), vv, chisqr[nv - NP(opt) - 1]);
            stat = 0;
        }
        else {
            RTKtrace(3,
                     "valpos : validation ok (%s nv=%d np=%d vv=%.2f cs=%.2f)\n",
                     rtk->tstr, nv, NP(opt), vv, chisqr[nv - NP(opt) - 1]);
        }
    }
#endif
    return stat;
}

/* relative positioning ------------------------------------------------------
 using IFLC combination to eliminate ionospheric delay*/
static int relpos(rtk_t *rtk, const obsd_t *obs, int nu, int nr,
                  const nav_t *nav) {
    //TODO:relpos
    prcopt_t *opt = &rtk->opt;
    gtime_t time = obs[0].time;
    double *rs, *dts, *var, *y, *e, *azel, *v, *H, *R, *xp, *Pp, *xa, *bias, dt;
    int i, j, f, n = nu + nr, ns, ny, nv, sat[MAXSAT], iu[MAXSAT],
    ir[MAXSAT], niter;
    int info, vflg[MAXOBS * NFREQ * 2 + 1], svh[MAXOBS * 2];
    int stat = rtk->opt.mode <= PMODE_DGPS ? SOLQ_DGPS : SOLQ_FLOAT;
    int nf = opt->ionoopt == IONOOPT_IFLC ? 1 : opt->nf;
    static int icount=0;
    double gpst;
    int week;
    icount++;
    RTKtrace(3, "relpos  : nx=%d nu=%d nr=%d\n", rtk->nx, nu, nr);
    
    dt = timediff(time, obs[nu].time);
    
    rs = mat(6, n);
    dts = mat(2, n);
    var = mat(1, n);
    y = mat(nf * 2, n);
    e = mat(3, n);
    azel = zeros(2, n);
    
    for (i = 0; i < MAXSAT; i++) {
        rtk->ssat[i].sys = satsys(i + 1, NULL);
        for (j = 0; j < NFREQ; j++)
            rtk->ssat[i].vsat[j] = rtk->ssat[i].snr[j] = 0;
    }
    /* satellite positions/clocks */
    satposs(time, obs, n, nav, opt->sateph, rs, dts, var, svh);
    //    gpst=time2gpst(time,&week);
    
    /* undifferenced residuals for base station */
    if (!zdres(1, obs + nu, nr, rs + nu * 6, dts + nu * 2, svh + nu, nav,
               rtk->rb, opt, 1, y + nu * nf * 2, e + nu * 3, azel + nu * 2)) {
        errmsg(rtk, "initial base station position error\n");
        
        free(rs);
        free(dts);
        free(var);
        free(y);
        free(e);
        free(azel);
        return 0;
    }
    //          RTKtrace(1, "basdts="); tracemat(1, dts + nu * 2, 2, 12, 13, 7);
    //       RTKtrace(1, "basrs="); tracemat(1, rs + nu * 6, 6, 12, 13, 7);
    //       RTKtrace(1, "basel="); tracemat(1, azel + nu * 2, 2, 12, 13, 7);
    /* time-interpolation of residuals (for post-processing) */
    if (opt->intpref) {
        dt = intpres(time, obs + nu, nr, nav, rtk, y + nu * nf * 2);
    }
    /* select common satellites between rover and base-station */
    if ((ns = selsat(obs, azel, nu, nr, opt, sat, iu, ir)) <= 0) {
        errmsg(rtk, "no common satellite\n");
        
        free(rs);
        free(dts);
        free(var);
        free(y);
        free(e);
        free(azel);
        return 0;
    }
    /* temporal update of states */
    udstate(rtk, obs, sat, iu, ir, ns, nav);
    
    RTKtrace(4, "x(0)=");
    tracemat(4, rtk->x, 1, NR(opt), 13, 4);
    
    xp = mat(rtk->nx, 1);
    Pp = zeros(rtk->nx, rtk->nx);
    xa = mat(rtk->nx, 1);
    matcpy(xp, rtk->x, rtk->nx, 1);
    
    ny = ns * nf * 2 + 2;
    v = mat(ny, 1);
    H = zeros(rtk->nx, ny);
    R = mat(ny, ny);
    bias = mat(rtk->nx, 1);
    
    /* add 2 iterations for baseline-constraint moving-base */
    niter = opt->niter + (opt->mode == PMODE_MOVEB && opt->baseline[0] >
                          0.0 ? 2 : 0);
    //      RTKtrace(1, "dts="); tracemat(1, dts, 2, ns, 13, 7);
    //        RTKtrace(1, "rs="); tracemat(1, rs, 6, ns, 13, 7);
    for (i = 0; i < niter; i++) {
        /* undifferenced residuals for rover */
        if (!zdres(0, obs, nu, rs, dts, svh, nav, xp, opt, 0, y, e, azel)) {
            errmsg(rtk, "rover initial position error\n");
            stat = SOLQ_NONE;
            break;
        }
        //        RTKtrace(1, "SDAmat="); tracemat(1, e, 3,ns+1, 13, 4);
        /* double-differenced residuals and partial derivatives */
        if ((nv = ddres(rtk, nav, dt, xp, Pp, sat, y, e, azel, iu, ir, ns, v, H,
                        R, vflg)) < 1) {
            errmsg(rtk, "no double-differenced residual\n");
            stat = SOLQ_NONE;
            break;
        }
        
        //        RTKtrace(1, "H=\n");
        //        tracemat(1, H, rtk->nx, nv, 7, 4);
        /* kalman RTKfilter measurement update */
        matcpy(Pp, rtk->P, rtk->nx, rtk->nx);
        if ((info = RTKfilter(xp, Pp, H, v, R, rtk->nx, nv))) {
            errmsg(rtk, "RTKfilter error (info=%d)\n", info);
            stat = SOLQ_NONE;
            break;
        }
        RTKtrace(4, "x(%d)=", i + 1);
        tracemat(4, xp, 1, NR(opt), 13, 4);
    }
    //	RTKtrace(1,"BFE x=");tracemat(1,xp,1,rtk->nx,13,4);
    //	RTKtrace(1,"BFE Qx=");tracemat(1,Pp,rtk->nx,rtk->nx,13,4);
    
    //	RTKtrace(1,"BFE R=");tracemat(1,R,nv,nv,13,4);
    if (stat != SOLQ_NONE && zdres(0, obs, nu, rs, dts, svh, nav, xp, opt, 0, y,
                                   e, azel)) {
        
        /* post-fit residuals for float solution */
        nv = ddres(rtk, nav, dt, xp, Pp, sat, y, e, azel, iu, ir, ns, v, NULL,
                   R, vflg);
        
        /* validation of float solution */
        if (valpos(rtk, v, R, vflg, nv, 4.0)) {
            
            /* update state and covariance matrix */
            matcpy(rtk->x, xp, rtk->nx, 1);
            matcpy(rtk->P, Pp, rtk->nx, rtk->nx);
            //            RTKtrace(1,"AFT x=");tracemat(1,xp,1,NR(opt),13,4);
            //            RTKtrace(1,"AFT Qx=");tracemat(1,Pp,NR(opt),NR(opt),13,4);
            /* update ambiguity control struct */
            rtk->sol.ns = 0;
            for (i = 0; i < ns; i++)
                for (f = 0; f < nf; f++) {
                    if (!rtk->ssat[sat[i] - 1].vsat[f])
                        continue;
                    rtk->ssat[sat[i] - 1].lock[f]++;
                    rtk->ssat[sat[i] - 1].outc[f] = 0;
                    if (f == 0)
                        rtk->sol.ns++; /* valid satellite count by L1 */
                }
            /* lack of valid satellites */
            if (rtk->sol.ns < 4)
                stat = SOLQ_NONE;
        }
        else
            stat = SOLQ_NONE;
    }
    /*savefloatsolsdualfreq(rtk,bias,xa);    */
    savefloatsols(rtk,bias,xa);/* add by yanqing to save float solutions*/
    /* resolve integer ambiguity by WL-NL */
    /*rtk->opt.modear=ARMODE_WLNL;for debug*/
    if (stat != SOLQ_NONE && rtk->opt.modear == ARMODE_WLNL) {
        
        if (resamb_WLNL(rtk, obs, sat, iu, ir, ns, nav, azel)) {
            stat = SOLQ_FIX;
        }
    }
    /* resolve integer ambiguity by TCAR */
    else if (stat != SOLQ_NONE && rtk->opt.modear == ARMODE_TCAR) {
        
        if (resamb_TCAR(rtk, obs, sat, iu, ir, ns, nav, azel)) {
            stat = SOLQ_FIX;
        }
    }
    /* resolve integer ambiguity by LAMBDA */
    // else if (stat != SOLQ_NONE && resamb_TSRC(rtk, bias, xa) > 1) {
    else if (stat != SOLQ_NONE && resamb_LAMBDA(rtk, bias, xa) > 1) {
        if (zdres(0, obs, nu, rs, dts, svh, nav, xa, opt, 0, y, e, azel)) {
            
            /* post-fit reisiduals for fixed solution */
            nv = ddres(rtk, nav, dt, xa, NULL, sat, y, e, azel, iu, ir, ns, v,
                       NULL, R, vflg);
            
            /* validation of fixed solution */
            if (valpos(rtk, v, R, vflg, nv, 4.0)) {
                
                /* hold integer ambiguity */
                if (++rtk->nfix >=
                    rtk->opt.minfix && rtk->opt.modear == ARMODE_FIXHOLD) {
                    holdamb(rtk, xa);
                }
                stat = SOLQ_FIX;
            }
            else{
                printf("fixed pos not valid, icount=%d\n",icount);
            }
        }
    }
    /* save solution status*/
    if (stat == SOLQ_FIX) {
        for (i = 0; i < 3; i++) {
            rtk->sol.rr[i] = rtk->xa[i];
            rtk->sol.qr[i] = (float)rtk->Pa[i + i * rtk->na];
        }
        rtk->sol.qr[3] = (float)rtk->Pa[1];
        rtk->sol.qr[4] = (float)rtk->Pa[1 + 2 * rtk->na];
        rtk->sol.qr[5] = (float)rtk->Pa[2];
    }
    else {
        for (i = 0; i < 3; i++) {
            rtk->sol.rr[i] = rtk->x[i];
            rtk->sol.qr[i] = (float)rtk->P[i + i * rtk->nx];
        }
        rtk->sol.qr[3] = (float)rtk->P[1];
        rtk->sol.qr[4] = (float)rtk->P[1 + 2 * rtk->nx];
        rtk->sol.qr[5] = (float)rtk->P[2];
        rtk->nfix = 0;
    }
    for (i = 0; i < n; i++)
        for (j = 0; j < nf; j++) {
            if (obs[i].L[j] == 0.0)
                continue;
            rtk->ssat[obs[i].sat - 1].pt[obs[i].rcv - 1][j] = obs[i].time;
            rtk->ssat[obs[i].sat - 1].ph[obs[i].rcv - 1][j] = obs[i].L[j];
        }
    for (i = 0; i < ns; i++)
        for (j = 0; j < nf; j++) {
            
            /* output snr of rover receiver */
            rtk->ssat[sat[i] - 1].snr[j] = obs[iu[i]].SNR[j];
        }
    for (i = 0; i < MAXSAT; i++)
        for (j = 0; j < nf; j++) {
            if (rtk->ssat[i].fix[j] == 2 && stat != SOLQ_FIX)
                rtk->ssat[i].fix[j] = 1;
            if (rtk->ssat[i].slip[j]&1)
                rtk->ssat[i].slipc[j]++;
        }
    free(rs);
    free(dts);
    free(var);
    free(y);
    free(e);
    free(azel);
    free(xp);
    free(Pp);
    free(xa);
    free(v);
    free(H);
    free(R);
    free(bias);
    
    if (stat != SOLQ_NONE)
        rtk->sol.stat = stat;
    
    return stat != SOLQ_NONE;
}

/* initialize rtk control ------------------------------------------------------
 * initialize rtk control struct
 * args   : rtk_t    *rtk    IO  rtk control/result struct
 *          prcopt_t *opt    I   positioning options (see rtklib.h)
 * return : none
 *----------------------------------------------------------------------------- */
extern void rtkinit(rtk_t *rtk, const prcopt_t *opt) {
    sol_t sol0 = { {0}};
    ambc_t ambc0 = { { {0}}};
    ssat_t ssat0 = {0};
    int i;
    
    RTKtrace(3, "rtkinit :\n");
    
    rtk->sol = sol0;
    for (i = 0; i < 6; i++)
        rtk->rb[i] = 0.0;
    rtk->nx = opt->mode <= PMODE_FIXED ? NX(opt) : pppnx(opt);
    rtk->na = opt->mode <= PMODE_FIXED ? NR(opt) : 0;
    rtk->tt = 0.0;
    rtk->x = zeros(rtk->nx, 1);
    rtk->P = zeros(rtk->nx, rtk->nx);
    rtk->xa = zeros(rtk->na, 1);
    rtk->Pa = zeros(rtk->na, rtk->na);
    rtk->nfix = rtk->neb = 0;
    for (i = 0; i < MAXSAT; i++) {
        rtk->ambc[i] = ambc0;
        rtk->ssat[i] = ssat0;
    }
    for (i = 0; i < MAXERRMSG; i++)
        rtk->errbuf[i] = 0;
    rtk->opt = *opt;
}

/* free rtk control ------------------------------------------------------------
 * free memory for rtk control struct
 * args   : rtk_t    *rtk    IO  rtk control/result struct
 * return : none
 *----------------------------------------------------------------------------- */
extern void rtkfree(rtk_t *rtk) {
    RTKtrace(3, "rtkfree :\n");
    
    rtk->nx = rtk->na = 0;
    free(rtk->x);
    rtk->x = NULL;
    free(rtk->P);
    rtk->P = NULL;
    free(rtk->xa);
    rtk->xa = NULL;
    free(rtk->Pa);
    rtk->Pa = NULL;
}

/* precise positioning ---------------------------------------------------------
 * input observation data and navigation message, compute rover position by
 * precise positioning
 * args   : rtk_t *rtk       IO  rtk control/result struct
 *            rtk->sol       IO  solution
 *                .time      O   solution time
 *                .rr[]      IO  rover position/velocity
 *                               (I:fixed mode,O:single mode)
 *                .dtr[0]    O   receiver clock bias (s)
 *                .dtr[1]    O   receiver glonass-gps time offset (s)
 *                .Qr[]      O   rover position covarinace
 *                .stat      O   solution status (SOLQ_???)
 *                .ns        O   number of valid satellites
 *                .age       O   age of differential (s)
 *                .ratio     O   ratio factor for ambiguity validation
 *            rtk->rb[]      IO  base station position/velocity
 *                               (I:relative mode,O:moving-base mode)
 *            rtk->nx        I   number of all states
 *            rtk->na        I   number of integer states
 *            rtk->ns        O   number of valid satellite
 *            rtk->tt        O   time difference between current and previous (s)
 *            rtk->x[]       IO  float states pre-RTKfilter and post-RTKfilter
 *            rtk->P[]       IO  float covariance pre-RTKfilter and post-RTKfilter
 *            rtk->xa[]      O   fixed states after AR
 *            rtk->Pa[]      O   fixed covariance after AR
 *            rtk->ssat[s]   IO  sat(s+1) status
 *                .sys       O   system (SYS_???)
 *                .az   [r]  O   azimuth angle   (rad) (r=0:rover,1:base)
 *                .el   [r]  O   elevation angle (rad) (r=0:rover,1:base)
 *                .vs   [r]  O   data valid single     (r=0:rover,1:base)
 *                .resp [f]  O   freq(f+1) pseudorange residual (m)
 *                .resc [f]  O   freq(f+1) carrier-phase residual (m)
 *                .vsat [f]  O   freq(f+1) data vaild (0:invalid,1:valid)
 *                .fix  [f]  O   freq(f+1) ambiguity flag
 *                               (0:nodata,1:float,2:fix,3:hold)
 *                .slip [f]  O   freq(f+1) slip flag
 *                               (bit8-7:rcv1 LLI, bit6-5:rcv2 LLI,
 *                                bit2:parity unknown, bit1:slip)
 *                .lock [f]  IO  freq(f+1) carrier lock count
 *                .outc [f]  IO  freq(f+1) carrier outage count
 *                .slipc[f]  IO  freq(f+1) cycle slip count
 *                .rejc [f]  IO  freq(f+1) data reject count
 *                .gf        IO  geometry-free phase (L1-L2) (m)
 *                .gf2       IO  geometry-free phase (L1-L5) (m)
 *            rtk->nfix      IO  number of continuous fixes of ambiguity
 *            rtk->neb       IO  bytes of error message buffer
 *            rtk->errbuf    IO  error message buffer
 *            rtk->tstr      O   time string for debug
 *            rtk->opt       I   processing options
 *          obsd_t *obs      I   observation data for an epoch
 *                               obs[i].rcv=1:rover,2:reference
 *                               sorted by receiver and satellte
 *          int    n         I   number of observation data
 *          nav_t  *nav      I   navigation messages
 * return : status (0:no solution,1:valid solution)
 * notes  : before calling function, base station position rtk->sol.rb[] should
 *          be properly set for relative mode except for moving-baseline
 *----------------------------------------------------------------------------- */
extern int rtkpos(rtk_t *rtk, const obsd_t *obs, int n, const nav_t *nav) {
    prcopt_t *opt = &rtk->opt;
    sol_t solb = { {0}};
    gtime_t time;
    int i, nu, nr;
    char msg[128] = "";
    
    RTKtrace(3, "rtkpos  : time=%s n=%d\n", time_str(obs[0].time, 3), n);
    RTKtrace(4, "obs=\n");
    traceobs(4, obs, n);
    /* RTKtrace(5,"nav=\n"); tracenav(5,nav); */
    
    /* set base staion position */
    if (opt->refpos <= 3 && opt->mode != PMODE_SINGLE && opt->mode !=
        PMODE_MOVEB) {
        for (i = 0; i < 6; i++)
            rtk->rb[i] = i < 3 ? opt->rb[i] : 0.0;
    }
    /* count rover/base station observations */
    for (nu = 0; nu < n && obs[nu].rcv == 1; nu++);
    for (nr = 0; nu + nr < n && obs[nu + nr].rcv == 2; nr++);
    
    time = rtk->sol.time; /* previous epoch */
    
    /* rover position by single point positioning */
    if (!pntpos(obs, nu, nav, &rtk->opt, &rtk->sol, NULL, rtk->ssat, msg)) {
        errmsg(rtk, "point pos error (%s)\n", msg);
        
        if (!rtk->opt.dynamics) {
            outsolstat(rtk);
            return 0;
        }
    }
    if (time.time != 0)
        rtk->tt = timediff(rtk->sol.time, time);
    
    /* single point positioning */
    if (opt->mode == PMODE_SINGLE) {
        outsolstat(rtk);
        return 1;
    }
    /* precise point positioning */
    if (opt->mode >= PMODE_PPP_KINEMA) {
        pppos(rtk, obs, nu, nav);
        pppoutsolstat(rtk, statlevel, fp_stat);
        return 1;
    }
    /* check number of data of base station and age of differential */
    if (nr == 0) {
        errmsg(rtk, "no base station observation data for rtk\n");
        outsolstat(rtk);
        return 1;
    }
    if (opt->mode == PMODE_MOVEB) { /* moving baseline */
        
        /* estimate position/velocity of base station */
        if (!pntpos(obs + nu, nr, nav, &rtk->opt, &solb, NULL, NULL, msg)) {
            errmsg(rtk, "base station position error (%s)\n", msg);
            return 0;
        }
        rtk->sol.age = (float)timediff(rtk->sol.time, solb.time);
        
        if (fabs(rtk->sol.age) > TTOL_MOVEB) {
            errmsg(rtk, "time sync error for moving-base (age=%.1f)\n",
                   rtk->sol.age);
            return 0;
        }
        for (i = 0; i < 6; i++)
            rtk->rb[i] = solb.rr[i];
        
        /* time-synchronized position of base station */
        for (i = 0; i < 3; i++)
            rtk->rb[i] += rtk->rb[i + 3] * rtk->sol.age;
    }
    else {
        rtk->sol.age = (float)timediff(obs[0].time, obs[nu].time);
        
        if (fabs(rtk->sol.age) > opt->maxtdiff) {
            errmsg(rtk, "age of differential error (age=%.1f)\n", rtk->sol.age);
            outsolstat(rtk);
            return 1;
        }
    }
    /* relative potitioning */
    relpos(rtk, obs, nu, nr, nav);
    outsolstat(rtk);
    
    return 1;
}
