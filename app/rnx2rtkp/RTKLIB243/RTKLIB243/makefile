# RTKLIB makefile for rnx2rtkp in MacOs

BINDIR  = /usr/local/bin
SRC     = ../../../../src
DATADIR = /Users/YanqingHou/Documents/Work/datacenter/BJFSAOE/expAOE/MMDD0823
DavideData=/Users/YanqingHou/Documents/Work/datacenter/DataCurtin
DYData=/Users/YanqingHou/Documents/Work/datacenter/dongyihangkong

CUTDATA=/Users/YanqingHou/Documents/Work/datacenter/DataCurtin/CUT


OPTS    = -DTRACE -DENAGLO -DENAQZS -DENAGAL -DENACMP -DNFREQ=3


# for no lapack
# CFLAGS  = -Wall -O3 -ansi -pedantic -Wno-unused-variable -I$(SRC) $(OPTS) -g
# LDLIBS  = -lm -lrt
# for lapack
CFLAGS  = -Wall -O3 -ansi -pedantic -Wno-empty-body -Wno-unused-function -Wno-unused-variable -I$(SRC) -DLAPACK $(OPTS)
LDLIBS  = -lm -llapack -lblas

# for gprof 
#CFLAGS  = -Wall -O3 -ansi -pedantic -Wno-unused-but-set-variable -I$(SRC) -DLAPACK $(OPTS) -pg
#LDLIBS  = -lm -lrt -llapack -lblas -pg

# for mkl
#MKLDIR  = /opt/intel/mkl
#CFLAGS  = -O3 -ansi -pedantic -Wno-unused-but-set-variable -I$(SRC) $(OPTS) -DMKL
#LDLIBS  = -L$(MKLDIR)/lib/intel64 -lm -lrt -lmkl_core -lmkl_intel_lp64 -lmkl_gnu_thread -liomp5 -lpthread

all        : rnx2rtkp
rnx2rtkp   : rnx2rtkp.o rtkcmn.o rinex.o rtkpos.o postpos.o solution.o
rnx2rtkp   : lambda.o geoid.o sbas.o preceph.o pntpos.o ephemeris.o options.o
rnx2rtkp   : ppp.o ppp_ar.o rtcm.o rtcm2.o rtcm3.o rtcm3e.o ionex.o qzslex.o

rnx2rtkp.o : rnx2rtkp.c
	$(CC) -c $(CFLAGS) rnx2rtkp.c
rtkcmn.o   : $(SRC)/rtkcmn.c
	$(CC) -c $(CFLAGS) $(SRC)/rtkcmn.c
rinex.o    : $(SRC)/rinex.c
	$(CC) -c $(CFLAGS) $(SRC)/rinex.c
rtkpos.o   : $(SRC)/rtkpos.c
	$(CC) -c $(CFLAGS) $(SRC)/rtkpos.c
postpos.o  : $(SRC)/postpos.c
	$(CC) -c $(CFLAGS) $(SRC)/postpos.c
solution.o : $(SRC)/solution.c
	$(CC) -c $(CFLAGS) $(SRC)/solution.c
lambda.o   : $(SRC)/lambda.c
	$(CC) -c $(CFLAGS) $(SRC)/lambda.c
geoid.o    : $(SRC)/geoid.c
	$(CC) -c $(CFLAGS) $(SRC)/geoid.c
sbas.o     : $(SRC)/sbas.c
	$(CC) -c $(CFLAGS) $(SRC)/sbas.c
preceph.o  : $(SRC)/preceph.c
	$(CC) -c $(CFLAGS) $(SRC)/preceph.c
pntpos.o   : $(SRC)/pntpos.c
	$(CC) -c $(CFLAGS) $(SRC)/pntpos.c
ephemeris.o: $(SRC)/ephemeris.c
	$(CC) -c $(CFLAGS) $(SRC)/ephemeris.c
options.o  : $(SRC)/options.c
	$(CC) -c $(CFLAGS) $(SRC)/options.c
ppp.o      : $(SRC)/ppp.c
	$(CC) -c $(CFLAGS) $(SRC)/ppp.c
ppp_ar.o   : $(SRC)/ppp_ar.c
	$(CC) -c $(CFLAGS) $(SRC)/ppp_ar.c
rtcm.o     : $(SRC)/rtcm.c
	$(CC) -c $(CFLAGS) $(SRC)/rtcm.c
rtcm2.o    : $(SRC)/rtcm2.c
	$(CC) -c $(CFLAGS) $(SRC)/rtcm2.c
rtcm3.o    : $(SRC)/rtcm3.c
	$(CC) -c $(CFLAGS) $(SRC)/rtcm3.c
rtcm3e.o   : $(SRC)/rtcm3e.c
	$(CC) -c $(CFLAGS) $(SRC)/rtcm3e.c
ionex.o    : $(SRC)/ionex.c
	$(CC) -c $(CFLAGS) $(SRC)/ionex.c
qzslex.o   : $(SRC)/qzslex.c
	$(CC) -c $(CFLAGS) $(SRC)/qzslex.c

rnx2rtkp.o : $(SRC)/rtklib.h
rtkcmn.o   : $(SRC)/rtklib.h
rinex.o    : $(SRC)/rtklib.h
rtkpos.o   : $(SRC)/rtklib.h
postpos.o  : $(SRC)/rtklib.h
solution.o : $(SRC)/rtklib.h
lambda.o   : $(SRC)/rtklib.h
geoid.o    : $(SRC)/rtklib.h
sbas.o     : $(SRC)/rtklib.h
preceph.o  : $(SRC)/rtklib.h
pntpos.o   : $(SRC)/rtklib.h
ephemeris.o: $(SRC)/rtklib.h
options.o  : $(SRC)/rtklib.h
ppp.o      : $(SRC)/rtklib.h
ppp_ar.o   : $(SRC)/rtklib.h
rtcm.o     : $(SRC)/rtklib.h
rtcm2.o    : $(SRC)/rtklib.h
rtcm3.o    : $(SRC)/rtklib.h
rtcm3e.o   : $(SRC)/rtklib.h
ionex.o    : $(SRC)/rtklib.h
qzslex.o   : $(SRC)/rtklib.h

clean :
	rm *.o

CMD    = ./rnx2rtkp
INPUT11 = $(DATADIR)/BJFS823.16O
INPUT12 = $(DATADIR)/AOE0823.16O
INPUT13 = $(DATADIR)/gfzs2360.16p 
INPUT14 = $(DATADIR)/gbm19112.sp3 
INPUT15 = $(DATADIR)/gbm19112.clk

CONF1    = $(DATADIR)/pppaoe0823.conf
CONF2    = $(DATADIR)/pppaoe0823gpsonly.conf
CONF3    = $(DATADIR)/pppaoe0823bdsonly.conf
CONF12    = $(DATADIR)/pppaoe0823nobdsgeo.conf

#Davide data
INPUT21=$(DavideData)/hksc022m.18o
INPUT22=$(DavideData)/PERT.18p
INPUT23=$(DavideData)/static_ublox.obs
CONF4=$(DavideData)/hkscspp.conf
CONF5=$(DavideData)/ublox.conf

# dongyi aviation data
INPUT31=$(DYData)/Rover.obs
INPUT32=$(DYData)/Base.obs
INPUT33=$(DYData)/Rover.nav
CONF34=$(DYData)/aviation.conf
CONF35=$(DYData)/basespp.conf
TS=2018/07/27 12:02:32.1

#cutab data
INPUT41=$(CUTDATA)/CUTA.18o
INPUT42=$(CUTDATA)/CUTB.18o
INPUT43=$(CUTDATA)/CUTA.18p
INPUT44=$(CUTDATA)/gbm20210.sp3
INPUT45=$(CUTDATA)/gbm20210.clk
CONF45=$(CUTDATA)/cutab.conf
testgb:
	$(CMD) -k $(CONF1)  $(INPUT11) $(INPUT12) $(INPUT13) $(INPUT14) $(INPUT15) -y 0 -x 0 -o $(DATADIR)/bjfs823mkgb.pos
testgbnogeo:
	$(CMD) -k $(CONF12)  $(INPUT11) $(INPUT12) $(INPUT13) $(INPUT14) $(INPUT15) -y 0 -x 0 -o $(DATADIR)/bjfs823mkgbnogeo.pos
testg:
	$(CMD) -k $(CONF2)  $(INPUT11) $(INPUT12) $(INPUT13) $(INPUT14) $(INPUT15) -y 0 -x 0 -o $(DATADIR)/bjfs823mkg.pos
testb:
	$(CMD) -k $(CONF3)  $(INPUT11) $(INPUT12) $(INPUT13) $(INPUT14) $(INPUT15) -y 0 -x 0 -o $(DATADIR)/bjfs823mkb.pos
testdvd:
	$(CMD) -k $(CONF5)  $(INPUT23) $(INPUT22) -y 0 -x 0 -o $(DavideData)/davideublox.pos
testdy:
	$(CMD) -k $(CONF34) $(INPUT31) $(INPUT32) $(INPUT33) -ts $(TS) -x 1 -o $(DYData)/DYaviation.pos
testdyspp:
	$(CMD) -k $(CONF35) $(INPUT32) $(INPUT33) -y 0 -x 0 -o $(DYData)/DYbasespp.pos
testcutab:
	$(CMD) -k $(CONF45) $(INPUT42) $(INPUT41) $(INPUT43) $(INPUT44) $(INPUT45) -y 0 -x 0 -o $(CUTDATA)/cutab.pos
install :
	cp rnx2rtkp $(BINDIR)