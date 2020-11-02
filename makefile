  FFLAGS=  -g -C  +es +Obb1000 -K +FPZOU
  ABSOFT=/apps/absoft/absoft-8.2/opt/absoft/
  FABSFLAGS=-O -V -W -f -s -N1 -B108 -B100 -N90 -N22 -N2 -N113
  INCLUDES=-I.
  EXTRAFLAGS=-DABSOFTFORTRAN
  FFLAGS= $(INCLUDES) $(FABSFLAGS) $(EXTRAFLAGS)
  FFLAG1=$(FFLAGS) -c
#  GFORTRAN :=$(ABSOFT)/bin/gfortran
  GFORTRAN :=gfortran
  direct=/group/c-xem2/cmorean/SCGSR/cs_model/
#CERN_ROOT=/usr/local/cern/pro
#CERN_ROOT=/usr/lib/cernlib/2006
CERN_ROOT = /apps/cernlib/i386_rhel3/2003
CERNLIBS = -L$(CERN_ROOT)/lib -lpawlib -lgraflib -lgrafX11 -lpacklib -lmathlib


xemmodel:  xem_model.o nform.o bdisnew4he3.o f1f220.o ./src/constants.inc
	gfortran -g -o XEM_model xem_model.o nform.o bdisnew4he3.o f1f220.o
	rm *.o
xem_model.o: ./src/xem_model.f
	gfortran -c -ffixed-line-length-none ./src/xem_model.f
nform.o: ./src/nform.f
	gfortran -c -ffixed-line-length-none ./src/nform.f
bdisnew4he3.o: ./src/bdisnew4he3.f
	gfortran -c -ffixed-line-length-none ./src/bdisnew4he3.f
f1f220.o: ./src/f1f220.f
	gfortran -c -ffixed-line-length-none ./src/f1f220.f


