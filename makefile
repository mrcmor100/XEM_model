ifeq ($(MYOS),Linux)
#  CERN_ROOT = /apps/cernlib/i386_fc8/2005
  FFLAGSA=-O -ffixed-line-length-132 -ff2c -fno-automatic -fdefault-real-8
#  FFLAGSA=-O -V -W -f -s -N1 -B108 -B100 -N90 -N22 -N2 -N113
#  INCLUDES=-I.,..,./sos,$(SH),./hrsr,./hrsl,./shms,$(Csoft)/SRC/INCLUDE
#  INCLUDES=-I.,..,./sos,$(SH),./hrsr,./hrsl,./shms,$(Csoft)/INCLUDE
  FFLAGS= $(INCLUDES) $(FFLAGSA)
  FFLAG1=$(FFLAGS) -c
  OTHERLIBS =  -L$(CERN_ROOT)/lib $(CERNLIBS) -L/usr/lib
  FC  := gfortran
  F77 :=gfortran
endif

xemmodel:  xem_model.o nform.o smear4all.o f1f221.o ./src/constants.inc
	gfortran -g -o XEM_model xem_model.o nform.o smear4all.o f1f221.o 
	$(OTHERLIBS)
	rm *.o
xem_model.o: ./src/xem_model.f
	gfortran -c -ffixed-line-length-none ./src/xem_model.f
nform.o: ./src/nform.f
	gfortran -c -ffixed-line-length-none ./src/nform.f
smear4all.o: ./src/smear4all.f
	gfortran -c -ffixed-line-length-none ./src/smear4all.f
f1f221.o: ./src/f1f221.f
	gfortran -c -ffixed-line-length-none ./src/f1f221.f


