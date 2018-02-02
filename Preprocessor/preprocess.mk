SHELL=/bin/sh

PROGRAM=preprocessor_1km
TOPO=wrf_to_topo2

FC=gfortran
LD=gfortran

FFLAGS   = -c -O $(APINCL) $(STEMINCL) -g -fopenmp -Wtabs
LDFLAGS  = -O0 -g77libs -g -fopenmp -Wtabs
APINCL=-I/usr/include -I/software/co2flux/LIBRARIES/netcdf-4.1.3/include -I/software/co2flux/LIBRARIES/ioapi-3.2/Linux2_x86_64gfort/
APILIB=-L/usr/lib64 -L/software/co2flux/LIBRARIES/ioapi-3.2/Linux2_x86_64gfort/ -lioapi -L/software/co2flux/LIBRARIES/netcdf-4.1.3/lib -lnetcdff

CMD=$(PROGRAM).x
TOPOCMD=$(TOPO).x
SRCDRV=$(PROGRAM).F90
OBJDRV=$(PROGRAM).o

all:  	$(CMD) $(TOPOCMD)

$(CMD):	$(OBJDRV)
	$(LD) $(LDFLAGS) -o $(CMD) $(OBJDRV) \
	$(APILIB)

$(TOPOCMD):	$(TOPO).o
	$(LD) $(LDFLAGS) -o $(TOPOCMD) $(TOPO).o \
	$(APILIB)

f.o:
	$(FC) $(FFLAGS) $(APINCL) $<

$(TOPO).o:	$(TOPO).f90
	$(FC) $(FFLAGS) $(APINCL) $<

clean:
	rm -f  $(OBJDRV) *.o *.mod $(CMD)
