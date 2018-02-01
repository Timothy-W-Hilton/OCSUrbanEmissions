SHELL=/bin/sh

PROGRAM=preprocessor_1km

FC=gfortran
LD=gfortran

FFLAGS   = -c -O $(APINCL) $(STEMINCL) -g -fopenmp -Wtabs
LDFLAGS  = -O -g77libs -g -fopenmp -Wtabs
APINCL=-I/usr/include -I/software/co2flux/LIBRARIES/netcdf-3.6.3/include -I/nfs/pic.es/user/t/thilton/Software/ioapi-3.2_tim/Linux2_x86_64gfort
APILIB=-L/usr/lib64 -L/nfs/pic.es/user/t/thilton/Software/ioapi-3.2_tim/Linux2_x86_64gfort -L$(PWD) -lioapi -L/software/co2flux/LIBRARIES/netcdf-3.6.3/lib -lnetcdf

CMD=$(PROGRAM).x
SRCDRV=$(PROGRAM).F90
OBJDRV=$(PROGRAM).o

all:  	$(CMD)

$(CMD):	$(OBJDRV)
	$(LD) $(LDFLAGS) -o $(CMD) $(OBJDRV) \
	$(APILIB)

f.o:
	$(FC) $(FFLAGS) $(APINCL) $<

f90.o:
	$(FC) $(FFLAGS) $(APINCL) $<

clean:
	rm -f  $(OBJDRV) *.o *.mod $(CMD)
