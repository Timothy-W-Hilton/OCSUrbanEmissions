SHELL=/bin/sh

PROGRAM=SiB_2_IOAPI

FC=gfortran
LD=gfortran

FFLAGS   = -c -O $(APINCL) $(STEMINCL) -g -fopenmp -Wtabs
LDFLAGS  = -O -g77libs -g -fopenmp -Wtabs
APINCL=-I/usr/include -I/software/co2flux/LIBRARIES/netcdf-3.6.3/include -I/nfs/pic.es/user/t/thilton/Software/ioapi-3.2_tim/Linux2_x86_64gfort
APILIB=-L/usr/lib64 -L/nfs/pic.es/user/t/thilton/Software/ioapi-3.2_tim/Linux2_x86_64gfort -L$(PWD) -lioapi_regrid_tools -lioapi -L/software/co2flux/LIBRARIES/netcdf-3.6.3/lib -lnetcdf

CMD=$(PROGRAM).x
SRCDRV=$(PROGRAM).F90
OBJDRV=$(PROGRAM).o array_funcs.o

all:  	$(CMD)

$(CMD):	lib $(OBJDRV)
	$(LD) $(LDFLAGS) -o $(CMD) $(OBJDRV) \
	$(APILIB)

lib:	ioapi_regrid_tools.o
	ar -rcs libioapi_regrid_tools.a ioapi_regrid_tools.o

ioapi_regrid_tools.o:	ioapi_regrid_tools.F90
	$(FC) $(FFLAGS) ioapi_regrid_tools.F90 $(APILIB)

array_funcs.o:array_funcs.F90
	$(FC) $(FFLAGS) $(APINCL) array_funcs.F90

%.o: %.F90 array_funcs.o
	$(FC) $(FFLAGS) $(APINCL) $<

clean:
	rm -f  $(OBJDRV) *.o *.mod
