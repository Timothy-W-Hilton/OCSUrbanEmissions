SHELL=/bin/sh

PROGRAM=SiB_2_IOAPI

FC=ftn
LD=ftn

FFLAGS= -O2 -c -C -I$(PROJ)/local/include
LDFLAGS= -openmp -lpthread
APINCL=-I$(PROJ)/local/include
APILIB= -L$(PROJ)/local/lib -lioapi_regrid_tools -lioapi -lnetcdf

CMD=$(PROGRAM).x
SRCDRV=$(PROGRAM).F90
OBJDRV=$(PROGRAM).o array_funcs.o

all:  	$(CMD)

$(CMD):  $(OBJDRV)
	$(LD) $(LDFLAGS) -o $(CMD) $(OBJDRV) \
	$(APILIB)

array_funcs.o:array_funcs.F90
	$(FC) $(FFLAGS) $(APINCL) array_funcs.F90

%.o: %.F90 array_funcs.o
	$(FC) $(FFLAGS) $(APINCL) $<

clean:
	rm -f  $(OBJDRV) *.o *.mod
