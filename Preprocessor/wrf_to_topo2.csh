#!/bin/csh -f
#PBS -q debug
#PBS -l mppwidth=24
#PBS -l walltime=00:20:00
#PBS -S /usr/bin/csh
setenv JOBNAME topo_gen
#PBS -N $JOBNAME
#PBS -o $JOBNAME.$PBS_JOBID.out
#PBS -V
cd $PBS_O_WORKDIR

#PATH=/home/thilton/bin:/usr/local/bin:/opt/pgi/linux86-64/2013/bin:/home/thilton/bin:/usr/local/bin:/opt/pgi/linux86-64/2013/bin:/home/thilton/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:/usr/games
echo "creating object file"
ftn -O2 -c -I$PROJ/local/include -c -o wrf_to_topo2.o wrf_to_topo2.f90

echo "creating executable file"
ftn -o wrf_to_topo2.x wrf_to_topo2.o -L$PROJ/local/lib -lioapi -lnetcdf -ldatetime -openmp

echo "setenv commands"
setenv WRF_FILE /global/homes/g/gara/WRFV3/run/wrfout_d01_2015-05-24_00:00:00
rm wrf_topo2.nc
setenv OUTPUT wrf_topo2-9km.nc
#setenv OUTPUT /mnt/home10/azumkehr/kettle_extraction/ioapi_conversion/kettle_soil_1x1.nc
#setenv EXECUTION_ID 1001B
setenv EXECUTION_ID 1001B
#if (-e $OUTPUT) rm -f $OUTPUT

aprun -n 1 ./wrf_to_topo2.x


