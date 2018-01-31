#!/bin/csh -f
#PBS -q debug
#PBS -l mppwidth=24
#PBS -l walltime=00:20:00
#PBS -S /usr/bin/csh
setenv JOBNAME wrf_2_stem_input
#PBS -N $JOBNAME
#PBS -o $JOBNAME.$PBS_JOBID.out
#PBS -V
cd $PBS_O_WORKDIR





#PATH=/home/thilton/bin:/usr/local/bin:/opt/pgi/linux86-64/2013/bin:/home/thilton/bin:/usr/local/bin:/opt/pgi/linux86-64/2013/bin:/home/thilton/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin:/usr/games
echo "creating object file"
ftn -O2 -c -I$PROJ/local/include -c -o preprocessor_1km.o preprocessor_1km.f

echo "creating executable file"
ftn -o preprocessor_1km.x preprocessor_1km.o -L$PROJ/local/lib -lioapi -lnetcdf -ldatetime -openmp

setenv EXECUTION_ID 1001B
#read -n1 -r -p "Press ctrl to abort, press y to continue" key 
set req
echo "Continue?"
set req = $<
#if [[ "$key" = 'y' ]]; then
#    printf "\ncontinue?\n"
#   ./preprocessor_1km.x
#   echo 'Program ended\n'
#fi



rm tmp.mm53nc.ini

cat > tmp.mm53nc.ini <<EOF
&control
 wrffile='/global/homes/g/gara/WRFV3/run/wrfout_d01_2015-05-24_00:00:00'
 begyear= 2015    
 begdate= 144 
 begtime= 00 
 dtstep = 6     
 wdiog =.false.   
 cdiog =.false.   
 numts=  5
 initialfile=.false.
&end
 
USGS-category     RADM_1    FACTOR_1   RADM_2   FACTOR_2
  1                 1          1.         0        0.       ! URBAN LAND
  2                 3          1.         0        0.       ! Dryland Cropland and Pasture
  3                 3          1.         0        0.       ! Irrigated Cropland and Pasture
  4                 3         0.5         8       0.5       ! Mixed Dryland/Irr. Cropl. and Pasture
  5                10          1.         0        0.       ! Cropland/Grassland Mosaic
  6                10          1.         0        0.       ! Cropland/Woodland Mosaic
  7                 3          1.         0        0.       ! Grassland
  8                11          1.         0        0.       ! Shrubland
  9                11          1.         0        0.       ! Mixed Shrubland/Grassland
 10                 3          1.         0        0.       ! Savanna
 11                 4          1.         0        0.       ! Deciduous Broadleaf Forest
 12                 5          1.         0        0.       ! Deciduous Needleleaf Forest
 13                 4          1.         0        0.       ! Evergreen Broadleaf Forest
 14                 5          1.         0        0.       ! Evergreen Needleleaf Forest
 15                 6          1.         0        0.       ! Mixed Forest
 16                 7          1.         0        0.       ! Water Bodies
 17                 9          1.         0        0.       ! Herbaceous Wetland
 18                 6          1.         0        0.       ! Wooded Wetland
 19                 8          1.         0        0.       ! Barren or Sparsely Vegetated
 20                 8         0.5        11       0.5       ! Herbaceous Tundra
 21                 6         0.5         9       0.5       ! Wooded Tundra
 22                 9          1.         0        0.       ! Mixed Tundra
 23                 8          1.         0        0.       ! Bare Ground Tundra
 24                 7          1.         0        0.       ! Snow or Ice
EOF
# Change below
#/usr/local/bin/namecut tmp.mm52nc.ini brazil.ini 
#rm tmp.mm52nc.ini

#setenv TOPO TOPO-preprocessor_1km.nc
#setenv TOPO "/home/ecampbell_lab/ftpanon.al.noaa.gov/wrfout_d03_2010-05-25_00_R23_fpcut.nc -v"
setenv METEO3D "meteo3d-preprocessor_9km.nc"
setenv TOPO "wrf_topo2-9km.nc -v"

setenv METEO2D "meteo2d-preprocessor_9km.nc -v"
setenv DRYDEP  "drydep-preprocessor_9km.nc -v"
setenv HEIGHT3D "wrfheight-preprocessor_9km.nc"


    rm $METEO2D
    rm $DRYDEP
    rm $HEIGHT3D
    rm $METEO3D

aprun -n 1 ./preprocessor_1km.x #> print-1.out


exit()


