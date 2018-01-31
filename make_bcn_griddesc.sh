#!/bin/sh

export TMP_WRFOUT=$PWD/wrfout_d01_subset.nc
export GRIDDESC_BCN=$PWD/GRIDDESC_BCN

# wrfgriddesc has a hard-coded limit of 12 for the number of WRF dimensions it can handle.  These wrfout files have 14 dimensions.  Therefore ran ncks first to extract a subset with fewer than MXWRFDIMS==12
# dimensions.  Needed to include all these variables and dimensions to
# get everything wrfgriddesc needed.
# as of 30 Jan 2018 the PIC installation of Models-3 I/O API requires
# netCDF version 3 input.  Therefore use "-3" flag to ncks.
ncks -3 -v Times,XLAT,XLONG,T,P,U,V,W,ZS /software/co2flux/ilavrador/WRFOUT_STEM/wrfout_d04_2016-02-01_00:00:00 $TMP_WRFOUT

# WRFFILE <path name for input WRF-netcdf file>
export WRFFILE=$TMP_WRFOUT
export OUTDESC=$GRIDDESC_BCN # <path name for output GRIDDESC file>
export CRDNAME=CRD_BCND04  # <GRIDDESC name for WRF coordinate system>
export CROGRID=CRO_BCND04  # <GRIDDESC name for     cross-point grid>
export DOTGRID=DOT_BCND04  # <GRIDDESC name for       dot-point grid>
export STXGRID=STX_BCND04  # <GRIDDESC name for X-stagger-point grid>
export STYGRID=STY_BCND04  # <GRIDDESC name for Y-stagger-point grid>

# remove GRIDDESC if it already exists
if [ -e "$OUTDESC" ]
then
   rm -fv $OUTDESC
fi
# now run wrfgriddesc, answering "Y" at the "continue" prompt
wrfgriddesc <<EOF
Y
EOF

rm -fv $TMP_WRFOUT

# now add description of the native SiB grid into the GRIDDESC file
awk -v n=2 -v s=" 'LAT_LON'\n     1,0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0" 'NR == n {print s} {print}' $GRIDDESC_BCN > tmp_out
awk -v n=17 -v s=" 'SIB_NATIVE'\n 'LAT_LON' -180, -90, 1.25, 1., 288, 181, 1" 'NR == n {print s} {print}' tmp_out > tmp_out2

# clean up
rm -f tmp_out
mv tmp_out2 $GRIDDESC_BCN
