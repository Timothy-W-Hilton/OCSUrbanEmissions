#!/bin/sh

export TMP_WRFOUT=$PWD/wrfout_d01_subset.nc
export GRIDDESC_BCN=$PWD/GRIDDESC_BCN

# wrfgriddesc has a hard-coded limit of 12 for the number of WRF dimensions it can handle.  These wrfout files have 14 dimensions.  Therefore ran ncks first to extract a subset with fewer than MXWRFDIMS==12
# dimensions.  Needed to include all these variables and dimensions to
# get everything wrfgriddesc needed.
# as of 30 Jan 2018 the PIC installation of Models-3 I/O API requires
# netCDF version 3 input.  Therefore use "-3" flag to ncks.
ncks -3 -o -v Times,XLAT,XLONG,T,P,U,V,W,ZS /software/co2flux/ilavrador/WRFOUT_STEM/wrfout_d01_2016-02-01_00:00:00 $TMP_WRFOUT

# WRFFILE <path name for input WRF-netcdf file>
export WRFFILE=$TMP_WRFOUT
export OUTDESC=$GRIDDESC_BCN # <path name for output GRIDDESC file>
export CRDNAME=CRD_BCND01  # <GRIDDESC name for WRF coordinate system>
export CROGRID=CRO_BCND01  # <GRIDDESC name for     cross-point grid>
export DOTGRID=DOT_BCND01  # <GRIDDESC name for       dot-point grid>
export STXGRID=STX_BCND01  # <GRIDDESC name for X-stagger-point grid>
export STYGRID=STY_BCND01  # <GRIDDESC name for Y-stagger-point grid>

# now run wrfgriddesc
wrfgriddesc

rm -fv $TMP_WRFOUT
