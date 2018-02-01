#!/bin/sh

##################################################
# script to use
# [wrfgriddesc](https://www.cmascenter.org/ioapi/documentation/all_versions/html/WRFGRIDDESC.html)
# to create a
# [GRIDDESC](https://www.cmascenter.org/ioapi/documentation/all_versions/html/GRIDDESC.html)
# file describing the grids for (1) the innermost domain (D04) from
# WRF simulations for Barcelona and (2) the 1.0 by 1.25 degree global
# grid used by SiB.  wrfgriddesc creates the GRIDDESC file directly
# from the grid description in a WRF output file. The GRIDDESC file
# completely describes each coordinate system and grid.  The GRIDDESC
# file is useful for several things, including: (1) regridding SiB
# data to the Barcelona domain using
# [mtxcalc](https://www.cmascenter.org/ioapi/documentation/all_versions/html/MTXCALC.html)
# and
# [mtxcple](https://www.cmascenter.org/ioapi/documentation/all_versions/html/MTXCPLE.html).
# (2) extracting the latitude and longitude coordinates of the
# Barcelona grid using
# [latlon](https://www.cmascenter.org/ioapi/documentation/all_versions/html/LATLON.html).
# (3) If creating data files using the I/O API fortran interface, the
# subroutine
# [DSCGRID](https://www.cmascenter.org/ioapi/documentation/all_versions/html/DSCGRID.html)
# can read the GRIDDESC file and populate the grid description
# parameters P_ALP, P_BET, P_GAM, etc.  This avoids needing to
# hard-code their values in the Fortran code.
#
# REQUIREMENTS
# (1) a working [EDSS/Models-3 I/O
# API](https://www.cmascenter.org/ioapi/documentation/all_versions/html/index.html)
# installation
# (2) wrfgriddesc must be on the search path ($PATH environment variable)
##################################################

# The result: the GRIDDESC file to be created
export GRIDDESC_BCN=$PWD/GRIDDESC_BCN
# the wrfout file whose grid will be described
export WRFOUTPUT=/software/co2flux/ilavrador/WRFOUT_STEM/wrfout_d04_2016-02-01_00:00:00

# wrfgriddesc has a hard-coded limit of 12 for the number of WRF
# dimensions it can handle.  These wrfout files have 14 dimensions.
# Therefore run ncks first to extract a subset of variables that
# together have fewer than MXWRFDIMS==12 dimensions.  Needed to
# include all the below variables and dimensions to get everything
# wrfgriddesc needed.
#
# as of 30 Jan 2018 the PIC installation of Models-3 I/O API requires
# netCDF version 3 input.  Therefore use "-3" flag to ncks.
export TMP_WRFOUT=$PWD/wrfout_d01_subset.nc
ncks -3 -v Times,XLAT,XLONG,T,P,U,V,W,ZS $WRFOUTPUT $TMP_WRFOUT

# set up "logical names"
# (https://www.cmascenter.org/ioapi/documentation/all_versions/html/LOGICALS.html)
# needed by wrfgriddesc.
export WRFFILE=$TMP_WRFOUT  # <path name for input WRF-netcdf file>
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
# now run wrfgriddesc, answering "Y" at the "continue?" prompt
wrfgriddesc <<EOF
Y
EOF

# we're done with the temporary WRFOUT file, so remove it
rm -fv $TMP_WRFOUT

# now add description of the native SiB grid into the GRIDDESC file
# define the SiB lat/lon coordinate system after line 2 of GRIDDESC
# awk cannot rewrite a file in place so use a temporary file tmp_out
awk -v n=2 -v s=" 'LAT_LON'\n     1,0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0" 'NR == n {print s} {print}' $GRIDDESC_BCN > tmp_out
# define the global 1.0 by 1.25 degree grid after line 17 of GRIDDESC
# awk cannot rewrite a file in place so use a temporary file tmp_out2
awk -v n=17 -v s=" 'SIB_NATIVE'\n 'LAT_LON' -180, -90, 1.25, 1., 288, 181, 1" 'NR == n {print s} {print}' tmp_out > tmp_out2

# clean up
rm -f tmp_out
mv tmp_out2 $GRIDDESC_BCN && echo "created $GRIDDESC_BCN"
