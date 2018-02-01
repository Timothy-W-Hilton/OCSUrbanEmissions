## make_GRIDDESC/ ##

The subdirectory make_GRIDDESC/ contains code to create a [GRIDDESC](https://www.cmascenter.org/ioapi/documentation/all_versions/html/GRIDDESC.html) file describing the grids for (1) the innermost domain (D04) from WRF simulations for Barcelona and (2) the 1.0 by 1.25 degree global grid used by SiB.

make_bcn_griddesc.sh is the executable script that creates the GRIDDESC file.  It contains a lot of comments describing how it works. GRIDDESC_BCN is the output.

### Requirements ###
1. a working [EDSS/Models-3 I/O API](https://www.cmascenter.org/ioapi/documentation/all_versions/html/index.html) installation
2. wrfgriddesc must be on the search path ($PATH environment variable)
