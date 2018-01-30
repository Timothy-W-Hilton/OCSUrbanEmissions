MODULE ioapi_regrid_tools

!======================================================================>
!
! MODULE: ioapi_regrid_tools
!
! VERSION: 0.0.0
!
! LAST UPDATE: 2014-06-11
!
! AUTHOR: Timothy W. Hilton
!         University of California-Merced
!         e-mail: thilton@ucmerced.edu
!
! DESCRIPTION: A Fortran module that provides several subroutines
!              useful for regridding geophysical data using the I/O API
!              (http://www.cmascenter.org/ioapi/documentation/3.1/html/).
!
! DISTRIBUTING:
!     Includes code written by Sarika Kulkarni
!     <sarika-kulkarni@uiowa.edu> and Y. F. Cheng
!     <yafang.cheng@gmail.com>.  Sarika must be notified before this
!     code is shared outside of the Campbell group.
!
! CONTAINS:
!
!     PUBLIC SUBROUTINES & FUNCTIONS:
!
!         SUBROUTINE describe_netcdf_dims
!         SUBROUTINE calc_grid_area
!         FUNCTION fixpole
!         FUNCTION units_convert
!         FUNCTION calc_fcos
!         SUBROUTINE per_cell_2_per_m2
!         SUBROUTINE open3_and_desc3
!
! REQUIRES
!      - I/O API libraries version 3.1:
!        http://www.cmascenter.org/ioapi/documentation/3.1/html/AVAIL.html
!      - netcdf libraries
!
!======================================================================>

CONTAINS
  !******************************************************************************
  SUBROUTINE describe_netcdf_dims(fname, varname)
    !----------------------------------------------------------------------
    ! DESC: return dimension names and lengths for a variable in a netCDF file
    !----------------------------------------------------------------------
    ! INPUT:
    ! fname: full path to the netCDF file
    ! varname: name of the variable to be described
    !----------------------------------------------------------------------
    ! OUTPUT:
    ! garea - grid area in m2
    !----------------------------------------------------------------------
    ! HISTORY:
    ! By: Timothy W. Hilton (thilton@ucmerced.edu)
    ! On: 23 Sep 2015
    !----------------------------------------------------------------------
    USE netcdf
    IMPLICIT NONE

    CHARACTER(len=512), INTENT(in)::fname
    CHARACTER(len=*), INTENT(in)::varname
    INTEGER :: status, ncid, dimlen, i, &
         COSVarId,       &
         numDims, numAtts
    INTEGER, DIMENSION(nf90_max_var_dims) :: COSDimIds
    CHARACTER(len = nf90_max_name) :: RecordDimName

    print *, "nf90_max_name:", nf90_max_name

    status = nf90_open(fname, nf90_NoWrite, ncid)
    status = nf90_inq_varid(ncid, varname, COSVarId)
    status = nf90_inquire_variable(ncid, COSVarId, ndims = numDims, natts = numAtts)
    status = nf90_inquire_variable(ncid, COSVarId, dimids = COSDimIds(:numDims))
    PRINT *, 'dimIDs: ', COSDimIds(1:numDims)
    DO i = 1, numDims
       status = nf90_inquire_dimension(ncid, COSDimIds(i), &
            name = RecordDimName, len = dimlen)
       PRINT*, 'dim id: ', COSDimIds(i), &
            'name:', trim(RecordDimName), &
            '    size:', dimlen
    ENDDO
  END SUBROUTINE describe_netcdf_dims

  !******************************************************************************
  SUBROUTINE calc_grid_area( garea, xlon1, xlon2, ylat1, ylat2 )
    !----------------------------------------------------------------------
    ! DESC: calculate grid area as a function of longitude and latitude
    !----------------------------------------------------------------------
    ! INPUT:
    ! xlon1 - starting longitude
    ! xlon2 - ending longitude
    ! xlat1 - starting latitude
    ! xlat2 - ending latitude
    !----------------------------------------------------------------------
    ! OUTPUT:
    ! garea - grid area in m2
    !----------------------------------------------------------------------
    ! HISTORY:
    ! By: Y. F. Cheng (yafang.cheng@gmail.com)
    ! On: April 23, 2010
    !----------------------------------------------------------------------

    REAL, PARAMETER :: EARTH_RADIUS = 6370000.
    REAL, PARAMETER :: PI = 3.141592653589793
    REAL, PARAMETER :: Rad_per_Deg = PI/180.

    REAL     :: garea
    REAL     :: xlon1, xlon2
    REAL     :: ylat1, ylat2
    REAL     :: cap1_area, cap2_area
    REAL     :: cap_ht

    !******************************************************************************
    ylat1 = fixpole(ylat1)
    ylat2 = fixpole(ylat2)

    cap_ht = EARTH_RADIUS*( 1- SIN( ylat1*Rad_per_Deg ) )
    cap1_area = 2*PI*EARTH_RADIUS*cap_ht

    cap_ht = EARTH_RADIUS*( 1- SIN( ylat2*Rad_per_Deg ) )
    cap2_area = 2*PI*EARTH_RADIUS*cap_ht

    garea = ABS( cap1_area-cap2_area )*ABS(xlon1-xlon2 )/REAL(360)

    !     garea =(EARTH_RADIUS**2)*(xlon2-xlon1)*Rad_per_Deg*
    !     1(sin(ylat2*Rad_per_Deg)-sin(ylat1*Rad_per_Deg))

    IF(garea.EQ.0.)THEN
       PRINT*, 'calc_land_area: calculated area = 0'
       PRINT*,xlon2,xlon1,(xlon2-xlon1)
       PRINT*,ylat2,ylat1,(SIN(ylat2*Rad_per_Deg)-  &
            SIN(ylat1*Rad_per_Deg))

    ENDIF
    RETURN

  END SUBROUTINE calc_grid_area

  !============================================================

  SUBROUTINE calc_land_area(pct, cell_EW, cell_NS, lon, lat, area)
    !--------------------------------------------------
    !DESC: calculate the land area of a specified percentage of each
    !cell in a grid defined by longitude-latitude cell dimensions
    !(e.g. 1.25 degree by 1 degree, 5 minute by 5 minute, etc.)
    !--------------------------------------------------
    !INPUT:
    !   pct: nlat by nlon real array (in); percentages to be converted to land area
    !   cell_EW: real; E-W length of cells in degrees
    !   cell_NS: real; N-S length of cells in degrees
    !--------------------------------------------------
    !OUTPUT
    !   area: nlat by nlon real array (out); land areas calculated from pct
    !--------------------------------------------------
    !AUTHOR: Timothy W. Hilton, UC Merced, 18 Feb 2015
    !--------------------------------------------------

    USE M3UTILIO

    IMPLICIT NONE

    REAL, DIMENSION(:, :, :, :), ALLOCATABLE, INTENT(in) :: pct
    REAL, DIMENSION(:, :, :, :), ALLOCATABLE, INTENT(out) :: area
    REAL, INTENT(in) :: cell_EW, cell_NS
    REAL, DIMENSION(:), ALLOCATABLE, INTENT(in) :: lon, lat
    INTEGER i, j, nlon, nlat, ierr
    REAL ll_lon, ll_lat  !lower left corner lon, lat

    nlat = SIZE(lat)
    nlon = SIZE(lon)

    ALLOCATE(area(1, 1, nlat, nlon), STAT=ierr)

    WRITE(*,*) "nlon: ", nlon, "   nlat: ", nlat
    WRITE(*,*) "cell_EW: ", cell_EW, "   cell_NW: ", cell_NS
    PRINT *,'default kind for real is', KIND(area)
    DO j=1,nlon
       ll_lon = lon(j) - (cell_EW / 2.0)
       DO i=1,nlat
          ll_lat = lat(i) - (cell_NS / 2.0)

          CALL calc_grid_area(area(1, 1, i, j), &
               & ll_lon, ll_lon + cell_EW, &
               & ll_lat, ll_lat + cell_NS)
          ! write(*,"(9999(G12.5,:,','))"), 'area', area(1, 1, i, j), &
          !      & 'lon LL, UR', ll_lon, ll_lon + cell_EW, &
          !      & 'lat LL, UR', ll_lat, ll_lat + cell_NS, &
          !      & 'center', lat(i), lon(j), i, j
          ! the preferred MODELS-3 I/O API method of testing for
          ! missing values, according to CHANGES.html in the
          ! documentation, is: x < AMISS3.  Here set missing values to
          ! zero because mtxcple does not seem to handle them.
          IF (pct(1, 1, i, j) .LT. AMISS3) THEN
             area(1, 1, i, j) = 0.0
          ELSE
             area(1, 1, i, j) = area(1, 1, i, j) * pct(1, 1, i, j)
          ENDIF
       ENDDO
    ENDDO

  END SUBROUTINE calc_land_area


  !******************************************************************************
  REAL FUNCTION fixpole(lat_in)
    !----------------------------------------------------------------------
    ! DESC: correct latitudes larger than 90.0 or smaller than -90.0
    !----------------------------------------------------------------------
    ! INPUT:
    !     latin - a latitude value
    !----------------------------------------------------------------------
    ! OUTPUT:
    !     if lat is > 90.0 return 90.0
    !     if lat is < -90.0 return -90.0
    !----------------------------------------------------------------------
    ! HISTORY:
    ! By: Timothy W. Hilton
    ! On: 21 May 2014
    !----------------------------------------------------------------------

    REAL lat_in, n_pole, s_pole

    n_pole = 90.0
    s_pole = -90.0

    fixpole = lat_in

    IF (lat_in .GT. n_pole) THEN
       fixpole = n_pole
    ELSEIF (lat_in .LT. s_pole) THEN
       fixpole = s_pole
    ENDIF

  END FUNCTION fixpole

  !******************************************************************************
  REAL FUNCTION units_convert(gpp_in, cell_area, fillval)
    !     convert GPP from "g*m-2*d-1" (MPI raw) to kg*s-1*gridcell-1.  gpp
    !     values equal to fillval are set to zero because mtxcple appears
    !     not to handle missing values; floating point equality is respected
    !     for this comparison.
    ! author: Timothy W. Hilton

    REAL gpp_in, cell_area, fillval
    LOGICAL gpp_is_nan
    !     define some constants
    REAL KG_PER_G, DAYS_PER_SEC, EPSILON
    INTEGER HOURS_PER_DAY

    KG_PER_G = 1 / 1000.0     !kilograms per gram
    DAYS_PER_SEC = 1.0 / (60.0 * 60.0 * 24.0) !days per second
    HOURS_PER_DAY = 24
    EPSILON = 1.0D-8   !small value for floating point equality test.
    !Values within epsilon of one another are
    !considered equal.

    gpp_is_nan = (ABS(gpp_in - fillval)) .LT. EPSILON
    IF (gpp_is_nan) THEN
       units_convert = 0.0
    ELSE
       units_convert = ABS(gpp_in) * cell_area * &
            KG_PER_G * DAYS_PER_SEC
       !         write(*,'(ES10.3,  ES10.3,  ES10.3,  ES10.3,  ES10.3)'),
       !     &        abs(gpp_in), cell_area, KG_PER_G,
       !     &        SECS_PER_DAY, units_convert
    ENDIF

  END FUNCTION units_convert

  !******************************************************************************
  REAL FUNCTION calc_fcos(gpp, LRU, cos_co2_ratio, cell_per_m2, debug)
    !     calculate fCOS from GPP, leaf relative uptake, and COS/CO2
    ! author: Timothy W. Hilton
    REAL gpp, LRU, cos_co2_ratio, cell_per_m2
    REAL g_per_kg, molC_per_gC, umol_per_mol, mol_per_pmol
    LOGICAL debug

    g_per_kg = 1e3            !grams per kilogram
    molC_per_gC = 1.0 / 12.011 !moles carbon per gram carbon
    umol_per_mol = 1e6  !micromoles per mole
    mol_per_pmol = 1e-12 !moles per picomole

    calc_fcos = gpp * LRU * cos_co2_ratio * &
         g_per_kg * molC_per_gC * &
         umol_per_mol * mol_per_pmol / cell_per_m2
    calc_fcos = -1.0 * calc_fcos
    RETURN
  END FUNCTION calc_fcos

  !******************************************************************************
  SUBROUTINE handle_err(rmarker,nf_status)
    ! writes a netcdf error message to stdout.
    ! author: Sarika Kulkarni
    INCLUDE "netcdf.inc"
    INTEGER nf_status
    CHARACTER*(*)        :: rmarker
    IF (nf_status .NE. nf_noerr) THEN
       WRITE(*,*)  'NetCDF error : ',rmarker
       WRITE(*,*)  '  ',nf_strerror(nf_status)
       STOP
    ENDIF
  END SUBROUTINE handle_err

  !******************************************************************************
  SUBROUTINE per_cell_2_per_m2(vname_arg, units_arg, vdesc_arg, m2_per_cell)
    ! creates an I/O API file containing fluxes per m2 from a file
    ! containing fluxes per STEM 124x124 North American grid cell.
    !
    ! pre-conditions: input and output netcdf files must be defined in
    ! environment variables INPUT and OUTPUT
    !
    ! author: Timothy W. Hilton

    USE M3UTILIO

    IMPLICIT NONE

    INCLUDE 'netcdf.inc'

    LOGICAL debug1
    DATA debug1/.FALSE./

    REAL, ALLOCATABLE,DIMENSION (:,:,:) :: this_t_val
    INTEGER ntimes,i,j,t,year,dec31_doy
    INTEGER jdate,jtime,tstep,tstep3d_in
    INTEGER cdate, ctime
    CHARACTER*16 inspec,outspec,sector
    CHARACTER*100 input
    CHARACTER*100 invars
    CHARACTER*16 outvname
    CHARACTER*16 outvunits
    CHARACTER*80 outvdesc
    CHARACTER*80 vname_arg
    CHARACTER*80 units_arg
    CHARACTER*80 vdesc_arg
    REAL m2_per_cell, LRU, cos_co2_ratio
    INTEGER LOGDEV  !unit for I/O API log file
    INTEGER ierr    !error status for array allocation
    INTEGER LAY1

    !namelist /control/input

    LOGDEV = INIT3()

    !     ##################################################
    !     get file description of input
    IF(.NOT.OPEN3('INPUT',FSREAD3,'gpp_per_m2')) THEN !input file does not exist!
       PRINT*, 'Error opening input file'
       STOP
    ELSE
       IF (.NOT. DESC3('INPUT') ) THEN ! if exit, get information
          PRINT*, 'Error getting info from input IOAPI'
          STOP
       ELSE
          PRINT*, 'successfully got description of input IOAPI'
       ENDIF
    ENDIF

    PRINT*, 'vars: ', VNAME3D(1), UNITS3D(1), VDESC3D(1)
    VNAME3D(1) = vname_arg
    UNITS3D(1) = units_arg
    VDESC3D(1) = vdesc_arg
    !     ##################################################
    !     create output file with same file description as input but cos, not gpp
    IF(.NOT.OPEN3('OUTPUT',FSCREA3,'gpp_per_m2')) THEN
       !     output file does not exist!
       PRINT*, 'Error opening output file'
       STOP
    ENDIF

    ALLOCATE(this_t_val(NLAYS3D,NROWS3D,NCOLS3D), STAT=ierr)

    jdate = sdate3d
    jtime = stime3d
    LAY1 = 1
    PRINT*, 'cell size:', m2_per_cell
    DO t = 1, MXREC3D
       ! read gpp from input
       IF (READ3('INPUT', vname_arg, ALLAYS3, jdate, jtime, this_t_val)) THEN
          ! calculate value for every grid cell
          DO i=1, NCOLS3D
             DO j=1, NROWS3D
                this_t_val(LAY1, i, j) = &
                     this_t_val(LAY1, i, j) / m2_per_cell
             ENDDO
          ENDDO
          ! write timestemp to output
          IF (.NOT. WRITE3('OUTPUT', vname_arg, jdate, jtime, this_t_val)) THEN
             PRINT*, 'error writing', vname_arg,  'for ', jdate, jtime
             STOP
          ENDIF
          CALL nextime(jdate,jtime,TSTEP3D)
       ELSE
          PRINT*, 'error reading input'
          STOP
       ENDIF
    ENDDO

  END SUBROUTINE per_cell_2_per_m2

  !******************************************************************************
  SUBROUTINE gpp_2_fcos(LRU, cos_co2_ratio, vdesc_arg)
    ! creates an I/O API file containing COS fluxes per m2 calculated
    ! using equation 1 of Campbell et al (2008).
    !
    ! pre-conditions: input and output netcdf files must be defined in
    ! environment variables INPUT and OUTPUT
    !
    ! REFERENCES
    !
    ! Campbell, J. E., Carmichael, G. R., Chai, T., Mena-Carrasco, M.,
    ! Tang, Y., Blake, D. R., Blake, N. J., Vay, S. A., Collatz,
    ! G. J., Baker, I., Berry, J. A., Montzka, S. A., Sweeney, C.,
    ! Schnoor, J. L., and Stanier, C. O.: Photosynthetic Control of
    ! Atmospheric Carbonyl Sulfide During the Growing Season, Science,
    ! 322, 1085–1088, doi:10.1126/science.1164015, 2008.
    !
    ! author: Timothy W. Hilton

    USE M3UTILIO

    IMPLICIT NONE
    INCLUDE 'netcdf.inc'

    LOGICAL debug1
    DATA debug1/.FALSE./
    INTEGER LOGDEV  !unit for I/O API log file
    INTEGER ierr    !error status for array allocation

    REAL, ALLOCATABLE,DIMENSION (:,:,:) :: this_t_gpp
    REAL, ALLOCATABLE,DIMENSION (:,:,:) :: this_t_fcos
    INTEGER ntimes,i,j,t,year,dec31_doy
    INTEGER jdate,jtime,tstep,tstep3d_in
    INTEGER cdate, ctime
    CHARACTER*16 inspec,outspec,sector
    CHARACTER*100 input
    CHARACTER*100 invars
    CHARACTER*16 outvname
    CHARACTER*16 outvunits
    CHARACTER*80 outvdesc
    REAL cell_per_m2
    CHARACTER*80, INTENT(in):: vdesc_arg
    REAL, INTENT(in):: LRU, cos_co2_ratio

    NAMELIST /control/input

    LOGDEV = INIT3()

    !     ##################################################
    !     get file description of input
    IF(.NOT.OPEN3('INPUT',FSREAD3,'gpp_2_fcos')) THEN !output file does not exist!
       PRINT*, 'Error opening input file'
       STOP
    ELSE
       IF (.NOT. DESC3('INPUT') ) THEN ! if exit, get information
          PRINT*, 'Error getting info from input IOAPI'
          STOP
       ELSE
          PRINT*, 'successfully got description of input IOAPI'
       ENDIF
    ENDIF

    PRINT*, 'vars: ', VNAME3D(1), UNITS3D(1), VDESC3D(1)
    VNAME3D(1) = 'cos'
    UNITS3D(1) = 'mol m-2 s-1'
    VDESC3D(1) = vdesc_arg
    !     ##################################################
    !     create output file with same file description as input but cos, not gpp
    IF(.NOT.OPEN3('OUTPUT',FSCREA3,'gpp_2_fcos')) THEN
       !     output file does not exist!
       PRINT*, 'Error opening output file'
       STOP
    ENDIF

    ALLOCATE(this_t_gpp(NLAYS3D,NROWS3D,NCOLS3D), STAT=ierr)
    ALLOCATE(this_t_fcos(NLAYS3D,NROWS3D,NCOLS3D), STAT=ierr)

    jdate = sdate3d
    jtime = stime3d
    cell_per_m2 = 1.0  ! gpp is already in units of flux m-2
    PRINT*, 'LRU:', LRU, 'ratio:', cos_co2_ratio, 'cell size:', cell_per_m2
    DO t = 1, MXREC3D
       ! read gpp from input
       IF (READ3('INPUT', 'GPP', ALLAYS3, jdate, jtime, this_t_gpp)) THEN
          ! calculate fCOS for every grid cell
          DO i=1, NCOLS3D
             DO j=1, NROWS3D
                this_t_fcos(1, i, j) = calc_fcos( &
                     this_t_gpp(1, i, j), &
                     LRU, &
                     cos_co2_ratio, &
                     cell_per_m2, &
                     debug1)
             ENDDO
          ENDDO
          ! write timestemp to output
          IF (.NOT. WRITE3('OUTPUT', 'cos', jdate, jtime, this_t_fcos)) THEN
             PRINT*, 'error writing fCOS for ', jdate, jtime
             STOP
          ENDIF
          CALL nextime(jdate,jtime,TSTEP3D)
       ELSE
          PRINT*, 'error reading input'
          STOP
       ENDIF
    ENDDO

  END SUBROUTINE gpp_2_fcos

  !******************************************************************************


  SUBROUTINE open3_and_desc3(file_arg, status_arg, prog_name_arg)
    !----------------------------------------------------------------------
    ! DESC: call OPEN3 and DESC3 on a specified I/O API file, with
    !       error handling.  If either OPEN3 or DESC3 fails writes an
    !       error message and issues a stop command.
    !----------------------------------------------------------------------
    ! INPUT:
    !     file_arg: string; logical name of file to be described (see
    !        I/O API documentation for logical name discussion.  Must
    !        be 16 characters or less.
    !     status_arg: integer; one of the file status indicators
    !        defined in I/O API's PARMS3.EXT
    !     prog_name_arg: string; the name of the program calling the
    !        subroutine (to be recorded by I/O API).  Must be 16
    !        characters or less.
    !----------------------------------------------------------------------
    ! OUTPUT:
    !     No explicit output; has the side effect of populating the
    !     common block variables defined in I/O API's FDESC3.EXT
    !----------------------------------------------------------------------
    ! HISTORY:
    ! By: Timothy W. Hilton
    ! On: 14 Jan 2015
    !----------------------------------------------------------------------

    USE M3UTILIO
    IMPLICIT NONE

    CHARACTER(len=*), INTENT(in) :: file_arg, prog_name_arg
    INTEGER, INTENT(in) :: status_arg

    IF(.NOT.OPEN3(file_arg, status_arg, prog_name_arg)) THEN
       !output file does not exist!
       PRINT*, 'Error opening input file'
       STOP
    ELSE
       IF (.NOT. DESC3(file_arg) ) THEN ! if exit, get information
          PRINT*, 'Error getting info from input IOAPI'
          STOP
       ELSE
          PRINT*, 'successfully got description of input IOAPI'
       ENDIF
    ENDIF

    !print*, 'in open3_and_desc3, vars: ', VNAME3D(1), UNITS3D(1), VDESC3D(1)

  END SUBROUTINE open3_and_desc3


  ! SUBROUTINE get_arctas_NA_grid_parmeters(gdname)

  !   !----------------------------------------------------------------------
  !   ! DESC: define the I/O API grid parmameters describing the 124x124x22
  !   !     STEM grid and place them in the FDESC3 commons structures.
  !   !----------------------------------------------------------------------
  !   ! INPUT:
  !   !     gdname: character string; the name of the grid in the GRIDDESC file
  !   !----------------------------------------------------------------------
  !   ! OUTPUT:
  !   !     no explicit outputs.
  !   !----------------------------------------------------------------------
  !   ! PRECONDITIONS: INCLUDE 'FDESC3.EXT', include 'PARMS3.EXT', file
  !   ! with logical name 'GRIDDESC' must exist and be a valid I/O API
  !   ! grid description file.
  !   !----------------------------------------------------------------------
  !   ! SEE ALSO:
  !   !    See I/O API documentation for DESC3 function and for FDESC3.EXT.
  !   !----------------------------------------------------------------------
  !   ! HISTORY:
  !   ! By: Timothy W. Hilton
  !   ! On: 19 Jan 2015
  !   !----------------------------------------------------------------------

  !   IMPLICIT NONE

  !   include 'PARMS3.EXT'	! I/O API
  !   include 'FDESC3.EXT'	! I/O API
  !   include 'IODECL3.EXT'	! I/O API

  !   character(len=80) :: coord_sys_name
  !   character(len=16), intent(in) :: gdname
  !   integer :: coord_sys_type

  !   ! read grid parameters from GRIDDESC
  !   CALL DSCGRID( 'ARCNAGRID', GDNAM3D, &
  !        GDTYP3D, P_ALP3D, P_BET3D, P_GAM3D, XCENT3D, YCENT3D, &
  !        XORIG3D, YORIG3D, XCELL3D, YCELL3D, NCOLS3D, NROWS3D, NTHIK3D )

  !   FTYPE3D = GRDDED3

  !   !define vertical coordinates
  !   VGTYP3D = VGSGPN3   ! non-hydrostatic sigma-P coords
  !   VGTOP3D = 1.0
  !   !VGLVS3D = (/1.0, 0.0/)  !TODO: figure this one out...
  !   ! GDNAM3D = "ARCNAGRID"

  ! end SUBROUTINE get_arctas_NA_grid_parmeters


END MODULE ioapi_regrid_tools
