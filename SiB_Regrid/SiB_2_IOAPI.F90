MODULE HELPER_ROUTINES

  USE M3UTILIO
  USE ioapi_regrid_tools
  IMPLICIT NONE

  INCLUDE 'netcdf.inc'

  TYPE SiB_real_variable
     CHARACTER(LEN=16) :: name
     CHARACTER(LEN=80) :: units
     CHARACTER(LEN=80) :: desc
     REAL, DIMENSION(:, :, :, :), ALLOCATABLE :: values
     INTEGER :: m3type
  END TYPE SiB_real_variable

CONTAINS

  SUBROUTINE get_ncdf_dimlen(ncid, dimname, dimlen)
    !--------------------------------------------------
    !DESC: get the size of a netcdf dimension given a netcdf ID
    !      and the dimension name
    !--------------------------------------------------
    !INPUT:
    !ncid: the netCDF ID
    !dimname: string; the dimension requested
    !--------------------------------------------------
    !OUTPUT
    !dimlen: integer; the length of the requested dimension
    !--------------------------------------------------
    !AUTHOR: Timothy W. Hilton, UC Merced, 1 Apr 2015
    !--------------------------------------------------
    USE ioapi_regrid_tools  ! for handle_err
    IMPLICIT NONE

    INTEGER, INTENT(in) :: ncid
    CHARACTER(len=*), INTENT(in) :: dimname
    CHARACTER(len=160) :: msg
    INTEGER, INTENT(out) :: dimlen
    INTEGER :: dimid, nf_status

    nf_status = nf_inq_dimid (ncid, dimname, dimid)
    CALL handle_err(dimname, nf_status)
    nf_status = nf_inq_dimlen(ncid, dimid, dimlen)
    CALL handle_err('dimension length', nf_status)

  END SUBROUTINE get_ncdf_dimlen

  SUBROUTINE read_SiB_latlonidx(ncid, lonidx, latidx)
    !--------------------------------------------------
    !DESC: read values the latindex and lonindex variable of a SiB
    !netcdf file
    !--------------------------------------------------
    !INPUT:
    !ncid: the netCDF ID
    !--------------------------------------------------
    !OUTPUT
    !latidx, lonidx: integer arrays; the values of lon_index and lat_index
    !--------------------------------------------------
    !AUTHOR: Timothy W. Hilton, UC Merced, 1 Apr 2015
    !--------------------------------------------------
    USE ioapi_regrid_tools  ! for handle_err
    IMPLICIT NONE

    INTEGER, INTENT(in) :: ncid
    INTEGER, ALLOCATABLE, DIMENSION(:), INTENT(out) :: lonidx, latidx
    INTEGER :: ierr, vid, nf_status, nvals

    CALL get_ncdf_dimlen(ncid, 'land_points', nvals)

    ALLOCATE(lonidx(nvals), STAT=ierr)
    ALLOCATE(latidx(nvals), STAT=ierr)

    nf_status = nf_inq_varid(ncid, 'lonindex', vid)
    nf_status = NF_GET_VAR_INT(ncid, vid, lonidx)
    CALL handle_err('read lon index', nf_status)

    nf_status = nf_inq_varid(ncid, 'latindex', vid)
    nf_status = NF_GET_VAR_INT(ncid, vid, latidx)
    CALL handle_err('read lat index', nf_status)

    ! print *, 'latidx:', latidx
    ! print *, 'lonidx:', lonidx


  END SUBROUTINE read_SiB_latlonidx

  SUBROUTINE read_SiB_var2D(ncid, varname, ntime, nlon, nlat, &
       & lonidx, latidx, sib_var)
    !--------------------------------------------------
    !DESC: read surface flux values from a netcdf file of SiB output into an 4-D
    !   allocatable array.  Reads values for CO2 GPP, CO2 RE, COS GPP,
    !   COS soil flux, and COS flux.  The dimensions of the array are
    !   time, layer, latitude, and longitude.  Layer is set to 1
    !   because these are all surface fluxes; having the layer
    !   dimension is useful to write out the values to a Models-3 I/O
    !   API file.
    !--------------------------------------------------
    !INPUT:
    !ncid; integer: the netCDF ID of the file to be read
    !varname; string: the variable requested
    !ntime; integer: the length of the time dimension
    !nlat; integer: the length of the latitude dimension
    !nlon; integer: the length of the longitude dimension
    !lonidx; integer array: the lon_index of the SiB file
    !latidx; integer array: the lat_index of the SiB file
    !--------------------------------------------------
    !OUTPUT
    !arr4D: four dimensional array containing the values from the
    !   requested variable
    !--------------------------------------------------
    !AUTHOR: Timothy W. Hilton, UC Merced, 1 Apr 2015
    !--------------------------------------------------
    USE m3utilio
    USE ioapi_regrid_tools  ! for handle_err
    USE ARRAY_FUNCS
    IMPLICIT NONE

    INTEGER, INTENT(in) :: ncid, ntime, nlon, nlat
    INTEGER, DIMENSION(:) :: lonidx, latidx
    INTEGER :: vid, nlayer, nvals, ierr, nf_status, i, t
    CHARACTER(len=*), INTENT(in) :: varname
    CHARACTER(len=160) :: msg
    REAL, DIMENSION(:, :), ALLOCATABLE :: arr2D
    TYPE(SiB_real_variable), INTENT(out) :: sib_var

    nlayer = 1
    nvals = SIZE(lonidx)

    ALLOCATE(arr2D(nvals, ntime), STAT=ierr)
    nf_status = nf_inq_varid(ncid, varname, vid)
    nf_status = NF_GET_VAR_REAL(ncid, vid, arr2D)
    CALL handle_err('read SiB var', nf_status)

    ALLOCATE(sib_var%values(ntime, 1, nlat, nlon), STAT=ierr)
    sib_var%values(:, :, :, :) = 0.0
    DO t = 1, ntime
       DO i = 1, SIZE(lonidx)
          sib_var%values(t, 1, latidx(i), lonidx(i)) = arr2D(i, t)
       END DO
    END DO

    DEALLOCATE(arr2d)
    sib_var%desc=''   ! gfortran fills with junk if not initialized
    sib_var%name = varname
    nf_status = NF_get_att_text(ncid, vid, 'units', sib_var%units)
    nf_status = NF_get_att_text(ncid, vid, 'long_name', sib_var%desc)
    sib_var%m3type = M3REAL  ! Models-3 I/O API type for real
                               ! variable, defined in PARMS3.EXT
  END SUBROUTINE read_SiB_var2D

  SUBROUTINE read_SiB_file(fname, GPP_CO2, RE_CO2, GPP_COS, soil_COS, &
       & flux_COS, lon, lat, debug)
    !--------------------------------------------------
    !DESC: read values from a netcdf file of SiB output into an 2-D
    !   allocatable array.  Reads values for CO2 GPP, CO2 RE, COS GPP,
    !   COS soil flux, and COS flux
    !--------------------------------------------------
    !INPUT:
    !fname: the name of the file to be read
    !debug: logical; if true, print debugging output
    !--------------------------------------------------
    !OUTPUT
    !GPP_CO2: array containing the CO2 GPP
    !RE_CO2: real; array containing the CO2 respiration
    !GPP_COS: real; array containing the COS plant flux
    !soil_COS: real: array containing the COS soil flux
    !flux_COS: real; array containing the COS flux
    !mask: integer; array containing 1 where a missing value is
    !   present in the input
    !--------------------------------------------------
    !AUTHOR: Timothy W. Hilton, UC Merced, 1 Apr 2015
    !--------------------------------------------------
    USE ioapi_regrid_tools
    USE ARRAY_FUNCS
    IMPLICIT NONE

    LOGICAL, INTENT(in) :: debug
    TYPE(SiB_real_variable), INTENT(out) :: &
         & GPP_CO2, RE_CO2, GPP_COS, soil_COS, flux_COS
    REAL, DIMENSION(:), ALLOCATABLE, INTENT(out) :: lon, lat
    INTEGER, DIMENSION(:), ALLOCATABLE :: lonidx, latidx
    REAL :: NC_BADVAL, epsilon
    INTEGER i, j, ierr, nlon, nlat, ntime, nland_points
    INTEGER ncid, dimid, nf_status, vid_farea, vid
    CHARACTER(len=1024)::format_str, discard
    CHARACTER(len=512), INTENT(in)::fname

    ! open netcdf file and read dimensions
    nf_status = nf_open(TRIM(fname), nf_nowrite, ncid)
    CALL handle_err('Error opening file', nf_status)

    ! read dimensions
    CALL get_ncdf_dimlen(ncid, 'time', ntime)
    CALL get_ncdf_dimlen(ncid, 'longitude', nlon)
    CALL get_ncdf_dimlen(ncid, 'latitude', nlat)
    CALL get_ncdf_dimlen(ncid, 'land_points', nland_points)
    PRINT*, 'nlat: ', nlat, 'nlon: ', nlon, 'ntime: ', ntime, &
         & 'nland_points: ', nland_points

    ! read longitude
    ALLOCATE(lon(nlon), stat=ierr)
    nf_status = nf_inq_varid(ncid, 'longitude', vid)
    nf_status = NF_GET_VAR_REAL(ncid, vid, lon)
    CALL handle_err('read longitude', nf_status)

    ! read latitude
    ALLOCATE(lat(nlat), stat=ierr)
    nf_status = nf_inq_varid(ncid, 'latitude', vid)
    nf_status = NF_GET_VAR_REAL(ncid, vid, lat)
    CALL handle_err('read latitude', nf_status)

    CALL read_SiB_latlonidx(ncid, lonidx, latidx)

    ! read in the fluxes
    CALL read_SiB_var2D(ncid, 'GPP', ntime, nlon, nlat, &
         & lonidx, latidx, GPP_CO2)
    CALL read_SiB_var2D(ncid, 'Rtotal', ntime, nlon, nlat, &
         & lonidx, latidx, RE_CO2)
    CALL read_SiB_var2D(ncid, 'OCS_gpp', ntime, nlon, nlat, &
         & lonidx, latidx, GPP_COS)
    CALL read_SiB_var2D(ncid, 'OCS_soil', ntime, nlon, nlat, &
         & lonidx, latidx, soil_COS)
    CALL read_SiB_var2D(ncid, 'OCS_flux', ntime, nlon, nlat, &
         & lonidx, latidx, flux_COS)

    nf_status = nf_close(ncid)

    DEALLOCATE(lonidx)
    DEALLOCATE(latidx)

  END SUBROUTINE read_SiB_file

!============================================================

  SUBROUTINE write_sib_vars_to_ioapi(variables, lon, lat, year, month)
    !--------------------------------------------------
    !DESC: write SiB variables to Models-3 I/O API file with logical
    !   name OUTPUT.
    !--------------------------------------------------
    !INPUT:
    !   variables: 5-element array of type SiB_real_variable
    !   lon, lat: real arrays; the latitude and longitude values of
    !      SiB grid points
    !   year: the year of the SiB data
    !   month: the month of the SiB data
    !--------------------------------------------------
    !OUTPUT
    !
    !--------------------------------------------------
    !AUTHOR: Timothy W. Hilton, UC Merced, 3 Apr 2015
    !--------------------------------------------------
    USE M3UTILIO
    USE ioapi_regrid_tools   ! for calc_grid_area
    USE ARRAY_FUNCS
    IMPLICIT NONE

    TYPE(SiB_real_variable) :: variables(5)
    REAL, ALLOCATABLE, DIMENSION(:), INTENT(in) :: lon, lat
    INTEGER, INTENT(in) :: year, month
    REAL, ALLOCATABLE, DIMENSION(:, :, :, :) :: area, pct
    REAL, ALLOCATABLE, DIMENSION(:, :, :) :: tmp3d
    REAL :: cell_EW, cell_NS
    LOGICAL status
    INTEGER ierr, dimid, nlon, nlat, ntimes
    INTEGER i, j, LOGDEV, this_var, this_t
    INTEGER jdate, jtime, this_date, this_time

    nlon = SIZE(variables(1)%values, 4)
    nlat = SIZE(variables(1)%values, 3)

    !--------------------------------------------------
    ! convert fluxes from umol m-2 s-1 to umol gridcell-1 s-1
    !--------------------------------------------------
    cell_EW = 1.25  ! cell E-W size in degrees longitude
    cell_NS = 1.0   ! cell N-S size in degrees latitude
    ALLOCATE(area(1, 1, nlat, nlon), stat=ierr)
    ALLOCATE(pct(1, 1, nlat, nlon), stat=ierr)
    pct(:, :, :, :) = 1.0  ! no scaling here
    CALL calc_land_area(pct, cell_EW, cell_NS, lon, lat, area)
    ntimes = SIZE(variables(1)%values, 1)
    PRINT *,'shape of array before: ', SHAPE(variables(1)%values)
    DO this_var = 1, SIZE(variables)
       do i = 1, ntimes
          variables(this_var)%values(i, :, :, :) = &
               & variables(this_var)%values(i, :, :, :) * area(1, :, :, :)
       end do
       variables(this_var)%units = 'umol grid-1 s-1'

    END DO
    PRINT *, 'shape of array after: ', SHAPE(variables(1)%values)

    ! initialize the I/O API data structures
    LOGDEV = INIT3()

    ! read grid parameters from GRIDDESC
    status = DSCGRID( 'SiB_grid', GDNAM3D, &
         & GDTYP3D, P_ALP3D, P_BET3D, P_GAM3D, XCENT3D, YCENT3D, &
         & XORIG3D, YORIG3D, XCELL3D, YCELL3D, NCOLS3D, NROWS3D, NTHIK3D )

    ! populate dimensions, date & time, file & variable descriptions, etc.
    FDESC3D(1) = "SiB global output received from Ian Baker 13 Mar 2015"
    nvars3d=SIZE(variables)   ! how many variables in file?
    ftype3d=GRDDED3           ! file is in gridded
    ncols3d=nlon
    nrows3d=nlat
    nlays3d=1
    vgtyp3d=VGZVAL3           ! vertical levels are height above ground, m
    vgtop3d=1.                ! domain top in meter
    vglvs3d(1)=0.             ! levels in meter

    SDATE3D = 1000 * year + JULIAN(year, month, 1)
    STIME3D = 0
    TSTEP3D = 10000   ! one hour expressed as HMMSS

    DO i=1, SIZE(variables)
       VNAME3D(i) = TRIM(variables(i)%name)
       UNITS3D(i) = TRIM(variables(i)%units)
       VDESC3D(i) = TRIM(variables(i)%desc)
       VTYPE3D(i) = variables(i)%m3type
    END DO

    ! create new I/O API file
    IF(.NOT.OPEN3('OUTPUT',FSNEW3,'SiB_2_IOAPI')) THEN ! output file does not exit
       PRINT*, 'Error opening output file'
       STOP
    ENDIF

    !--------------------------------------------------
    ! write the flux data to the I/O API file
    !--------------------------------------------------
    PRINT*, 'writing SiB to I/O API'

    ALLOCATE(tmp3d(nlon, nlat, 1), STAT=ierr)
    DO this_var = 1, SIZE(variables)
       this_date = SDATE3D
       this_time = STIME3D

       tmp3d(:, :, :) = 0.0
       DO this_t = 1, ntimes
          DO i = 1, nlon
             DO j = 1, nlat
                !get the data dimensions into the order that WRITE3
                !expects (but not, oddly, the order that they appear
                !in the in the I/O API netcdf file.  The only place I
                !can find that this behavior is "documented" is in the
                !examples section of the WRITE3 documentation.
                tmp3d(i, j, 1) = variables(this_var)%values(this_t, 1, j, i)
             ENDDO
          ENDDO

          IF (.NOT.write3('OUTPUT', trim(variables(this_var)%name), &
               & this_date, this_time, tmp3d)) STOP
          CALL nextime(this_date, this_time, TSTEP3D)
       END DO
    END DO
    DEALLOCATE(tmp3d)

    PRINT *, '-----------------------------------------------'

    RETURN

  END SUBROUTINE write_sib_vars_to_ioapi

END MODULE HELPER_ROUTINES

!============================================================

PROGRAM SiB_to_IOAPI

  USE M3UTILIO
  USE HELPER_ROUTINES
  USE ARRAY_FUNCS
  IMPLICIT NONE

  LOGICAL debug_flag
  DATA debug_flag/.FALSE./

  TYPE(SiB_real_variable) :: GPP_CO2, RE_CO2, GPP_COS, &
       soil_COS, flux_COS
  TYPE(SiB_real_variable) all_vars(5)
  REAL, DIMENSION(:), ALLOCATABLE :: lon, lat
  CHARACTER*512 :: fname, buf
  INTEGER nlon, nlat, success, ierr, year, month

  CALL GETARG(1, buf)
  READ(buf, *) year

  CALL GETARG(2, buf)
  READ(buf, *) month

  print *, year, month

  ! fname = '/home/ecampbell_lab/COS/Additional_Flux_Models/SiB_From_Ian_2015-03/flux_hourly_200808p001.nc'
      CALL ENVSTR('RAW_SIB_FILE', 'infile: the SiB file to process', &
     & './default.nc', fname, ierr)
  IF (ierr .NE. 0) THEN
     print *, 'unable to obtain RAW_SIB_FILE physical file'
     STOP
  ENDIF
  print *, 'processing RAW_SIB_FILE: ', fname


  CALL read_SiB_file(fname, GPP_CO2, RE_CO2, GPP_COS, soil_COS, &
       & flux_COS, lon, lat, debug_flag)

  all_vars(1) = GPP_CO2
  all_vars(2) = RE_CO2
  all_vars(3) = GPP_COS
  all_vars(4) = soil_COS
  all_vars(5) = flux_COS

  PRINT *, 'GPP shape:', shape(GPP_CO2%values)
  PRINT *, 'writing I/O API file now'
  CALL write_sib_vars_to_ioapi(all_vars, lon, lat, year, month)

  ! indicates successful completion (Models-3 I/O API docs for
  ! m3exit)
  success = 0
  CALL M3EXIT('SiB_to_IOAPI', sdate3d, stime3d, 'program completed', success)

  DEALLOCATE(GPP_CO2%values, stat=ierr)
  DEALLOCATE(RE_CO2%values, stat=ierr)
  DEALLOCATE(GPP_COS%values, stat=ierr)
  DEALLOCATE(soil_COS%values, stat=ierr)
  DEALLOCATE(flux_COS%values, stat=ierr)

END PROGRAM SiB_to_IOAPI
