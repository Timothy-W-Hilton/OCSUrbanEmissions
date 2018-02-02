!--------------------------------------------------
!     REQUIREMENTS
!     EDSS/Models-3 I/O API logical names WRF_FILE, OUTPUT, and GRIDDESC must
!     be set.  WRF_FILE should point to the WRF output file to be converted
!     to a TOPO file.  OUTPUT should point to the TOPO file to create.
!--------------------------------------------------
      program wrf_to_ioapi
      USE M3UTILIO
      USE NETCDF
      IMPLICIT NONE


! variables list
      CHARACTER(len=255) :: infile
      CHARACTER(len=255) :: outfile
      CHARACTER(len=255) :: griddesc
      integer::LOGDEV
      integer::jdate,jtime      !integers used for TFLAG during write
      integer::ncid, xlatid, xlonid, hgtid, MFid, LUid !netcdf variable id's
      integer::nrows, ncols     !data dimensions
      integer::x,y
      real, DIMENSION(:,:), ALLOCATABLE :: xlat
      real, DIMENSION(:,:), ALLOCATABLE :: xlon
      real, DIMENSION(:,:), ALLOCATABLE :: hgt
      real, DIMENSION(:,:), ALLOCATABLE :: MF
      INTEGER, DIMENSION(:,:), ALLOCATABLE :: LU
      logical log_stat

      print*, '---------------checking environment:   '
      CALL get_environment_variable("WRF_FILE", infile)
      CALL get_environment_variable("OUTPUT",outfile)
      CALL get_environment_variable("GRIDDESC",griddesc)
      WRITE (*,*) TRIM(infile)
      WRITE (*,*) TRIM(outfile)
      WRITE (*,*) TRIM(griddesc)
      print*,"infile###############",infile

      log_stat = DSCGRID( 'CRO_BCND04', GDNAM3D, GDTYP3D, P_ALP3D, &
     &     P_BET3D, P_GAM3D, XCENT3D, YCENT3D, XORIG3D, YORIG3D, &
     &     XCELL3D , YCELL3D, NCOLS3D, NROWS3D, NTHIK3D)

      print *, 'nrows3d, ncols3d', NROWS3D, NCOLS3D

      allocate(xlat(NCOLS3D, NROWS3D))
      allocate(xlon(NCOLS3D, NROWS3D))
      allocate(hgt(NCOLS3D, NROWS3D))
      allocate(MF(NCOLS3D, NROWS3D))
      allocate(LU(NCOLS3D, NROWS3D))

! Open the net cdf file. NF90_NOWRITE Just means read only.
      print*,"Opening file..."
      call check( nf90_open(infile, NF90_NOWRITE, ncid))

! get the varid of the data and read the data
      print*,"Reading data variable XLAT"
      call check( nf90_inq_varid(ncid, "XLAT", xlatid))
      call check( nf90_get_var(ncid,xlatid,xlat))

      print*,"Reading data variable XLONG"
      call check( nf90_inq_varid(ncid, "XLONG", xlonid))
      call check( nf90_get_var(ncid,xlonid,xlon))

      print*,"Reading data variable HGT"
      call check( nf90_inq_varid(ncid, "HGT", hgtid))
      call check( nf90_get_var(ncid,hgtid,hgt))

      print*,"Reading data variable MAPFAC_M"
      call check( nf90_inq_varid(ncid, "MAPFAC_M", MFid))
      call check( nf90_get_var(ncid,MFid,MF))

      print*,"Reading data variable LU_INDEX"
      call check( nf90_inq_varid(ncid, "LU_INDEX", LUid))
      call check( nf90_get_var(ncid,LUid,LU))

      print*, "Data imported succsessfully!"

! Check the data.
!     print*, "Checking the data:"
!     do x = 1, ncols
!     do y = 1, nrows
!     print *, "lat:", xlat(x,y), ", lon:", xlon(x,y),", hgt:", hgt(y, x)
!
!     end do
!     end do

!     print*, "lower left", "lat",xlat(209,1), "lon", xlon(209,1)
!     print*, "lower left", "lat",xlat(1,1), "lon", xlon(1,1)
! Close the netcdf file
      call check( nf90_close(ncid))

! START I/O API CREATION OF TOPO
      LOGDEV = INIT3()
!     setup the variables used to write the I/O api header.  Many of
!     these have been read from the GRIDDESC file; setup a few remaining
!     ones.
      NVARS3D=5                 ! emission species number
      ftype3d=GRDDED3           ! file is in grided, Global dobson file header
      NLAYS3D=1                 ! documentation is vague on this
!     coordinate I/O API documentation section "vertical grids"
!     (https://www.cmascenter.org/ioapi/documentation/all_versions/html/GRIDS.html#vert)
!     suggests setting VGTYP3D to IMISS3 and VGTOP3D to BADVAL3 and when
!     there is only one layer
      VGTYP3D=IMISS3           !  single-layer data
      VGTOP3D=BADVAL3          ! domain top not applicable for single layer

      VGLVS3D(1)=1. ! vertical levels
      units3d(1)='Degree'
      vname3d(1) = 'LAT'
      vdesc3d(1) = 'Latitude'
      vtype3d(1)=m3real

      vglvs3d(2)=1.             ! vertical levels
      units3d(2)='Degree'
      vname3d(2) = 'LON'
      vdesc3d(2) = 'Longitude'
      vtype3d(2)=m3real

      vglvs3d(3)=1.             ! vertical levels
      units3d(3)='meter'
      vname3d(3) = 'TOPO'
      vdesc3d(3) = 'Topograph'
      vtype3d(3)=m3real

      vglvs3d(4)=1.             ! vertical levels
      units3d(4)='RATIO'
      vname3d(4) = 'MFACTOR'
      vdesc3d(4) = 'MASS MAP FACTOR'
      vtype3d(4)=m3real

      vglvs3d(5) = 1.           ! vertical levels
      units3d(5) = 'Integer CATEGORY'
      vname3d(5) = 'LANDUSE'
      vdesc3d(5) = 'LANDUSE CATEGORIES'
      vtype3d(5)=m3int

      sdate3d = 2015064
      stime3d = 0
      tstep3d = 0

! attempt to open the I/O api file. Create one if it does not exist
      print*, 'Attempting to open file...'
      if(.not.OPEN3('OUTPUT',FSRDWR3,'wrf_to_ioapi')) then ! output file does not exit
         print*, 'File does not exist, attempting file creation...'
         if(.not.OPEN3('OUTPUT',FSCREA3,'wrf_to_ioapi')) then ! FSCREA3 FSUNKN3
            print*, 'Error opening output file'
            stop
         endif
      else
         print*,'Reading from previously created file!'
         if (.not. DESC3('OUTPUT') ) then ! if exit, get information
            print*, 'Error getting info from OUTPUT nc'
            stop
         endif
      endif


! Setup time and iteration information.
      jdate = 0
      jtime = 0


! Attempt an I/O api data write. For whatever reason the convention is to
!     loop through time and write one 2d grid at a time.
      if(.not.write3('OUTPUT',vname3d(1),jdate,jtime,xlat)) then
         print*,'writing error'
         stop
      endif

      if(.not.write3('OUTPUT',vname3d(2),jdate,jtime,xlon)) then
         print*,'writing error'
         stop
      endif

      if(.not.write3('OUTPUT',vname3d(3),jdate,jtime,hgt)) then
         print*,'writing error'
         stop
      endif

      if(.not.write3('OUTPUT',vname3d(4),jdate,jtime,MF)) then
         print*,'writing error'
         stop
      endif

      if(.not.write3('OUTPUT',vname3d(5),jdate,jtime,LU)) then
         print*,'writing error'
         stop
      endif
! Close I/O api
      print*, 'SHUT3()=',SHUT3()


      contains
      subroutine check(status)
      integer, intent ( in) :: status

      if(status /= nf90_noerr) then
         print *, trim(nf90_strerror(status))
         stop "Stopped"
      end if
      end subroutine check
      end program wrf_to_ioapi
