
program wrf_to_ioapi
   use netcdf
   implicit none
   include 'netcdf.inc'
   include 'PARMS3.EXT'	! i/o API
   include 'FDESC3.EXT'	! i/o API
   include 'IODECL3.EXT'	! i/o API

   ! variables list
   CHARACTER(len=255) :: infile
   CHARACTER(len=255) :: outfile 
   integer::LOGDEV
   integer::jdate,jtime !integers used for TFLAG during write
   integer::ncid, xlatid, xlonid, hgtid, MFid, LUid !netcdf variable id's
   integer::nrows, ncols!data dimensions
   integer::x,y
   real, DIMENSION(:,:), ALLOCATABLE :: xlat
   real, DIMENSION(:,:), ALLOCATABLE :: xlon
   real, DIMENSION(:,:), ALLOCATABLE :: hgt 
   real, DIMENSION(:,:), ALLOCATABLE :: MF 
   INTEGER, DIMENSION(:,:), ALLOCATABLE :: LU 

   print*, '---------------checking invironemnt:   '
   CALL get_environment_variable("WRF_FILE", infile)
   CALL get_environment_variable("OUTPUT",outfile)
   WRITE (*,*) TRIM(infile)
   WRITE (*,*) TRIM(outfile)
   print*,"infile###############",infile


   ! allocate the input data variable
   !nrows = 310
   !ncols = 202
   nrows = 50 
   ncols = 50 
!   zlevs = 59
   allocate(xlat(ncols,nrows))
   allocate(xlon(ncols,nrows))
   allocate(hgt(ncols,nrows))
   allocate(MF(ncols,nrows))
   allocate(LU(ncols,nrows))


   ! Open the net cdf file. NF90_NOWRITE Just means read only.
   print*,"Opening file..."
   call check( nf90_open(infile, NF90_NOWRITE, ncid))
   
   ! get the varid of the data
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
   ! read the data

   print*, "Data imported succsessfully!"
   
   ! Check the data.
!   print*, "Checking the data:"
!   do x = 1, ncols
!      do y = 1, nrows
!            print *, "lat:", xlat(x,y), ", lon:", xlon(x,y),", hgt:", hgt(y, x)
!      
!      end do
!   end do

!   print*, "lower left", "lat",xlat(209,1), "lon", xlon(209,1)
!   print*, "lower left", "lat",xlat(1,1), "lon", xlon(1,1)
   ! Close the netcdf file
   call check( nf90_close(ncid))

   ! START I/O API CREATION OF TOPO
   LOGDEV = INIT3()


      ! setup the variables used to write the I/O api header
      nvars3d=5                 ! emission species number
      ftype3d=GRDDED3           ! file is in grided, Global dobson file header
      gdtyp3d= 2!Lambert                !
!      p_alp3d=90.              !unused in lat-lon
!      p_bet3d=90.              !unused in lat-lon
      xcent3d=37.4136                !unused in lat-lon
      ycent3d=-121.886                !unused in lat-lon
      xorig3d=-122.0788
      yorig3d=30.65546
      xcell3d=9000.
      ycell3d=9000.
      ncols3d=ncols
      nrows3d=nrows
      nlays3d=1                 ! documentation is vague on this
      vgtyp3d=VGSGPN3           !  non-hydrostatic sigma-p vertical coordinate
      vgtop3d=1.                ! domain top in meter
      nthik3d=1
      vglvs3d(1)=1.             ! levels in meter
      units3d(1)='Degree'
      vname3d(1) = 'LAT'
      vdesc3d(1) = 'Latitude'
      vtype3d(1)=m3real

      vglvs3d(2)=1.             ! levels in meter
      units3d(2)='Degree'
      vname3d(2) = 'LON'
      vdesc3d(2) = 'Longitude'
      vtype3d(2)=m3real

      vglvs3d(3)=1.             ! levels in meter
      units3d(3)='meter'
      vname3d(3) = 'TOPO'
      vdesc3d(3) = 'Topograph'
      vtype3d(3)=m3real
      
      vglvs3d(4)=1.             ! levels in meter
      units3d(4)='RATIO'
      vname3d(4) = 'MFACTOR'
      vdesc3d(4) = 'MASS MAP FACTOR'
      vtype3d(4)=m3real
      
      vglvs3d(5) = 1.             ! levels in meter
      units3d(5) = 'Integer CATEGORY'
      vname3d(5) = 'LANDUSE'
      vdesc3d(5) = 'LANDUSE CATEGORIES'
      vtype3d(5)=m3int
      !sdate3d = 2010185
      !stime3d =  120000  
      !tstep3d = 003000
      !sdate3d = 1999365
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
