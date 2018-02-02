      program wrf2nc
c-------------------------------------------------------------------------------------
c convert WRF output to netCDF format, and perform interpolation
C if necessary
C
C  use Ustar and PBL from WRF directly except for the first timestep
C-------------------------------------------------------------------------------------
      USE M3UTILIO
      IMPLICIT NONE

      include 'netcdf.inc'

      logical iflag_log, nflag
      integer imax, jmax, kmax, invar2d, invar3d, logunit,
     &     imaxout, jmaxout,
     &     kmaxout, ndspecies, nusgs,
     &     imaxstag, jmaxstag, kmaxstag, LDUC, ioff, joff, koff,
     &     nratio, dimid, i, j, k, id_ntime, id_times, id_znu, id_znw,
     &     idate, ierr, ihour, ihour1, ii, iminute, iminute1,
     &     imonth, isecond, iyear, jj, mtime, n, ncid, nf_status,
     &     nland, nowdate, nowtime, nowtstring , nowyear , nseason ,
     &     ntimewrf , ntmp , numts , nx , ny , nz , rdry , sumz ,
     &     ust1 , ust2, ksoil, L, iflag, in_rain1, in_rain2,
     &     in_rwater, in_smois, in_swdown, in_t, in_tg,
     &     in_tsoil, in_u10, in_ustar, in_v10, in_w, in_glw, in_pl1,
     &     in_p1, in_p2, in_pbl, in_q, indep_div, in_cldfra,
     &     in_cwater, kdiv, ktop, lb, lkh, lkv, lp, lrh, lt,
     &     lu, lv, lw
      real temk, frac, g, esat1, esat, qsat, pcb

      parameter(imax=50,jmax=50,imaxstag=imax+1,jmaxstag=jmax+1,
     1 kmax=59,kmaxstag=kmax+1,ioff=0,
     1 joff=0,nratio=1,imaxout=50,jmaxout=50,kmaxout=59,
     2 ndspecies=16)

      parameter(invar3d=8,invar2d=12,LDUC=11,nusgs=24)         ! RADM landuse category
      character  bhic(50,20)*80,bhrc(20,20)*80,wrffile*180, name*9,
     1  start_date*24,current_date*24,staggering*4,ordering*4,units*25,
     2  description*46

      integer  begyear,begdate,begtime,dtstep,mapusgs(nusgs,2),
     1 start_index(4), end_index(4), bhi(50,20), modate(300),
     1 wrftstep,hour,d_off
      real s,xlat(imaxout,jmaxout),xlon(imaxout,jmaxout),
     a topo(imaxout,jmaxout),XLUC(0:LDUC),
     1 xland(imaxout,jmaxout),zfull(imaxout,jmaxout,kmax+1),
     2 z(imaxout,jmaxout,kmax),t(imaxout,jmaxout,kmax),
     3 u(imaxout,jmaxout,kmax),v(imaxout,jmaxout,kmax),
     4 p(imaxout,jmaxout,kmax),rh(imaxout,jmaxout,kmax),
     5 qv(imaxout,jmaxout,kmax),pfull(imaxout,jmaxout,kmax+1),
     6 rwater(imaxout,jmaxout,kmax),cwater(imaxout,jmaxout,kmax),
     7 prate(imaxout,jmaxout), air(imaxout,jmaxout,kmax), ! air density
     8 w(imaxout,jmaxout,kmax),kh(imaxout,jmaxout,kmax),
     8 kv(imaxout,jmaxout,kmax),kz(kmax),tg(imaxout,jmaxout),
     9 drydep(imaxout,jmaxout,ndspecies),vd(ndspecies),
     a ustar(imaxout,jmaxout),pblht(imaxout,jmaxout),
     b wstar(imaxout,jmaxout),amol(imaxout,jmaxout),
     c pv(imaxout,jmaxout,kmax),dust(imaxout,jmaxout),
     d mfactor(imaxout,jmaxout),pstar(imaxout,jmaxout),
     e u10(imaxout,jmaxout),v10(imaxout,jmaxout),
     f coriolis(imaxout,jmaxout),pstate(imaxout,jmaxout,kmax),
     g rainold(imaxout,jmaxout),ss1(imaxout,jmaxout),
     h ss2(imaxout,jmaxout),ss3(imaxout,jmaxout),ss4(imaxout,jmaxout),
     i ztmp(imax,jmax,kmax+1),ztmp2(imax,jmax,kmax+1),
     j p1(imaxout,jmaxout),t1(imaxout,jmaxout),ccover(imaxout,jmaxout),
     k cldfra(imaxout,jmaxout,kmax),tsoil(imaxout,jmaxout),
     l smois(imaxout,jmaxout),swdown(imaxout,jmaxout),
     m glw(imaxout,jmaxout)

      real rwlocal(kmax),cwlocal(kmax),rhlocal(kmax),zlocal(kmax),
     & plocal(kmax),tlocal(kmax),qlocal(kmax),cldlocal(kmax),
     &     zsigmaf(kmax+1)

      real  data3din(imax,jmax,kmax+1,invar3d),   ! the lowest layer is used for computing surfacr property
     1 sigmaf(kmax+1),data2din(imax,jmax,invar2d),fland(nusgs,2),
     2 sigmah(kmax),datau(imax+1,jmax,kmax),datav(imax,jmax+1,kmax)

      real, allocatable, dimension (:,:,:) :: tsoilin, smoisin

      common /surface/zfull,z,t,u,v,p,qv,prate,zsigmaf,sigmaf   ! common block for bdrydep     

      logical wdiog,cdiog,first,timematch,debug1,initialfile
      data debug1/.false./

      character dprefix*3,dryname(ndspecies)*14
      data dryname/'SO2','SO4','NO2','NO','O3',
     1 'HNO3','H2O2','ALD2', 'HCHO','OP',
     2 'PAA','ORA','NH3','PAN','HNO2','CO'/

      character m3dnam(invar3d)*9,wrftimes(300)*19

c      data m3dnam/'U     ','V      ','T      ','QVAPOR',
c     1 'P      ','QCLOUD','QRAIN','CLDFRA'/
c      character m2dnam(invar2d)*9
c      data m2dnam/'RAINC','RAINNC','U10   ','V10     ',
c     1 'TSK','UST','PBLH','CCOVER','SWDOWN','GLW','SMOIS','TSLB'/

      data dprefix/'vd_'/

      data ust1/0.6/, ust2/0.4/        ! ustar threshold value for dust emission


      data lu/1/,lv/2/,lw/3/,lt/4/,lp/5/,lrh/6/,lkv/7/,lkh/8/
      data in_t/1/,in_q/2/,in_p1/3/,in_p2/4/,in_cwater/5/,in_rwater/6/,
     1  in_cldfra/7/,in_w/8/,
     2  in_rain1/1/,in_rain2/2/,in_u10/3/,in_v10/4/,in_tg/5/,
     3  in_ustar/6/,in_pbl/7/,in_swdown/8/,in_glw/9/,in_smois/10/,
     4  in_tsoil/11/

      integer istart3d(4),istart2d(3),istart1d(2),
     1 icount3d(4),icount2d(3),icount1d(2)

      data istart3d/1,1,1,1/
      data istart2d/1,1,1/
      data istart1d/1,1/
      namelist /control/wrffile,begyear,begdate,
     1 begtime,dtstep,wdiog,cdiog,numts,initialfile

      real c303,c302
      parameter(C303=19.83,C302=5417.4)
c --- end WORK_AREAS declarations
      ESAT(TEMK)=.611*EXP(C303-C302/TEMK)       ! for calculating saturated water vapor pressure  
      QSAT(ESAT1,PCB)=ESAT1*.622/(PCB-ESAT1)    ! TEMK is ambient temperature in K, PCB is the pressue in KPa
                                                ! QSAT is the saturated humidity in kg/kg
      rdry=287.
      g=9.8
      open(7,file='tmp.mm53nc.ini')
      read(7,control)
      call aq_locate(7,'USGS-',iflag)
      if(iflag.ne.0) then
       print*,' Can not find LANDUSE MAP factor'
       stop
      endif
      do n=1,nusgs
       read(7,*)nland,(mapusgs(nland,k),fland(nland,k),k=1,2)
       if(nland.lt.1.or.nland.gt.nusgs.or.mapusgs(nland,1).gt.lduc.
     1   or.mapusgs(nland,2).gt.lduc.or.fland(nland,1).gt.1.or.
     2	fland(nland,2).gt.1) then
       print*,' landuse input parameters error !'
       stop
       endif
      end do
      close(7)
      print*,'####################',wrffile


c------read in latitude, longitude, topographj and landuse of grids
      logunit = INIT3()
      if(.not.open3('TOPO',FSREAD3,'preprocessor')) then
       print*, 'failed to open TOPO'
       stop
      endif
      if(.not.READ3('TOPO','LAT',ALLAYS3, 0, 0, xlat)) then
        print*, 'Error in reading Latitude'
	stop
      endif
      if(.not.READ3('TOPO','LON',ALLAYS3, 0, 0, xlon)) then
        print*, 'Error in reading Latitude'
	stop
      endif
      if(.not.READ3('TOPO','TOPO',ALLAYS3, 0, 0, topo)) then
        print*, 'Error in reading topograph'
	stop
      endif
      if(.not.READ3('TOPO','LANDUSE',ALLAYS3, 0, 0, xland)) then ! USGS 25 categrories
        print*, 'Error in reading Landuse'
	stop
      endif
      if (.not. DESC3('TOPO') ) then   ! get grid information from meteorological 3d file to fill the description of 3d chemical output
        print*, 'Error getting info from TOPO'
        stop
      endif
      if(ncols3d.ne.imaxout.or.nrows3d.ne.jmaxout) then
       print*,'dimension does not match'
       print*,'nrows3d,ncols3d',ncols3d,nrows3d
       print*,'imaxout,jmaxout',imaxout,jmaxout
       stop
      endif
      iflag_log = close3('TOPO')



c      print*, 'wrffile',trim(wrffile)

      nf_status = nf_open(trim(wrffile), nf_nowrite, ncid)
      call handle_err('Error opening file',nf_status)
      nf_status = nf_inq_dimid (ncid, 'west_east_stag', dimid)
      call handle_err('west_east_stag',nf_status)
      nf_status = nf_inq_dimlen (ncid, dimid, nx)
      call handle_err('Get NX',nf_status)
      nx = nx-1
      print*, nx,'testnx'

      nf_status = nf_inq_dimid (ncid, 'south_north_stag', dimid)
      call handle_err('south_north_stag',nf_status)
      nf_status = nf_inq_dimlen (ncid, dimid, ny)
      call handle_err('Get NY',nf_status)
      ny = ny-1

      nf_status = nf_inq_dimid (ncid, 'bottom_top_stag', dimid)
      call handle_err('bottom_top_stag',nf_status)
      nf_status = nf_inq_dimlen (ncid, dimid, nz)
      call handle_err('Get NZ',nf_status)
      nz = nz-1

      if(nx.ne.imax.or.ny.ne.jmax.or.nz.ne.kmax) then
       print*,'wrf dimension does not match dog',nx,imax,ny,jmax,nz,kmax
       stop
      endif

      nf_status = nf_inq_dimid (ncid, 'soil_layers_stag', dimid)
      call handle_err('soil_layers_stag',nf_status)
      nf_status = nf_inq_dimlen (ncid, dimid, ksoil)
      call handle_err('Get ksoil',nf_status)

      print*,'find soil layer ',ksoil

      allocate(tsoilin(imax,jmax,ksoil),STAT=ierr)
      allocate(smoisin(imax,jmax,ksoil),STAT=ierr)

      nf_status=nf_inq_dimid(ncid,'Time',id_ntime)       ! total time steps
      call handle_err('Time',nf_status)
      nf_status=nf_inq_dimlen(ncid,id_ntime,ntimewrf)
      call handle_err('get Time',nf_status)

      nf_status=nf_inq_varid(ncid,'Times',id_times)
      call handle_err('Times',nf_status)
      nf_status=nf_get_var_text(ncid,id_times,wrftimes)
      call handle_err('get Times',nf_status)
      do ntmp=1,ntimewrf
       read(wrftimes(ntmp),'(i4,1x,i2,1x,i2,1x,i2,1x,i2,1x,i2)')
     1   iyear,imonth,idate,ihour,iminute,isecond
       modate(ntmp)=(iyear*1000+julian(iyear,imonth,idate))*100+ihour  ! YYYYDDDHH
c      assume wrftstep less then 24
       if(ntmp.eq.1) then
        ihour1=ihour
        iminute1=iminute
       elseif(ntmp.eq.2) then
        wrftstep= (ihour-ihour1)+(iminute-iminute1)/60.
        if (wrftstep .lt. 0) then
          wrftstep=24+wrftstep
        endif
       endif
      enddo

      wrftstep=6
      print*,'WRF wrftstep ',wrftstep
      print*,'WRF file time ',modate(1:ntimewrf)


c----- time-series start

      first=.true.
      mtime=1
      hour = 0
      d_off = 1
      print*, numts
      do while(mtime.le.numts)

       if (hour.gt.23)then
         hour = 0
         d_off = d_off+1
       endif
       nowtime=mod(begtime+(hour)*dtstep,24)  ! GMT time
       nowdate=begdate+(begtime+(mtime-1)*dtstep)/24
c       nowdate=begdate+(begtime+(mtime-1)*dtstep*d_off)/24
       nowyear=begyear

       nowtstring=(begyear*1000+nowdate)*100+nowtime
       print*,'hour',hour
       print*,'d_off',d_off
       print*,'mtime',mtime
       print*,'dtstep',dtstep
       print*,'begtime',begtime
       print*,'nowtime',nowtime
       print*,'begyear',begyear
       print*,'nowyear',nowyear
       print*,'begdate',begdate
       print*,'nowdate',nowdate
       print*,'nowstring ', nowtstring
c       read(*,*)
C***********************************************************************
c    utilize season values for NSEASON, determined from Julian date
c     WINTER = 1 = DEC, JAN, FEB
c     SPRING = 2 = MAR, APR, MAY
c     SUMMER = 3 = JUNE, JULY, AUG, SEP
c     FALL   = 4 = OCT, NOV.
C     compute leapyear bias (1 or 0), and season.
C     first, extract year from JUDATE
C***********************************************************************
      IF(MOD(nowyear,4).EQ.0) THEN
        LB=1
      ELSE
        LB=0
      END IF

C     compute NSEASON to summer, early fall, later fall, winter and transitional spring.
      IF(nowdate.GT.(59+LB).AND.nowdate.LE.(151+LB)) THEN               ! transitional spring
        nseason = 5
      ELSE IF(nowdate.GT.(151+LB).AND.nowdate.LE.(226+LB)) THEN         ! summer
        nseason = 1
      ELSE IF(nowdate.GT.(226+LB).AND.nowdate.LE.(288+LB)) THEN         ! early fall
        nseason = 2
      ELSE IF(nowdate.GT.(288+LB).AND.nowdate.LE.(335+LB)) THEN         ! later fall
        nseason = 3
      ELSE IF(nowdate.GT.(335+LB).OR.nowdate.LE.(59+LB)) THEN          ! winter
        nseason = 4
      END IF

      do ntmp=1,ntimewrf
       if(modate(ntmp).eq.nowtstring) exit
      enddo
      if(ntmp.gt.ntimewrf) then
       print*,'can not find suitable time in WRF', nowtstring,
     1  modate(1:ntimewrf)
       stop
      endif

      istart3d(4)=ntmp
      istart2d(3)=ntmp
      istart1d(2)=ntmp

      if(first) then                    ! check coordinate information

       icount1d(1)=kmax
       icount1d(2)=1
       nf_status=nf_inq_varid(ncid,'ZNU',id_znu)
       call handle_err('ZNU',nf_status)
       nf_status=nf_get_vara_real(ncid,id_znu,istart1d,icount1d,sigmah)  ! get WRF sigma level 
       call handle_err('get ZNU',nf_status)

       icount1d(1)=kmax+1
       icount1d(2)=1
       nf_status=nf_inq_varid(ncid,'ZNW',id_znw)
       call handle_err('ZNW',nf_status)
       nf_status=nf_get_vara_real(ncid,id_znw,istart1d,icount1d,sigmaf)  ! get WRF sigma level 
       call handle_err('get ZNW',nf_status)

       call get_var_3d_real_cdf(trim(wrffile),'PH',ztmp,         ! perturbation geopotential
     1                      imax,jmax,kmax+1,ntmp,debug1)
       call get_var_3d_real_cdf(trim(wrffile),'PHB',ztmp2,       ! base geopotential
     1                      imax,jmax,kmax+1,ntmp,debug1)

       do i=1,imaxout
        ii=i+ioff
        do j=1,jmaxout
	 jj=j+joff
         do k=1,kmax+1
          zfull(i,j,k)=(ztmp(ii,jj,k)+ztmp2(ii,jj,k))/9.81-topo(i,j)
	 enddo

	 do k=1,kmax
	  z(i,j,k)=0.5*(zfull(i,j,k)+zfull(i,j,k+1))
	 enddo
	enddo
       enddo

       if(.not.OPEN3('METEO3D',FSRDWR3,'wrf52nc')) then ! if OUTPUT file does not exit       
                                                       ! create it
        vglvs3d(1:kmax)=sigmah(1:kmax)

        vgtop3d=1.
        nlays3d=kmaxout

        vname3d(1)='U       '
        vname3d(2)='V       '
        vname3d(3)='W       '
        vname3d(4)='T       '
        vname3d(5)='P       '
        vname3d(6)='CLDFRA  '
        vname3d(7)='KV      '
        vname3d(8)='KH      '
        vname3d(9)='RWATER  '
        vname3d(10)='CWATER  '
        vname3d(11)='VAPOR'

        units3d(1)=' m/s    '
        units3d(2)=' m/s    '
        units3d(3)=' m/s    '
        units3d(4)=' K      '
        units3d(5)=' Pa     '
        units3d(6)=' none   '
        units3d(7)=' m2/s   '
        units3d(8)=' m2/s   '
        units3d(9)=' Kg/Kg  '
        units3d(10)=' Kg/Kg  '
        units3d(11)=' Kg/Kg  '

        nvars3d=11
       do L=1,nvars3d
        vtype3d(L)=M3REAL
        vdesc3d(L)='WRF '//vname3d(L)
       enddo
       vdesc3d(6)='Cloud Fraction'
       vdesc3d(9)='Rain Water Content'
       vdesc3d(10)='Cloud Water Content'
       vdesc3d(11)='Water Vapor Content'

       sdate3d=begyear*1000+begdate
       stime3d=begtime*10000
       tstep3d=dtstep*10000
       if(.not.OPEN3('METEO3D',FSUNKN3,'wrf52nc'))
     1   then	! FSCREA3 FSUNKN3
        print*, 'Error opening output file METEO3D'
        stop
       endif

       nvars3d=2
       vname3d(1)='AGL'
       units3d(1)='meter'
       vdesc3d(1)='WRF altitude above ground (half sigmaf level)'

       vname3d(2)='AGLF'
       units3d(2)='meter'
       vdesc3d(2)='WRF altitude above ground (full sigmaf level)'

       tstep3d=0
       if(.not.OPEN3('HEIGHT3D',FSUNKN3,'wrf52nc'))   ! output height
     1   then	! FSCREA3 FSUNKN3
        print*, 'Error opening output file METEO3D'
        stop
       endif
       if(.not.WRITE3('HEIGHT3D','AGL', sdate3d,
     1	stime3d, z)) then
        print*,'write error in HEIGHT3D'
	stop
       endif
       if(.not.WRITE3('HEIGHT3D','AGLF', sdate3d,  ! start from the second layer 
     1	stime3d, zfull(1,1,2))) then
        print*,'write error in HEIGHT3D'
	stop
       endif
       nflag=close3('HEIGHT3D')

       print*, 'sigmaf =',sigmaf

C---2D meteorological file head, include dust and sea salt emission
      nlays3d=1
      tstep3d=dtstep*10000
      vname3d(1)='PRATE '
      vname3d(2)='USTAR '
      vname3d(3)='PBLHT '
      vname3d(4)='WSTAR '
      vname3d(5)='MOL   '
      vname3d(6)='TSKIN '
      vname3d(7)='CCOVER'
      vname3d(8)='P     '
      vname3d(9)='T     '
      vname3d(10)='TSOIL'
      vname3d(11)='SMOIS'
      vname3d(12)='SWDOWN'
      vname3d(13)='GLW'
      vname3d(14)='DUST'
      vname3d(15)='U10'
      vname3d(16)='V10'
      vname3d(17)='SS1'
      vname3d(18)='SS2'
      vname3d(19)='SS3'
      vname3d(20)='SS4'

      units3d(1)=' mm/hr    '
      units3d(2)=' m/s    '
      units3d(3)=' m    '
      units3d(4)=' m/s  '
      units3d(5)=' m  '
      units3d(6)=' K  '
      units3d(7)=' none  '
      units3d(8)='Pa '
      units3d(9)='K '
      units3d(10)='K'
      units3d(11)='m3 m-3'
      units3d(12)='W/m2'
      units3d(13)='W/m2'
      units3d(14)=' gram/cm2/sec  '
      units3d(15)='m/s  '
      units3d(16)='m/s  '
      units3d(17)=' molecules/cm2/sec  '
      units3d(18)=' molecules/cm2/sec  '
      units3d(19)=' molecules/cm2/sec  '
      units3d(20)=' molecules/cm2/sec  '

      nvars3d=20
      do L=1,nvars3d
       vtype3d(L)=M3REAL
       vdesc3d(L)='WRF '//vname3d(L)
      enddo
      vdesc3d(1)='Precipication Rate'
      vdesc3d(2)='Surface Frictional Velocity'
      vdesc3d(3)='Boundary Layer Height'
      vdesc3d(5)='Monin-Obukhov Length'
      vdesc3d(6)='Ground Skin Temperature'
      vdesc3d(7)='Cloud Coverage'
      vdesc3d(8)='Pressure'
      vdesc3d(9)='Air Temperature'
      vdesc3d(10)='Top Soil Temperature'
      vdesc3d(11)='Top Soil Moisture'
      vdesc3d(12)='DOWNWARD SHORT WAVE FLUX AT GROUND SURFACE'
      vdesc3d(13)='DOWNWARD LONG WAVE FLUX AT GROUND SURFACE'
      vdesc3d(14)='Dust Emission'
      vdesc3d(15)='10m U'
      vdesc3d(16)='10m V'
      vdesc3d(17)='Sea Salt Emission (diameter 0.1-0.3 um)'
      vdesc3d(18)='Sea Salt Emission (diameter 0.3-1.0 um)'
      vdesc3d(19)='Sea Salt Emission (diameter 1.0-2.5 um)'
      vdesc3d(20)='Sea Salt Emission (diameter 2.5-10. um)'

!      if(.not.OPEN3('METEO2D',FSCREA3,'rams2nc'))
      if(.not.OPEN3('METEO2D',FSUNKN3,'rams2nc'))
     1   then	! FSCREA3 FSUNKN3
        print*, 'Error opening output file METEO2D'
        stop
      endif

C---dry Deposition file head
      nvars3d=ndspecies
      do L=1,nvars3d
       vname3d(L)=dprefix//dryname(L)
       units3d(L)=' m/s '
       vtype3d(L)=M3REAL
       vdesc3d(L)='Dry Deposition Velocity of '//dryname(L)
      enddo
      if(.not.OPEN3('DRYDEP',FSUNKN3,'wrf2nc'))
     1   then	! FSCREA3 FSUNKN3
        print*, 'Error opening output file DRYDEP'
        stop
      endif                        ! end of constructing output files

      else

       if(.not.OPEN3('METEO2D',FSRDWR3,'wrf2nc'))  ! open existed files
     1   then	! FSCREA3 FSUNKN3
        print*, 'Error opening output file METEO2D'
        stop
       endif
       if(.not.OPEN3('DRYDEP',FSRDWR3,'wrf2nc'))
     1   then	! FSCREA3 FSUNKN3
        print*, 'Error opening output file DRYDEP'
        stop
       endif

      endif

      call get_var_2d_real_cdf(trim(wrffile),'RAINC',
     1  data2din(1,1,in_rain1),imax,jmax,ntmp,debug1)
      call get_var_2d_real_cdf(trim(wrffile),'RAINNC',
     1  data2din(1,1,in_rain2),imax,jmax,ntmp,debug1)
      do i=1,imaxout
       ii=i+ioff
       do j=1,jmaxout
        jj=j+joff
	rainold(i,j)=data2din(ii,jj,in_rain1)+data2din(ii,jj,in_rain2)
       enddo
      enddo

      endif                   ! end of first time

C----process time-series data

      call get_var_2d_real_cdf(trim(wrffile),'U10',
     1 data2din(1,1,in_u10), imax,jmax,ntmp,debug1)
      call get_var_2d_real_cdf(trim(wrffile),'V10',
     1 data2din(1,1,in_v10), imax,jmax,ntmp,debug1)
      call get_var_2d_real_cdf(trim(wrffile),'TSK',
     1 data2din(1,1,in_tg), imax,jmax,ntmp,debug1)
c      print*,'in_ustar,imax,jmax,ntmp=',in_ustar,imax,jmax,ntmp
      call get_var_2d_real_cdf(trim(wrffile),'UST',
     1  data2din(1,1,in_ustar),imax,jmax,ntmp,debug1)
      call get_var_2d_real_cdf(trim(wrffile),'PBLH',
     1  data2din(1,1,in_pbl),imax,jmax,ntmp,debug1)
      call get_var_2d_real_cdf(trim(wrffile),'RAINC',
     1  data2din(1,1,in_rain1),imax,jmax,ntmp,debug1)
      call get_var_2d_real_cdf(trim(wrffile),'RAINNC',
     1  data2din(1,1,in_rain2),imax,jmax,ntmp,debug1)

      call get_var_2d_real_cdf(trim(wrffile),'SWDOWN',
     1  data2din(1,1,in_swdown),imax,jmax,ntmp,debug1)
      call get_var_2d_real_cdf(trim(wrffile),'GLW',
     1  data2din(1,1,in_glw),imax,jmax,ntmp,debug1)

      call get_var_3d_real_cdf(trim(wrffile),'U',datau,
     1  imax+1,jmax,kmax,ntmp,debug1)
      call get_var_3d_real_cdf(trim(wrffile),'V',datav,
     1  imax,jmax+1,kmax,ntmp,debug1)
      call get_var_3d_real_cdf(trim(wrffile),'P',data3din(1,1,1,in_p1),
     1  imax,jmax,kmax,ntmp,debug1)
      call get_var_3d_real_cdf(trim(wrffile),'PB',
     1 data3din(1,1,1,in_p2),imax,jmax,kmax,ntmp,debug1)
      call get_var_3d_real_cdf(trim(wrffile),'T',
     1 data3din(1,1,1,in_t),imax,jmax,kmax,ntmp,debug1)
      call get_var_3d_real_cdf(trim(wrffile),'QVAPOR',
     1 data3din(1,1,1,in_q),imax,jmax,kmax,ntmp,debug1)
      call get_var_3d_real_cdf(trim(wrffile),'QCLOUD',
     1 data3din(1,1,1,in_cwater),imax,jmax,kmax,ntmp,debug1)
      call get_var_3d_real_cdf(trim(wrffile),'QRAIN',
     1 data3din(1,1,1,in_rwater),imax,jmax,kmax,ntmp,debug1)
      call get_var_3d_real_cdf(trim(wrffile),'CLDFRA',
     1 data3din(1,1,1,in_cldfra),imax,jmax,kmax,ntmp,debug1)
      call get_var_3d_real_cdf(trim(wrffile),'W',
     1 data3din(1,1,1,in_w),imax,jmax,kmax+1,ntmp,debug1)

      call get_var_3d_real_cdf(trim(wrffile),'SMOIS',
     1 smoisin,imax,jmax,ksoil,ntmp,debug1)
      call get_var_3d_real_cdf(trim(wrffile),'TSLB',
     1 tsoilin,imax,jmax,ksoil,ntmp,debug1)
cAZedit
c        AD update read previous prec assume the file exist
c      if (first .and. (.not.(initialfile))) then
c         OPEN(10,FILE='preclasth.dat',FORM='UNFORMATTED',
c     1 ACCESS='DIRECT',RECL=imaxout*jmaxout*4,STATUS='OLD')
c         read(10,rec=1) rainold(1:imaxout,1:jmaxout)
c         print*, 'READ rainold:',rainold(1,1)
c         print*, 'READ rainold:',rainold(imaxout,jmaxout)
c         close(10)
c      endif
cAZedit
      do i=1,imaxout
       ii=i+ioff
       do j=1,jmaxout
        rainold(i,j) = 0.0!AZedit
        jj=j+joff
	tg(i,j)=data2din(ii,jj,in_tg)         ! ground temperature
	u10(i,j)=data2din(ii,jj,in_u10)       ! 10 meter wind
	v10(i,j)=data2din(ii,jj,in_v10)
c	ustar(i,j)=data2in(ii,jj,in_ustar)
c	pblht(i,j)=data2in(ii,jj,in_pbl)

        swdown(i,j)=data2din(ii,jj,in_swdown)
	glw(i,j)=data2din(ii,jj,in_glw)
	tsoil(i,j)=tsoilin(ii,jj,1)         ! top layer
	smois(i,j)=smoisin(ii,jj,1)

c  AD update
        if(first) then
         if(.not.(initialfile)) then
           prate(i,j)=(data2din(ii,jj,in_rain1)+data2din(ii,jj,in_rain2)
     1       -rainold(i,j))/wrftstep
         else
	   prate(i,j)= data2din(ii,jj,in_rain1)+data2din(ii,jj,in_rain2)
c precipitation rate in mm/hr
         endif
        else
	 prate(i,j)=(data2din(ii,jj,in_rain1)+data2din(ii,jj,in_rain2)
     1     -rainold(i,j))/wrftstep                 ! precipitation rate in mm/hr
!        AD print
!         if(mtime .eq. 4) then
!            print*, 'prate on',i,j,' :',prate(i,j)
!         end if
         prate(i,j)=amax1(prate(i,j),0.)
        endif

        do k=1,kmax
	 p(i,j,k)=data3din(ii,jj,k,in_p1)+data3din(ii,jj,k,in_p2)    ! Pa
	 t(i,j,k)=data3din(ii,jj,k,in_t)                             ! K
	 t(i,j,k)=(t(i,j,k)+300.)*(p(i,j,k)/100000.)**0.28571

	 air(i,j,k)=2.687e+19*p(i,j,k)/101300.*273.15/t(i,j,k)
         	 ! 2.687e+19=avo/22400=6.02e+23/22400 molecules/cm3
		 ! air density in molecules/cm3

         cldfra(i,j,k)=data3din(ii,jj,k,in_cldfra)
	 qv(i,j,k)=amax1(data3din(ii,jj,k,in_q),0.)                   ! water vapor mixing Kg/Kg
	 rwater(i,j,k)=amax1(data3din(ii,jj,k,in_rwater),0.)          ! rain water mixing ratio, Kg/Kg
	 cwater(i,j,k)=amax1(data3din(ii,jj,k,in_cwater),0.) ! Cloud water mixing ratio, Kg/Kg
         ! computing relative humidity, convert to KPa
	 rh(i,j,k) = 100. * qv(i,j,k) /
     &        QSAT(ESAT(t(i,j,k)), p(i,j,k) / 1000)


	 u(i,j,k)=0.5*(datau(ii,jj,k)+datau(ii+1,jj,k))            ! convert to cross point
	 v(i,j,k)=0.5*(datav(ii,jj,k)+datav(ii,jj+1,k))
	 w(i,j,k)=0.5*(data3din(ii,jj,k,in_w)+data3din(ii,jj,k+1,in_w))
	end do

        p1(i,j)=p(i,j,1)
	t1(i,j)=t(i,j,1)

C---computing surface characteristics and dry deposition
cAZedit
c       do nland=1,lduc
c        xluc(nland)=0.
c       end do
c
c       nland=int(xland(i,j))
c
c       xluc(mapusgs(nland,1))=fland(nland,1)
c       xluc(mapusgs(nland,2))=fland(nland,2)
c       call BDRYDEP(I,J,nowdate,nowtime,nseason,tg(i,j), !ground  temprature
c     1  pblht(i,j),                           !   boundary layer height
c     2  ustar(i,j),wstar(i,j),amol(i,j),vd,           !  dry deposition
c     3  XLUC,                                 ! RADM 11-category landuse
c     4  xlat(i,j),xlon(i,j))
c
c
c       if(.not.(abs(pblht(i,j)).lt.1e13.and.abs(ustar(i,j)).lt.1e13
c     1	.and.abs(wstar(i,j)).lt.1e13.and.
c     2  abs(amol(i,j)).lt.1e13)) then
c         print*,' data Overflow after calling BDRYDEP'
c	 print*, i,j,pblht(i,j),ustar(i,j),wstar(i,j),amol(i,j)
c	 stop
c	endif
c	do L=1,ndspecies
c	 drydep(i,j,L)=vd(L)
c	enddo

C----computing KV

	zlocal(1:kmax)=zfull(i,j,1:kmax)
	rhlocal(1:kmax)=rh(i,j,1:kmax)
	rwlocal(1:kmax)=rwater(i,j,1:kmax)
	cwlocal(1:kmax)=cwater(i,j,1:kmax)

	call kz_integ(kmax,kmax-1,zlocal,ustar(i,j),amol(i,j),
     1	pblht(i,j),wstar(i,j),kz)
        do k=1,kmax-1
	 kh(i,j,k)=amax1(4.*kz(k),100.)
	 kz(k)=amax1(kz(k),5.)
	 kv(i,j,k)=kz(k)
	end do
	call kz_cloud(kmax,kmax-1,kz,rhlocal,cwlocal,prate(i,j)) ! adjust cloud KZ
	kv(i,j,1:kmax-1)=kz(1:kmax-1)

C------computing cloud coverage
        if(cdiog) then

	 zlocal(1:kmax)=z(i,j,1:kmax)
	 plocal(1:kmax)=p(i,j,1:kmax)
	 tlocal(1:kmax)=t(i,j,1:kmax)
	 qlocal(1:kmax)=qv(i,j,1:kmax)
	 call jclouds(kmax,zlocal,pblht(i,j),plocal,prate(i,j),tlocal,
     1	   qlocal,cldlocal,ccover(i,j),ktop)
         cldfra(i,j,1:kmax)=cldlocal(1:kmax)
        else                 ! use WRF's cloud fraction
	 frac=0
	 sumz=0
	 do k = 1, kmax-1
          frac  = frac + cldfra(i,j,k) * (z(i,j,k+1)-z(i,j,k))
          sumz  = sumz + (z(i,j,k+1)-z(i,j,k))
         enddo
	 ccover(i,j)=frac/sumz
	endif

       end do
      end do

c      call pvs(data3din(1,1,1,in_u),data3din(1,1,1,in_v),           ! compute PV
c     1  data3din(1,1,1,in_t),data2din(1,1,in_dmf),data2din(1,1,in_xmf),
c     2  data2din(1,1,in_f),p,sngl(xcell3d),imm5,jmm5,kmax,ioff,joff,
c     3  imaxout,jmaxout,pv)

      print*,'Print out example of 2D-data'
      print*,pblht(40,40),ustar(40,40),wstar(40,40),amol(40,40),
     1 tg(40,40),t(40,40,1),qv(40,40,1),xlat(40,40),xlon(40,40)
      print*,'sample of Dry deposition'
      print*,(drydep(40,40,l),l=1,ndspecies)
c      write(*,'(100E9.2)')((drydep(i,j,1),i=1,imaxout),j=1,jmaxout)
c      if(pblht(40,40).gt.1e13) print*,' PBL40 overflow'
c      if(.not.(pblht(40,40).gt.4)) print*,' PBL40 not overflow'
c      if(pblht(40,40).lt.4) print*, ' PBL40 looks strange'

      if((.not.first).or.(.not.initialfile) ) then    ! load ustar and pblht from WRF
       do i=1,imaxout
        ii=i+ioff
        do j=1,jmaxout
         jj=j+joff
c	 pblht(i,j)=data2din(ii,jj,in_pbl)        ! WRF's PBL definition has some thing strange
	 ustar(i,j)=data2din(ii,jj,in_ustar)
	enddo
       enddo
      endif

      do i=1,imaxout
       do j=1,jmaxout
        if(.not.(abs(pblht(i,j)).lt.1e13.and.abs(ustar(i,j)).lt.1e13
     1	.and.abs(wstar(i,j)).lt.1e13.and.
     2  abs(amol(i,j)).lt.1e13)) then
         print*,' 2d data Overflow'
	 print*, i,j,pblht(i,j),ustar(i,j),wstar(i,j),amol(i,j)
	 stop
	endif
       end do
      end do

C
C     Generate dust and sea salt emissions


c      if(first) then                    ! MM5 have not output U10, V10 at the initial time
c       roughness=0.001                ! roughtness over water
c
c       do i=1,imaxout
c        do j=1,jmaxout
c	  u10(i,j)=u(i,j,1)*log(10/roughness)/log(z(i,j,1)/roughness)
c	  v10(i,j)=v(i,j,1)*log(10/roughness)/log(z(i,j,1)/roughness)
c	              ! convert to 10m wind based on power law in surface layer
c	enddo
c       enddo
c       endif

       call calc_seasalt(imaxout,jmaxout,u10,v10,rh(1,1,1),xland,ss1,
     1  ss2,ss3,ss4)

       do i=1,imaxout
        do j=1,jmaxout

	if(xlat(i,j).ge.40.and.xlat(i,j).le.50.and.xlon(i,j).ge.-95
     1 	  .and.xlon(i,j).le.-77) then               !  great lake is with fresh water
         ss1(i,j)=0.
	 ss2(i,j)=0.
         ss3(i,j)=0.
	 ss4(i,j)=0.
	endif

        if(prate(i,j).lt.1..and.nint(xland(i,j)).eq.19) then

	 if(xlat(i,j).ge.36.and.xlat(i,j).le.46.and.         ! Taklamagan Desert
     1	  xlon(i,j).ge.75.and.xlon(i,j).le.96.and.
     2    nint(xland(i,j)).eq.19.and.ustar(i,j).gt.ust2) then

          dust(i,j)=1.42e-6*ustar(i,j)**3*(ustar(i,j)-ust2)/125.  ! gram/cm^2/sec, convert 0-PM50 

	 else if((.not.(xlat(i,j).ge.26.and.xlat(i,j).le.36.and.          ! exclude Tibet Plateau
     1    xlon(i,j).ge.80.and.xlon(i,j).le.100.))
     2    .and.xlon(i,j).le.115.and.xlon(i,j).ge.0
     3    .and.ustar(i,j).gt.ust2) then
	  dust(i,j)=1.42e-6*ustar(i,j)**3*(ustar(i,j)-ust2)/125.   ! gram/cm^2/sec, convert 0-PM50 to 0-PM10
         else if(xlat(i,j).ge.50.0.and.xlat(i,j).le.90.0) then    ! exclude anything above 50N latitude
          dust(i,j)=0.0   ! gram/cm^2/sec, convert 0-PM50 to 0-PM10
	 else if(ustar(i,j).gt.ust1) then     ! define dust domain
          dust(i,j)=1.42e-6*ustar(i,j)**3*(ustar(i,j)-ust1)/125.   ! gram/cm^2/sec, convert 0-PM50 to 0-PM10
         else
	  dust(i,j)=0.
	 endif
        else
         dust(i,j)=0.
        endif
	enddo
       enddo


c	do i=1,imaxout  	! couple u,v with map factor
c	 ii=i+ioff
c	 do j=1,jmaxout
c	  jj=j+joff
c	  do k=1,kmax
c	  u(i,j,k)=u(i,j,k)/data2din(ii,jj,in_xmf)
c	  v(i,j,k)=v(i,j,k)/data2din(ii,jj,in_xmf)
c	  enddo
c	 enddo
c	enddo

C---if there is not vertical velocity, calculate it
c       if(wdiog) call divergence(6, imaxout, jmaxout, kmax,
c     1   imaxout, jmaxout, kmaxout,
c     2   sngl(xcell3d), sngl(ycell3d), topo, z, u, v, w, .true., 1.e-05)
       if(wdiog) call wind_sub(imaxout,jmaxout,kmax,u,v,w,
     1  sngl(xcell3d),sngl(ycell3d),z,air)

       indep_div=2
       kdiv=7
       do i=1,imaxout
        do k=kdiv,kmaxout
         do j=1,indep_div
	 w(i,j,k)=0.        ! see w to zero for south and north boundary
	 enddo
	 do j=jmaxout-indep_div,jmaxout
	 w(i,j,k)=0.
         enddo
	enddo
       enddo


C----output 3D Meterological data
       if(.not.WRITE3('METEO3D','U', nowyear*1000+nowdate,
     1	nowtime*10000, u)) then
        print*,'write error in METEO3D'
	stop
       endif
       if(.not.WRITE3('METEO3D','V', nowyear*1000+nowdate,
     1	nowtime*10000, v)) then
        print*,'write error in METEO3D'
	stop
       endif
       if(.not.WRITE3('METEO3D','T', nowyear*1000+nowdate,
     1	nowtime*10000, t)) then
        print*,'write error in METEO3D'
	stop
       endif
       if(.not.WRITE3('METEO3D','W', nowyear*1000+nowdate,
     1	nowtime*10000, w)) then
        print*,'write error in METEO3D'
	stop
       endif
       if(.not.WRITE3('METEO3D','P', nowyear*1000+nowdate,
     1	nowtime*10000, p)) then
        print*,'write error in METEO3D'
	stop
       endif
       if(.not.WRITE3('METEO3D','CLDFRA', nowyear*1000+nowdate,
     1	nowtime*10000, cldfra)) then
        print*,'write error in METEO3D'
	stop
       endif
       if(.not.WRITE3('METEO3D','KV', nowyear*1000+nowdate,
     1	nowtime*10000, kv)) then
        print*,'write error in METEO3D'
	stop
       endif
       if(.not.WRITE3('METEO3D','KH', nowyear*1000+nowdate,
     1	nowtime*10000, kh)) then
        print*,'write error in METEO3D'
	stop
       endif
       if(.not.WRITE3('METEO3D','CWATER', nowyear*1000+nowdate,
     1	nowtime*10000, cwater)) then
        print*,'write error in METEO3D'
	stop
       endif
       if(.not.WRITE3('METEO3D','RWATER', nowyear*1000+nowdate,
     1	nowtime*10000, rwater)) then
        print*,'write error in METEO3D'
	stop
       endif
       if(.not.WRITE3('METEO3D','VAPOR', nowyear*1000+nowdate,
     1	nowtime*10000, qv)) then
        print*,'write error in METEO3D'
	stop
       endif

C----output 2D Meterological data
       if(.not.WRITE3('METEO2D','PRATE', nowyear*1000+nowdate,
     1	nowtime*10000, prate)) then
        print*,'write error in METEO2D'
	stop
       endif
       if(.not.WRITE3('METEO2D','USTAR', nowyear*1000+nowdate,
     1	nowtime*10000, ustar)) then
        print*,'write error in METEO2D'
	stop
       endif
       if(.not.WRITE3('METEO2D','WSTAR', nowyear*1000+nowdate,
     1	nowtime*10000, wstar)) then
        print*,'write error in METEO2D'
	stop
       endif
       if(.not.WRITE3('METEO2D','MOL', nowyear*1000+nowdate,
     1	nowtime*10000, amol)) then
        print*,'write error in METEO2D'
	stop
       endif
       if(.not.WRITE3('METEO2D','PBLHT', nowyear*1000+nowdate,
     1	nowtime*10000, pblht)) then
        print*,'write error in METEO2D'
	stop
       endif
       if(.not.WRITE3('METEO2D','TSKIN', nowyear*1000+nowdate,
     1	nowtime*10000, tg)) then
        print*,'write error in METEO2D'
	stop
       endif
       if(.not.WRITE3('METEO2D','CCOVER', nowyear*1000+nowdate,
     1	nowtime*10000, ccover)) then
        print*,'write error in METEO2D'
	stop
       endif
       if(.not.WRITE3('METEO2D','P', nowyear*1000+nowdate,
     1	nowtime*10000, p1)) then
        print*,'write error in METEO2D'
	stop
       endif
       if(.not.WRITE3('METEO2D','T', nowyear*1000+nowdate,
     1	nowtime*10000, t1)) then
        print*,'write error in METEO2D'
	stop
       endif
c       do i=1,imaxout
c        do j=1,jmaxout
c	 if(t1(i,j).lt.200) then
c	  print*,'surface air temperature wrong ',i,j,t1(i,j)
c	 stop
c	endif
c	enddo
c       enddo
       if(.not.WRITE3('METEO2D','SWDOWN', nowyear*1000+nowdate,
     1	nowtime*10000, swdown)) stop
       if(.not.WRITE3('METEO2D','GLW', nowyear*1000+nowdate,
     1	nowtime*10000, glw)) stop
       if(.not.WRITE3('METEO2D','TSOIL', nowyear*1000+nowdate,
     1	nowtime*10000, tsoil)) stop
       if(.not.WRITE3('METEO2D','SMOIS', nowyear*1000+nowdate,
     1	nowtime*10000, smois)) stop

       if(.not.WRITE3('METEO2D','DUST', nowyear*1000+nowdate,
     1	nowtime*10000, dust)) then
        print*,'write error in METEO2D'
	stop
       endif
        if(.not.WRITE3('METEO2D','U10', nowyear*1000+nowdate,
     1	nowtime*10000, u10)) then
        print*,'write error in METEO2D'
	stop
       endif
        if(.not.WRITE3('METEO2D','V10', nowyear*1000+nowdate,
     1	nowtime*10000, v10)) then
        print*,'write error in METEO2D'
	stop
       endif
       if(.not.WRITE3('METEO2D','SS1', nowyear*1000+nowdate,
     1	nowtime*10000, ss1)) then
        print*,'write error in METEO2D'
	stop
       endif
        if(.not.WRITE3('METEO2D','SS2', nowyear*1000+nowdate,
     1	nowtime*10000, ss2)) then
        print*,'write error in METEO2D'
	stop
       endif
       if(.not.WRITE3('METEO2D','SS3', nowyear*1000+nowdate,
     1	nowtime*10000, ss3)) then
        print*,'write error in METEO2D'
	stop
       endif
        if(.not.WRITE3('METEO2D','SS4', nowyear*1000+nowdate,
     1	nowtime*10000, ss4)) then
        print*,'write error in METEO2D'
	stop
       endif

C----output 2D dry deposition
       if(.not.WRITE3('DRYDEP',ALLVAR3, nowyear*1000+nowdate,
     1	nowtime*10000, drydep)) then
        print*,'write error in DRYDEP'
	stop
       endif
c       do n=1,ndspecies
c        if(.not.WRITE3('DRYDEP',vname3d(n), nowyear*1000+nowdate,
c     1	nowtime*10000, drydep(1,1,n))) then
c        print*,'write error in DRYDEP'
c 	stop
c        endif
c       enddo

       mtime=mtime+1
       hour = hour+1

      do i=1,imaxout
       ii=i+ioff
       do j=1,jmaxout
        jj=j+joff
	rainold(i,j)=(data2din(ii,jj,in_rain1)+data2din(ii,jj,in_rain2))
     1                                           ! accumuluated rain fall in mm	
       enddo
      enddo

      first=.false.

      enddo

c     AD update
      print*, 'pippoprecold',imaxout,jmaxout, rainold(imaxout,jmaxout)
      print*, 'PRINT rainold:',rainold(1,1)
      print*, 'PRINT rainold:',rainold(imaxout,jmaxout)

      if(initialfile) then
        OPEN(10,FILE='preclasth.dat',FORM='UNFORMATTED',
     1 ACCESS='DIRECT',RECL=imaxout*jmaxout*4,STATUS='REPLACE')
        write(10,REC=1) ((rainold(i,j),i=1,imaxout),j=1,jmaxout)
        close(10)
      endif

 199  iflag_log=shut3()
      print*,'model conversion complete'
      end


C## cblcalz.f for Met Preprocessor (p54)
C##
      SUBROUTINE CBLCALZ(THETAV,DDSGH,ZZF,KPBLHT,PBLHT,KKMAX,
     1                   KCBLHT,CBLHT,DZZ)
C**********************************************************************
C SUBROUTINE TO COMPUTE BLACKADAR'S CONVECTIVE PBL HEIGHT
C THIS VERSION IS MODIFIED FROM CBLCALC BY JON PLEIM TO BE USED
C IN MET PREPROCESSOR. COMPUTES CBLHT WHEN THETAV(1)>THETAV(2).
C
C INPUT:
C   THETAV(KKMAX), VIRTUAL POTENTIAL TEMPERATURE PROFILE
C   DDSGH(KKMAX),  DELTA SIGMA OF EACH K LAYER
C   ZZF(0:KKMAX),  SIGMA SURFACE LEVEL HEIGHT
C   KPBLHT,        VERTICAL LAYER INDEX FOR PBLHT
C   PBLHT,         PLANETARY BOUNDARY LAYER HEIGHT
C   KKMAX,         NUMBER OF VERTICAL LAYERS
C
C OUTPUT:
C   KCBLHT,        VERTICAL LAYER INDEX FOR CBLHT
C   CBLHT,         ESTIMATED CBL HEIGHT
C
C**********************************************************************
      DIMENSION THETAV(KKMAX),DDSGH(KKMAX),ZZF(0:KKMAX),DZZ(KKMAX)
      DATA ENTR/0.2/
C
      REAL NSUM,NPSUM
      CBLHT = PBLHT
      KCBLHT= KPBLHT
C
C     ** CHECK IF THE CONVECTIVE ASSUMPTION IS SATISFIED
C
      IF( THETAV(1).LE. THETAV(2) ) THEN
        RETURN
      ELSE
        PSUM   = 0.0
        NSUM   = 0.0
        FRACTN = 1.0
C
        DO 10 K=2,KKMAX
           DTHV   = THETAV(1)-THETAV(K)
           DZZ(K) = ZZF(K)-ZZF(K-1)
           IF(DTHV.GE.0.0) THEN
             PSUM=PSUM-DTHV*DDSGH(K)   ! NOTE THAT DDSGH IS NEGATIVE
           ELSE
             EPSUM = ENTR*PSUM
             NPSUM = -DTHV*DDSGH(K)

C            ** CBL HT IS DEFINED AT NSUM=PSUM*ENTR (ENTR=0.2)
             IF(EPSUM+NSUM+NPSUM.LE.0.0) THEN
               FRACTN=-(EPSUM+NSUM)/(NPSUM-1.E-9)
               KCBLHT=K
               CBLHT  = ZZF(K-1)+FRACTN*DZZ(K)
               RETURN
             END IF
             NSUM=NSUM+NPSUM
           END IF
   10  CONTINUE
C
      END IF
C
      RETURN
      END

C######################################################################
C##
C## bdrydep.f for Met Preprocess (p54)
C##
      SUBROUTINE BDRYDEP(I,J,nowdate,nowtime,nseason,tgd, !ground  temprature
     1  pblht,    !   boundary layer height
     2  ustar,wstar, amol,  !  Monin-Obukhov length
     3  vd,xluc,  !  dry deposition
     4  XLAT,XLON)
C*********************************************************************
C*DWB M54 91110  APRIL-1991   BY  DAEWON BYUN
C*
C* FOR KMAXL-LAYER METEOROLOGY
C     NEW METHOD FOR THE ESTIMATION OF SURFACE FLUXES NOV 1990
C     NOW PBL PARAMETERS ARE ESTIMATED WITH 1ST MM4 LAYER (OUT OF 15)
C     METETEOROLOGICAL VARIABLES
C
C     Z0B IS ESTIMATED DIFFERENTLY FROM PREVIOUS VERSION
C     FOR NEW Z0B FORMULA, REFER TO BYUN AND BINKOWSKI (1991)
C
C     PBL PARAMETERIZATION USED - BYUN'S PBL SIMILARITY THEORY (1991)
C     AERODYNAMIC RESISTANCE    - ACCORDING TO BYUN (1990)
C
C-- NLEVDEP = NUMBER OF 15-LEVEL MM4 LAYERS AVERAGED INTO
C         LOWEST RADM LAYER TO COMPUTE DEPOSITION VELOCITIES
C         NLEVMAX IS SET TO 3 TO ALLOW 1, 2, OR THREE LAYERS
C         FOR THE DEPOSITION DEPENDING ON THE VERTICAL RESOLUTION
C         KMAXL.
C IF KMAXL = 15, NLEVDEP = 1
C    KMAXL = 6,  NLEVDEP = 2
C         *NOTE: STILL KMAX = 15 AS WE ARE USING 15-LAYER MM4 OUTPUT
C
C   NLEV   = NUMBER OF 15-LEVEL MM4 LAYERS USED TO ESTIMATE
C         PBL PARAMETERS. PREVIOUSLY, THERE WAS NO DISTINCTION
C         BETWEEN NLEV AND NLEVDEP. NOW WE USE LOWEST LAYER MM4
C         INFORMATION TO COMPUTE USTAR, RIB, WSTAT, MOL SINCE THE
C         FORMULATIONS ARE VALID ONLY IN THE SURFACE LAYER.
C         TO GET MOST OUT OF MM-4'S RESOLUTION, NLEV IS SET TO 1.
C-- PBLMIN: FLOOR OF PBL IS LIMITED BY THE TOP OF 1ST LAYER
C-- PBLMAX: CEILING OF PBL HEIGHT IS KEPT AT THE TOP OF 8-TH LAYER
COF THE 15-LAYER MM-4 LAYERS
C*********************************************************************
      PARAMETER (NLEV=1,NLEVMAX=3,imax=50,jmax=50,KMAX=59, !imax and jmax is equal to  imaxout and jmaxout at main model
     1 lduc=11,numsea=5,  ! lduc is landuse category
     2 LTOTG=16,LDDEP=16)

C
C
C     *** JP 4/91 - FOR CBL CALC
      REAL VWORK0(0:KMAX),VWORK1(KMAX),VWORK2(KMAX),VJI(kmax),
     1 pdum(kmax),xluc(0:lduc)
ctyh
      real tyh1(kmax),mol
C
C Z0B1 IS TEMPORARY ARRAY USED FOR ESTIMATING CELL-AVERAGED ROUGHNESS
C ZREF IS AN ARBITRARY REFERENCE HEIGHT OF WIND WITH LOG-LINEAR PROFILE
C      WE CHOSE ZREF TO BE 10 M
      DIMENSION Z0(LDUC),USTS(LDUC),DZZ(kmax+1),PK(kmax+1)
C
      REAL PBLMAXS,PBLMINS,PBLMIN,ZREF,ZSL,ZU1,ZU2,X1,X2,BETAB,A1,B1,
     1     CKWSTA,PP,RR,PPSAVE,RRSAVE,ZETA,BETAM,BETAH,PRO,ETA,VKAR,
     2     F0(LTOTG),RDIF(LTOTG),HSTR(LTOTG), RI(NUMSEA,LDUC),
     4               RLU(NUMSEA,LDUC),          RAC(NUMSEA,LDUC),
     5               RGSS(NUMSEA,LDUC),         RGSO(NUMSEA,LDUC),
     6               RCLS(NUMSEA,LDUC),         RCLO(NUMSEA,LDUC),
     7     VDSMAX(LDUC,NUMSEA),VD(LDDEP),Z00(LDUC,NUMSEA)
C
C ADD RAA(JMAX),RBA(JMAX),RCA(JMAX) FOR RESISTENCE COMPARISON (SJ)

      LOGICAL FIRSTT
      DATA FIRSTT/.TRUE./
C
       real zfull(imax,jmax,kmax+1),
     2 z(imax,jmax,kmax),t(imax,jmax,kmax),
     3 u(imax,jmax,kmax),v(imax,jmax,kmax),
     4 p(imax,jmax,kmax),qv(imax,jmax,kmax),
     5 prate(imax,jmax),zsigmaf(kmax+1),sigmaf(kmax+1)

      common /surface/zfull,z,t,u,v,p,qv,prate,zsigmaf,sigmaf

      DATA PBLMIN /50.0/, KPBLMAX /10/, KPBLMIN /1/, ZREF /10./
     1    ,BETAM /4.7/, BETAH /6.35/, G/9.8/, EOMEGA /7.29E-5/
     2    ,PRO /0.74/, VKAR /0.4/, GAMMAM /15.0/, GAMMAH /9.0/
     3    ,SMALL /1.E-04/

      data ktgd/1/,kps/2/,kwst/3/,krt/4/,kpbl/5/,kmol/6/,khfx/7/,
     1 kqfx/8/,kust/9/,kfrac/10/,kcldb/11/,kcldt/12/,kwbar/13/
C
C     ** SPECIES INDEX FOR RADM DEPOSITION SPECIES
C
      DATA LSO2 /1/, LSO4 /2/, LNO2 /3/, LNO  /4/, LO3  /5/,
     1     LHNO3/6/, LH2O2/7/, LALD /8/, LHCHO/9/, LOP /10/,
     2     LPAA/11/, LORA/12/, LNH3/13/, LPAN/14/, LHNO2/15/,
     3     LCO/16/
C
      DATA HSTR / 1.E5,    0.,   .01,  .002,  .01,
     1            1.E14, 1.E5,   15., 6000., 240.,
     2            540.,  4.E6,  2.E4, 3.6, 1.E5, 0.001/
      DATA F0   / 0.0, 0.0, 0.1, 0.0, 1.0,
     1            0.0, 1.0, 0.0, 0.0, 0.1,
     2            0.1, 0.0, 0.0, 0.1, 0.1, 0./
      DATA RDIF / 1.9, 0.0, 1.6, 1.3, 1.6,
     1            1.9, 1.4, 1.6, 1.3, 1.6,
     2            2.0, 1.6, 1.0, 2.6, 1.6, 1.2/
C
C********** Values of surface resistance components (s/m)
C       #1, Midsummer with lush vegetation
      DATA (RI(1,L),L=1,LDUC) /1.E5, 60., 120., 70., 130.,
     &    100., 1.E5, 1.E5, 80., 100., 150./
      DATA (RLU(1,L),L=1,LDUC) /1.E5, 2000., 2000., 2000., 2000.,
     &   2000., 1.E5, 1.E5, 2500., 2000., 4000./
      DATA (RAC(1,L),L=1,LDUC) /100., 200., 100., 2000., 2000.,
     &   2000., 0.1, 0.1, 300., 150., 200./
      DATA (RGSS(1,L),L=1,LDUC) /400., 150., 350., 500., 500.,
     &   100., 0.1, 1000., 0.1, 220., 400./
      DATA (RGSO(1,L),L=1,LDUC) /300., 150., 200., 200., 200.,
     &   300., 2000., 400., 1000., 180., 200./
      DATA (RCLS(1,L),L=1,LDUC) /1.E5, 2000., 2000., 2000., 2000.,
     &  2000., 1.E5, 1.E5, 2500., 2000., 4000./
      DATA (RCLO(1,L),L=1,LDUC) /1.E5, 1000., 1000., 1000., 1000.,
     &  1000., 1.E5, 1.E5, 1000., 1000., 1000./
C       #2, Autumn with unharvested cropland
      DATA (RI(2,L),L=1,LDUC) /1.E5, 1.E5, 1.E5, 1.E5, 250.,
     &  500., 1.E5, 1.E5, 1.E5, 1.E5, 1.E5/
      DATA (RLU(2,L),L=1,LDUC) /1.E5, 9000., 9000., 9000., 4000.,
     &  8000., 1.E5, 1.E5,  9000.,  9000., 9000./
      DATA (RAC(2,L),L=1,LDUC) /100., 150., 100., 1500., 2000.,
     &  1700., 0.1, 0.1, 200., 120., 140./
      DATA (RGSS(2,L),L=1,LDUC) /400., 200., 350., 500., 500.,
     &  100., 0.1, 1000., 0.1, 300., 400./
      DATA (RGSO(2,L),L=1,LDUC) /300., 150., 200., 200., 200.,
     &   300., 2000., 400., 800., 180., 200./
      DATA (RCLS(2,L),L=1,LDUC) /1.E5, 9000., 9000., 9000., 2000.,
     &   4000., 1.E5, 1.E5, 9000., 9000., 9000./
      DATA (RCLO(2,L),L=1,LDUC) /1.E5, 400., 400., 400., 1000.,
     &   600., 1.E5, 1.E5 , 400., 400., 400./
C       #3, Late autumn after frost, no snow
      DATA (RI(3,L),L=1,LDUC) /1.E5, 1.E5, 1.E5, 1.E5, 250.,
     &   500., 1.E5, 1.E5, 1.E5, 1.E5, 1.E5/
      DATA (RLU(3,L),L=1,LDUC) /1.E5, 1.E5, 9000., 9000., 4000.,
     &  8000., 1.E5, 1.E5, 9000., 9000,  9000./
      DATA (RAC(3,L),L=1,LDUC) /100., 10., 100., 1000., 2000.,
     &   1500., 0.1, 0.1, 100., 50., 120./
      DATA (RGSS(3,L),L=1,LDUC) /400., 150., 350., 500., 500.,
     & 200., 0.1, 1000., 0.1, 200., 400./
      DATA (RGSO(3,L),L=1,LDUC) /300., 150., 200., 200., 200.,
     &  300., 2000., 400., 1000., 180., 200./
      DATA (RCLS(3,L),L=1,LDUC) /1.E5, 1.E5, 9000., 9000., 3000.,
     &  6000., 1.E5, 1.E5, 9000., 9000., 9000./
      DATA (RCLO(3,L),L=1,LDUC) /1.E5, 1000., 400., 400., 1000.,
     & 600., 1.E5, 1.E5, 800., 600., 600./
C
      DATA (RI(4,L),L=1,LDUC) /1.E5, 1.E5, 1.E5, 1.E5, 400.,
     & 800., 1.E5, 1.E5, 1.E5, 1.E5, 1.E5/
      DATA (RLU(4,L),L=1,LDUC) /1.E5, 1.E5, 1.E5, 1.E5, 6000.,
     & 9000., 1.E5, 1.E5, 9000., 9000., 9000./
      DATA (RAC(4,L),L=1,LDUC) /100., 10., 10., 1000., 2000.,
     & 1500., 0.1, 0.1, 50., 10., 50./
      DATA (RGSS(4,L),L=1,LDUC) /100., 100., 100., 100., 100.,
     & 100., 0.1, 1000., 100., 100., 50./
      DATA (RGSO(4,L),L=1,LDUC) /600., 3500., 3500., 3500., 3500.,
     & 3500., 2000., 400., 3500., 3500., 3500./
      DATA (RCLS(4,L),L=1,LDUC) /1.E5, 1.E5, 1.E5, 9000., 200.,
     & 400., 1.E5, 1.E5, 9000., 1.E5, 9000./
      DATA (RCLO(4,L),L=1,LDUC) /1.E5, 1000., 1000., 400., 1500.,
     & 600., 1.E5, 1.E5, 800., 1000., 800./
C       #5, Transitional spring with partial green coverage
      DATA (RI(5,L),L=1,LDUC) /1.E5, 120., 240., 140., 250.,
     & 190., 1.E5, 1.E5 , 160., 200., 300./
      DATA (RLU(5,L),L=1,LDUC) /1.E5, 4000., 4000., 4000., 2000.,
     & 3000., 1.E5, 1.E5, 4000., 4000., 8000./
      DATA (RAC(5,L),L=1,LDUC) /100., 50., 80., 1200., 2000.,
     & 1500., 0.1, 0.1, 200., 60., 120./
      DATA (RGSS(5,L),L=1,LDUC) /500., 150., 350., 500., 500.,
     & 200., 0.1, 1000., 0.1, 250., 400./
      DATA (RGSO(5,L),L=1,LDUC) /300., 150., 200., 200., 200.,
     & 300., 2000., 400., 1000., 180., 200./
      DATA (RCLS(5,L),L=1,LDUC) /1.E5, 4000., 4000., 4000., 2000.,
     & 3000., 1.E5, 1.E5, 4000., 4000., 8000./
      DATA (RCLO(5,L),L=1,LDUC) /1.E5, 1000., 500., 500., 1500.,
     & 700., 1.E5, 1.E5, 600., 800., 800./
C
c&&&C*** #1 MIDSUMMER, MAX SURF SO4 DEP VEL(M/S), LANDTYPES 1-10
      DATA (VDSMAX(L,1),L=1,LDUC)
     1 /.001,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,.01/
C*** #2 AUTUM UNHARVESTED CROPLAND, MAX SO4 DEP VEL FOR TYPES 1-10
      DATA (VDSMAX(L,2),L=1,LDUC)
     1 /.001,0.01,0.01,0.001,0.008,0.004,0.01,0.01 ,0.01,0.01,0.01/
C*** #3 LATE AUTUMN AFTER FROST, MAX DEP VEL FOR LANDTYPES 1-10
      DATA (VDSMAX(L,3),L=1,LDUC)
     1 /.001,0.01,0.01,0.001,0.008,0.004,0.01,0.01,0.01,0.01,0.01/
C*** #4 WINDTER-SNOW ON GROUND, MAX VDS4 FOR LANDTYPES 1-10
      DATA (VDSMAX(L,4),L=1,LDUC)
     1 /.001,0.01,0.01,0.001,0.008,0.004,0.01,0.01,0.01,0.01,0.01/
C*** #5 TRANSITIONAL, SPRING & PARTIALLY GREEN MAX VDS4
      DATA (VDSMAX(L,5),L=1,LDUC)
     1 /.001,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01/
C
C****************************************
C*** #1 MIDSUMMER, LUSH VEGETATION Z0(M) FOR TYPE 1-10
      DATA (Z00(L,1),L=1,LDUC)
     1   /1.,.25,.05,1.,1.,1.,.0001,.002,.15,.1,.1/
C*** #2 AUTUMN WITH UNHARVESTED CROPLAND
      DATA (Z00(L,2),L=1,LDUC)
     1   /1.,.1,.05,1.,1.,1.,.0001,.002,.1,.08,.08/
C*** #3 LATE AUTUMN AFTER FROST Z0 (  METER; 10I6)
      DATA (Z00(L,3),L=1,LDUC)
     1   /1.,.005,.05,1.,1.,1.,.0001,.002,.1,.02,.06/
C*** #4 WINTER-SNOW ON GROUND, SUBFREEZING Z0 IN   METER
      DATA (Z00(L,4),L=1,LDUC)
     1   /1.,.001,.001,1.,1.,1.,.0001,.002,.001,.001,.04/
C*** #5 TRANSITIONAL, SPRING & PARTIALLY GREEN SHORT ANNUALS Z0
      DATA (Z00(L,5),L=1,LDUC)
     1   /1.,.03,.02,1.,1.,1.,.0001,.002,.1,.03,.06/
C****************************************
C       LAND-USE TYPES
C
C     1 URBAN LAND
C     2 AGRICULTURE
C     3 RANGE
C     4 DECIDUOUS FOREST
C     5 CONIFEROUS FOREST
C     6 MIXED FOREST WETLAND
C     7 WATER
C     8 BARREN LAND
C     9 NON-FORESTED WETLAND
C    10 MIXED AGRICULTURE/RANGELAND
C    11 ROCKY OPEN AREAS WITH LOW SHRUBS
C****************************************
C    V89223 change. References to CZJOFI replaced all ISESN, below. JKV

C    V89223 CHANGE. CZLABEL SETS UP DESCRIPTIONS FOR CZONE. JKV

      LOGICAL*1 FIRST
      DATA FIRST/.TRUE./
C
C     ** DWB FUNCTION DEFINITION FOR AERODYNAMIC RESISTANCE
C     ** RA FOR STABLE SURFACE LAYER
      RASS(D,ZU,AMOLK,CKUST) = (.74*ALOG(ZU/D)+4.7*(ZU-D)/AMOLK)/CKUST
C
C     ** RA FOR UNSTABLE SURFACE LAYER
      RAUS(A1,B1,CKUST) = 0.74*ALOG(
     1    ((2.+A1)-2.*SQRT(1.+A1))*((2.+B1)+2.*SQRT(1.+B1)) /
     2    (A1*B1) ) / CKUST
C
C     ** RA FOR STABLE PBL
      RASP(X1,X2,BETAB,CKUST)=0.74*( 2.*(BETAB+1.)*
     1  ( 1./SQRT(1.-X2) - 1./SQRT(1.-X1) ) +
     2  ALOG( ABS( (-1.+SQRT(1.-X2))*(1.+SQRT(1.-X1))/
     3 ((1.+SQRT(1.-X2))*(-1.+SQRT(1.-X1))) ) ) )/CKUST
C
C     ** RA FOR MIXED LAYER
      RAUM(X1,X2,CKWST) = ALOG( X2*(1.-X1)/(X1*(1.-X2)) )/CKWST
C
C     ** FUNCTIONS NEEDED FOR SURFACE LAYER SIMILARITY
C     ** ZETA = Z/AMOL
      FX(ZETA) = CVMGP(0.0,ABS(1.-GAMMAM*ZETA)**.25,ZETA)
      FY(ZETA) = CVMGP(0.0,ABS(1.-GAMMAH*ZETA)**.5,ZETA)
      FPSIM(ZETA,X) = CVMGP(-BETAM*ZETA,
     1  2.*ALOG(X+1.)+ALOG(1.+X*X)-2.*ATAN(X), ZETA)
      FPSIH(ZETA,Y) = CVMGP(-BETAH*ZETA, 2.*ALOG(Y+1.), ZETA)
C
C******************************************************************C
C* DEFINITION OF VARIABLES USED IN THE FUNCTION DEFINITION
C*   D     : DISPLACEMENT HEIGHT. HERE SURFACE ROUTNESS Z0 IS USED.
C*   ZU    : UPPER LIMIT OF INTEGRATION
C*   AMOL  : MONIN-OBUKHOV LENGTH
C*   CKUST : VKAR*USTAR
C*   A1    : 9.*ZU/ABS(AMOL)
C*   B1    : 9.*D /ABS(AMOL)
C*   Z1,Z2 : LOW AND UPPER BOUNDARY OF INTEGRATION
C*   X1,X2 : Z1/PBLHT, Z2/PBLHT
C*   BETAB : (4.7/0.74)*PBLHT/AMOL
C*   CKWST : VKAR*WSTAR
C******************************************************************C
C
C##################################################################
C MODIFICATION OF THE TABULATED RASISTENCE VALUES ACCORDING TO
C Bill Massman's RECOMMENDATION (B) OF 2-16-93, BY S. JIN 2/18/93
C
      IF(FIRSTT) THEN
C      FIRSTT=.FALSE.
C
C 1) FOR AGRICULTURAL LAND (LDUC=2)
      DO NS = 1,NUMSEA
C GRAPE IS TAKEN AS AGRICULTURE
 	IF(RI(NS,2).LT.9999.) RI(NS,2)= RI(NS,2)* 0.333
 	RLU(NS,2)  = RLU(NS,2)  * 1.5
        RAC(NS,2)  = RAC(NS,2)  * 1.5
	RGSS(NS,2) = RGSS(NS,2) * 1.5
	RGSO(NS,2) = RGSO(NS,2) * 1.5
	RCLS(NS,2) = RCLS(NS,2) * 1.5
	RCLO(NS,2) = RCLO(NS,2) * 1.5
       ENDDO
C
C 2) FOR GRASS (RANGE) LAND (LDUC=3)
      DO NS = 1,NUMSEA
 	RLU(NS,3)  = RLU(NS,3)  * 2.5
        RAC(NS,3)  = RAC(NS,3)  * 2.5
	RGSS(NS,3) = RGSS(NS,3) * 2.5
	RGSO(NS,3) = RGSO(NS,3) * 2.5
	RCLS(NS,3) = RCLS(NS,3) * 2.5
	RCLO(NS,3) = RCLO(NS,3) * 2.5
      ENDDO
C
C 3) FOR MIXED AGRICULTURAL AND GRASS (RANGE) LAND (LDUC=10)
      DO NS = 1,NUMSEA
 	IF(RI(NS,10).LT.9999.) RI(NS,10)= RI(NS,10)*0.666 !AGR. AND RANGE
 	RLU(NS,10)  = RLU(NS,10)  * 2.
        RAC(NS,10)  = RAC(NS,10)  * 2.
	RGSS(NS,10) = RGSS(NS,10) * 2.
	RGSO(NS,10) = RGSO(NS,10) * 2.
	RCLS(NS,10) = RCLS(NS,10) * 2.
	RCLO(NS,10) = RCLO(NS,10) * 2.
      ENDDO
C
      ENDIF

      do k=1,kmax  ! initalize vji
       vji(k)=0.
      enddo
C
C
C*** CZJOFI REPORTED OUT FOR I=23, ONCE ONLY.  89223 JKV***********C
C      IF (FIRST.AND.(I.EQ.23)) THEN
C         FIRST = .FALSE.
C        DO 101 J = 1, JMAX
C          WRITE(6,*) ' CLIMATE AT I=23,  J=',J,' IS ',
C    &      CZLABEL(CZJOFI(J))
C  101    CONTINUE
C      ENDIF
      RAINT0 = .01
      EXPON=1./3.
C
      G     =  9.81
      RHOO  =  1.225E-3
      CPAIR =  0.24
      H2OMW =  18.016
      AIRMW =  28.9644
C
C PARAMETERS NEEDED TO FIND MAX POSSIBLE SENSIBLE HEAT FLUX FROM
C ENERGY BALANCE
      STEFA   = 5.7E-8  ! W/M^2/K^4 STEFAN-BOLTZMANN CONSTANT
      EMSVTYG = 0.97  ! EMISSIVITY OF GRASSY SOIL (FOR ALL LANDUSE)
      EMSVTYA = 0.90  ! EMISSIVITY OF SKY
      ALBEDO  = 0.2    ! ROUGH ESTIMATION OF EARTH'S ALBEDO
C
      TOG   = 273.15/G
      TOGK  = TOG/.4
      GAMMA = 9.8E-3
      ROG   = 1000./34.1
C
      CON3  = 0.16*9.4
      CON4  = 0.16/0.74
C
C CON5 IS REPLACED WITH RECENT VALUES 0.0185/9.81
C     CON5  = .032/G
C
      PI    = 4.*ATAN(1.)
      ANGRAD= PI/180.
      ANG1  = ANGRAD*90./91.3125
C
C     I=IPASS
C
C Z AND ZF IS CONSTANT IN NON-HYDROSTATIC MODEL
C     DO 110 J=1,JMAX
C        PDUM(J,0)= SIGMAF(0)*PSTAR(J,I)+PTOP
C        ZF(J,0)  = 0.0
C 110 CONTINUE
C
      DO K=1,KMAX
       PK (K)   = P(i,j,k)
       DZZ(k) = zfull(i,j,k+1)-zfull(i,j,k)
      end do
C
      PBLMINS = ZFULL(i,J,KPBLMIN)
      PBLMAXS = ZFULL(i,J,KPBLMAX)
      PN      = P(i,j,NLEV)/1000.  ! Convert to KPa
      ZN      = Z(i,J,NLEV)
      ZNDEP   = Z(i,J,NLEV)
      UN      = u(i,j,1)
      VN      = v(i,j,1)
      FCORIS  = 2.*EOMEGA*SIN(ANGRAD*XLAT)
      TN=T(i,j,1)
      QVN=qv(i,j,1)

      K=1


C
C SATURATION VAPOR PRESSURE (MB) OVER WATER
      ES  = 6.1078*EXP(5384.21*(1./273.15-1./TGD))  ! TGD is the ground temperature in K
C PRESSURE IN hPa
      PSURF=P(i,j,1)/100.
      QSS = ES*.622/(PSURF-ES)
C  COMPUTE VIRTUAL TEMPERATURE OF GROUND AND AIR ABOVE GROUND
C      (assume QVground=QV(k=1) but le QVsat(Tg))
      TO = TGD*(1.+.6077*AMIN1(QSS,QVN))
      TH0 = TO*(1000./PSURF)**.286
      TV1 = TN*(1.+.6077*QVN)
      TH = TV1*(100./PN)**.286
      DTMP = TH - TH0
      DTMP = CVMGZ(1.E-10,DTMP,DTMP)


      WS  = AMAX1(UN*UN + VN*VN,1.0)
C
C  COMPUTE RICHARDSON NUMBER
      RIB = ZN *DTMP /(TOG*WS)
      RIB = CVMGZ(1.E-10,RIB,RIB)
      WS  = SQRT(WS)

C
      DO K=1,KMAX
       print*, "################AZedit",P(I,J,K)
       !AZedit begin
       IF (P(I,J,K) == 0) THEN
          PDUM = 1.0
       ELSE
       !END AZedit
       PDUM(K) = T(I,J,K)*(1.+.6077*QV(I,J,K)) !qv is the water vapor in Kg/Kg
     &            *(100000./P(I,J,K))**.286
       if(.not.(PDUM(K).le.1e13)) then
        print*,'PDUM overflow ',i,j,k,PDUM(K),T(I,J,K),P(I,J,K),
     1	 QV(I,J,K)
	stop
       end if
       ENDIF !AZedit
      END DO
C

      DO K=2,KMAX
        KLEVEL=K-1
        IF(PDUM(K).GT.TH0) GOTO 1
      ENDDO
 1    IF(KLEVEL.EQ.1)THEN
        VJI(KPBL) = ZFULL(I,J,2)
        PBLHT = VJI(KPBL)   ! STORE OLD PBLHT CALC FOR COMPARISON
      ELSE
C
C* NOW, PBL HEIGHT IS DETERMINED BY THE TOP OF THE LAYER
C* WHERE V.P.TEMP FIRST BECOMES LARGER THAN THE SFC VALUE.
C* ALSO, IMPOSE CEILING ON PBL HT FOR LONG RANGE TRANSPORT
C* CELINING IS THE TOP OF KPBLMAX(10)-TH LAYER OF 15-LAYER SAQM.
C
        VJI(KPBL) = AMIN1(PBLMAXS,ZFULL(I,J,KLEVEL+1))
C*DWB V. M53 COMPARE WITH FREE CONVECTIVE MIXING HT, CHOOSE MINIMUM
C*DWB  IF(PDUM(J,2).LE.PDUM(J,1)) THEN

        VWORK0(0)=0.
        DO KC=1,KMAX
         VWORK1(KC)=PDUM(KC)
         VWORK2(KC)=(zsigmaf(KC)-zsigmaf(KC+1))/zsigmaf(kmax)
         VWORK0(KC)=ZFULL(I,J,KC+1)
        ENDDO
        PBLHT=VJI(KPBL)
C**
        CALL CBLCALZ(VWORK1,VWORK2,VWORK0,KLEVEL,PBLHT,KMAX,
     &               KCBLHT,CBLHT,tyh1)
        VJI(KPBL)=AMIN1(PBLHT,CBLHT)
C*
      ENDIF

      PBLHT=VJI(KPBL)
      if(.not.(PBLHT.le.1e13)) then
        print*,'PBLHT overflow ',i,j,PBLHT,VWORK1,VWORK2,VWORK0,
     1	 KLEVEL,CBLHT
	stop
       end if

C
C  *** DETERMINE INSOLATION AT GIVEN GRID POINT FROM SOLAR ANGLE FORMULA
C
C        LOCAL HOUR ANGLE
C  THIS IS SUPPOSED TO BE THE TIME IN GMT AND LONG. SHOULD BE -
C  FOR NORTH AMERICA
      PDUM(1) = (NOWTIME + XLON/15. - 12)*15.*ANGRAD     ! nowtime is in GMT hour
C
C        INCLINATION ANGLE
      PDUM(2) = ANGRAD*23.5*SIN((NOWDATE+nowtime/24.-81.1875)*ANG1) ! nowdate is in Julian date
C
C        LATITUDE ANGLE (RAD)
      PDUM(3) = ANGRAD*XLAT
C *** RADIATION CALCULATION
C        IN WATT/M^2
      RADIAT = 1370.*(1.-ALBEDO)*(SIN(PDUM(3))*SIN(PDUM(2))+
     &            COS(PDUM(3))*COS(PDUM(2))*COS(PDUM(1)))
      RADIAT = AMAX1(RADIAT,0.)  ! SHOR-WAVE RADIATION
C*DWB ADDED HFXMAX TO CONTROL RUN-AWAY SIMILARITY PREDICTION
      RADNET = RADIAT - (EMSVTYG*STEFA*TGD**4
     &           - EMSVTYA*STEFA*T(i,j,1)**4)        ! NET RADIATION
C* REFERENCE KASAHARA AND WASHINGTON(1971)
C*     G=LAMDA*RADNET, LAMDA=0.19 (RADNET>0) OR 0.32 (RADNET<0)
      HFXMAXS = (1.-CVMGP(0.19,0.32,RADNET))*RADNET   ! IN WATT/M^2
      HFXMAXS = AMAX1(SMALL,HFXMAXS/1255.8)  ! IN (K M/S)
C
C   Compute if ground temperature is cooler than dew point of air above
C
      EMAX   = AMAX1(QVN*PN/(.622+QVN),1.E-20)
      TDPT   = 5417.4/(19.83-ALOG(EMAX/.611))
      DEW = CVMGP(0.0, 1.0, (TGD - TDPT))
      DEW = CVMGP(DEW, 0.0, (TGD - 273.15))
      DRAIN  = prate(i,j)/10                       ! precipitation rate in cm/hr
      RAIN = CVMGP(1.0, 0.0, (DRAIN - RAINT0))
      RAIN = CVMGP(RAIN, 0.0, (TGD - 273.15))


C
C  INITIALIZE ALL DEPOSITION VELOCITIES TO ZERO
C     WRITE (6,*) ' DO 440'
C
      DO JDUM=1,LTOTG
         VD(JDUM)=0.0
      ENDDO
C
	 RAA = 0.
	 RBA = 0.
	 RCA = 0.
C
C     WRITE(6,*) ' DO 450'
      DO ILU=1,LDUC
        Z0(ILU) = Z00(ILU,nseason)
      ENDDO
C
      ICHARN=0
C
      NWAT=4
C  COMPUTE THE LOGARITHMIC AVERAGE SURFACE ROUGHNESS FOR EACH GRID AREA
C
      DO 500 IX=1,NWAT

        Z0B = 0.
        Z0B1 = 0.
C
C     WRITE(6,*) ' DO 520'
      DO 520 ILU=1,LDUC
C
       Z0B  = ALOG(Z0(ILU))*XLUC(ILU)+Z0B
       Z0B1 = Z0B1 + XLUC(ILU)*SQRT(ALOG(ZREF/Z0(ILU)))
  520 CONTINUE

      Z0B = ZREF*EXP(-Z0B1*Z0B1)

C NOW Z0B(J) CONTAINS REGULAR ROUGHNESS LENGTH IN m
C +++ CEL AVERAGE VALUE +++
C RECOMPUTE SEA-SURFACE ROUGNESS USING CHARNOCK'S RELATION
C RECOMPUTE CELL-AVERAGE SURFACE FLUXES
C  --- THIS PART IS NOT LAND USE DEPENDENT EXCEPT FOR ADJUSTING FOR WATER
C     WRITE(6,*) ' DO 530'

      ALATDEG= XLAT
C
      WM     = WS
      THM    = TH
      THREF  = TH0
      ZF1    = ZFULL(I,J,2)
      ZMEAN1 = ZN
      Z01    = Z0B
      PBL1   = VJI(KPBL)
      PBLMIN1 = PBLMINS ! for now do not impose PBLMINS
      HFXMAX1 = HFXMAXS
C
C  CALCULATE (U*, TH*, W*, AMOL) FOLLOWING BYUN (1991)
C     WRITE(6,*) ' ZN(J)= ',ZN(J),' Z0B(J)= ',Z0B(J),' J= ',J

C
C************************************************************************
C************************************************************************
C
      CALL SFCFLUX(ALATDEG,WM,THM,THREF,ZF1,ZMEAN1,Z01,PBL1,PBLMIN,
     1             HFXMAX1,      USTAR, THSTAR, WSTAR, AMOL, PP,RR)

C
C************************************************************************
C************************************************************************
C
      IF(.not.(abs(USTAR).LT.1e13.and.abs(THSTAR).LT.1e13.and.
     1  abs(WSTAR).LT.1e13.and.abs(AMOL).LT.1e13).or.
     2  (i.eq.40.and.j.eq.40)) THEN
      WRITE(6,*) '*** DEBUG ****',i,j
      write(6,*) 'TN, WS, QVN, PN, HFXMAX1 =',TN,WS,QVN,PN,HFXMAX1
      write(6,*) 'XLUC= ',XLUC
      WRITE(6,*) 'ALAT,WM,THM,THREF input=',ALATDEG,WM,THM,THREF
      WRITE(6,*) 'ZF1,ZMEAN1,Z01,PBL1 = ',ZF1,ZMEAN1,Z01,PBL1
      WRITE(6,*) 'USTAR,THSTAR,WSTAR output= ',USTAR,THSTAR,WSTAR
      WRITE(6,*) 'AMOL, PP, RR,  PBL1 = ',AMOL, PP, RR,  PBL1

      END IF
C
      UST  = USTAR
      TST  = THSTAR
      WST  = AMAX1(1.E-05,WSTAR) ! just to prevent divide by zero
      MOL  = AMOL
      UUSTC = UST*WS
C PBL HEIGHT FOR STABLE CONDITIONS UPDATED
      PBLS = CVMGP(PBL1,PBLS,RIB)
C SAVE POWERS FOR LATER USE (SUCH AS COMPUTING PSIH FUNCTION)
      PPSAVE = PP
      RRSAVE = RR

C
C     WRITE(6,*) ' DO 540, RRSAVE = ',(RRSAVE(J),J=1,JMAX)

C *** Ocean/water roughness computed from friction velocity
C CHARNOCK'S RELATION A = 0.0185 (WU, 1982)
C     CHARNOCK, H. (1955) WIND STRESS ON A WATER SURFACE. Q.J.R.M.
C SOC. 81, 639-640.
C     WU, J. (1982) WIND-STRESS COEFFICIENTS OVER SEA SURFACE
C FROM BREEZE TO HURICANE. J. GEOPHYS. RES. 87,9704-9706.
      Z0(7)=0.0185*UST**2/G + .0001
C MODIFY PBL HEIGHT TO HAVE REASONABLE TIME SERIES
C  -- STABLE BOUNDARY LAYER HEIGHT IS LIMITED TO THE TOP OF 4TH OF 15 LAYER
C     -- MINIMUM UNSTABLE PBL HEIGHT IS MODIFIED BY 0.07 USTAR / FCORI
      VJI(KPBL)=CVMGP(AMIN1(PBLS,ZFULL(I,J,4)),
     &              VJI(KPBL),RIB )

C
C
  500 CONTINUE
C     WRITE (6,*) ' DO 500 ENDED'
C
C DO LOOP 500 ESTIMATED CELL AVERAGED VALUES OF
C        UUSTC(J) = WS(J)*UST(J)
C ASSUMING UUSTC(J) IS CONSTANT FOR A CELL INDEPENDENT OF LANDUSE,
C RECOMPUTE USTAR FOR EACH LANDUSE IN DO LOOP 600
C
C     WRITE (6,*) ' DO 600'
      DO 600 ILU=1,LDUC
C IF WATER SURFACE, ITERATE 6 TIMES TO HAVE CONVERGENCE ON USTAR, Z0 AND
C    MONIN-OBUKHOV LENGTH ESTIMATION. OTHERWISE, NO ITERATION - DWB
C
C     WRITE(6,*) '   ILU = ',ILU
      IF(ILU.EQ.7) THEN
        IWAT = 06
C       * WE NEED TO UPDATE WM,Z0 FOR WATER SFC. - ITERATE SIX  TIMES
      ELSE
        IWAT = 01
      END IF
C
C     WRITE(6,*) ' DO 610'
      DO 610 IW=1,IWAT
C

C
C% 3/8/91: SINCE WE JUST CARE HOW THE FRICTION VELOCITY WILL BE MODIFIED
C%    FOR THE DIFFERENT SURFACE ROUGHNESS FOR EACH LANDUSE,
C%    JUST USE NEUTRAL FORMULA TO UPDATE LANDUSE DEPENDENT FRICTION
C%    VELOCITY: THIS IS THE SAME AS ASSUMING THAT STABILITY AND PBL HEIGHT
C%    OVER DIFFERENT LANDUSE IS THE SAME AS FOR THE CELL AVERAGE VALUES.
C%    (since sfc temp, profiles are cell averaged value, we cannot do
C%     any better)  Next Eq. is from Byun and Wesley (1991)
      UST = SQRT(VKAR*UUSTC/ALOG(ZREF/Z0(ILU)))

C
      IF(ILU.EQ.7) Z0(7)=0.0185*UST**2/G + .0001

C
  610 CONTINUE

      USTS(ILU) = UST

C
C
C  Determine stable and unstable stability correction functions
C     WRITE (6,*) ' DO 640'
C
C WE ASSUME DEPOSITION LAYER BELONGS TO SURFACE LAYER + OUTER LAYER
C WHERE POLYNOMIAL CORRECTION IS ADDED FOR SCALAR PROFILES
C
C     ETA IS NONDIMENSIONAL HEIGHT OF THE CENTER OF DEPOSITION LAYER
C         IT IS LIMITED BY (10.*Z0/PBL <= ETA <= 0.9 )
      ETA = AMAX1( 10.*Z0(ILU)/VJI(KPBL), ZNDEP/VJI(KPBL) )
      ETA = AMIN1( ETA, 0.9 )
      ZETA= ZNDEP/MOL
      Y   = FY(ZETA)
C     AT IS COEFFICIENT FOR THE POLYNOMIAL CORRECTION
      AT  = (- FPSIH(ZETA,Y)/PRO - 1.)/RRSAVE

C     WRITE (6,*) ' i,j, ETA, ZETA, AT = ',I,J,ETA, ZETA, AT
C     RRSAVE IS POWER FOR THE POLYNOMIAL
C--------------- junk it
C      PDUM(J,6) = FPSIH(ZETA,Y)  -  AT*ETA**RRSAVE(J)
C FOR  SFC. LAYER   PLUS OUTER LAYER
C
C NOW PDUM(J,6) CONTAINS PSIH(Z/L) NEEDED FOR SO2 DRY DEPOSITION
C     COMPUTATION
C
C AERODYNAMIC RESISTANCE
C*DWB
C     IF(J.EQ.J1-1) THEN
C      WRITE(6,*) ' COMPUTE RA AND RB AT J,AMOL,ZNDEP,Z0 = ',J,MOL(J),
C    &  ZNDEP(J),Z0(J,ILU)
C      WRITE(6,*) ' PBLS, UST, WST, ZETA = ',VJI(J,KPBL),UST(J),WST(J),ZETA
C     END IF
      ZSL   = 0.1*VJI(KPBL)
      ZTEMP = AMIN1(ZNDEP,ZSL)
      ZU1   = AMAX1( Z0(ILU), ZTEMP)

C     IF(J.EQ.J1) WRITE(6,*) 'ZU1 =',ZU1
C* WHEN SETTING INTEGRATION LIMITS, WE NEED TO AVOID SINGULAR POINTS:
C* MINIMUM HT IS SET TO Z0 AND MAXIMUM HT IS SET TO 0.9*PBLHT+0.05(M)
      ZU2   = AMAX1(AMIN1(ZNDEP,VJI(KPBL)),ZU1)+0.05

      X1    = AMIN1(0.9,ZU1/VJI(KPBL))
      X2    = AMIN1(0.9+0.05/VJI(KPBL),ZU2/VJI(KPBL))

      BETAB = 6.345*VJI(KPBL)/MOL
      A1    = 9.*ZU1/ABS(MOL)
      B1    = 9.*Z0(ILU)/ABS(MOL)

      RA  = CVMGP( ( RASS(Z0(ILU),ZU1,MOL,VKAR*UST)+
     &        RASP(X1,X2,BETAB,VKAR*UST) ),
     &      ( RAUS(A1,B1,0.4*UST)+RAUM(X1,X2,VKAR*WST) ),
     &        MOL)
C
C SUB-LAYER RESISTANCE
      RB=2.83/(VKAR*UST)

C
C*** FOLLOWING CODE HAS SAME TECHNIQUE FOR VD ESTIMATION
C*** REPORTED BY WALCEK, ET AL. (1986) AE VOL. 30 NO 5, 946-964, 1986
C*** SOME OF COMPUTATION PROCESSES ARE NOT UNDERSTOOD FULLY BY DWB.
C*** FOR DETAILED BACKGROUND, CONSULT THE PAPER OR CONTACT JULIUS
C*** CHANG'S GROUP AT SUNY ALBANY. ( DEC 01, 1990 BY DAEWON BYUN)
C
C*** C********************* DEPOSITION OF SO2 ************************
C      WRITE(6,*) ' DO 650 '
C      SO2 GAS DIFFUSION IN AIR
      PDUM(2) = 0.116E-4/PN*101.325*(TN/273.15)**1.75
C     PDUM(2) IS DG
      PDUM(4) = CVMGP(AMIN1(VJI(KPBL)*.5/ZNDEP,1.),1.,MOL)
C     PDUM(4) IS PBL HT DIVIDED BY THE DEPOSITION LAYER THICKNESS
      PDUM(5) = 1./(VKAR*UST)
C     PDUM(5) IS INVERSE OF VKAR*USTAR
      ARGSO2    = PDUM(2)*PDUM(5)
      PDUM(3) =  RA+RB
C Surface resistance
C     RTRM  = 200./(RADIAT + .1)
      RTRM  = 3300./(RADIAT + .0001)
      TC    = TGD - 273.15
C     TTRM1 = 100.
C     IF(TC.GT.0.1.AND.TC.LT.39.9) TTRM1= 400./TC/(40. - TC)
      TTRM1 = 1.0
      TTRM2 = 1000.*EXP(-TC-4.)
      WET   = 1.+ 2.*AMIN1(1.,(DEW + RAIN))
      RI0   = RI(nseason,ILU)*WET*(1. + RTRM)

C
C Stomatal and mesophylic resistance
      RS   = RI0*RDIF(LSO2) + 1./(HSTR(LSO2)/3000. + F0(LSO2)*100.)
C
      RL1  = RLU(nseason,ILU) + TTRM2
      RL   = RL1
      IF(RL.LT.9999..AND.DEW.GT..9) RL= 100.
      IF(RL.LT.9999..AND.RAIN.GT..9) RL= 1./(.0002 + .333/RL1)
      IF(ILU.EQ.1.AND.WET.GT.1.5) RL= 50.
C
C     RDCL = 100.*(1.+1000./(RADIAT+10.))+RCLS(nseason,ILU)+TTRM2
      RDCL = RCLS(nseason,ILU)+TTRM2
      RG   = RGSS(nseason,ILU) + TTRM2 + RAC(nseason,ILU)
      RCAN = 1./(1./RS + 1./RL + 1./RDCL +  1./RG)
      V1   = PDUM(4)/(PDUM(3) + RCAN)
      VD(LSO2) = VD(LSO2) + XLUC(ILU)*V1

C
C********************* DEPOSITION OF SO4
C
C     WRITE(6,*) ' DO 660'

C        SCALE MAXIMUM DEPOSITION VELOCITY ACCORDING TO STABILITY
      PDUM(1) = CVMGP(1.,(1.+(ABS(-300./MOL))**0.6667),MOL)
      PDUM(1) = CVMGP(.002*PDUM(1),
     & .0009*(ABS(-VJI(KPBL)/MOL))**0.6667 ,VJI(KPBL)/MOL + 30.)
C
C*** SET VDS TO BE LESS THAN VDSMAX
      PDUM(1) = AMIN1(VDSMAX(ILU,nseason),UST*PDUM(1))
      VD(LSO4)= VD(LSO4) + XLUC(ILU)/(RA+1./PDUM(1))
C
C
C ********************* DEPOSITION OF O3

      TC    = TGD - 273.15
      TTRM2 = 1000.*EXP(-TC-4.)
C
      RS    = RI0*RDIF(LO3) + 1./(HSTR(LO3)/3000. + F0(LO3)*100.)
      RL1   = RLU(nseason,ILU) + TTRM2
      RLO3  = RL1
      IF(RL.LT.9999..AND.DEW.GT..9) RLO3 = 1./(.333E-3 + .333/RL1)
      IF(RL.LT.9999..AND.RAIN.GT..9) RLO3 = 1./(.001 + .333/RL1)
C
C
      RDCL=RCLO(nseason,ILU)+TTRM2
      RG  = RGSO(nseason,ILU) + TTRM2 + RAC(nseason,ILU)
      RCAN= 1./(1./RS + 1./RLO3 + 1./RDCL +  1./RG)
      VO3 = PDUM(4)/(RA + RB + RCAN)
C
      RAA  = RAA  + XLUC(ILU) * RA
      RBA = RBA + XLUC(ILU) * RB
      RCA = RCA + XLUC(ILU) * RCAN
C
      VD(LO3) = VD(LO3) + XLUC(ILU)*VO3

C********************** DEPOSITION OF ALL OTHER GASES
C     WRITE(6,*) ' DO 680'
      DO 680 L=3,LTOTG
       IF(L.EQ.5) GOTO 680

        TC    = TGD - 273.15
        TTRM2 = 1000.*EXP(-TC-4.)

        WET   = 1.+ 2.*AMIN1(1.,(DEW + RAIN))
C
        RS    = RI0*RDIF(L) + 1./(HSTR(L)/3000. + F0(L)*100.)
C
        RL    = (RLU(nseason,ILU) + TTRM2)/(HSTR(L)/1.E5 + F0(L))
        IF(RL.LT.9999..AND.WET.GT.1.5)
     &       RL = 1./(HSTR(L)/1.E7 + F0(L)/RLO3 + .333/RL)
C
C       RDC   = 100.*(1.+ 1000./(RADIAT +10.))
        RDC   = 0.
        RCL1  = 1./(HSTR(L)*1.E-5/(RCLS(nseason,ILU)+TTRM2)+
     &          F0(L)/(RCLO(nseason,ILU)+TTRM2))
        RGS1  = 1./(HSTR(L)*1.E-5/(RGSS(nseason,ILU)+TTRM2)+
     &          F0(L)/(RGSO(nseason,ILU)+TTRM2))
        RCAN  = 1./(1./RS + 1./RL + 1./(RDC+RCL1) +
     &          1./(RAC(nseason,ILU) + RGS1))
        V1    = PDUM(4)/(RA + RB + RCAN)
C
        VD(L) = VD(L) + XLUC(ILU)*V1
C

  680 CONTINUE
      if(i.eq.40.and.j.eq.40) then
      print*,'***********debug for dry deposition'
      print*,'WET,RS,RL,RCL1,RGS1,RCAN,V1= ',WET,RS,RL,RCL1,RGS1,RCAN,
     1 V1
      endif
C
C*** END OF CURRENT DEPOSITION CALCULATIONS ************************
  600 CONTINUE
      if(i.eq.40.and.j.eq.40) then
      print*,'***********debug for dry deposition'
      print*,'VD= ',VD
      endif

C*** END OF LANDUSE LOOP
C
C WRITE OUT THE THREE RESISTENCE HERE
C      WRITE(54,*) RAA,RBA,RCA
C
C IMPOSE CEILING AND FLOOR OF PBL HEIGHT
C now, impose PBLMINS
      VJI(KPBL) = AMAX1(VJI(KPBL),PBLMINS)
      VJI(KPBL) = AMIN1(VJI(KPBL),PBLMAXS)
      PBLHT=VJI(KPBL)

C
C CELL-AVERAGED SURFACE ROUGHNESS

      Z0B = 0.0
      Z0B1 = 0.0

      DO  ILU = 1, LDUC
        Z0B1 = Z0B1 + XLUC(ILU)*SQRT(ALOG(ZREF/Z0(ILU)))
      ENDDO

      ZRFB = ZREF*EXP(-Z0B1*Z0B1)

C
C Compute cell average friction velocity following Byun and Wesley (1991)
      DO ILU = 1, LDUC

        VJI(KUST) = VJI(KUST)+XLUC(ILU)*USTS(ILU)**2
     &          *( ALOG(ZREF/Z0(ILU)) / ALOG(ZREF/ZRFB) )

      ENDDO


      VJI(KUST) = SQRT(VJI(KUST))
      ustar= VJI(KUST)

      VJI(KHFX) = -VJI(KUST)*TST
      VJI(KQFX) = TST*(QSS-qv(I,J,1))/(TH0-TH)
      VJI(KMOL) = -TGD*VJI(KUST)**3
     &   /CVMGZ(1.0E-6,(0.4*9.8*VJI(KHFX)),VJI(KHFX) )
      AMOL=VJI(KMOL)

      VJI(KWST) = VJI(KUST)*( VJI(KPBL)/
     &   (VKAR*ABS(VJI(KMOL)) ))**0.3333333
      VJI(KWST) = CVMGP(0.0,VJI(KWST),VJI(KMOL))
      WSTAR=VJI(KWST)

      IF(.not.(abs(USTAR).LT.1e13.and.abs(pblht).lt.1e13.and.abs(WSTAR)
     1  .LT.1e13.and.abs(AMOL).LT.1e13)) THEN
      print*,' data overflow in the end of BDRYDEP !'
      print*, 'ustar',ustar
      print*, 'pblht',pblht
      print*, 'amol',amol
      print*, 'wstar', wstar
      print*, 'xluc',xluc
      print*, 'vd',vd
       stop
      endif
C
C%   IT HAS BEEN CONVENTIONAL FOR RADM PREPROCESSOR TO PROVIDE
C%   NEGATIVE HEAT FLUX = USTAR*THSTAR = - HFLUX FOR RADM
C%   ALTHOUGH I DON'T LIKE IT, I HAVE TO FOLLOW THE CONVENTION
C%   THUS, VJI(J,KHFX) CONTAINS -HFLUX IN (M K/S)
C     DO 750 J=1,J1
C      DO 750 J=1,JMAX
C      VJI(KHFX) = -VJI(J,KHFX)
C  750 CONTINUE

      RETURN
      END

C######################################################################
C######################################################################
C######################################################################
C######################################################################
C##
C## sfcflux.f for Met Preprocessor (p54)
C##
      SUBROUTINE SFCFLUX(ALATDEG,WM,THM,THREF,ZF1,ZMEAN1,ZO, HPBL,
     1              PBLMIN,HFXMAX, USTAR, THSTAR,WSTAR, AMOL, P, R)
C***********************************************************************
C SUBROUTINE TO FIND SURFACE FLUXES FROM MM-4 IST LAYER MEAN WIND
C AND TEMPERATURE USING THE POWER-LAW SIMILARITY OF PROFILE FUNCTIONS
C
C REFER TO PAPER BY BYUN(1991)
C
C   " DETERMINATION OF SIMILARITY FUNCTIONS OF
C     THE RESISTANCE LAWS FOR THE PLANETARY BOUNDARY LAYER
C     USING SURFACE LAYER SIMILARITY FUNCTIONS "
C
C  - HERE WE NEGLECT TURNING OF WIND TO SIMPLIFY OUR CALCULATION
C  - MINIMUM BOUNDARY LAYER HEIGHT = TOP OF THE FIRST LAYER
C
C  INPUTS
C      ALATDEG: LATITUDE IN DEGREE
C      WM    : 1ST LAYER U(E-W) WIND SPEED FROM MM-4 15-LAYER    (M/S)
C      THM   : 1ST LAYER P.TEMP FROM MM-4 15-LAYER                (K)
C      THREF : GROUND TEMPERATURE FROM MM-4 15-LAYER              (K)
C      ZF1   : HEIGHT OF TOP OF FIRST LAYER                       (M)
C      ZMEAN1: SIGMA MID-HEIGHT OF 1ST LAYER
C      ZO    : SURFACE ROUGHNESS                                  (M)
C      HPBL  : BOUNDARY LAYER HEIGHT                              (M)
C            : FOR STABLE   - INPUT MODEL MINIMUM VALUE
C            : FOR UNSTABLE - INPUT INVERSION HEIGHT
C      HFXMAX: MAX. SENSIBLE HEAT FLUX FROM ENERGY BALANCE        (K M/S)
C
C  OUTPUTS
C      USTAR : FRICTION VELOCITY                                  (M/S)
C      THSTAR: TEMPERATURE SCALE                                  (K)
C      WSTAR : CONVECTIVE VELOCITY SCALE                          (M/S)
C      AMOL  : MONIN-OBUKHOV LENGTH                               (M)
C      HPBL  : BOUNDARY LAYER HEIGHT                              (M)
C            : FOR STABLE - INTERPOLATED BETWEEN NEUTRAL/STABLE
C      P     : POWER FOR MOMENTUM PROFILE
C      R     : POWER FOR TEMPERATURE PROFILE
C
C  VERSION: 04/12/91
C          BY DAEWON BYUN
C***********************************************************************
      CHARACTER*80 FILENAME
C
C DATA
C
C** DOES NOT USE NEW PARAMETERIC VALUES A LA HOGSTROM(1987) YET
C DATA FOR CONSTANTS - SURFACE LAYER SIMILARITY THEORY
C          SIMILIM - LIMIT AT WHICH SIMILARITY THEORY IS ALLOWED
      DATA SIMILIM /0.008/
C*OLD
      DATA BETAM/4.7/, BETAH/6.3513514/, PRO/0.74/, ALAMDAO /0.07/
C      DATA BETAM/6.0/, BETAH/8.2105263/, PRO/0.95/, ALAMDAO /0.3/
C*OLD
      DATA GAMMAM/15.0/,GAMMAH/9.0/, VKAR/0.40/,RICR/0.21/,RIMAX/0.20/
      DATA RIMIN /-4.75/
C      DATA GAMMAM/19.3/,GAMMAH/11.6/,VKAR/0.4/,
C     &     RICR/0.1666/,RIMAX/0.1666/
      DATA G /9.8/, EOMEGA /7.292E-5/, DEGRAD /0.0174532/, TOL /1.E-3/
      DATA BIG /1.0E04/, SMALL /1.0E-4/, ITMAX/20/, CH /0.8/
C ERROR LIMIT NEEDED TO FIND P AND R IS BASED ON A(0), B(0)
      DATA ANEUT/1.7/, BNEUT/4.5/
C A1 = -1./15., USED FOR FINDING ZETA FOR UNSTABLE ATMOSPHERES
      DATA A1 /-0.06666667/
C FLOOR AND CEILING OF POWERS (REFFER TO BYUN, 1991)
      DATA PMIN /1.0/, PMAX /3.0/, RMIN /1.0/, RMAX /2.5/
C
C234567890123456789012345678901234567890123456789012345678901234567890
C
C***
C FUNCTION DEFINITION
C
C*** FUNCTIONS FOR STABLE AND UNSTABLE CASES
C
      FX(ZETA) = CVMGP(0.0,ABS(1.-GAMMAM*ZETA)**.25,ZETA)
      FY(ZETA) = CVMGP(0.0,ABS(1.-GAMMAH*ZETA)**.5,ZETA)
C
      FPHIM(ZETA)=CVMGP(1.0+BETAM*ZETA,ABS(1.-GAMMAM*ZETA)**(-.25),ZETA)
      FPHIH(ZETA)=CVMGP( PRO*(1.0+BETAH*ZETA),
     &                   ABS(1.-GAMMAH*ZETA)**(-.5), ZETA)
C
      FPSIM(ZETA,X) = CVMGP(-BETAM*ZETA,
     &  2.*ALOG(X+1.)+ALOG(1.+X*X)-2.*ATAN(X), ZETA)
      FPSIH(ZETA,Y) = CVMGP(-BETAH*ZETA, 2.*ALOG(Y+1.), ZETA)
C
      FGAMU(ZETA,X) = CVMGP(-.5*BETAM*ZETA*ZETA,
     & -( 4.*X**3/3.-X**4 +
     &  (X**4-1.)*(ALOG(1.+X*X)+2.*ALOG(1.+X)-2.*ATAN(X)) )/GAMMAM
     & ,ZETA)
      FGAMH(ZETA,Y) = CVMGP(-.5*BETAH*ZETA*ZETA,
     & -( 2.*(Y*Y-1.)*ALOG(1.+Y)+2.*Y-Y*Y )/GAMMAH
     & ,ZETA)
C
C FUNCTIONS FOR SURFACE POWER
      FPS(AMU,ETAS,ETAO,XS,XO) = FPHIM(ETAS*AMU)/
     & (ALOG(ETAS/ETAO)-FPSIM(ETAS*AMU,XS)+FPSIM(ETAO*AMU,XO))
      FRS(AMU,ETAS,ETAO,YS,YO) = FPHIH(ETAS*AMU)/
     & (PRO*(ALOG(ETAS/ETAO)-FPSIH(ETAS*AMU,YS)+FPSIH(ETAO*AMU,YO)))
C
C FUNCTIONS FOR COEFFICIENTS FOR PBL PROFILES
C
      FAU(AMU,ETAO,P,XA,XO)=-((P+1.)/(P*(1.-ETAO)**P))*
     &       ( ETAO*ALOG(ETAO)/(1.-ETAO) + 1. - FPSIM(AMU,XA)
     &     + (FGAMU(AMU,XA) - FGAMU(AMU*ETAO,XO))/(AMU*(1.-ETAO)) )
      FAV(AMU,ETAO,ALAMDA,Q)=((Q+1.)/(Q*(1.-ETAO)**(Q+1)))*VKAR/ALAMDA
      FAT(AMU,ETAO,AKAPA,R)= (1./(R*(1-ETAO)**(R-1)))*
     &                       ( VKAR*AKAPA/PRO-FPHIH(AMU)/PRO-1. )
C
C FUNCTIONS FOR NORMALIZED PBL PROFILES
C
C     UPBLN = VKAR*UPBL/USTAR
      UPBLN(AMU,ETA,ETAO,P,AU,X,XO) =
     &                FPSIM(AMU*ETAO,XO)+ ALOG(ETA/ETAO)
     &              - FPSIM(AMU*ETA,X) + AU*(ETA-ETAO)**P
C     VPBLN = VKAR*VPBL/USTAR
      VPBLN(AMU,ETA,ETAO,Q,AV) = -AV*(ETA-ETAO)**Q
C     DTPBLN = VKAR*(TPBL-THTETO)/(PRO*THSTAR)
      DTPBLN(AMU,ETA,ETAO,R,AT,Y,YO) =
     &                 FPSIH(AMU*ETAO,YO) + ALOG(ETA/ETAO)
     &               - FPSIH(AMU*ETA,Y) + AT*(ETA-ETAO)**R
C
      USFCN(AMU,ETA,ETAO,X,XO) =
     &                FPSIM(AMU*ETAO,XO)+ ALOG(ETA/ETAO)
     &              - FPSIM(AMU*ETA,X)
      DTSFCN(AMU,ETA,ETAO,Y,YO) =
     &                 FPSIH(AMU*ETAO,YO) + ALOG(ETA/ETAO)
     &               - FPSIH(AMU*ETA,Y)
C
C FUNCTIONS NEEDED TO FIND P AND R
C
      FHM(AMU,ETAS,ETAO,XA,XS,XO) =
     &  ABS( ALOG(ETAS/ETAO) - FPSIM(AMU*ETAS,XS) + FPSIM(AMU*ETAO,XO))
     &/ ABS( ETAO*ALOG(ETAO)/(1.-ETAO) + 1. - FPSIM(AMU,XA)
     &+ (FGAMU(AMU,XA) - FGAMU(AMU*ETAO,XO))/(AMU*(1.-ETAO)) )
      FHH(AMU,ETAS,ETAO,YS,YO,AKAPA) =
     &  ABS( ALOG(ETAS/ETAO) - FPSIH(AMU*ETAS,YS) + FPSIH(AMU*ETAO,YO))
     &/ (ABS( VKAR*AKAPA/PRO-FPHIH(AMU)/PRO-1. )*(1.-ETAO))
C
C   SURFACE LAYER  ANALYTICAL SOLUTION FOR RI
C   STABLE CASE
      FZS(RI,BM,BH,PRO)
     &  =(-2.*BH*RI+1.-SQRT(1.+4.*(BH-BM)*RI/PRO))
     &  /(2.*BH*(BM*RI-1.))
C
C   UNSTABLE CASE
      FZU(QI,THTI,A1)=-2.*SQRT(QI)*COS(THTI/3.)-A1/3.
      FQI(GM,GH,SI) =(1./GM**2+3.*GH*SI*SI/GM)/9.
      FPI(GM,GH,SI) = ( -2./GM**3+9.*(-GH/GM+3.)*SI*SI/GM )/54.
      FTHTI(PI,QI)=ACOS(PI/QI**1.5)
C
C
C*** END OF FUNCTION DEFINITION
C
      PHI = 4.*ATAN(1.0)
      HPBL=AMAX1(HPBL,ZF1,ZO/SIMILIM)
C
C COMPUTE SURFACE PARAMETERS
C
C
      ALAT = DEGRAD*ALATDEG
      FCORI = 2.*EOMEGA*abs(SIN(ALAT))
      AKAPA = 0.0
      SPEED = AMAX1(1.0,WM)
      ZZO    = AMAX1(ZMEAN1/ZO, 10.0)
C LIMIT Z/ZO 10.0
      ZINV   = 1./ZZO
      ALNZZO = ALOG(ZZO)
C      ZETAMAX = ALNZZO*FZS(RIMAX,BETAM,BETAH,PRO)/(1.-ZINV)
C      WRITE(*,*) ' ZETAMAX,RIMAX =', ZETAMAX,RIMAX
C
C ESTIMATE FLUXES USING SURFACE LAYER SIMILARITY - ANALYTICAL SOLUTION
C
C 1. ESTIMATE BULK RICHARDSON NUMBER
C
      RIBEST = G*(ZMEAN1-ZO)*(THM-THREF)/(THREF*SPEED*SPEED)
C
C LIMIT THE TEMPERATURE DIFFERENCE (THM-THREF) TO FOLLOW RIMAX ALLOWED
C
      IF(RIBEST.GT.RIMAX) THEN
       THREF1 = G*(ZMEAN1-ZO)*THM/
     &          (SPEED*SPEED*RIMAX+G*(ZMEAN1-ZO))
C     WRITE(*,*) ' INPUT, MODIFIED SFC TEMP IN K = ', THREF, THREF1
      ELSE IF(RIBEST.LT.RIMIN) THEN  ! RIMIN=-4.75
       THREF1 = G*(ZMEAN1-ZO)*THM/
     &          (SPEED*SPEED*RIMIN+G*(ZMEAN1-ZO))
      ELSE
       THREF1 = THREF
      END IF
C
      IF(RIBEST.GT.0.005) THEN
       ZETA1 = ALNZZO*FZS(RIBEST,BETAM,BETAH,PRO)/(1.-ZINV)
       ZETA1 = AMIN1(2.0, ZETA1) ! LIMITATION OF CURRENT THEORY
      ELSE IF (RIBEST.LT.-SMALL) THEN
       SB   = RIBEST/PRO
       QB   = FQI(GAMMAM,GAMMAH,SB)
       PB   = FPI(GAMMAM,GAMMAH,SB)
C TEST DETERMINANT
       DET  = QB**3 - PB*PB
       IF(DET.GE.0.0) THEN
        THTB = FTHTI(PB,QB)
        ZETA1 = FZU(QB,THTB,A1)
       ELSE
C*        WRITE(*,*)'L,K,DET,PB,QB = ',L,K,DET,PB,QB
        ARG2      = (SQRT(-DET) + ABS(PB))**(1./3.)
        ZETA1 = -( ARG2 + QB/ARG2 ) + 1./(3.*GAMMAM)
       END IF
        ZETA1 = ALNZZO*ZETA1/(1.-ZINV)
        ZETA1 = AMAX1(-5.0, ZETA1)  ! LIMITATION OF CURRENT THEORY
      ELSE
        ZETA1=SIGN(1.,RIBEST)/BIG  ! NEUTRAL CASE
      END IF
C
      ZETAO = ZETA1*ZO/ZMEAN1
      XO    = FX(ZETAO)
      YO    = FY(ZETAO)
      X1    = FX(ZETA1)
      Y1    = FY(ZETA1)
      USTAR1 = VKAR*SPEED/(ALNZZO-FPSIM(ZETA1,X1)+FPSIM(ZETAO,XO))
      THSTAR1= VKAR*(THM-THREF1)/
     &        (  PRO*( ALNZZO-FPSIH(ZETA1,Y1)+FPSIH(ZETAO,YO) )  )
C
C EVERYTHING ABOVE IS BASED ON THE SURFACE LAYER SIMILARITY
C NOW, TRY TO FIND FLUXES USING PBL SIMILARITY THEORY
C
      THREF = THREF1
      USTAR = USTAR1
      THSTAR= THSTAR1
      HPBL = CVMGP(HPBL,AMAX1(ALAMDAO*USTAR/FCORI,HPBL),ZETA1)
C      WRITE(*,*) ' ============================== '
C      WRITE(*,*) ' SURFACE SIMILARITY RESULT      '
C      WRITE(*,*) ' HPBL, RIBEST, ZETA = ', HPBL, RIBEST,ZETA1
C      WRITE(*,*) ' THREF1, THSTAR1, USTAR1 = ',THREF1,THSTAR1,USTAR1
C
C                        SIMILIM IS MAXIMUM ALLOWED ETAO FOR SIMILARITY
      IF( (ZO/HPBL) .GT. SIMILIM ) THEN
C RAISE HPBL AT LEAST ZO*(1/SIMILIM). CURRENTLY SIMILIM = 0.008
        HPBL = ZO/SIMILIM
      END IF
C
C
C  ASSUME NEUTRAL CONDITION
      IF(ABS(RIBEST).LT.SMALL) THEN
C% NEUTRAL CASE - WE KNOW PNEUT ALREADY
      DO 100 IT=1,2  ! ITERATE ONCE TO UPDATE
      ETAO = AMIN1(ZO/HPBL, 0.01)
      ETAO = AMAX1(1.0E-6,ETAO)
      ETAS = AMAX1(10.*ETAO, 0.1)
      ETAF = AMIN1(1.0, ZF1/HPBL)
      ETAF = AMAX1(1.1*ETAS,ETAF)
      ETERM1 = ETAO*ALOG(ETAO)/(1.- ETAO)
      AMOL1 = SIGN(1.,(THM-THREF1))*BIG
      PNEUT = (ETERM1+1.)/(ANEUT - ETERM1 - 1.0)
      QNEUT = 1./(ALAMDAO*BNEUT*(1-ETAO)/VKAR-1.0)
      RNEUT = PNEUT
      AUNEUT= -((PNEUT+1.)/(PNEUT*(1.-ETAO)**PNEUT))*
     &       ( ETAO*ALOG(ETAO)/(1.-ETAO) + 1. )
      ATNEUT= FAT(0.0,ETAO,AKAPA,RNEUT)
C     WRITE(*,*) ' %% AU, AT =',AU,AT
C
      TERMO = ETAF*ALOG(ETAF/ETAO) + ETAO - ETAF
C
      TERM3 = AUNEUT*(ETAF-ETAO)**(P+1.)/(P+1.)
C
      TERM6 = ATNEUT*(ETAF-ETAO)**(R+1.)/(R+1.)
C
      USNEW = VKAR*SPEED*(ETAF-ETAO)/(TERMO+TERM3)
      THSNEW= VKAR*(THM-THREF)*(ETAF-ETAO)/PRO
C
      THSTAR = THSNEW/(TERMO+TERM6)
      USTAR  = USNEW
      AMOL   = AMOL1  ! ESSENTIALLY A BIG NUMBER
      WSTAR  = 0.0
      HPBL   = AMAX1(ALAMDAO*USTAR/FCORI,HPBL)
  100 CONTINUE
C% END OF NEUTRAL COMPUTATION
      GOTO 2
C%
      ELSE
C% FOR DIABATIC CASES
      AMOL1 = THREF1*USTAR1*USTAR1/(VKAR*G*THSTAR1)
      AMOL  = SIGN(1.,THSTAR1)*AMIN1(10000.,AMOL1)
      ITER = 0
C     COMPUTE POWERS FOR SURFACE LAYER
C*** FOR VARIABLE VALUE OF P ( P,R DEPENDENT ON AMU)
      ETAO = AMIN1(ZO/HPBL, 0.01)
      ETAO = AMAX1(1.0E-6,ETAO)
      ETAS = AMAX1(10.*ETAO, 0.1)
      ETAF = AMIN1(1.0, ZF1/HPBL)
      ETAF = AMAX1(1.1*ETAS,ETAF)
      ETERM1 = ETAO*ALOG(ETAO)/(1.- ETAO)
      AMU  = HPBL/ AMOL
      PNEUT = (ETERM1+1.)/(ANEUT - ETERM1 - 1.0)
      QNEUT = 1./(ALAMDAO*BNEUT*(1-ETAO)/VKAR-1.0)
      AMOL1 = SIGN(1.,(THM-THREF1))*BIG
      RNEUT = PNEUT
      AUNEUT= -((PNEUT+1.)/(PNEUT*(1.-ETAO)**PNEUT))*
     &       ( ETAO*ALOG(ETAO)/(1.-ETAO) + 1. )
C
      DELTAM = ABS(AUNEUT*(0.1)**PNEUT/ALOG(.1/ETAO))
      DELTAH = DELTAM
C
      XO = FX(ETAO*AMU)
      XS = FX(ETAS*AMU)
C
      YO = FY(ETAO*AMU)
      YS = FY(ETAS*AMU)
C
      PS = FPS(AMU,ETAS,ETAO,XS,XO)
      RS = FRS(AMU,ETAS,ETAO,YS,YO)
      P  = PS
      R  = RS
C      WRITE(*,*) 'PS, RS =', PS,RS
    1 CONTINUE
      ITER = ITER + 1
      IF(ITER.GT.ITMAX) THEN
        WRITE(*,*) ' TOO MUCH ITERATION REQUIRED: '
        WRITE(*,*) ' ESTIMATION OF PBL FLUXES FAILED '
        GOTO 999
      END IF
C
C FOR STABLE CASE
      IF(AMOL.GT.0.0) THEN
C** ESTIMATE HPBL USING ZILITINKEVICH'S INTERPOLATION
       ARGCH = VKAR*CH*CH*AMOL/ALAMDAO
       HPBL = .5*( - ARGCH +
     &   SQRT(ARGCH*ARGCH + 4.*VKAR*USTAR*CH*CH*AMOL/FCORI) )
C
C SIMILARITY THEORY REQUIRES ( HPBL > 100*ZO) AT LEAST
C FOR THIS, USE SIMILIM = 0.008
       HPBL = AMAX1(HPBL,ZO/SIMILIM)
C
       IF((HPBL.LT.ZF1).AND.(ITER.GE.ITMAX) ) THEN
C IF PBL METHOD FAILS, USE SURFACE LAYER SIMILARITY ESTIMATED VALUES
        THREF = THREF1
        AMOL  = AMOL1
        USTAR = USTAR1
        THSTAR= THSTAR1
        HPBL = AMAX1(HPBL,ZF1,PBLMIN)
        GOTO 2
       ELSE
        HPBL = AMAX1(HPBL,ZF1,PBLMIN)
       END IF
      END IF

C      WRITE(*,*)'At 1797 of SFCFLUX      '
C      WRITE(*,*)'HPBL,ETAO,MOL,ARGCH,CH = ', HPBL, ETAO,AMOL,ARGCH,CH
C      WRITE(*,*)'FCORI,THREF, THSTAR, USTAR = ',FCORI,THREF,THSTAR,USTAR

C
      AMU  = HPBL/ AMOL
      ALAMDA = FCORI*HPBL/USTAR
      ZETAO= ZO  / AMOL
      ETAO = AMIN1(ZO/HPBL, 0.01)
      ETAO = AMAX1(1.0E-6,ETAO)
      ETAS = AMAX1(10.*ETAO, 0.1)
      ETAF = AMIN1(1.0, ZF1/HPBL)
      ETAF = AMAX1(1.1*ETAS, ETAF)
C
      XO = FX(ETAO*AMU)
      XA = FX(1.*  AMU)
      XS = FX(ETAS*AMU)
      XF = FX(ETAF*AMU)
C
      YO = FY(ETAO*AMU)
      YA = FY(1.*  AMU)
      YS = FY(ETAS*AMU)
      YF = FY(ETAF*AMU)
C
C FIND POWER P AND R
C
C     WRITE(*,*) 'XO,YO,XS,YS,XA,YA=',
C    &   XO,YO,XS,YS,XA,YA
C     WRITE(*,*) 'AMU=',AMU
      ANS = (ETAS-ETAO)/(1.-ETAO)
      ALOGNS = ALOG(ANS)
C
C** ESTIMATE P AND R USING THE ASYMPTOTIC MAGNITUDE COMPARISONS
C
      P1 = P
      FHM1 = FHM(AMU,ETAS,ETAO,XA,XS,XO)
      DO 520 IT=1,ITMAX
C     WRITE(*,*) 'P1 = ',P1
C*      PNEW=ALOG(ABS(P1*DELTAM*FHM1/(1.+P1)))/ALOGNS
C* NEWTON-RAPSON METHOD
      GP = ALOG((P1+1.)/P1) + P1*ALOGNS - ALOG(DELTAM*FHM1)
      G1P= ALOGNS - 1./(P1*(1.+P1))
      PNEW = P1 - GP/G1P
C SOMETIMES, THE PROFILES PROVIDED BY THE MM4-DOES NOT MAKE SENSE
C COMPARED TO THE IDEALIZED HORIZONTALLY HOMOGENEOUS STEADY STATE PBL
C --- TO PREVENT ABNORMAL POWERS, LIMIT FLOOR AND CEILING OF THE POWER
      PNEW = AMIN1(PMAX,PNEW)
      PNEW = AMAX1(PMIN,PNEW)
C
      IF(ABS((PNEW-P1)/PNEW).LT.TOL) THEN
        P1  = PNEW
        ITP = IT
        GOTO 521
      ELSE
        P1 = PNEW
      END IF
  520 CONTINUE
C
C NOT CONVERGED WITHIN ITMAX
      WRITE(6,*) '*************************************************'
      WRITE(6,*) ' CONVERGENCE FAILED FOR POWER P '
      WRITE(6,*) ' FAILED AT ALAMDA,AKAPA,AMU = ',ALAMDA,AKAPA,AMU
      WRITE(6,*) '*************************************************'
      STOP
  521 CONTINUE
C     WRITE(*,*) 'P1 IS DETERMINED = ',P1
C
C** USE ITERATION FOR R STARTING WITH INITIAL R = RS
C
      R1 = R
      FHH1 = FHH(AMU,ETAS,ETAO,YS,YO,AKAPA)
      DO 540 IT=1,ITMAX
C     WRITE(*,*) 'R1 = ',R1
C      RNEW = ALOG(R1*DELTAH*FHH1)/ALOGNS
C* USE NEWTON-RAPSON METHOD
      GR   = -ALOG(R1) + R1*ALOGNS - ALOG(DELTAH*FHH1)
      G1R  = -1./R1 + ALOGNS
      RNEW = R1 - GR/G1R
C --- TO PREVENT ABNORMAL POWERS, LIMIT FLOOR AND CEILING OF THE POWER
      RNEW = AMIN1(RMAX,RNEW)
      RNEW = AMAX1(RMIN,RNEW)
C
      IF(ABS((RNEW-R1)/RNEW).LT.TOL) THEN
        R1  = RNEW
        ITR = IT
        GOTO 541
      ELSE
        R1 = RNEW
      END IF
  540 CONTINUE
C
C NOT CONVERGED WITHIN ITMAX
      WRITE(6,*) '*************************************************'
      WRITE(6,*) ' CONVERGENCE FAILED FOR POWER R '
      WRITE(6,*) ' FAILED AT ALAMDA,AKAPA,AMU = ',ALAMDA,AKAPA,AMU
      WRITE(6,*) '*************************************************'
      STOP
  541 CONTINUE
C     WRITE(*,*)  ' R1 IS DETERMINED  ',R1
      P = AMIN1(P1,PMAX)
      P = AMAX1(P, PMIN)
      R = AMIN1(R1,RMAX)
      R = AMAX1(R, RMIN)
C
C COMPUTE NEW USTAR, THSTAR AND AMOL
C
      ALNETAO = ALOG(ETAO)
      AU    = FAU(AMU,ETAO,P,XA,XO)
      AT    = FAT(AMU,ETAO,AKAPA,R)
C
      TERMO = ETAF*ALOG(ETAF/ETAO) + ETAO - ETAF
C
      TERM1 = (ETAF-ETAO)*FPSIM(AMU*ETAO,XO)
      TERM2 = -(FGAMU(AMU*ETAF,XF)-FGAMU(AMU*ETAO,XO))/AMU
      TERM3 = AU*(ETAF-ETAO)**(P+1.)/(P+1.)
C
      TERM4 = (ETAF-ETAO)*FPSIH(AMU*ETAO,YO)
      TERM5 = -(FGAMH(AMU*ETAF,YF)-FGAMH(AMU*ETAO,YO))/AMU
      TERM6 = AT*(ETAF-ETAO)**(R+1.)/(R+1.)
C
      USNEW = VKAR*SPEED*(ETAF-ETAO)/(TERMO+TERM1+TERM2+TERM3)
      THSNEW= VKAR*(THM-THREF)*(ETAF-ETAO)/PRO
      THSNEW= THSNEW/(TERMO+TERM4+TERM5+TERM6)
C
      ERR0  = ABS( (AMOL-AMOLNEW)/AMOL )
      ERR1  = ABS( (USTAR-USNEW)/USTAR    )
      IF(ABS(THSTAR).GT.SMALL*0.01) THEN
       AMOLNEW = THREF*USNEW*USNEW/(VKAR*G*THSNEW)
       ERR2  = ABS( (THSTAR-THSNEW)/THSTAR )
      ELSE
       AMOLNEW = SIGN(1.,(THM-THREF))*BIG
       ERR2 = 0.0
      END IF
      AMOLNEW = SIGN(1.,AMOLNEW)*AMIN1(ABS(AMOLNEW),BIG)
      THSTAR = THSNEW
      USTAR  = USNEW
      AMOL   = AMOLNEW
C* LIMIT USTAR*THSTAR WITH HFXMAX
      IF( AMU .LT. -2.0) THEN
         HRATIO = - USTAR*THSTAR/HFXMAX
C* ! UNSTABLE AND TOO HIGH HFLUX
C* LIMIT SFC-AIR TEMPERATURE DIFFERENCE WHEN TOO LARGE HEAT FLUX
C* BUT KEEP THE AMOL THE SAME
       IF (HRATIO.GT.1.0) THEN
         USTAR = USTAR/HRATIO**.333333
	 THSTAR = THSTAR/HRATIO**.666666
C         THREF = THM - PRO*THSTAR*(TERMO+TERM4+TERM5+TERM6)
C     &   /(VKAR*(ETAF-ETAO))
C         AMOL = THREF*USTAR*USTAR/(VKAR*G*THSTAR)
C       write(6,*) ' new routine test '
       END IF
      END IF
C
      WSTAR  = USTAR*(HPBL/(VKAR*ABS(AMOL)) )**0.333333
      WSTAR  = CVMGP(0.0, WSTAR, AMOL)
C      IF( AMAX1(ERR1,ERR2).LT.TOL ) THEN
      IF (ERR0.LT.TOL)    GOTO 2
C
      GOTO 1
C% END OF DIABATIC CASES
      END IF
C
    2 CONTINUE
      RETURN
C
  999 CONTINUE
C
C SINCE CONVERSION FOR PBL SIMILARITY IS FAILED, USE PREDICTION BY
C SURFACE LAYER SIMILARITY
      THREF = THREF1
      AMOL  = AMOL1
      USTAR = USTAR1
      THSTAR= THSTAR1
      ITER = 0
C
      RETURN
      END


      function cvmgp(r1,r2,t)
      cvmgp=r1
      if(t.ge.0.0) return
      cvmgp=r2
      return
      end
c
      function cvmgz(r1,r2,t)
      cvmgz=r1
      if(t.eq.0.0) return
      cvmgz=r2
      return
      end


      subroutine divergence
     + (io_msg, maxx, maxy, maxz, nx, ny, nz, dx,dy, orog, grid,
     +  u,v,w,lobrien,divlimit)
c------------------------------------------------------------------------------
c C. Silibello, 30.10.98
c
c LAST REV.:
c
c PURPOSE: ITERATIVE SCHEME TO MINIMIZE DIVERGENCE (adapted from original CALMET model)
c
c INPUT:
c   io_msg                      I       message unit
c   maxx, maxy, maxz		I	arrays dimensions
c   nx, ny, nz			I	actual # of points in each direction
c   dx,dy                       R       deltax, deltay [m]
c   grid(maxx,maxy,maxz)        R       3-D grid [m above terrain]
c   U,V,W(maxx,maxy,maxz)	R	gridded U,V,W wind field components
c   lobrien                     L       O'Brien flag
c   divlimit                    R       maximum divergence allowed
c
c OUTPUT:
c   U,V,W(maxx,maxy,maxz)	R	U,V,W wind field components
c
c CALLS: DIV_CEL
c------------------------------------------------------------------------------
      implicit none

      integer io_msg
      integer maxx, maxy, maxz
      integer nx, ny, nz
      real dx, dy, grid(maxx,maxy,maxz)
      real U(maxx,maxy,maxz)
      real V(maxx,maxy,maxz)
      real W(maxx,maxy,maxz)
      logical lobrien
      real divlimit

      integer i, j, k
      real dym, dxm(maxy), dz(maxx,maxy,maxz),
     +     dxi, dyi, dzi(maxx,maxy,maxz), cellzb(maxz)
      integer itnum
      real DIV(maxx,maxy,maxz)
      real OROG(maxx,maxy)
      real divmax, ofac

c
c     calculate grid spacing of cell centers in km
c
      do j=1,ny
        dxm(j)=dx / 1000.
      enddo
      dym=dy / 1000.
      dyi=1.0/(2.0*dy)
      dxi=1.0/(2.0*dx)

      do i=1,nx
        do j = 1,ny
          dz(i,j,1) = grid(i,j,1)
	  dzi(i,j,1)= 1.0 / dz(i,j,1)
          do k=2,nz-1
            dz (i,j,k) = grid(i,j,k+1)-grid(i,j,k)
            dzi(i,j,k) = 1.0 / dz(i,j,k)
          enddo
          dz(i,j,nz)= dz(i,j,nz-1)
          dzi(i,j,nz) = 1.0 / dz(i,j,nz)
        enddo
      enddo

c     write (*,*) maxx, maxy, maxz, nx, ny, nz, dx,dy,
c    +  dxm(10), dym, lobrien,divlimit
c     do k = 1,nz
c     write (*,*) k,grid(10,10,k),dz (10,10,k),
c    +   u (10,10,k),v (10,10,k),w (10,10,k)
c     enddo

      do k = 1,nz
        divmax=-1.0e+9
        call div_cell
     + (maxx, maxy, maxz, nx, ny, k, dxi, dyi, dzi,
     +  orog, u, v, w, div, divmax)
      enddo

      do i=1,nx
        do j = 1,ny
c         w(i,j,1) = -dz(i,j,1)*div(i,j,1)
          w(i,j,1) = 0.
          do k=2,nz
            w(i,j,k)=-dz(i,j,k)*div(i,j,k)+w(i,j,k-1)
          enddo
        enddo
      enddo
c
      if (lobrien) then
c
        do i = 1,nx
          do j = 1,ny
            do k=1,nz
              cellzb(k)= grid(i,j,k)-grid(i,j,1)
            enddo
            ofac = w(i,j,nz)/cellzb(nz)
            do k =1,nz
              w(i,j,k)=w(i,j,k)-cellzb(k) * ofac
            enddo
          enddo
        enddo
c
c     minimize divergence
c
C        write(*,*)'entering min_div...'
        itnum = 50
        call min_div
     + (io_msg, maxx, maxy, maxz, nx, ny, nz,
     +  dym, dxm, dxi, dyi, dzi,
     +  orog, u, v, w, itnum, divlimit, div)
c
        write(*,*)'finished adjustment of vertical velocities...'
c
c
c    recalculate vertical velocities
c
        do k = 1,nz
          do i = 1,nx
            do j = 1,ny
              w(i,j,k) = 0.0
            enddo
          enddo
          call div_cell
     + (maxx, maxy, maxz, nx, ny, k, dxi, dyi, dzi,
     +  orog, u, v, w, div, divmax)
        enddo
        do i=1,nx
          do j = 1,ny
            w(i,j,1) = -dz(i,j,1)*div(i,j,1)
            do k=2,nz
              w(i,j,k)=-dz(i,j,k)*div(i,j,k)+w(i,j,k-1)
            enddo
          enddo
        enddo
      endif

      return
      end


      subroutine div_cell
     + (maxx, maxy, maxz, nx, ny, k, dxi, dyi, dzi,
     +  orog, u, v, w, div, divmax)
c------------------------------------------------------------------------------
c C. Silibello, 18.01.99
c
c LAST REV.:
c
c PURPOSE: COMPUTES 3-D DIVERGENCE IN EACH GRID CELL (adapted from original CALMET model)
c
c INPUT:
c   maxx, maxy, maxz		I	arrays dimensions
c   nx, ny     			I	actual # of points in each direction
c   k          			I	vertical level
c   dxi,dyi,dzi                 R       inverse of deltax, deltay and deltaz [m]
c   OROG (maxx,maxy)		R	orography
c   U,V,W(maxx,maxy,maxz)	R	gridded U,V,W wind field components
c
c OUTPUT:
c   DIV(maxx,maxy,maxz)		R	divergence
c   divmax                      R       computed maximum divergence
c
c CALLS: none
c------------------------------------------------------------------------------
      implicit none

      integer maxx, maxy, maxz
      integer nx, ny, k
      real dym, dxm(maxy), dxi, dyi, dzi(maxx,maxy,maxz)
      real U(maxx,maxy,maxz), V(maxx,maxy,maxz), W(maxx,maxy,maxz)
      real DIV(maxx,maxy,maxz), divmax
      real OROG(maxx,maxy)

      integer i,j
      real wmh, wph, umh, uph, vmh, vph,
     +     dudx, dvdy, dwdz, divabs
      real corr, ups, ums, vps, vms, huph, humh, hvph, hvmh


c
c  computes the 3-d divergence in each grid cell
c
      do j=1,ny
        do i=1,nx
          div(i,j,k)=0.

          wmh = 0.
          if (k .gt. 1) wmh = w(i,j,k-1)
          wph = w(i,j,k)

          umh=u(1,j,k)
          if(i.gt.1) umh=u(i-1,j,k)
          uph=u(nx,j,k)
          if(i.lt.nx) uph=u(i+1,j,k)

          vmh=v(i,1,k)
          if(j.gt.1) vmh=v(i,j-1,k)
          vph=v(i,ny,k)
          if(j.lt.ny) vph=v(i,j+1,k)

c
c orography correction
c
          ums = 0.
          if (k .gt. 1) ums = u(i,j,k-1)
          ups = u(i,j,k)
          vms = 0.
          if (k .gt. 1) vms = v(i,j,k-1)
          vps = v(i,j,k)

          humh=orog(1,j)
          if(i.gt.1) humh=orog(i-1,j)
          huph=orog(nx,j)
          if(i.lt.nx) huph=orog(i+1,j)

          hvmh=orog(i,1)
          if(j.gt.1) hvmh=orog(i,j-1)
          hvph=orog(i,ny)
          if(j.lt.ny) hvph=orog(i,j+1)
c
c     divergence is evaluated using central differences
c
          corr = (ups-ums)*dzi(i,j,k) * (- (huph - humh) * dxi)
          dudx = (uph - umh) * dxi + corr
          corr = (vps-vms)*dzi(i,j,k) * (- (hvph - hvmh) * dyi)
          dvdy = (vph - vmh) * dyi + corr
          dwdz = (wph-wmh)*dzi(i,j,k)
          div(i,j,k) = dwdz + dudx + dvdy
          divabs=abs(div(i,j,k))
          divmax=amax1(divmax,divabs)
        enddo
      enddo

      return
      end


      subroutine min_div
     + (io_msg, maxx, maxy, maxz, nx, ny, nz,
     +  dym, dxm, dxi, dyi, dzi,
     +  orog, u, v, w, niter, divlim, div)
c------------------------------------------------------------------------------
c C. Silibello, 18.01.99
c
c LAST REV.:
c
c PURPOSE: ITERATIVE SCHEME TO MINIMIZE DIVERGENCE (adapted from original CALMET model)
c
c INPUT:
c   io_msg                      I       message unit
c   maxx, maxy, maxz		I	arrays dimensions
c   nx, ny, nz     		I	actual # of points in each direction
c   dym, dxm			R	grid spacing of cell centers in km
c   dxi,dyi,dzi                 R       inverse of deltax, deltay and deltaz [m]
c   OROG (maxx,maxy)		R	orography
c   U,V,W(maxx,maxy,maxz)	R	gridded U,V,W wind field components
c   niter                       I       maximum number of iterations
c   divlim                      R       maximum divergence
c   DIV(maxx,maxy,maxz)		R	divergence to be minimized
c
c OUTPUT:
c   DIV(maxx,maxy,maxz)		R	minimized divergence
c
c CALLS: div_cell
c------------------------------------------------------------------------------
      implicit none

      integer io_msg
      integer maxx, maxy, maxz
      integer nx, ny, nz
      real dym, dxm(maxy), dxi, dyi, dzi(maxx,maxy,maxz)
      real OROG(maxx,maxy)
      real U(maxx,maxy,maxz), V(maxx,maxy,maxz), W(maxx,maxy,maxz)
      integer niter
      real divlim, DIV(maxx,maxy,maxz)


      integer i, j, k, ii, jj, iter, idir
      integer ip1, im1, jp1, jm1
      real uim1, uip1, vjm1, vjp1, ut, vt
      real alpha1, alpha2, alpha3, alpha4, al1234
      real divmax

      write(io_msg,2809)

      do k=1,nz
        iter=0
c
c     compute divergence
c
        divmax=-1.0e+09
        call div_cell
     + (maxx, maxy, maxz, nx, ny, k, dxi, dyi, dzi,
     +  orog, u, v, w, div, divmax)
        if (divmax .le. divlim) cycle
c
c     adjust horizontal wind fields
c
        do iter=1,niter
          do idir=1,4
            do jj=1,ny
              do ii=1,nx
                if (idir .eq. 1) then
                  i=ii
                  j=jj
                elseif (idir .eq. 2) then
                  i=nx-ii+1
                  j=jj
                elseif (idir .eq. 3) then
                  i=ii
                  j=ny-jj+1
                elseif (idir .eq. 4) then
                  i=nx-ii+1
                  j=ny-jj+1
                endif
                if (div(i,j,k) .eq. 0.) cycle
                ip1=i+1
                im1=i-1
                jp1=j+1
                jm1=j-1
                uim1=u(1,j,k)
                if (i .gt. 1) uim1=u(im1,j,k)
                uip1=u(nx,j,k)
                if (i .lt. nx) uip1=u(ip1,j,k)
                vjm1=v(i,1,k)
                if (j .gt. 1) vjm1=v(i,jm1,k)
                vjp1=v(i,ny,k)
                if (j .lt. ny) vjp1=v(i,jp1,k)
                alpha1=0.5
                alpha2=0.5
                alpha3=0.5
                alpha4=0.5
                al1234=alpha1+alpha2+alpha3+alpha4
                if (al1234 .lt. 1.e-6) cycle
                ut=-2.*(div(i,j,k)*dxm(j))/al1234
                vt=-2.*(div(i,j,k)*dym)/al1234
                uip1=uip1+alpha1*ut
                uim1=uim1-alpha2*ut
                vjp1=vjp1+alpha3*vt
                vjm1=vjm1-alpha4*vt
                if (i .gt. 1) u(im1,j,k)=uim1
                if (i .lt. nx) u(ip1,j,k)=uip1
                if (j .gt. 1) v(i,jm1,k)=vjm1
                if (j .lt. ny) v(i,jp1,k)=vjp1
              enddo			! loop over ii
            enddo			! loop over jj
            divmax=-1.0e+09
            call div_cell
     + (maxx, maxy, maxz, nx, ny, k, dxi, dyi, dzi,
     +  orog, u, v, w, div, divmax)
          enddo				! loop over idir
c
c     convergence test for divergence magnitude
c
c ...if minim is not converging to a solution, then abort program...

          if(divmax.gt.10.) then
            write(*,*)'error in minim: divmax= ',divmax
            write(*,*)'level=                  ',k
            write(*,*)'iter=                   ',iter
            write(io_msg,*)'error in minim--failed to converge'
            write(io_msg,2810)k,iter,divmax
            write(io_msg,*)'div= ',(div(18,j,k),j=1,ny)
          endif
          if (divmax .le. divlim) exit
        enddo				! loop over iter
        write(io_msg,2810) k,iter,divmax
      enddo				! loop over k
c
 2809 format(//,5x,'summary of divergence minimization'
     1,//,'  level   iterations   maximum divergence (/sec)')
 2810 format(3x,i2,i11,12x,e10.3)
      return
      end

      subroutine KZ_INTEG
     +  (maxz, nz, zlev, ustar, lmonin, zpbl, wstar, kz)
c------------------------------------------------------------------------------
c G. Calori, 1.11.99
c
c LAST REV.:
c
c PURPOSE: Compute a Kz profile according to Byun & Dennis
c          (Businger et al. inside the SL; Brost & Wyngaard inside the PBL;
c           constant inside the FT)
c
c INPUT:
c maxz       	I  max # of vertical levels
c nz		I  actual # of levels
c zlev(maxz) 	R  vertical levels (faces heights) (m above topography)
c ustar      	R  friction velocity at surface (m/s)
c lmonin     	R  Monin-Obukhov lenght at surface (m)
c zpbl       	R  PBL height (m)
c wstar      	R  convective velocity scale at surface (m/s)
c
c OUTPUT:
c kz(nz)     	R  vertical diffusivities (m2/sec)
c
c CALLS: -
c-----------------------------------------------------------------------------
      implicit none

      real vk, pr, vkpr, b0, b1
      parameter (vk = 0.4, pr = 0.74, vkpr = vk / pr)
      parameter (b0 = 4.7, b1 = b0 / pr)
      real lmin, lmax
      parameter (lmin = 1., lmax = 999.)
      real freek, kzmin
      parameter (freek = 1.)	! free troposphere kz
      parameter (kzmin = 1.)

      integer maxz, nz
      real zlev(maxz), lmonin, ustar, zpbl, wstar, kz(maxz)

      integer k
      real obuk, zsl, z1, z2, dz, pre, pre1,
     +     alpha, beta, beta1, a, a2, r1, r2
      real stsl, unsl, stpbl, unpbl

      stsl(z1,z2,beta1) =
     +  (z2 - z1 - alog((beta1 * z2 + 1.) / (beta1 * z1 + 1.)) / beta1)
     +  / beta1

      unsl(z1,z2,alpha) =
     +  ((3. * alpha * z2 - 2.) * sqrt(1. + alpha * z2)**3 -
     +   (3. * alpha * z1 - 2.) * sqrt(1. + alpha * z1)**3
     +  ) / (7.5 * alpha*alpha)

      stpbl(r1,r2,a,a2) =
     +  0.2 * (r1**5 - r2**5) + (a2 - 1.) *
     +  ( (r1**3 - r2**3) / 3. +
     +    a2 * (r1 - r2 - 0.5 * a *
     +          alog(((a-r2) * (a+r1)) / ((a+r2) * (a-r1))) ) )

      unpbl(z1,z2,wstar,zpbl) =
     +  vk * wstar *
     +  (z2**2 * (0.5 - z2 / (3. * zpbl)) -
     +   z1**2 * (0.5 - z1 / (3. * zpbl)))



      obuk = lmonin
      if (obuk .gt. 0.) then
        obuk = amax1(amin1(obuk, lmax), lmin)
      else
        obuk = amin1(amax1(obuk, -lmax), -lmin)
      endif

      zsl = amin1(0.1 * zpbl, 50.)
      pre = vkpr * ustar

      if (obuk .lt. 0.) then        ! unstable

        do k = 1,nz-1

          z1 = zlev(k)
          z2 = zlev(k+1)
          dz = z2 - z1

          if (z1 .lt. zsl) then
            alpha = -9. / obuk
            if (z2 .lt. zsl) then			! inside SL
              kz(k) = pre * unsl(z1,z2,alpha) / dz
            else if (z2.lt.zpbl.and.z2.ge.zsl) then	! across SLtop
              kz(k) = (pre * unsl(z1,zsl,alpha) +
     +                 unpbl(zsl,z2,wstar,zpbl)) / dz
            else					! across BLtop
              kz(k) = (pre * unsl(z1,zsl,alpha) +
     +                 unpbl(zsl,zpbl,wstar,zpbl) +
     +                 (z2 - zpbl) * freek) / dz
            endif
          else if(z1.lt.zpbl.and.z1.ge.zsl) then
            if (z2 .lt. zpbl) then			! inside BL
              kz(k) = unpbl(z1,z2,wstar,zpbl) / dz
            else					! across BLtop
              kz(k) = (unpbl(z1,zpbl,wstar,zpbl) +
     +                 (z2 - zpbl) * freek) / dz
            endif
	  else						! FT
	    kz(k) = freek
          endif

        enddo

      else                            ! stable and neutral

        do k = 1,nz-1

          z1 = zlev(k)
          z2 = zlev(k+1)
          dz = z2 - z1

          beta1 = b1 / obuk
          beta = zpbl * beta1
          a2 = (1 + beta) / beta
          a = sqrt(a2)
          pre1 = 2. * pre * zpbl * zpbl / beta

          if (z1 .lt. zsl) then
            if (z2 .lt. zsl) then			! inside SL
              kz(k) = pre * stsl(z1,z2,beta1) / dz
            else if (z2 .lt. zpbl) then			! across SLtop
              r1 = amax1(1.e-5, sqrt(1. - zsl/zpbl))
              r2 = sqrt(1. - z2/zpbl)
              kz(k) = (pre * stsl(z1,zsl,beta1) +
     +                 pre1 * stpbl(r1,r2,a,a2)) / dz
            else					! across BLtop
              r1 = amax1(1.e-5, sqrt(1. - zsl/zpbl))
              r2 = 0.
              kz(k) = (pre * stsl(z1,zsl,beta1) +
     +                 pre1 * stpbl(r1,r2,a,a2) +
     +                 (z2 - zpbl) * freek) / dz
            endif
          else if(z1 .lt. zpbl) then
            if (z2 .lt. zpbl) then			! inside BL
              r1 = sqrt(1. - z1/zpbl)
              r2 = sqrt(1. - z2/zpbl)
              kz(k) = pre1 * stpbl(r1,r2,a,a2) / dz
            else					! across BLtop
              r1 = sqrt(1. - z1/zpbl)
              r2 = 0.
              kz(k) = (pre1 * stpbl(r1,r2,a,a2) +
     +                 (z2 - zpbl) * freek) / dz
            endif
	  else						! FT
	    kz(k) = freek
          endif

        enddo

      endif

      kz(nz) = kz(nz-1)

      do k = 1,nz
        kz(k) = amax1(kz(k), kzmin)
      enddo

      return
      end

C*****
      subroutine kz_cloud(kmax,nz,kz,rh,cw,prate)
C-----------------------------------------------------------------
C     compute the KZ value due to cloud
C     Input:
C      KZ       original vertical diffusion coefficient (m2/s)
C      RH       relative humdity (%)
C      CW       Cloud Water Content (Kg/Kg)
C      kmax     vertical dimension of all these arrays
C      nz       actual # of layers
C      prate    precipitation rate in mm/hr
C
C     Output:
C      kz      adjusted vertical diffusion coefficient (m2/s)
C
C   Author: Youhua Tang
C-----------------------------------------------------------------
      real kz(kmax),rh(kmax),cw(kmax)
      parameter ( RHC=98.,     ! cloud relative humdity
     1  CWMIN=5e-4, CWMAX=2e-3,  ! define if there is cloud
     2  CKMIN=10.,CKMAX=100., PAI=3.14159265)

      CK(cwater)=ckmin+sin(cwater/cwmax*pai/2)*ckmax  ! kz profile

      cover_max=0.
      do k=1,nz
       if(cover_max.lt.cw(k)) then
        cover_max=cw(k)
	kcover=k
       endif
      enddo

      if(cover_max.ge.cwmin) then  ! cloud exist
        do k=kcover,nz
	 if(cw(k).lt.cwmin) goto 20
	enddo
 20     ktop=k
       do k=kcover,1,-1
	 if(cw(k).lt.cwmin) goto 30
	enddo
 30    kbottom=k

c       if(kbottom.gt.2) then      !set kz below cloud
c        do k=2,kbottom-1
c	 kz(k)=amax1(kz(k),ckmin)
c	enddo
c       endif

       if(prate.gt.1) kbottom=1    ! cumulus cloud from surface to cloud top

       do k=kbottom,ktop           ! set kz within cloud
        cw(k)=amin1(cw(k),cwmax)
c        kz(k)=amax1(kz(k),ckmin+sin(cw(k)/cwmax*pai/2)*ckmax)
        kz(k)=amax1(kz(k),ck(cw(k)))
       enddo

      else
       kcover=1
       do k=1,nz
        if(rh(k).ge.rhc) then   ! thin cloud
	 kcover=k
	endif
       enddo

       if(kcover.ne.1) then
        do k=1,kcover
         kz(k)=amax1(kz(k),ckmin)
        enddo
       endif
      endif
      end


      subroutine wind_sub(ix,iy,iz,u2,v2,w,dx,dy,z,air)
c----------------------------------------------------------------
      dimension u2(ix,iy,iz),v2(ix,iy,iz),w(ix,iy,iz),dz(100),
     1  air(ix,iy,iz)
      dimension z(ix,iy,iz)
      real u(ix,iy,iz),v(ix,iy,iz)
      write(6,*) 'insde wind_sub',ix,iy,iz,htop,dx,dy

      iz1=iz-1
      do 20 i=1,ix
      do 20 j=1,iy
      do 20 k=1,iz1
c	u(i,j,k)=0.5*(u2(i,j,k)+u2(i,j,k+1))*air(i,j,k) ! encode air density
c	v(i,j,k)=0.5*(v2(i,j,k)+v2(i,j,k+1))*air(i,j,k)
	 u(i,j,k)=0.5*(u2(i,j,k)+u2(i,j,k+1))
	 v(i,j,k)=0.5*(v2(i,j,k)+v2(i,j,k+1))
 20    continue
cc

      do 30 i=1,ix
      do 30 j=1,iy
        do k=1,iz1
         dz(k)=z(i,j,k+1)-z(i,j,k)
        enddo
	 w(i,j,1)=0.e0
      do 30 k=1,iz1
	 if (i.eq.1) then
	    a=(-3.*u(i,j,k)+4.*u(i+1,j,k)-u(i+2,j,k))/(2.*dx)
            if (j.eq.1) then
	       b=(-3.*v(i,j,k)+4.*v(i,j+1,k)-v(i,j+2,k))/(2.*dy)
	    else if(j.lt.iy) then
               b=(v(i,j+1,k)-v(i,j-1,k))/(2.*dy)
            else
               b=(3.*v(i,j,k)-4.*v(i,j-1,k)+v(i,j-2,k))/(2.*dy)
            endif
         else if(i.lt.ix) then
            a=(u(i+1,j,k)-u(i-1,j,k))/(2.*dx)
	    if(j.eq.1) then
	       b=(-3.*v(i,j,k)+4.*v(i,j+1,k)+v(i,j+2,k))/(2.*dy)
            else if(j.lt.iy) then
               b=(v(i,j+1,k)-v(i,j-1,k))/(2.*dy)
	    else
               b=(3.*v(i,j,k)-4.*v(i,j-1,k)+v(i,j-2,k))/(2.*dy)
            endif
	 else
            a=(3.*u(i,j,k)-4.*u(i-1,j,k)+u(i-2,j,k))/(2.*dx)
	    if(j.eq.1) then
	       b=(-3.*v(i,j,k)+4.*v(i,j+1,k)+v(i,j+2,k))/(2.*dy)
	    else if(j.lt.iy) then
               b=(v(i,j+1,k)-v(i,j-1,k))/(2.*dy)
	    else
               b=(3.*v(i,j,k)-4.*v(i,j-1,k)+v(i,j-2,k))/(2.*dy)
            endif
         endif
c         w(i,j,k+1)=w(i,j,k)+dz(k)*(-(a+b))/air(i,j,k+1) !decode from air density
         w(i,j,k+1)=w(i,j,k)+dz(k)*(-(a+b))
30    continue

c
      return
      end

      subroutine calc_seasalt(nx,ny,u10,v10,rh,xland,ss1,ss2,ss3,ss4)
c------------------------------------------------------------------------------
c by. Youhua Tang adopt from Gong S. L. (2003), Global Biogeochemical Cycles, 17(4), 1097, doi:10.1029/2003GB002079
c
c PURPOSE: To Caluclate the seasalt emissions
c
c INPUT:
c nx, ny 			I  actual dimension of 3D grid
C U10, V10                      10 meter wind
c
c OUTPUT:
c ss1			     Sea Salt emission (molecular/cm^2/s) (diameter 0.1-0.3 um)
c ss2                        sea Salt emission (molecular/cm^2/s) (diameter 0.3-1.0 um)
c ss3                        sea Salt emission (molecular/cm^2/s) (diameter 1.0-2.5 um)
c ss4                        sea Salt emission (molecular/cm^2/s) (diameter 2.5-10. um)
c------------------------------------------------------------------------------

       integer nx,ny,nz,i,j
       real u10(nx,ny),v10(nx,ny),rh(nx,ny)
       real ss1(nx,ny),ss2(nx,ny),ss3(nx,ny),ss4(nx,ny),xland(nx,ny)

       parameter (pai=3.14159, roughness=0.001,       ! roughtness over water
     1  sswt=58.44, afconst=6.023e23)        ! NaCl molecular weight

c      sea salt 20.0032 part/cm3, 39.5ug/m3
C      number density: acc 20.(99.98%), coa 3.2E-3(0.02%)
C      mass : acc 38.6 (97.72%), coa 0.9 (2.28%)
c      Nacl dry mass density is 2.16g/cm3
       parameter(c25=1.93e-6 ,c10=1.544e-5, c75=9.8816e-4,
     1   c100=7.585672e-2,c1k=1.93, theta=30., density=2.16)
          ! convert factor ug/particle  for each size, assuming same mass density
	  ! average size of each bin: 2.5um, 5um, 20um, 85um, 250um


c       do 10 i=1,nx
c        do 10 j=1,ny
c         ss1(i,j) = 0.
c	 ss2(i,j) = 0.
c         ss3(i,j) = 0.
c	 ss4(i,j) = 0.
c10     continue


       do i=1,nx
       do j=1,ny

	 if(abs(xland(i,j)-16.).lt.1e-20.and.rh(i,j).ge.70.) then

         wind10=sqrt(u10(i,j)**2+v10(i,j)**2)

	 ss1(i,j)=0.                  ! unit #particles/m2/s
	 dr=(0.3-0.1)/2/50
         do nr=1,50                     ! integrate the emission of particles with diamter 0.1-0.3um
	  r=dr*nr+0.05                  ! sea salt radii in um
	  a = 4.7*(1+ theta*r)**(-0.017*r**(-1.44))
	  b = (0.433-log(r))/0.433
	  term1=1.373*(wind10**3.41)*(r**(-a))
	  term2 = 1+.057*(r**3.45)
	  term3=10**(1.607*exp(-(b**2)))
	  volume=4./3.*pai*(r/1.e4)**3                 ! volume in cm3
	  ss1(i,j)=ss1(i,j)+term1*term2*term3*dr
     1     *volume*density/sswt*afconst/10000.           ! convert to molecules/cm2/s	  
	 enddo

	 ss2(i,j)=0.                  ! unit #particles/m2/s
	 dr=(1.0-0.3)/2/50
         do nr=1,50                ! integrate the emission of particles with diamter 0.3-1.0um
	  r=dr*nr+0.15                  ! sea salt radii in um
	  a = 4.7*(1+ theta*r)**(-0.017*r**(-1.44))
	  b = (0.433-log(r))/0.433
	  term1=1.373*(wind10**3.41)*(r**(-a))
	  term2 = 1+.057*(r**3.45)
	  term3=10**(1.607*exp(-(b**2)))
	  volume=4./3.*pai*(r/1.e4)**3                 ! volume in cm3
	  ss2(i,j)=ss2(i,j)+term1*term2*term3*dr
     1     *volume*density/sswt*afconst/10000.           ! convert to molecules/cm2/s
	 enddo

	 ss3(i,j)=0.                  ! unit #particles/m2/s
	 dr=(2.5-1.0)/2/50
         do nr=1,50                ! integrate the emission of particles with diamter 1.0-2.5um
	  r=dr*nr+0.5                  ! sea salt radii in um
	  a = 4.7*(1+ theta*r)**(-0.017*r**(-1.44))
	  b = (0.433-log(r))/0.433
	  term1=1.373*(wind10**3.41)*(r**(-a))
	  term2 = 1+.057*(r**3.45)
	  term3=10**(1.607*exp(-(b**2)))
	  volume=4./3.*pai*(r/1.e4)**3                 ! volume in cm3
	  ss3(i,j)=ss3(i,j)+term1*term2*term3*dr
     1     *volume*density/sswt*afconst/10000.           ! convert to molecules/cm2/s
	 enddo

	 ss4(i,j)=0.                  ! unit #particles/m2/s
	 dr=(10.-2.5)/2/50
         do nr=1,50                   ! integrate the emission of particles with diamter 2.5-10. um
	  r=dr*nr+1.25                  ! sea salt radii in um
	  a = 4.7*(1+ theta*r)**(-0.017*r**(-1.44))
	  b = (0.433-log(r))/0.433
	  term1=1.373*(wind10**3.41)*(r**(-a))
	  term2 = 1+.057*(r**3.45)
	  term3=10**(1.607*exp(-(b**2)))
	  volume=4./3.*pai*(r/1.e4)**3                 ! volume in cm3
	  ss4(i,j)=ss4(i,j)+term1*term2*term3*dr
     1     *volume*density/sswt*afconst/10000.           ! convert to molecules/cm2/s
	 enddo

	 endif
        enddo
       enddo

       return
      end

      subroutine aq_locate(iunit,char,iflag)
c***********************************************************************
      character*(*) char
      character*80 dum1
      nchar=len(char)
      iflag=0
      do iter=1,10000
      read(iunit,'(a)',end=98) dum1(1:nchar)
c      print*,'dum1= ',dum1(1:nchar)
      if(dum1(1:nchar).eq.char) return
      enddo
98    iflag=1
c      print*,'dum1= ',dum1(1:nchar)
      return
      end


****subroutine calculate PV
*      DMF   dot point map factor
*      XMF   cross point map factor
*       F    coriolis term
      SUBROUTINE PVS(U,V,T,DMF,XMF,F,P,DS,I1,J1,K1,ioff,
     1               joff,iout,jout,PV)
C
      REAL U(I1,J1,K1), V(I1,J1,K1), T(I1,J1,K1),
     *          DMF(I1,J1), XMF(I1,J1), F(I1,J1),
     *          P(iout,jout,k1),
     *          PV(iout,jout,k1)
      REAL VOR (iout,jout),DTDS(iout,jout),
     *          DUDS(iout,jout), DVDS(iout,jout),
     *          DTDX(iout,jout), DTDY(iout,jout)

C
C     ... EVERYTHING IS BY SLABS
C
      G=9.81
C      SCALE=-1.E6
C
      DO 5000 K=1,K1
C
C        ... COMPUTE VERTICAL DERIVATIVES: DU/DP, DV/DP, DTHETA/DP
C
         IF(K.EQ.1) THEN
            K0=K
            K2=K+1
         ELSE IF(K.EQ.K1) THEN
            K0=K-1
            K2=K
         ELSE
            K0=K-1
            K2=K+1
         ENDIF
         DO 100 J=1,JOUT
	 JJ=J+JOFF
         DO 100 I=1,IOUT
	 II=I+IOFF
            DUDS(I,J)=(U(II+1,JJ+1,K0)-U(II+1,JJ+1,K2)  +
     *                 U(II+1,JJ ,K0)-U(II+1,JJ ,K2)  +
     *                 U(II ,JJ+1,K0)-U(II ,JJ+1,K2)  +
     *                 U(II ,JJ ,K0)-U(II ,JJ ,K2)) *
     *                 0.25 / (P(I,J,K0)-P(I,J,K2)) * XMF(II,JJ)
            DVDS(I,J)=(V(II+1,JJ+1,K0)-V(II+1,JJ+1,K2)  +
     *                 V(II+1,JJ ,K0)-V(II+1,JJ ,K2)  +
     *                 V(II ,JJ+1,K0)-V(II ,JJ+1,K2)  +
     *                 V(II ,JJ ,K0)-V(II ,JJ  ,K2)) *
     *                 0.25 / (P(I,J,K0)-P(I,J,K2)) * XMF(II,JJ)
100      CONTINUE
         DO 110 J=1,JOUT
	 JJ=J+JOFF
         DO 110 I=1,IOUT
	  II=I+IOFF
            T1=T(II,JJ,K0)*(1e5/P(I,J,K0))**0.286
            T2=T(II,JJ,K2)*(1e5/P(I,J,K2))**0.286
            DTDS(I,J)=(T1-T2)/(P(I,J,K0)-P(I,J,K2))
110      CONTINUE
C
C        ... COMPUTE HORIZONTAL DERIVATIVES: DTHETA/DX, DTHETA/DY
C
         DS2R=1./(DS * 2.)
         DO 400 J=2,JOUT-1
	 JJ=J+JOFF
         DO 400 I=2,IOUT-1
	  II=I+IOFF
            T1=T(II+1,JJ ,K)*(1e5/P(I,J,K))**0.286
            T2=T(II-1,JJ ,K)*(1e5/P(I,J,K))**0.286
            T3=T(II ,JJ+1,K)*(1e5/P(I,J,K))**0.286
            T4=T(II ,JJ-1,K)*(1e5/P(I,J,K))**0.286
            DTDY(I,J)=(T1 - T2)*DS2R
            DTDX(I,J)=(T3 - T4)*DS2R
400      CONTINUE
         CALL FILLIT(DTDX,IOUT,JOUT,1,IOUT,JOUT,2,IOUT-1,2,JOUT-1)
         CALL FILLIT(DTDY,IOUT,JOUT,1,IOUT,JOUT,2,IOUT-1,2,JOUT-1)
C
C        ... COMPUTE SLAB ABSOLUTE VORTICITY
C
         DO 500 J=1,JOUT
	 JJ=J+JOFF
         DO 500 I=1,IOUT
	  II=I+IOFF
            U1=U(II ,JJ ,K)/DMF(II ,JJ )
            U2=U(II+1,JJ ,K)/DMF(II+1,JJ )
            U3=U(II ,JJ+1,K)/DMF(II ,JJ+1)
            U4=U(II+1,JJ+1,K)/DMF(II+1,JJ+1)
            V1=V(II ,JJ ,K)/DMF(II ,JJ )
            V2=V(II+1,JJ ,K)/DMF(II+1,JJ )
            V3=V(II ,JJ+1,K)/DMF(II ,JJ+1)
            V4=V(II+1,JJ+1,K)/DMF(II+1,JJ+1)
            FX=(F(II,JJ)+F(II+1,JJ)+F(II,JJ+1)+F(II+1,JJ+1))*0.25
            VOR(I,J)=XMF(II,JJ)**2 * DS2R*((V4-V2+V3-V1)-(U2-U1+U4-U3))
     *               + FX
500      CONTINUE
C
C        ... GROUP TERMS
C
         DO 600 J=1,JOUT
         DO 600 I=1,IOUT
            PV(I,J,K)=-G *
     *       (  VOR (I,J) * DTDS(I,J) -
     *          DVDS(I,J) * DTDX(I,J) + DUDS(I,J) * DTDY(I,J)  )
         If(abs(PV(I,J,K)).gt.1) then
	   print*,'PV too high',PV(I,J,K),i,j,k
	   print*,'VOR, DTDS, DVDS, DTDX, DUDS, DTDY=',VOR(i,j),
     1	    DTDS(i,j), DVDS(i,j), DTDX(i,j), DUDS(i,j), DTDY(i,j)
           print*,'FX, XMF, DMF=',FX,XMF(i+ioff,j+joff),
     1	    DMF(i+ioff,j+joff)
           print*,'P(i,j,K0), P(i,j,K2)=',P(i,j,K0), P(i,j,K2)
	   print*,'T(ii,j,K0), T(ii,j,K2)=',T(i+ioff,j+joff,K0),
     1	    T(i+ioff,j+joff,K2)
           stop
	  endif
600      CONTINUE
5000  CONTINUE
C
      RETURN
      END
C
CCCCCCCCCCCCCCCCCCCCCCC
C
      SUBROUTINE FILLIT(F,IX,JX,KX,IMX,JMX,IFIRST,ILAST,JFIRST,JLAST)
      implicit none
      integer IX,JX,KX,IMX,JMX,IFIRST,ILAST,JFIRST,JLAST
C
C     SECTION  TOOLS
C     PURPOSE  FILL DATA OUT TO IMX,JMX FROM AN INTERIOR DOMAIN
C
      REAL F(IX,JX,KX)

      integer i,j,k
C
      DO 1000 K=1,KX
         DO 300 J=JFIRST,JLAST
            DO 100 I=1,IFIRST-1
               F(I,J,K)=F(IFIRST,J,K)
100         CONTINUE
            DO 200 I=ILAST+1,IMX
               F(I,J,K)=F(ILAST,J,K)
200         CONTINUE
300      CONTINUE
         DO 600 I=1,IMX
            DO 400 J=1,JFIRST-1
               F(I,J,K)=F(I,JFIRST,K)
400         CONTINUE
            DO 500 J=JLAST+1,JMX
               F(I,J,K)=F(I,JLAST,K)
500         CONTINUE
600      CONTINUE
1000  CONTINUE
      RETURN
      END


      subroutine handle_err(rmarker,nf_status)

      include "netcdf.inc"
      integer nf_status
      character*(*)        :: rmarker
      if (nf_status .ne. nf_noerr) then
         write(*,*)  'NetCDF error : ',rmarker
         write(*,*)  '  ',nf_strerror(nf_status)
         stop
      endif

      end subroutine handle_err

c-------------------------------------------------------------------------------

      subroutine get_gl_att_real_cdf( file, att_name, value, debug )

      implicit none
      include 'netcdf.inc'

      character (len=*), intent(in) :: file
      character (len=*), intent(in) :: att_name
      logical, intent(in ) :: debug
      real,  intent(out) :: value

      integer cdfid, rcode, id_time, length

      cdfid = ncopn(file, NCNOWRIT, rcode )
      length = max(1,index(file,' ')-1)
      if( rcode == 0) then
       if(debug) write(6,*) ' open netcdf file ',trim(file)
      else
       write(6,*) ' error openiing netcdf file ',trim(file)
       stop
      end if

      rcode = NF_GET_ATT_REAL(cdfid, nf_global, att_name, value )

      call ncclos(cdfid,rcode)

      if(debug) write(6,*) ' global attribute ',att_name,' is ',value

      end subroutine get_gl_att_real_cdf

c--------------------------------------------------------------------

       subroutine get_var_2d_real_cdf( file, var, data,
     1                              i1, i2, time, debug )

       implicit none

       include 'netcdf.inc'

       integer, intent(in)  ::  i1, i2, time
       character (len=*), intent(in) :: file
       logical, intent(in ) :: debug
       character (len=*), intent(in) :: var
       real, dimension(i1,i2), intent(out) :: data

       integer cdfid, rcode, id_data
       character (len=80) :: varnam, time1
       integer :: ndims, natts, idims(10), istart(10),iend(10),
     1  dimids(10)
       integer :: i, ivtype, length

       cdfid = ncopn(file, NCNOWRIT, rcode )
       length = max(1,index(file,' ')-1)
       if( rcode == 0) then
         if(debug) write(6,*) ' open netcdf file ',trim(file)
       else
         write(6,*) ' error openiing netcdf file ',trim(file)
         stop
       end if

       id_data = ncvid( cdfid, var, rcode )

       rcode = nf_inq_var( cdfid, id_data, varnam, ivtype, ndims,
     1   dimids, natts )
       if(debug) then
         write(6,*) ' number of dims for ',var,' ',ndims
       endif
       do i=1,ndims
         rcode = nf_inq_dimlen( cdfid, dimids(i), idims(i) )
         if(debug) write(6,*) ' dimension ',i,idims(i)
       enddo

!  check the dimensions

       if( (i1 /= idims(1)) .or.
     1  (i2 /= idims(2)) .or.
     2  (time > idims(3))     )  then

        write(6,*) ' error in 2d_var_real read, dimension problem'
        write(6,*) i1, idims(1)
        write(6,*) i2, idims(2)
        write(6,*) time, idims(4)
        write(6,*) ' error stop '
        stop

       end if

!  get the data

       istart(1) = 1
       iend(1) = i1
       istart(2) = 1
       iend(2) = i2
       istart(3) = time
       iend(3) = 1

       call ncvgt( cdfid,id_data,istart,iend,data,rcode)

       call ncclos(cdfid,rcode)

       end subroutine get_var_2d_real_cdf

!--------------------------------------------------------------------

       subroutine get_var_3d_real_cdf( file, var, data,
     1                             i1, i2, i3, time, debug )

       implicit none

       include 'netcdf.inc'

       integer, intent(in)  ::  i1, i2, i3, time
       character (len=*), intent(in) :: file
       logical, intent(in ) :: debug
       character (len=*), intent(in) :: var
       real, dimension(i1,i2,i3), intent(out) :: data

       integer cdfid, rcode, id_data
       character (len=80) :: varnam, time1
       integer :: ndims,natts,idims(10),istart(10),iend(10),dimids(10)
       integer :: i, ivtype, length

       cdfid = ncopn(file, NCNOWRIT, rcode )
       length = max(1,index(file,' ')-1)
       if( rcode == 0) then
       if(debug) write(6,*) ' open netcdf file ',trim(file)
        else
       write(6,*) ' error openiing netcdf file ',trim(file)
        stop
       end if

       id_data = ncvid( cdfid, var, rcode )

       rcode=nf_inq_var(cdfid,id_data,varnam,ivtype,ndims,dimids,natts)
       if(debug) then
        write(6,*) ' number of dims for ',var,' ',ndims
       endif
       do i=1,ndims
        rcode = nf_inq_dimlen( cdfid, dimids(i), idims(i) )
       if(debug) write(6,*) ' dimension ',i,idims(i)
       enddo

!  check the dimensions

       if( (i1 /= idims(1)) .or.
     1  (i2 /= idims(2)) .or.
     2   (i3 /= idims(3)) .or.
     3   (time > idims(4)) )  then

        write(6,*) ' error in 3d_var_real read, dimension problem '
        write(6,*) i1, idims(1)
        write(6,*) i2, idims(2)
        write(6,*) i3, idims(3)
        write(6,*) time, idims(4)
        write(6,*) ' error stop '
        stop

       end if

!  get the data

       istart(1) = 1
       iend(1) = i1
       istart(2) = 1
       iend(2) = i2
       istart(3) = 1
       iend(3) = i3
       istart(4) = time
       iend(4) = 1

       call ncvgt( cdfid,id_data,istart,iend,data,rcode)

       call ncclos(cdfid,rcode)

       end subroutine get_var_3d_real_cdf

!-------------------------------------------------------------------------------

       subroutine get_gl_att_int_cdf( file, att_name, value, debug )

       implicit none

       include 'netcdf.inc'

       character (len=*), intent(in) :: file
       character (len=*), intent(in) :: att_name
       logical, intent(in ) :: debug
       integer, intent(out) :: value

       integer cdfid, rcode, id_time, length

       cdfid = ncopn(file, NCNOWRIT, rcode )
       length = max(1,index(file,' ')-1)
       if( rcode == 0) then
         if(debug) write(6,*) ' open netcdf file ',trim(file)
       else
         write(6,*) ' error openiing netcdf file ',trim(file)
       stop
       end if

       rcode = NF_GET_ATT_INT(cdfid, nf_global, att_name, value )

       call ncclos(cdfid,rcode)

       if(debug) write(6,*) ' global attribute ',att_name,' is ',value

       end subroutine get_gl_att_int_cdf


      subroutine jclouds(iz,z,pbl,p,prate,t,vapor,cfrac,ccover,ktop)
C-----------------------------------------------------------------
C     compute the cloud coverage
C     Input:
C      iz       demension
C      Z        model grid height (m)
C      PBL      pbl height (m)
C      P        Pressure   (Pa)
C      PRATE    Precipitation rate (mm/hr)
C      T        Temperature (k)
c      vapor    water vapor (Kg/Kg)
Cccc   cw       Cloud Water Content (Kg/Kg)
C
C     Output:
c      cfrac    cloud fraction
C      ccover   cloud coverage
C      Ktop     index of clound top height
C
C   Author: Youhua Tang
C-----------------------------------------------------------------
      parameter(CWMIN=1e-5,   ! define if there is cloud
     2  lv0 = 2.501e6, prmin=0.1, nwet=4 )          ! define if there is precipitation
      parameter(rgasuniv=8.314510, mwair=28.9628,
     1  rdgas =1.0e3 *rgasuniv/mwair, cpd=7.0*rdgas/2.0,
     2  mwwat=18.0153, mvoma=mwwat/mwair)

      integer iz,ktop
      real z(iz),pbl,p(iz),prate,t(iz),vapor(iz),
     1 jcorrect(iz), wetk(iz,nwet)
      real cfrac(iz)

      real c303,c302,vp0
      parameter(C303=19.83,C302=5417.4,vp0 = 611.29)

      ESAT(TEMK)=.611*EXP(C303-C302/TEMK)       ! for calculating saturated water vapor pressure  
      QSAT(ESAT1,PCB)=ESAT1*.622/(PCB-ESAT1)    ! TEMK is ambient temperature in K, PCB is the pressue in KPa
                                                ! QSAT is the saturated humidity in kg/kg
      e_aerk(tempc) = vp0 * EXP( 17.625 * tempc / ( 243.04 + tempc ) )

      ktop=0
      do k=1,iz
	do L=1,4
 	 wetk(k,L)=0.
        enddo

         rh=vapor(k)/QSAT(ESAT(         ! computing relative humudity
     1		t(k)),p(k)/1000)         !convert to KPa
          if(z(k).lt. pbl ) then        ! within Convective Boundary
	  RHC=0.98                           ! Cloud relative humididty
	  if(rh.gt.rhc) then
	   cfrac(k)=0.34 *(RH - RHC)/(1. - RHC) ! SCHUMANN 89, AND WYNGAARD AND BROST 84
	  else
	   cfrac(k)=0.
	  endif
	 else
	  SG1= P(K)/P(1)
          RHC= 1. - 2.*SG1*(1.-SG1)*(1+1.732*(SG1-.5))
          IF(RH.GT.RHC) THEN
           cfrac(K)= ((RH - RHC)/(1. - RHC))**2  ! Geleyn et al 82
          ELSE
           cfrac(K) = 0.0
          END IF
	endif

      enddo

      cover_max=0.
      do k=2,iz-1
       if(cover_max.lt.cfrac(k)) then
        cover_max=cfrac(k)
	kcover=k
       endif
      enddo

      if(cover_max.gt.0) then  ! cloud exist

      ! Look for cloud top and base layer up and down from level of max RH.

      DO k = kcover, iz
        ktop = k - 1
        IF ( cfrac(k) < 0.5*cover_max ) exit
      ENDDO

      DO k = kcover, 1, -1
        kbase = k + 1
        IF ( cfrac(k) < 0.5*cover_max ) exit
      ENDDO

      cbase=0
      ctop=0
      DO k = 1, ktop
        IF ( k < kbase ) cbase = cbase + (z(k+1)-z(k))  ! cloud base height
        ctop = ctop + (z(k+1)-z(k))                     ! cloud top height
      ENDDO


      plcl = p(kbase-1)
      tlcl = t(kbase-1)
      pbase = p(kbase-1)
      tbase = t(kbase-1)

      iflag=0
      frac=0.
      wtbar=0.
      sumz=0.

      ! Follow moist adiabat up.

      DO k = kbase, ktop

        dp   = pbase - p(k)
        pbar = pbase - dp / 2.0
        tbar = tbase

        DO itr = 1, 5
          x1= lv0 * qsat( e_aerk( tbar-273.15), pbar/1000)/(rdgas*tbar)
          dtdp = rdgas * tbar / pbar / cpd * ( ( 1.0 + x1 ) /
     1            ( 1.0 + mvoma * lv0 / cpd / tbar * x1 ) )
          tad  = tbase - dp * dtdp
          tbar = ( tad + tbase ) * 0.5
        ENDDO

        ! Determine water content by fraction of adiabatic.

        tad   = MAX(tad, 150.0)
        IF ( tad > p(k) ) iflag = 1

        ! Pressure in Pascal = cb*1000

        wl    = 0.7 * EXP( ( p(k) - plcl ) / 8000.0 ) + 0.2
        qwsa  = qsat( e_aerk(tad - stdtemp), p(k)/1000  )

        qwat  = wl * ( qlcl - qwsa )
        qwat  = MAX(qwat, 0.0)

        twc   = qwat * p(k) * 1.0e3 / rdgas / t(k)

        wtbar = wtbar + twc * (z(k+1)-z(k))

        frac  = frac + cfrac(k) * (z(k+1)-z(k))
        sumz  = sumz + (z(k+1)-z(k))
        tbase = tad
        pbase = p(k)

      ENDDO

      if(sumz.lt.1e-4) then
       print*,'wrong sumz ',kbase,ktop,z
       stop
      endif

      ccover = frac / sumz

      else

      ccover = 0.

      endif
      print*,"end"
      end
