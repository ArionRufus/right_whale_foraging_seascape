program test_model_parallel

  !====================================================================72
  !
  ! This program advects particles into OPA currents fields.
  !
  ! MapsF 2011 & BenkortD 2018, from ChasseJ & LambertN
  !
  !====================================================================72

  use netcdf
  use partClass_mod
  use partFun_mod
  use lagrangian_motion_mod
#ifdef PAR
  use omp_lib
#endif

  implicit none



  ! Indices

  integer :: i, j, k, it, it2, ii, it3, itp1

  ! Time

  integer :: timer_start, timer_end, clock_rate
  real    :: timer

  integer :: nloop, daily_freq

  ! Advection & interpolation
  integer :: nb_layers  !the number of depth layers used for the model

  real*8, dimension(m,n) :: e1v, e2u

  real,   dimension(ilo+1) :: gdepw, gdept, e3w, e3t

  ! I/O

  integer :: ioerr

  ! I/O netcdf

  integer :: ncid, nc_fid, ndims, nvars, nglobalatts, unlimdimid

  integer :: timeid, lonid, latid, uid, vid, wid, tid, sid
  
  real    :: missing_u, missing_v

#ifdef YUMY
  integer :: diatom_id, flag_id, mesozoo_id, microzoo_id, pon_id, trn_id, parsurf_id

#endif

  integer :: mnc, nnc, ilonc, tnc, mnc_f, nnc_f, tnc_f, ilonc_f, group_f
  

  integer, dimension(2) :: timenc = (/ 0, 86401 /)
#ifdef CIOPS
  real       (kind = 8) :: timenc1, timenc2                                             !the times of first 2 files if data = CIOPS -> kind has to be 8
  
  real       (kind = 8) :: dtnc8                                                        !for CIOPS, temp dtnc of type 8
  integer    (kind = 8) :: istart8, iend8, kd_ini8                                      !for CIOPS, to set the time
  integer    (kind = 8) :: itp18, kade8                                                 !new dates in kind = 8
#endif
  integer               :: tday,tmon,tyear,kade, it4, hour, hh, bi6                     !find the new date
  character (len = 100) :: file_nc, char_tyear, char_tmon, char_tday, char_bi6, char_hh !write file's name from new date

  integer               :: dim

  character(len=20)     :: dim_name

  ! Particles

  type(partList)    :: zooList

  class(*), pointer :: curr

  integer :: npart_init
  real    :: tp, fd, ph
  !real    :: tp, diat, flag, micz, mesz, max_chlA, fd

  !====================================================================72
  ! Date & time

  call system_clock(timer_start)

  call idate(today) ! today(1)=month, (2)=day, (3)=year

  call itime(now)   ! now(1)=hour, (2)=minute, (3)=second

  write(*,'(/a,i2,1x,i2,1x,i4,a,i2,a,i2,a,i2,a)') ' Date: ',(today(i),i=1,3), ' : ', now(1),'h ',now(2),'m ',now(3),'s'

  dt = 1800 ! in sec

  ! First kstep of the run
  initial = 1

  kstep   = 0 ! used to update time

  ! Start at 0h00 -> nighttime
  light_day = .false.

  write(*,'(/a)') ' ! HAVE TO PROVIDE THE START & END DATES IN THE NAMELIST FILE "Run/run.list" !'
  write(*,'(/a)') ' ! HAVE TO PROVIDE SOME PARTICLES RELATED VARIABLES IN THE NAMELIST TOO      !'

  !====================================================================72
  ! Initialization

  ! Initialize arrays
  call zeros

  ! Initialize pseudo-random number generator
  call init_random_seed

  open(10,file='run.list',status='old',action='read',err=1,iostat=ioerr)
1 if(ioerr/=0) stop '!!! Problem opening file run.list !!!'
  read(10, nml=run_nml)
  close(10)

  open(10,file='run.list',status='old',action='read',err=2,iostat=ioerr)
2 if(ioerr/=0) stop '!!! Problem opening file run.list !!!'
  read(10, nml=part_nml)
  close(10)

  select case (part_spec)
  case(0)
     write(*,'(/a)')  ' *************************'
     write(*,'(a,/)') ' *** Generic particles ***'
  case(1)
     write(*,'(/a)')  ' **************************************************'
     write(*,'(a,/)') ' *** Krill species = Meganychtiphanes norvegica'
  case(2)
     write(*,'(/a)')  ' **************************************************'
     write(*,'(a,/)') ' *** Krill species = Thysanoessa raschii'
  case(3)
     write(*,'(/a)')  ' **************************************************'
     write(*,'(a,/)') ' *** Calanus species = Calanus hyperboreus'
  case(4)
     write(*,'(/a)')  ' ***********************************************************'
     write(*,'(a,/)') ' *** Calanus species = Calanus finmarchicus (PLACEHOLDER)'
  case(5)
     write(*,'(/a)')  ' ***********************************************************'
     write(*,'(a,/)') ' *** Late hyperboreus'
  case(6)
     write(*,'(/a)')  ' ***********************************************************'
     write(*,'(a,/)') ' *** Young hyperboreus'
  case(7)
     write(*,'(/a)')  ' ***********************************************************'
     write(*,'(a,/)') ' *** Late finmarchicus'
  case(8)
     write(*,'(/a)')  ' ***********************************************************'
     write(*,'(a,/)') ' *** Young finmarchicus'

  end select

  ! Output frequency; convert unit from h to dt (=kstep)
  output_freq2 = int(output_freq*3600/dt)

  !--------------------------------------------------------------------72
  ! Define the domain

  ! Get cells dx and dy

  call get_dx_dy(e1v,e2u,m,n,ilo)

  dlx = e1v ! Array used in trajectory

  dly = e2u ! Array used in trajectory

  ! Get layers thickness
!!! WARNING: hard-coded the 46 or 100 original layers of NEMO here, otherwise the
!!! layers thickness would be all wrong !!!
#ifdef CIOPS
  nb_layers = ilo+1
#else
  nb_layers = 46
#endif


! Compute vertical grid size

  call gdep(nb_layers,gdepw,gdept,e3w,e3t)

!dz gets rid of the 0 value, and goes to the first value deeper than 550m :
#ifdef CIOPS
  dz(1:ilo)  = gdepw(2:ilo+1)
#else
  dz(1:ilo-1)  = gdepw(2:ilo)

  dz(ilo)      = gdepw(ilo)+e3w(ilo)
#endif

!dd is the bin width :
  dd(1)        = dz(1)

  up_depth(1)  = 0.

  mid_depth(1) = dz(1)*.5 !Joel


  do j = 2,ilo

     dd(j)        = dz(j)-dz(j-1)

     up_depth(j)  = dz(j-1) ! Upper interface depth (positive)

     mid_depth(j) = (dz(j)+dz(j-1))*.5 !Joel

  enddo

#ifndef CIOPS
  ! Get number of layers
  ! for NEMO, just read nlayer
  open(10,file='data/bathy_level2.dat',status='old',action='read',err=3,iostat=ioerr)
3 if(ioerr/=0) stop ' !!! Problem opening file data/bathy_level2.dat !!!'

  do i = 1,m
     read(10,*) (nlayer(i,j),j=1,n)
  enddo

  close(10)

  where (nlayer>int2(ilo)) nlayer = int2(ilo)
#endif

  !--------------------------------------------------------------------72
  ! Open NetCDF physical input file

  ioerr = nf90_open(trim(infile),nf90_nowrite,ncid)
  if(ioerr/=nf90_noerr) then
     print*, trim(nf90_strerror(ioerr))
     stop ' !!! Problem opening NetCDF physical input file !!!'
  endif

  ioerr = nf90_inquire(ncid, ndims, nvars, nglobalatts, unlimdimid)
  if(ioerr/=nf90_noerr) then
     print*, trim(nf90_strerror(ioerr))
     stop ' Problem getting information from NetCDF physical input file !!!'
  endif

  ! Getting array dimensions and check against declaration in module part_def

  mnc   = 0 ! Initialisation of dimensions for 3D fields

  nnc   = 0

  ilonc = 0

  tnc   = 0


  do i = 1,ndims

     ioerr = nf90_inquire_dimension(ncid,i,dim_name,dim)
     if(ioerr/=nf90_noerr) then
        print*, trim(nf90_strerror(ioerr))
        stop ' !!! Problem getting dimension name from NetCDF physical input file !!!'
     endif

     if(dim_name=='x') mnc = dim ! dimensions are backward in netcdf

     if(dim_name=='y') nnc = dim

#ifndef CIOPS
     if(dim_name=='depthu') ilonc = dim

     if(dim_name=='time_counter') tnc = dim
#endif

  enddo

#ifdef DEBUG
  write(*,'(/a)')    ' *** Physical forcing domain dimensions '
  write(*,'(/a,i5)') ' mnc   = ',mnc
  write(*,'(a,i5)')  ' nnc   = ',nnc
  write(*,'(a,i5)')  ' ilonc = ',ilonc
  write(*,'(a,i5)')  ' tnc   = ',tnc
#endif

#ifdef CIOPS
  if (mnc==0.or.nnc==0)              stop ' !!! A dimension (m|n)     is missing !!!'
#else
  if (mnc==0.or.nnc==0.or.ilonc==0)  stop ' !!! A dimension (m|n|ilo) is missing !!!'

  if (mnc<m .or.nnc<n .or.ilonc<ilo) stop ' !!! A dimension (m|n|ilo) is wrong   !!!'

  if (tnc==0)                        stop ' !!! Problem with time_counter        !!!'
#endif

  ! Get lon/lat IDs

  ioerr = nf90_inq_varid(ncid,"nav_lon",lonid)
  if(ioerr/=nf90_noerr) then
     print*, trim(nf90_strerror(ioerr))
     stop ' !!! Problem getting lon ID from NetCDF physical input file !!!'
  endif

  ioerr = nf90_inq_varid(ncid,"nav_lat",latid)
  if(ioerr/=nf90_noerr) then
     print*, trim(nf90_strerror(ioerr))
     stop ' !!! Problem getting lat ID from NetCDF physical input file !!!'
  endif

  ! Get lon/lat values

  ioerr = nf90_get_var(ncid,lonid,lon)
  if(ioerr/=nf90_noerr) then
     print*, trim(nf90_strerror(ioerr))
     stop ' !!! Problem getting lon value from NetCDF physical input file !!!'
  endif

  ioerr = nf90_get_var(ncid,latid,lat)
  if(ioerr/=nf90_noerr) then
     print*, trim(nf90_strerror(ioerr))
     stop ' !!! Problem getting lat value from NetCDF physical input file !!!'
  endif

  ! Get u/v/w IDs

#ifdef CIOPS
  ioerr = nf90_inq_varid(ncid,"uo",uid)
#else
  ioerr = nf90_inq_varid(ncid,"vozocrtx",uid)
#endif
  if(ioerr/=nf90_noerr) then
     print*, trim(nf90_strerror(ioerr))
     stop ' !!! Problem getting u ID from NetCDF physical input file !!!'
  endif


#ifdef CIOPS
  ioerr = nf90_inq_varid(ncid,"vo",vid)
#else
  ioerr = nf90_inq_varid(ncid,"vomecrty",vid)
#endif
  if(ioerr/=nf90_noerr) then
     print*, trim(nf90_strerror(ioerr))
     stop ' !!! Problem getting v ID from NetCDF physical input file !!!'
  endif

#ifdef VERT
#ifndef CIOPS
  ioerr = nf90_inq_varid(ncid,"vovecrtz",wid)
  if(ioerr/=nf90_noerr) then
     print*, trim(nf90_strerror(ioerr))
     stop ' !!! Problem getting w ID from NetCDF physical input file !!!'
  endif
#endif
#endif

  ! Get salinity ID

#ifdef CIOPS
  ioerr = nf90_inq_varid(ncid,"so",sid)
#else
  ioerr = nf90_inq_varid(ncid,"vosaline",sid)
#endif
  if(ioerr/=nf90_noerr) then
     print*, trim(nf90_strerror(ioerr))
     stop ' !!! Problem getting S ID from NetCDF physical input file !!!'
  endif

  ! Get temperature ID

#ifdef CIOPS
  ioerr = nf90_inq_varid(ncid,"thetao ",tid)
#else
  ioerr = nf90_inq_varid(ncid,"votemper",tid)
#endif
  if(ioerr/=nf90_noerr) then
     print*, trim(nf90_strerror(ioerr))
     stop ' !!! Problem getting Temp ID from NetCDF physical input file !!!'
  endif

  !---------------------------------------------------------------------
  ! Define inner time loop according to input file record frequency

#ifdef CIOPS
  ! Get first time :
  ioerr = nf90_inq_varid(ncid,"time_counter",timeid)
  if(ioerr/=nf90_noerr) then
     print*, trim(nf90_strerror(ioerr))
     stop ' !!! Problem getting time ID from NetCDF input file !!!'
  endif

  ioerr = nf90_get_var(ncid,timeid,timenc1)
  if(ioerr/=nf90_noerr) then
     print*, trim(nf90_strerror(ioerr))
     stop ' !!! Problem getting first 2 time values from NetCDF input file !!!'
  endif

#ifdef DEBUG
  write(*,'(/a,f12.0)') ' First  timestamp of NetCDF input file : ', timenc1
#endif

  ! Close first file :
  ioerr = nf90_close(ncid)
  if(ioerr/=nf90_noerr) then
     print*, trim(nf90_strerror(ioerr))
     stop ' !!! Problem closing forcing files !!!'
  endif

  ! Get second time :
  ioerr = nf90_open(trim(infile2),nf90_nowrite,ncid)
  if(ioerr/=nf90_noerr) then
     print*, trim(nf90_strerror(ioerr))
     stop ' !!! Problem opening NetCDF physical input file !!! '
  endif

  ioerr = nf90_inquire(ncid, ndims)
  if(ioerr/=nf90_noerr) then
     print*, trim(nf90_strerror(ioerr))
     stop ' Problem getting information from NetCDF physical input file !!!'
  endif

  ioerr = nf90_inq_varid(ncid,"time_counter",timeid)
  if(ioerr/=nf90_noerr) then
     print*, trim(nf90_strerror(ioerr))
     stop ' !!! Problem getting time ID from second NetCDF input file !!!'
  endif

  ioerr = nf90_get_var(ncid,timeid,timenc2)
  if(ioerr/=nf90_noerr) then
     print*, trim(nf90_strerror(ioerr))
     stop ' !!! Problem getting first 2 time values from second NetCDF input file !!!'
  endif

#ifdef DEBUG
  write(*,'(a,f12.0)') ' Second timestamp of NetCDF input file : ', timenc2
#endif

  ! Close second file
  ioerr = nf90_close(ncid)
  if(ioerr/=nf90_noerr) then
     print*, trim(nf90_strerror(ioerr))
     stop ' !!! Problem closing second forcing files !!!'
  endif

#else
  !for NEMO inputs :

  ! Get first 2 time axis values to find input file timestep (must be constant!)
  ioerr = nf90_inq_varid(ncid,"time_counter",timeid)
  if(ioerr/=nf90_noerr) then
     print*, trim(nf90_strerror(ioerr))
     stop ' !!! Problem getting time ID from NetCDF input file !!!'
  endif

  ioerr = nf90_get_var(ncid,timeid,timenc,(/1/),(/2/))
  if(ioerr/=nf90_noerr) then
     print*, trim(nf90_strerror(ioerr))
     stop ' !!! Problem getting first 2 time values from NetCDF input file !!!'
  endif

#endif


! Determine the time step :
! Assuming time_counter unit is always in second

#ifdef CIOPS

#ifdef BACK
  dtnc8 = timenc1 - timenc2
#else
  dtnc8 = timenc2 - timenc1
#endif
  dtnc  = int(dtnc8)

#else
  dtnc  = timenc(2) - timenc(1)
#endif

  nloop = int(dtnc/dt)

  write(*,'(/a,i5,a)')          ' Input file timestep                   : ',dtnc,'s'

  write(*,'(/a,i5,a)')          ' Simulation timestep                   : ',dt,'s'

  write(*,'(/a,i5)')            ' Inner time loop                       : ',nloop


  !--------------------------------------------------------------------72
  ! Fix the problem with date in NEMO
  ! The date in nc files is in seconds since tnc_origin in the julian calendar.
  ! Change it into kd_ini : days since year 0 in the gregorian calendar.

  ! CAUTION ! HAVE TO CHANGE TIME ORIGIN IN THE NAMELIST ACCORDING TO INPUT FILE
  call cday2(tnc_origin(1),tnc_origin(2),tnc_origin(3),kd_ini)

  ! Get the first and last simulation day
  call cday2(dstart,mstart,ystart,istart)
  call cday2(dend,mend,yend,iend)

#ifdef CIOPS
  ! Detemine the time of start / end of the model since historical time kd_ini in time steps dtnc :
  istart8 = int(istart, 8) ! change from real4 for to 8 because the numbers are too large
  kd_ini8 = int(kd_ini, 8)

  istart8 = (istart8-kd_ini8)*86400/dtnc ! unit = input file record period

  iend8   = int(iend, 8)
  iend8   = (iend8-kd_ini8)*86400/dtnc

  ! Change back into int*4 to work with other variables :
  istart  = int(istart8, 4)
  iend    = int(iend8,   4)

#else
  istart  = (istart-kd_ini)*86400/dtnc ! unit = input file record period
  iend    = (iend-kd_ini)*86400/dtnc
#endif

#ifdef BACK
  ! Enables computation backward in time: currents read backward and sign changed
  if (iend>istart) stop ' !!! CHECK run.list FILE -> END DATE SHOULD BE < START DATE IN BACKWARD COMPUTATION !!!'
  final = (istart-iend+1)*nloop ! unit = dt (kstep)
#else
  if (iend<istart) stop ' !!! CHECK run.list FILE -> END DATE SHOULD BE > START DATE IN FORWARD COMPUTATION !!!'
  final = (iend-istart+1)*nloop
#endif

  write(*,'(/a,i5)')            ' Number of model time steps            : ', final


  write(*,'(/a)')               ' *****************************************************'
  write(*,'(a,i2,1x,i2,1x,i4)') ' *** Simulation starting day           : ', dstart, mstart, ystart

  !--------------------------------------------------------------------72
  ! Get initial u/v/w values

#ifdef CIOPS

  ! open again first file :
  ioerr = nf90_open(trim(infile),nf90_nowrite,ncid)
  if(ioerr/=nf90_noerr) then
     print*, trim(nf90_strerror(ioerr))
     stop ' !!! Problem opening NetCDF physical input file !!! '
  endif

  ! get "u"
  ioerr = nf90_get_var(ncid,uid,u2)
  if(ioerr/=nf90_noerr) then
     print*, trim(nf90_strerror(ioerr))
     stop ' !!! Problem getting u value from NetCDF input file !!!'
  endif

  ioerr = nf90_get_att(ncid,uid,"missing_value",missing_u)
  if(ioerr/=nf90_noerr) then
     print*, trim(nf90_strerror(ioerr))
     stop ' !!! Problem getting missing value from NetCDF input file !!!'
  endif

  ! get "v"
  ioerr = nf90_get_var(ncid,vid,v2)
  if(ioerr/=nf90_noerr) then
     print*, trim(nf90_strerror(ioerr))
     stop ' !!! Problem getting v value from NetCDF input file !!!'
  endif

  ioerr = nf90_get_att(ncid,vid,"missing_value",missing_v)
  if(ioerr/=nf90_noerr) then
     print*, trim(nf90_strerror(ioerr))
     stop ' !!! Problem getting missing value from NetCDF input file !!!'
  endif



  
  ! close first file
  ioerr = nf90_close(ncid)
  if(ioerr/=nf90_noerr) then
     print*, trim(nf90_strerror(ioerr))
     stop ' !!! Problem closing first file !!!'
  endif

  ! get numbers of layers (nlayer)
  bathy = 1
  
  where( u2 .eq. missing_u ) bathy = 0 ! use the flag for missing value in the netcdf file (~ nan)
  where( v2 .eq. missing_v ) bathy = 0
  
  nlayer = int2( sum(bathy,3) )
  where ( nlayer>int2(ilo) ) nlayer = int2(ilo)

  ! remove missing/fill values from currents
  where( u2 .eq. missing_u ) u2 = 0
  where( v2 .eq. missing_v ) v2 = 0

#else
  !get initial values for nemo imput files :
  ioerr = nf90_get_var(ncid,uid,u2,(/1,1,1,istart/))
  if(ioerr/=nf90_noerr) then
     print*, trim(nf90_strerror(ioerr))
     stop ' !!! Problem getting u value from NetCDF input file !!!'
  endif

  ioerr = nf90_get_var(ncid,vid,v2,(/1,1,1,istart/))
  if(ioerr/=nf90_noerr) then
     print*, trim(nf90_strerror(ioerr))
     stop ' !!! Problem getting v value from NetCDF input file !!!'
  endif

#ifdef VERT
  ioerr = nf90_get_var(ncid,wid,w2,(/1,1,1,istart/))
  if(ioerr/=nf90_noerr) then
     print*, trim(nf90_strerror(ioerr))
     stop ' !!! Problem getting v value from NetCDF input file !!!'
  endif
#endif

#endif

#ifdef BACK
  u2 = -1.*u2
  v2 = -1.*v2
#ifdef VERT
  w2 = -1.*w2
#endif
#endif


  !====================================================================72
  ! MAIN LOOP
#ifdef BACK
  do it3 = istart,iend+1,-1 ! Main loop in time, following input records frequency
     it  = it3              ! Do not apply nc_repeat for back tracking
#else
  do it3 = istart,iend-1

     ! repeat reading netcdf if nc_repeat is "on" and simulation time is greater than last time stamp of nc time
     if (nc_repeat .eqv. .true. ) then
        it = it3-int(it3/tnc)*tnc
      else
        it = it3
      endif

      if (it==0) then
         it   = tnc
         itp1 = 1
      else
         itp1 = it+1
      endif
#endif

#ifdef DEBUG
      if( it3 .eq. istart ) then
         write(*,'(/a)')  ' *****************************************************'
         write(*,'(a,/)') ' *** Read physical forcing values'
      endif
#endif

      !-----------------------------------------------------------------72
      ! Linear interpolation ! MapsF

      u1 = u2

      v1 = v2

#ifdef VERT
      w1 = w2
#endif

      u  = u1 ! for output purpose

      v  = v1

#ifdef VERT
      w  = w1
#endif

#ifdef CIOPS
      !--
      !--Find the name of the new input file :

#ifdef BACK
      it4 = it3-1
#else
      it4 = it3+1         ! the file to open is the one at time t+1
#endif

      itp18 = int(it4, 8) ! have to change the value in type 8 because of its size

      kade8 = itp18 * dtnc8 / 86400
      kade8 = kade8 + kd_ini8
      kade  = int(kade8, 4)

      ! find the hour of the day from the time step :
      hour = modulo( (it4-istart) * dtnc / 3600, 24 )

      ! the name of the file includes an hour modulo 6, we compute it here :
      hh  = modulo( hour, 6 )
      bi6 = hour - hh
      
      ! find the date in the gregorian calendar :
      call dmy2(tday,tmon,tyear,kade)
#ifdef DEBUG
      if( it3 .eq. istart ) then
         write(*,'(a,i7)') ' year  from Gregorian calendar : ', tyear
         write(*,'(a,i7)') ' month from Gregorian calendar : ', tmon
         write(*,'(a,i7)') ' day   from Gregorian calendar : ', tday
         write(*,'(a,i7)') ' hour of the day               : ', hour
      endif
#endif

      ! concatenate the differents times to create the name of the file :
      write( unit=char_tyear, fmt=* ) tyear
      write( unit=char_tmon,  fmt=* ) tmon
      write( unit=char_tday,  fmt=* ) tday
      write( unit=char_bi6,   fmt=* ) bi6
      write( unit=char_hh,    fmt=* ) hh

      if (tmon < 10) then ! padd with 0 to fit the name of the file
         char_tmon = trim( adjustl("0")) // trim(adjustl(char_tmon) )
      end if

      if (tday < 10) then
         char_tday = trim( adjustl("0")) // trim(adjustl(char_tday) )
      end if

      if (bi6 < 10) then
         char_bi6     = trim( adjustl("0")) // trim(adjustl(char_bi6) )
      end if

      char_hh = trim( adjustl("00")) // trim(adjustl(char_hh) )


      file_nc = trim(adjustl(char_tyear)) // trim(adjustl(char_tmon)) // trim(adjustl(char_tday))
      file_nc = trim(adjustl(file_nc))       // trim(adjustl(char_bi6))     // trim(adjustl("_"))
      file_nc = trim(adjustl(infile_loc))    // trim(adjustl(file_nc))      // trim(adjustl(char_hh)) // trim(adjustl(".nc"))

#ifdef DEBUG
      if( it3 .eq. istart ) then
         write(*,'(/a,a)') ' filename of the NetCDF input file : ', file_nc
      endif
#endif

      !--
      !--Open and read the new input file :

      ! Open the physical input file of the loop :
      ioerr = nf90_open(trim(file_nc),nf90_nowrite,ncid)
      if(ioerr/=nf90_noerr) then
         print*, trim(nf90_strerror(ioerr))
         stop ' !!! Problem opening NetCDF physical input file in the loop !!! '
      endif


      ! Get new file's currents :
      ioerr = nf90_inq_varid(ncid,"uo",uid)
      if(ioerr/=nf90_noerr) then
         print*, trim(nf90_strerror(ioerr))
         stop ' !!! Problem getting u ID from NetCDF physical input file !!!'
      endif

      ioerr = nf90_get_var(ncid,uid,u2)
      if(ioerr/=nf90_noerr) then
         print*, trim(nf90_strerror(ioerr))
         stop ' !!! Problem getting u value from NetCDF input file !!!'
      endif

      ioerr = nf90_get_att(ncid,uid,"missing_value",missing_u)
      if(ioerr/=nf90_noerr) then
         print*, trim(nf90_strerror(ioerr))
         stop ' !!! Problem getting missing value from NetCDF input file !!!'
      endif

      where( u2 .eq. missing_u ) u2 = 0.


      ioerr = nf90_inq_varid(ncid,"vo",vid)
      if(ioerr/=nf90_noerr) then
         print*, trim(nf90_strerror(ioerr))
         stop ' !!! Problem getting v ID from NetCDF physical input file !!!'
      endif

      ioerr = nf90_get_var(ncid,vid,v2)
      if(ioerr/=nf90_noerr) then
         print*, trim(nf90_strerror(ioerr))
         stop ' !!! Problem getting v value from NetCDF input file !!!'
      endif

      ioerr = nf90_get_att(ncid,vid,"missing_value",missing_v)
      if(ioerr/=nf90_noerr) then
         print*, trim(nf90_strerror(ioerr))
         stop ' !!! Problem getting missing value from NetCDF input file !!!'
      endif

      where( v2 .eq. missing_v ) v2 = 0.

#else
      ! get values from nemo input files :

#ifdef BACK
      ioerr = nf90_get_var(ncid,uid,u2,(/1,1,1,it-1/),(/m,n,ilo,1/))
#else
      ioerr = nf90_get_var(ncid,uid,u2,(/1,1,1,itp1/),(/m,n,ilo,1/))
#endif
      if(ioerr/=nf90_noerr) then
         print*, trim(nf90_strerror(ioerr))
         stop ' !!! NEMO Problem getting u2 value from NetCDF input file !!!'
      endif

#ifdef BACK
      ioerr = nf90_get_var(ncid,vid,v2,(/1,1,1,it-1/),(/m,n,ilo,1/))
#else
      ioerr = nf90_get_var(ncid,vid,v2,(/1,1,1,itp1/),(/m,n,ilo,1/))
#endif
      if(ioerr/=nf90_noerr) then
         print*, trim(nf90_strerror(ioerr))
         stop ' !!! NEMO Problem getting v2 value from NetCDF input file !!!'
      endif

#ifdef VERT
#ifdef BACK
      ioerr = nf90_get_var(ncid,wid,w2,(/1,1,1,it-1/),(/m,n,ilo,1/))
#else
      ioerr = nf90_get_var(ncid,wid,w2,(/1,1,1,itp1/),(/m,n,ilo,1/))
#endif
      if(ioerr/=nf90_noerr) then
         print*, trim(nf90_strerror(ioerr))
         stop ' !!! NEMO Problem getting w2 value from NetCDF input file !!!'
      endif
#endif
#endif


#ifdef BACK
      u2 = -1.*u2
      v2 = -1.*v2
#ifdef VERT
      w2 = -1.*w2
#endif
#endif

#ifdef DEBUG 
      if ( it3 .eq. istart ) then
         write(*,'(/a)') ' * Horizontal components U & V'
      endif
#endif


      ! Get salinity :
#ifdef CIOPS
      ioerr = nf90_inq_varid(ncid,"so",sid)
      if(ioerr/=nf90_noerr) then
         print*, trim(nf90_strerror(ioerr))
         stop ' !!! Problem getting S ID from NetCDF physical input file !!!'
      endif

      ioerr = nf90_get_var(ncid,sid,salt,count = (/m,n,2,1/)) ! Need only the first 2 layers
      if(ioerr/=nf90_noerr) then
         print*, trim(nf90_strerror(ioerr))
         stop ' !!! Problem getting S value from NetCDF input file !!!'
      endif
#else

      ioerr = nf90_get_var(ncid,sid,salt,(/1,1,1,it/),(/m,n,2,1/)) ! Need only the first 2 layers
      if(ioerr/=nf90_noerr) then
         print*, trim(nf90_strerror(ioerr))
         stop ' !!! Problem getting S value from NetCDF input file !!!'
      endif
#endif

      salt2 = sum(salt,3)*0.5 ! Average of the first 2 layers
#ifdef DEBUG
      if ( it3 .eq. istart ) then
         write(*,'(a)') ' * Salinity'
      endif
#endif
      where (nlayer==0) salt2 = nf90_fill_real ! For later use when recording values in netcdf file


! Get temperature :
#ifdef CIOPS
      ioerr = nf90_inq_varid(ncid,"thetao ",tid)
      if(ioerr/=nf90_noerr) then
         print*, trim(nf90_strerror(ioerr))
         stop ' !!! Problem getting Temp ID from NetCDF physical input file !!!'
      endif

      ioerr = nf90_get_var(ncid,tid,temp,count = (/m,n,4,1/))
      if(ioerr/=nf90_noerr) then
         print*, trim(nf90_strerror(ioerr))
         stop ' !!! Problem getting Temp value from NetCDF input file !!!'
      endif
#else

      ioerr = nf90_get_var(ncid,tid,temp,(/1,1,1,it/),(/m,n,ilo,1/))
      if(ioerr/=nf90_noerr) then
         print*, trim(nf90_strerror(ioerr))
         stop ' !!! Problem getting Temp value from NetCDF input file !!!'
      endif
#endif

      temp2=sum(temp(1:m,1:n,1:4),3)*0.25  ! average of top 4 layer for predation module

      where (nlayer==0) temp2 = nf90_fill_real

#ifdef DEBUG
      if ( it3 .eq. istart ) then
         write(*,'(a)') ' * Temperature'
      endif
#endif


      ! Get food :
#ifdef YUMY

      if(part_spec .eq. 3) then

         ioerr = nf90_get_var(ncid,diatom_id,diatom,(/1,1,1,it/),(/m,n,ilo,1/))
           if(ioerr/=nf90_noerr) then
            print*, trim(nf90_strerror(ioerr))
            stop ' !!! Problem getting diatom value from NetCDF input file !!!'
         endif

#ifdef DEBUG
         if ( it3 .eq. istart ) then
            write(*,'(a)') ' * Diatoms'
         endif

#endif
         ioerr = nf90_get_var(ncid,flag_id,flag,(/1,1,1,it/),(/m,n,ilo,1/))
         if(ioerr/=nf90_noerr) then
            print*, trim(nf90_strerror(ioerr))
            stop ' !!! Problem getting flag value from NetCDF input file !!!'
         endif

#ifdef DEBUG
         if ( it3 .eq. istart ) then
            write(*,'(a)') ' * Flagellates'
         endif

#endif
         ioerr = nf90_get_var(ncid,microzoo_id,microzoo,(/1,1,1,it/),(/m,n,ilo,1/))
         if(ioerr/=nf90_noerr) then
            print*, trim(nf90_strerror(ioerr))
            stop ' !!! Problem getting microzoo value from NetCDF input file !!!'
         endif

#ifdef DEBUG
         if ( it3 .eq. istart ) then
            write(*,'(a)') ' * Microzooplankton'
         endif

#endif
         food    = diatom + flag + microzoo            ! for ingestion computation
         phyto   = diatom(1:m,1:n,1) + flag(1:m,1:n,1) ! for predation computation
         phyto3d = diatom + flag                       ! 3D phytoplankton

         ioerr = nf90_get_var(ncid,parsurf_id,parsurf,(/1,1,it/),(/m,n,1/))
         if(ioerr/=nf90_noerr) then
            print*, trim(nf90_strerror(ioerr))
            stop ' !!! Problem getting parsurf value from NetCDF input file !!!'
         endif

#ifdef DEBUG
         if ( it3 .eq. istart ) then
            write(*,'(a)') ' * Parsurf'
         endif

#endif
      endif

#endif

      !--------------------------------------------------------------72

#ifdef DEBUG
      if ( it3 .eq. istart ) then
#ifdef BACK
         write(*,'(/a,i6,a,i6)') ' *** Time step : ', istart - it3 + 1, '   of ', istart - iend
#else
         write(*,'(/a,i6,a,i6)') ' *** Time step : ', it3 - istart + 1, '   of ', iend - istart
#endif
      endif
#endif

#ifdef CIOPS
      if (it==istart) then
#else
      if (it==istart .and. it3<tnc) then
#endif
         tzero    = 0

         tzero(4) = dstart

         tzero(5) = mstart

         tzero(6) = ystart

         time     = tzero

         call sun_time(time_sunrise,time_sunset,48.5,-63.) ! NL -69.3

         p_night   = (24 - (time_sunset - time_sunrise))/24

         !time_go_down = 4.0 ! KB 4am

         !time_go_up = 22.0  ! KB 10pm

  !!! ADD A case(part_spec)
#ifdef PHYSIO
         ratio_ing = 1 / p_night
#endif

         ! Initialize particles

         select case (part_spec)
         case (0)
            call generateZoo(                zooList,'run.list' )
         case (1)
            call generateKrill(             zooList,'run.list' )
         case (2)
            call generateKrill(             zooList,'run.list' )
         case (3)
            call generateCalanus(           zooList,'run.list' )
         case (5)
            call generateLate_hyperboreus(   zooList,'run.list', in_diapause )
         case (6)
            call generateYoung_hyperboreus(  zooList,'run.list' )
         case (7)
            call generateLate_finmarchicus( zooList,'run.list', in_diapause )
         case (8)
            call generateYoung_finmarchicus( zooList,'run.list' )
         end select

! Create netCDF file & record initial particles data
#ifdef FTLE
         call output_ftle_netcdf(zooList, part_spec, 0)
#else
         call netcdfPart(zooList, part_spec, 0)
#endif

         npart_init = npart

      endif


      daily_freq = 86400/dt
      if (mod(kstep,daily_freq)==0) then
         write(*,'(/a,3(i4))') '--- Time : ', (time(ii),ii=6,4,-1)
      endif

      !--------------------------------------------------------------72
      ! Inner time loop

      du = u2 - u1

      dv = v2 - v1

#ifdef VERT
      dw = w2 - w1
#endif

      do it2 = 0,nloop-1
         
         u = u1 + du * it2 / nloop

         v = v1 + dv * it2 / nloop

#ifdef VERT
         w = w1 + dw * it2 / nloop
#endif

         !call optimal_depth(zooList)
           
         call trajectory(zooList, part_spec)
     

         !-----------------------------------------------------------72
         ! Physiological routine

         select case (part_spec)
         case (1)
            call evolveKrill(zooList)
         case (2)
            call evolveKrill(zooList)
         case (3)
            call evolveCalanus(zooList,'run.list')
         case (5)  
            call evolveLate_hyperboreus(zooList, in_diapause)
         case (6)
            call evolveYoung_hyperboreus(zooList, in_diapause)
         case (7)
            call evolveLate_finmarchicus(zooList, in_diapause)
         case (8)
            call evolveYoung_finmarchicus(zooList, in_diapause)
         end select

         !-----------------------------------------------------------72
         ! Output
         if (kstep>0.and.mod(kstep,output_freq2)==0) then
#ifdef FTLE
          !  write(*,'(/a,3(i4))') 'call output FTLE'
            call output_ftle_netcdf(zooList, part_spec, 1)
#else
            call netcdfPart(zooList, part_spec, 1)
#endif
         endif

         kstep = kstep+1

         call update2(kstep)

      enddo

  ioerr = nf90_close(ncid)
  if(ioerr/=nf90_noerr) then
     print*, trim(nf90_strerror(ioerr))
     stop ' !!! Problem closing forcing files !!!'
  endif

  !-----------------------------------------------------------------72
  ! MAIN LOOP END
  enddo

! close file
#ifdef FTLE
  call output_ftle_netcdf(zooList, part_spec, 2)
#else
  call netcdfPart(zooList, part_spec, 2)
#endif

  ! Print computation time
  call system_clock(timer_end, clock_rate)
  timer = real(timer_end-timer_start) / real(clock_rate)

  write(*,'(/a,3i1)') '***   THE END   ',timer_end,timer_start,clock_rate
  write(*,'(/a,i10,a,i10,a,f10.2,a/)') 'Computation time = ', &
        floor(timer/3600.),'h', &
        floor(amod(timer/60.,60.)),'m', &
        amod(timer,60.),'s'

!=====================================================================72

end program test_model_parallel
