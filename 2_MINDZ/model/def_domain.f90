subroutine get_dx_dy(e1v,e2u,mpass,npass)

  use netcdf

  implicit none

  integer,intent(in) :: mpass, npass                        !resolution lon lat

  integer            :: e1vid, e2uid                        !advection vectors localisation ids

  real*8,intent(out) :: e1v(mpass,npass), e2u(mpass,npass)  !advection vectors localisation

  integer            :: ioerr, ncid                         !for NetCDF manipulations

  character(len=20)  :: dim_name                            !names of NetCDF dimensions

  integer            :: ndims, dim, m, n, ilo, t            !dimensions of NetCDF

  integer            :: i                                   !for the loop


 
  !--
  !-- Verify that the dimensions of NetCDF files match with the model :

  ! open domain NetCDF file :
#ifdef CIOPS

  ioerr = nf90_open(path = "data/GEOMETRY/coordinates_NWA36_crop.nc", mode = nf90_nowrite, ncid = ncid)
  if (ioerr /= nf90_noerr) then
    print*, trim(nf90_strerror(ioerr))
    stop ' !!! Problem opening data/GEOMETRY/coordinates_NWA36_crop.nc'
  endif
#else
  ioerr = nf90_open('data/coordinates.nc', nf90_share, ncid)
  if(ioerr/=nf90_noerr) then
    print*, trim(nf90_strerror(ioerr))
    stop ' !!! Problem opening data/coordinates.nc'
  endif
#endif

  ! get number of dimensions :
  ioerr = nf90_inquire(ncid, ndims)
  if(ioerr/=nf90_noerr) then
     print*, trim(nf90_strerror(ioerr))
     stop ' !!! Problem getting information from geometry file !!!'
  endif

  m   = 0

  n   = 0

  ilo = 0 !Initialisation of dimensions for 3D field

  t   = 0


  do i = 1,ndims

     ioerr = nf90_inquire_dimension(ncid,i,dim_name,dim)
     if(ioerr/=nf90_noerr) then
        print*, trim(nf90_strerror(ioerr))
        stop ' !!! Problem getting dimension name from geometry file !!!'
     endif

     if(dim_name=='x')   m = dim  ! dimensions are backward in netcdf

     if(dim_name=='y')   n = dim

#ifndef CIOPS
     if(dim_name=='z') ilo = dim

     if(dim_name=='time') t = dim
#endif

  enddo

  if (m/=mpass.or.n/=npass) stop '!!! A dimension (m|n) does not match in geometry file !!!'

  !--
  !-- Get the dx,dy (= e1v, e2u, vector locations with C Arakawa grid for particle advection) :

  ! get e1v id :
  ioerr = nf90_inq_varid(ncid,"e1v",e1vid)
  if(ioerr/=nf90_noerr) then
     print*, trim(nf90_strerror(ioerr))
     stop ' !!! Problem getting e1v ID from geometry file !!!'
  endif

  ! get e2u id :
  ioerr = nf90_inq_varid(ncid,"e2u",e2uid)
  if(ioerr/=nf90_noerr) then
     print*, trim(nf90_strerror(ioerr))
     stop ' !!! Problem getting e2u ID from geometry file !!!'
  endif

  ! get e1v values :
  ioerr = nf90_get_var(ncid,e1vid,e1v)
  if(ioerr/=nf90_noerr) then
     print*, trim(nf90_strerror(ioerr))
     stop ' !!! Problem getting e1v value from geometry file !!!'
  endif

  ! get e2u values :
  ioerr = nf90_get_var(ncid,e2uid,e2u)
  if(ioerr/=nf90_noerr) then
     print*, trim(nf90_strerror(ioerr))
     stop ' !!! Problem getting e2u value from geometry file !!!'
  endif

  !--
  !-- Close NetCDF file :
  ioerr = nf90_close(ncid)
  if(ioerr/=nf90_noerr) then
     print*, trim(nf90_strerror(ioerr))
     stop ' !!! Problem closing geometry file !!!'
  endif

end subroutine get_dx_dy

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!72

subroutine gdep(jpk,gdepw,gdept,e3w,e3t)

  use netcdf

  implicit none

  integer,intent(in) :: jpk                                        !number of depth bins in CIOPS data

  integer            :: ioerr, ncid                                !for NetCDF manipulations

  integer            :: jk

  integer            :: depth_id                                   !id of depth bins

  real               :: gdepw(jpk), gdept(jpk), e3w(jpk), e3t(jpk) !depth bins values

  real               :: ppkth, ppacr, ppdzmin, pphmax

  real               :: zkth, zacr, zdzmin, zhmax, za1, za0, zsur, zw, zt

  ! This subroutine is pererated in 2 parts, the first if the imput file is CIOPS-E and the second for NEMO imput file

#ifdef CIOPS
  !--
  !-- Get the depth values corresponding to the different depth bins :

  ! open NetCDF file "full_geometry_3D.nc":
  ioerr = nf90_open(path = "data/GEOMETRY/full_geometry_3D_crop.nc", mode = nf90_nowrite, ncid = ncid)
  if (ioerr /= nf90_noerr) then
     print*, trim(nf90_strerror(ioerr))
     stop ' !!! Problem opening data/GEOMETRY/full_geometry_3D_crop.nc !!!'
  endif

  ! get depth id :
  ioerr = nf90_inq_varid(ncid,"depth",depth_id)
  if(ioerr/=nf90_noerr) then
     print*, trim(nf90_strerror(ioerr))
     stop ' !!! Problem getting depth_id from data/GEOMETRY/full_geometry_3D.nc !!!'
  endif

  ! get depth values :
  ioerr = nf90_get_var(ncid,depth_id,gdepw,(/1/),(/jpk/))
  if(ioerr/=nf90_noerr) then
     print*, trim(nf90_strerror(ioerr))
     stop ' !!! Problem getting gdepw value from data/GEOMETRY/full_geometry_3D.nc !!!'
  endif

  ! close NetCDF file :
  ioerr = nf90_close(ncid)
  if (ioerr /= nf90_noerr) stop '!!! Problem closing file data/GEOMETRY/full_geometry_3D.nc !!!'


  ! give values to variables used with NEMO to prevent having empty variables
  gdept = gdepw
  e3w   = gdepw
  e3t   = gdepw

#else
  !--
  !-- Get the OPA parameters

  ppkth   = 23.56;

  ppacr   = 9.000

  ppdzmin = 6.0   !Sopa definition

  pphmax  = 5750

  zkth    = ppkth

  zacr    = ppacr

  zdzmin  = ppdzmin

  zhmax   = pphmax


  za1  =  ( ppdzmin - pphmax/(jpk-1) )                                   &
        / (  tanh((1-ppkth)/ppacr)                                       &
           - ppacr/(jpk-1)                                               &
           * ( log(cosh((jpk-ppkth)/ppacr))-log(cosh((1-ppkth)/ppacr)) ) &
          )

  za0  = ppdzmin - za1*tanh( (1-ppkth)/ppacr )

  zsur = -za0 - za1*ppacr*log( cosh( (1-ppkth)/ppacr ) )


  do jk = 1,jpk

     zw = jk

     zt = jk+0.5

     gdepw(jk) = ( zsur + za0*zw + za1*zacr*log( cosh( (zw-zkth)/zacr ) ) )

     if(jk==1) gdepw(jk) = 0.

     gdept(jk) = ( zsur + za0*zt + za1*zacr*log( cosh( (zt-zkth)/zacr) ) )

     e3w(jk)   = za0 + za1*tanh( (zw-zkth)/zacr )

     e3t(jk)   = za0 + za1*tanh( (zt-zkth)/zacr )

  enddo

#endif

  return

end subroutine gdep

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!72

subroutine sun_time(SunR,SunS,lat,lon)

  use partClass_mod, only : tzero, time, pi

  implicit none

  real, intent(in) :: lat, lon

  real, intent(out) :: SunR, SunS

  real :: r2d, d2r, tzone

  real :: gamma, eqtime, decl, ha

  integer :: D0, JD


  ! Angle conversion
  r2d = 180./pi ! radian to degree

  d2r = pi/180. ! degree to radian


  ! Time variables
  tzone = -5 ! Time Zone. - for West. Model input in UTC!

  call cday2(31,12,tzero(6)-1,D0)

  call cday2(time(4),time(5),time(6),JD)


  ! Equations from NOAA
  gamma  =  2*pi/365 * (JD-D0)

  eqtime =  229.18 * (  0.000075 + 0.001868*cos(gamma) &
          - 0.032077*sin(gamma) &
          - 0.014615*cos(2*gamma) &
          - 0.040849*sin(2*gamma) )

  decl   =  0.006918 - 0.399912*cos(gamma) + 0.070257*sin(gamma) &
          - 0.006758*cos(2*gamma) + 0.000907*sin(2*gamma) &
          - 0.002697*cos(3*gamma) + 0.00148*sin(3*gamma)

  ha     = acos( cos(90.833*d2r)/(cos(lat*d2r)*cos(decl)) - tan(lat*d2r)*tan(decl) )


  SunR   = mod( aint( 720+4.*(lon-ha*r2d)-eqtime )-tzone*60, 1440. )/60.

  SunS   = mod( aint( 720+4.*(lon+ha*r2d)-eqtime )-tzone*60, 1440. )/60.

  return

end subroutine sun_time

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!72

subroutine update2(jpass)

  !===================================================================72
  ! Subroutine to update time array
  !===================================================================72

  use partSubClass_mod, only : tzero, time, dt

  implicit none

  integer, intent(in) :: jpass

  integer :: j, isec, imin, ihour, iday, imon, iyear, kd, ishift

  real*8 :: isecnew, iminnew, ihournew ! CRITICAL: problem with floor function otherwise


  j     = jpass

  iday  = tzero(4)

  imon  = tzero(5)

  iyear = tzero(6)

  call cday2(iday,imon,iyear,kd)


  ! Enables computation backward in time
#ifdef BACK
  isecnew  = tzero(1) - dt*j
#else
  isecnew  = tzero(1) + dt*j
#endif

  isec     = modulo(int(isecnew),60)

  iminnew  = tzero(2) + floor(isecnew/60)

  imin     = modulo(int(iminnew),60)

  ihournew = tzero(3) + floor(iminnew/60)

  ihour    = modulo(int(ihournew),24)

  kd       = kd + floor(ihournew/24)

  call dmy2(iday,imon,iyear,kd)


  time(1) = isec

  time(2) = imin

  time(3) = ihour

  time(4) = iday

  time(5) = imon

  time(6) = iyear


  return

end subroutine update2


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!72

subroutine cday2(idd,imm,iyy_pass,kd)

   !===================================================================72
   !
   ! Given day, month (each 2 digits) and year (4 digits), 
   ! cday returns the day#, kd, based on the gregorian calendar.
   ! The gregorian calendar, currently 'universally' in use, 
   ! was initiated in europe in the sixteenth century. 
   ! Note that cday is valid only for gregorian calendar dates.
   !
   ! kd = 1 corresponds to january 1, 0000
   !       
   ! Note that the gregorian reform of the julian calendar omitted 10 days
   ! in 1582 in order to restore the date of the vernal equinox to march 21 
   ! (the day after oct 4, 1582 became oct 15, 1582), and revised the leap 
   ! year rule so that centurial years not divisible by 400 were not leap years.
   !
   ! This routine was written by Eugene Neufeld, at ios, in june 1990.
   ! Modified by Joel Chasse (January 20th, 1995)
   !
   !===================================================================72
 
   implicit none
 
   integer, intent(in out) :: idd, imm, iyy_pass
 
   integer, intent(out) :: kd
 
   integer :: i, icc, iyy, jcc, jfh, jfy, jyy, jyyy, kk, kkd, l, lp
 
   integer, dimension(13) :: ndp
 
   integer, dimension(12) :: ndm

   data ndp /0,31,59,90,120,151,181,212,243,273,304,334,365/

   data ndm /31,28,31,30,31,30,31,31,30,31,30,31/
 

   icc = int(iyy_pass/100.)

   iyy = iyy_pass-icc*100

   lp  = 6

   ! Test for invalid input:

1  format(' input error. icc = ',i7)

2  format(' input error. iyy = ',i7)

3  format(' input error. imm = ',i7)

4  format(' input error. idd = ',i7)
 
   if(icc<0)then
      write(lp,1) icc
      stop
   endif
 
   if(iyy<0.or.iyy>99)then
      write(lp,2) iyy
     stop
   endif
 
   if(imm<=0.or.imm>12)then
      write(lp,3) imm
      stop
   endif

   if(idd<=0)then
      write(lp,4) idd
      stop
   endif
 
   if(imm/=2.and.idd>ndm(imm))then
      write(lp,4) idd
      stop
   endif
 
   if(imm==2.and.idd>29)then
      write(lp,4) idd
      stop
   endif
 
   if(imm==2.and.idd>28.and.((iyy/4)*4-iyy/=0.or. &
     (iyy==0.and.(icc/4)*4-icc/=0))) then
      write(lp,4) idd
      stop
   endif
 
   ! Calculate day# of last day of last century:
   kd = icc*36524 + (icc+3)/4
 
   ! Calculate day# of last day of last year:
   kd = kd + iyy*365 + (iyy+3)/4
 
   ! Adjust for century rule:
   !  viz. no leap-years on centurys except when 
   !  the 2-digit century is divisible by 4.
   if(iyy>0.and.(icc-(icc/4)*4)/=0) kd = kd-1
 
   ! kd now truly represents the day# of the last day of last year.
 
   ! Calculate day# of last day of last month:
   kd = kd + ndp(imm)
 
   ! Adjust for leap years:
   if( imm>2           .and. &
     ((iyy/4)*4-iyy)==0.and. &
     ((iyy/=0)         .or.  &
     (((icc/4)*4-icc)==0)) ) kd = kd+1

   ! kd now truly represents the day# of the last day of the last month.

   ! Calculate the current day#:
   kd = kd + idd

   return

   !-------------------------------------------------------------------72

   entry dmy2(idd,imm,iyy_pass,kd)

   !===================================================================72
   !
   ! Given the (gregorian) day#, kd, as calculated above in this routine,
   ! entry dmy2 returns the (gregorian) day, month, year and century.
   !
   !===================================================================72

   ! Test for valid input:

5  format(' kd = ',i7,' invalid input. dmy2 stop.')

   if(kd<=0) write(lp,5) kd

   ! Save kd
   kkd = kd
   ! Calculate icc and subtract the number of days represented by icc from kkd
   ! jfh is the number of 400 year intervals up to kkd
   ! jcc is the number of additional centuries up to kkd

   jfh = kkd/146097

   kkd = kkd - jfh*146097

   if(kkd<36525)then

      jcc = 0

   else

      kkd = kkd - 36525

      jcc = 1 + kkd/36524

      kkd = kkd - (jcc-1)*36524

   endif

   icc = 4*jfh + jcc

   if(kkd==0)then

      icc = icc-1

      iyy = 99

      imm = 12

      idd = 31

      iyy_pass = icc*100+iyy

      return
   endif

   ! Calculate iyy. jfy is the number of four year intervals in the current century
   ! The first four year interval is short (1460 days rather than 1461) if 
   ! the current century is not divisible by 4, and in this case jcc.ne.0 as
   ! calculated above.

   ! Calculate jfy:

   jfy = 0

   if(jcc==0) goto 10

   if(kkd<1460) goto 10

   jfy = 1

   kkd = kkd - 1460

10 kk  = kkd/1461

   jfy = jfy + kk

   kkd = kkd - kk*1461

   ! Calculate jyy, the remaining years of the current century up to the current day:

   jyy = 0

   ! The next year is not a leap year if jfy=0 and jcc.ne.0.

   if(jfy==0.and.jcc/=0) goto 20

   if(kkd<366) goto 30

   jyy  = 1
   kkd  = kkd - 366

20 jyyy = kkd/365

   jyy  = jyy + jyyy

   kkd  = kkd - jyyy*365

30 iyy  = 4*jfy + jyy

   if(kkd==0) then

      iyy = iyy-1

      imm = 12

      idd = 31

      iyy_pass = icc*100+iyy

      return

   endif

   ! Set l=1 if we have a leap year.

   l = 0

   if(iyy-(iyy/4)*4/=0) goto 40

   if(iyy==0.and.(icc-(icc/4)*4)/=0) goto 40

   l = 1

   ! Calculate imm and idd
40 if(kkd>31) goto 50
   imm = 1

   idd = kkd

   iyy_pass = icc*100+iyy

   return

50 if(kkd>59) goto 60

   imm = 2

   idd = kkd-31

   iyy_pass = icc*100+iyy

   return

60 if(kkd>60) goto 70

   if(l==0) goto 70

   imm = 2

   idd = 29

   iyy_pass = icc*100+iyy

   return

70 if(l==1) kkd = kkd-1

   do i = 4,13

      if(kkd>ndp(i)) cycle

      imm = i-1
      idd = kkd - ndp(i-1)

      iyy_pass = icc*100+iyy

      return

   enddo

6  format(' error in dmy2')

   write(lp,6)

   stop

end subroutine cday2

!--------------------------------------------------------------------72
