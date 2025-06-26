MODULE OPA_GET_XY_LL

  USE netcdf

  IMPLICIT NONE

  PRIVATE

  ! Internal Stuff
  CHARACTER(LEN=100) :: GRID_FILE,BATHY_FILE,COORD_FILE 

  ! Local declarations
  INTEGER :: id, status, variable_id

  INTEGER,DIMENSION(2) :: dimension_ids, dimension_sizes

  INTEGER :: i,j, k, kk ! loop dummy index

  REAL :: dist,distm

  LOGICAL :: first = .true.

  ! Public stuff accessable to outside world

  INTEGER, PUBLIC :: ni, nj, nk, nt ! dimension of GRID_FILE

  REAL,DIMENSION(:,:),ALLOCATABLE, PUBLIC :: nav_lon,nav_lat

  REAL,DIMENSION(:,:),ALLOCATABLE, PUBLIC :: e1t,e2t 

  REAL, PUBLIC :: lon00,lat00,ri00,rj00   

  REAL :: delta_lat,delta_lon

  REAL :: range_lat,range_lon, slope,angle

  INTEGER :: ip,jp

  INTEGER, PUBLIC :: i00, j00

  PUBLIC ::  opa_proj, opa_proj_inv

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
CONTAINS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE opa_proj

    IMPLICIT NONE

    REAL :: lonBL,latBL,lonBR,latBR,lonTL,latTL

    REAL :: r1,r2,Dx,Dy,dx1,dy1       

    INTEGER :: i1,i2,j1,j2

    INTEGER,parameter :: ijr = 15 ! Search radius around guessed point: Might need to be increase
    ! if the grid is too irregular

    if (first) then 
       call opa_proj_init
       first = .false.

       !-------------------------------------------------------------------------------------------
       !JC: will be used to Guess postion on assumption that grid is nearly rectangular
       !JC: This will work only for small rotation of the grird: Solution: Build a general case
       range_lon = nav_lon(ni,nint(nj/2.))-nav_lon(1,nint(nj/2.))  !Axe in middle of grid.
       range_lat = nav_lat(ni,nint(nj/2.))-nav_lat(1,nint(nj/2.))
       slope = range_lat/range_lon
       angle = atan(slope)*180./3.14159
       print*,'Rotation of grid in the middle (degree) = ', angle
       delta_lon = (nav_lon(ni,nint(nj/2.))-nav_lon(1,nint(nj/2.)))/ni
       delta_lat = (nav_lat(nint(ni/2.),nj)-nav_lat(nint(ni/2.),1))/nj
       !       print*, 'delta_lon,delta_lat = ',delta_lon,delta_lat
       !-------------------------------------------------------------------------------------------
    end if

    !JC:!     Find nearest wet i,j to lon00,lat00 
    !JC:       distm = 1.e10  !Large initial distance
    !JC:       i00 = 0
    !JC:       j00 = 0
    !JC:       do i = 1,ni
    !JC:       do j = 1,nj           
    !JC:               dist = geproj_dist(nav_lon(i,j),nav_lat(i,j),lon00,lat00) 
    !JC:               if (dist < distm ) then
    !JC:                  distm = dist
    !JC:                  i00 = i
    !JC:                  j00 = j
    !JC:               endif
    !JC:        enddo
    !JC:        enddo
    !
    !   Guess position
    ip = nint((lon00-nav_lon(1,nint(nj/2.)))/delta_lon)
    jp = nint((lat00-nav_lat(nint(ni/2.),1))/delta_lat)
    !    print*,'lon00,lat00 = ',lon00,lat00
    !    print*,'ip,jp = ',ip,jp

    !   Find nearest wet i,j to lon00,lat00 
    distm = 1.e20  !Large initial distance
    !JC:       i00 = 0
    !JC:       j00 = 0
    !JC:       do i = 1,ni
    !JC:       do j = 1,nj           
    i1  = max(1,ip-ijr)
    i2  = min(ni,ip+ijr)
    j1  = max(1,jp-ijr)
    j2  = min(nj,jp+ijr)
    i00 = ni !Maximum value at top right corner
    j00 = nj !Maximum value at top right corner
    !
    do i = i1,i2
       do j = j1,j2          
          dist = geproj_dist(nav_lon(i,j),nav_lat(i,j),lon00,lat00) 
          if (dist < distm ) then
             distm = dist
             i00 = i
             j00 = j
          endif
       enddo
    enddo

    !Approximating relese location to nearest t-point
    !       ri00 = i00-0.5
    !       rj00 = j00-0.5
    !Exact calculation using bottom left and right locations and triangulation
    lonBL = (nav_lon(i00,j00)+nav_lon(i00-1,j00) & !Bot Left 
           + nav_lon(i00,j00-1) + nav_lon(i00-1,j00-1))/4.
    latBL = (nav_lat(i00,j00)+nav_lat(i00-1,j00) & !Bot Left 
           + nav_lat(i00,j00-1) + nav_lat(i00-1,j00-1))/4.

    lonBR = (nav_lon(i00,j00)+nav_lon(i00+1,j00)  & !Bot Right
           + nav_lon(i00,j00-1) + nav_lon(i00+1,j00-1))/4.

    latBR = (nav_lat(i00,j00)+nav_lat(i00+1,j00)  & !Bot Right
           + nav_lat(i00,j00-1) + nav_lat(i00+1,j00-1))/4.


    !    lonTL = (nav_lon(i00,j00)+nav_lon(i00-1,j00)  & !top Left 
    !           + nav_lon(i00,j00+1) + nav_lon(i00-1,j00+1))/4.
    !       
    !    latTL = (nav_lat(i00,j00)+nav_lat(i00-1,j00)  & !top Left 
    !           + nav_lat(i00,j00+1) + nav_lat(i00-1,j00+1))/4.

    r1 = geproj_dist(lonBL,latBL,lon00,lat00)   
    r2 = geproj_dist(lonBR,latBR,lon00,lat00)          
    Dx = e1t(i00,j00)
    !    Dx = geproj_dist(lonBR,latBR,lonBL,latBL)
    Dy = e2t(i00,j00)
    !    Dy = geproj_dist(lonTL,latTL,lonBL,latBL)

    dx1 = (Dx*Dx+r1*r1 - r2*r2)/2./Dx
    dy1 = sqrt(abs(r1*r1-dx1*dx1))          

    ri00 = i00-1. + dx1/Dx
    rj00 = j00-1. + dy1/Dy
    !JC: Bound value close to boundary
    ri00 = max(0.,ri00)
    ri00 = min(real(ni+1),ri00)
    rj00 = max(0.,rj00)
    rj00 = min(real(nj+1),rj00)

    !    print*,'lonBL,latBL,lonBR,latBR,lonTL,latTL'
    !    print*,lonBL,latBL,lonBR,latBR,lonTL,latTL
    !    print*,'r1,r2,Dx,Dy,dx1,dy1,ri00,rj00'
    !    print*,r1,r2,Dx,Dy,dx1,dy1,ri00,rj00

  END SUBROUTINE opa_proj

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE opa_proj_inv 

    IMPLICIT NONE

    REAL :: lonBL,latBL,lonBR,latBR,lonTL,latTL,lonTR,latTR

    REAL :: r1,r2,Dx,Dy,dx1,dy1,dri00,drj00       

    if (first) then 
       call opa_proj_init
       first = .false.
    end if

    !   Find nearest wet i,j to lon00,lat00 
    i00   = min(ni,max(1,int(ri00)))
    j00   = min(nj,max(1,int(rj00)))

    dri00 = ri00-1.0*(i00)
    drj00 = rj00-1.0*(j00)

    lonBL = nav_lon(i00,j00)
    latBL = nav_lat(i00,j00)
    lonTL = nav_lon(i00,min(nj,j00+1))
    latTL = nav_lat(i00,min(nj,j00+1))
    lonBR = nav_lon(min(ni,i00+1),j00)
    latBR = nav_lat(min(ni,i00+1),j00)
    lonTR = nav_lon(min(ni,i00+1),min(nj,j00+1))
    latTR = nav_lat(min(ni,i00+1),min(nj,j00+1))

    lon00 = lonBL*(1-dri00)*(1-drj00) + lonTL*(1-dri00)*drj00 + &
            lonBR*dri00*(1-drj00) + lonTR*dri00*drj00

    lat00 = latBL*(1-dri00)*(1-drj00) + latTL*(1-dri00)*drj00 + &
            latBR*dri00*(1-drj00) + latTR*dri00*drj00

    return

  END SUBROUTINE opa_proj_inv

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE opa_proj_init 

    IMPLICIT NONE      

#ifdef CIOPS
    COORD_FILE = 'data/GEOMETRY/coordinates_NWA36_crop.nc'

#ifdef DAY_STEP
    COORD_FILE = 'data/GEOMETRY/coords_2017.nc'
#endif

#else
    COORD_FILE = 'data/coordinates.nc'
#endif
    !WRITE(0,*) "COORD_FILE:  ", COORD_FILE
    !print*,'lon00,lat00 = ',lon00,lat00

    dimension_sizes = 0
    !JC:      ! Get dimension sizes for U-grid file
    !JC:      CALL NC_INFO(COORD_FILE, "nav_lon", dimension_sizes, status)
    !-----------------------------------------------------------------------

    ! Open the netCDF file
    status = NF90_OPEN(coord_file, NF90_WRITE, id)
    IF (status /= NF90_NOERR) THEN
       WRITE(*,*) "ERROR: Unable to open file:", coord_file
       WRITE(*,*) "      ", NF90_STRERROR(status),'Stopping'
       stop 
    END IF

    ! Get variable ID number
    status = NF90_INQ_VARID(id, "nav_lon", variable_id)
    IF (status /= NF90_NOERR) THEN
       WRITE(*,*) "ERROR: Unable to read variable ID:", coord_file
       WRITE(*,*) "      ", NF90_STRERROR(status),'Stopping'
       stop
    END IF

    ! Get dimension ID numbers
    status = NF90_INQUIRE_VARIABLE(id, variable_id, dimids = dimension_ids)
    IF (status /= NF90_NOERR) THEN
       WRITE(*,*) "ERROR: Unable to read dimension IDs:", coord_file
       WRITE(*,*) "      ", NF90_STRERROR(status),'Stopping'
       stop
    END IF

    ! Get dimension info
    DO i = 1, SIZE(dimension_ids)
       status = NF90_INQUIRE_DIMENSION(id, dimension_ids(i), len = dimension_sizes(i))
       IF (status /= NF90_NOERR) THEN
          WRITE(*,*) "ERROR: Unable to read dimension lengths:", coord_file
          WRITE(*,*) "      ", NF90_STRERROR(status)
       END IF
    END DO

    !-----------------------------------------------------------------------

    IF (status /= 0) THEN
       WRITE(*,*) "ERROR: Couldn't read grid file"
    END IF

    ni = dimension_sizes(1)
    nj = dimension_sizes(2)

    !write(0,*)'Dimensions: ',(dimension_sizes(kk),kk=1,2)

    !  A.D -- grid file nav_lat
    ALLOCATE(nav_lat(ni,nj))
    ! Get nav_lat  variable ID number
    status = NF90_INQ_VARID(id, "nav_lat", variable_id)
    IF (status /= NF90_NOERR) THEN
       WRITE(*,*) "ERROR: Unable to read variable ID:", coord_file
       WRITE(*,*) "      ", NF90_STRERROR(status),'Stopping'
       stop
    END IF
    status = NF90_GET_VAR(id, variable_id, nav_lat,(/ 1, 1/), (/ ni, nj/))
    !print*,nav_lat(20,20)

    !  A.D -- grid file nav_lon
    ALLOCATE(nav_lon(ni,nj))
    ! Get nav_lat  variable ID number
    status = NF90_INQ_VARID(id, "nav_lon", variable_id)
    IF (status /= NF90_NOERR) THEN
       WRITE(*,*) "ERROR: Unable to read variable ID:", coord_file
       WRITE(*,*) "      ", NF90_STRERROR(status),'Stopping'
       stop
    END IF
    status = NF90_GET_VAR(id, variable_id, nav_lon,(/ 1, 1/), (/ ni, nj/))
    !print*,nav_lon(20,20)

    !  A.D -- grid file e1t
    ALLOCATE(e1t(ni,nj)) !dx through t-point
    ! Get nav_lat  variable ID number
    status = NF90_INQ_VARID(id, "e1t", variable_id)
    IF (status /= NF90_NOERR) THEN
       WRITE(*,*) "ERROR: Unable to read variable ID:", coord_file
       WRITE(*,*) "      ", NF90_STRERROR(status),'Stopping'
       stop
    END IF
    status = NF90_GET_VAR(id, variable_id, e1t,(/ 1, 1/), (/ ni, nj/))
    !print*,e1t(20,20)

    !  A.D -- grid file e2t
    ALLOCATE(e2t(ni,nj)) !dx through t-point
    ! Get nav_lat  variable ID number
    status = NF90_INQ_VARID(id, "e2t", variable_id)
    IF (status /= NF90_NOERR) THEN
       WRITE(*,*) "ERROR: Unable to read variable ID:", coord_file
       WRITE(*,*) "      ", NF90_STRERROR(status),'Stopping'
       stop
    END IF
    status = NF90_GET_VAR(id, variable_id, e2t,(/ 1, 1/), (/ ni, nj/))
    !print*,e2t(20,20)
    return

  END SUBROUTINE opa_proj_init

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  REAL FUNCTION geproj_dist(lon,lat,lon0,lat0) 

    ! % result is in meters
    ! % REARTH = 6.3675E6; 
    ! % REARTH = 6367500; %Webdrgoue???

    REAL :: lon,lat,lon0,lat0

    REAL :: REARTH,pi

    parameter(REARTH = 6378137.0,pi = 3.1415)! m_ll2dist

    REAL ::  dLat,dLon,sdLat,sdLon,EE,cc,clat0

    !   % this is an approximation that assumes region is flat enough for
    !   % pathagoras to work. Use Haversine formula
    !   %y = (lat-lat0)*REARTH*pi/180;
    !   %x = (lon-lon0)*REARTH*(pi/180) .*  cos(lat*pi/180);    
    !Haversine formula
    dLat  =  (lat-lat0)*pi/180;
    dLon  =  (lon-lon0)*pi/180 ; 
    sdLat = sin(dLat/2)
    sdLon = sin(dLon/2)   
    clat0 = cos(lat0 *pi/180 )
    EE    =  sdLat*sdLat +  clat0*clat0*sdLon*sdLon 
    cc    =  2 * atan2(sqrt(EE),sqrt(1-EE)); 
    geproj_dist = REARTH * cc;

  END FUNCTION geproj_dist

END MODULE OPA_GET_XY_LL
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
