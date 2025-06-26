module lagrangian_motion_mod

  use partClass_mod
  use partSubClass_mod

  use netcdf      ! netCDF file format module
  use rng_par_zig ! random generator module

#ifdef PAR
  use omp_lib
#endif

  implicit none 

contains
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!72

  real function forcing_value(forcing_type, x_grid, y_grid, z_grid)

    character(len=15), intent(in) :: forcing_type
    integer,           intent(in) :: x_grid, y_grid, z_grid

    if (forcing_type == 'temperature') then

        forcing_value = temp(x_grid, y_grid, z_grid)

    endif

    return

  end function forcing_value
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!72
  
  real function OU(theta, mu, sigma, dt, var)

    ! OU(theta,mean,sigma,dt,var) returns a value following an Ornstein-Uhlenbeck
    ! stochastic process.
    !
    ! Input:
    !   theta(>0)  : rate of mean reversion ~ swimming speed   [time^-1] 
    !   mu         : long-term mean of the process             [unit]
    !   sigma(>0)  : volatility = average magnitude of the random 
    !                fluctuations modelled as Brownian motions [time^-0.5]
    !   dt         : time-step                                 [time]
    !   var        : value of the variable at time t           [unit]
    !
    ! Output:
    !   var        : value of the variable updated at t+dt     [unit]
    !
    ! http://planetmath.org/encyclopedia/OrnsteinUhlenbeckProcess.html
    !
    ! MapsF 2011

    integer, intent(in) :: dt
    real,    intent(in) :: theta, mu, sigma, var
    real                :: mean, sd, rn


    mean = mu + (var-mu) * exp(-theta*dt)

    sd   = sigma * sqrt( (1-exp(-2*theta*dt)) / 2 / theta )

    call rng_norm(rn)
    OU   = rn * sd + mean

    return

  end function OU

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!72
  
  real function vert_mig(z, zday, znight, bottom, theta, sigma, dvm)

    ! Return depth of particle according to a diel cycle
    ! Includes a random component (Ornstein-Uhlenbeck distribution)
    !
    ! Input:
    !   z        = depth of the particule at the previous timestep
    !   zday     = target depth of the particule during daytime
    !   znight   =       "                 "            nighttime
    !   bottom   = bottom depth
    !   theta    = vertical "swimming" speed coefficient
    !   sigma    = vertical spread coefficient
    !   dvm      = particle migrate (.true.) or not (.false.)
    !
    ! Output:
    !   vert_mig = vertical position of the particle
    !
    ! MapsF 2011


    real,    intent(in) :: z, zday, znight, bottom
    real,    intent(in) :: theta, sigma
    logical, intent(in) :: dvm

    real                :: sig, zopt, r1

    ! Nighttime
    if (.not.light_day .and. dvm) then
       ! Forage close to the surface
       zopt = min(bottom, znight)
    else
       zopt = min(bottom, zday)
    endif

    ! sigma (vertical spread) is function of zopt - bottom
    ! 0.2 < f(zopt-bottom) < 1.0
    ! with f(<0)  = 1
    !      f(50)  = 1/2
    !      f(100) = 1/3
    !      ...
    !      f(inf) = 1/5

    sig = sigma * ( 0.2 + 50. / ( 50 + max(0.,zopt-bottom) ) ) / 1.2

    ! no particles above the surface | below bottom

    !vert_mig = max(0.1, min(bottom, abs( OU(theta,zopt,sig,dt,z) ) ) )
    vert_mig = max(0.1, abs( OU(theta,zopt,sig,dt,z) ) )

    ! if the particle is deeper than the bottom, the it is moved randomly
    ! between bdep and 10% shallower than bdep :
    if (vert_mig > bottom) then
       call rng_uni(r1) !random number between 0 and 1
       vert_mig = bottom - (r1 * (0.2 * bottom) )
    endif

    return

  end function vert_mig

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!72
  
  real function update(A3d, nm, nn, nl, index, xp, yp, zp)

    !------------------------------------------------------------------72
    !
    ! This routine interpolates linearly 3D array at a given point
    ! (x,y,z) where (x,y) are in grid units (i,k) and z is in meter
    ! (positive downward to keep the depths positives).
    !
    !------------------------------------------------------------------72
    !
    ! index can be u,v,w or s(sigmat). One can use s for temperature and 
    ! salinity. 
    !
    !------------------------------------------------------------------72

    integer,intent(in) :: nm, nn, nl
    integer            :: i, k, j, j1, j2, ii, kk

    real,   intent(in) :: xp, yp, zp, A3d(nm,nn,nl)
    real               :: dx, dy, d_z, thick, val1, val2

    character(len=1),intent(in) :: index


    ! Find local coordinates (dx,dy)

    if (index=='w'.or.index=='s') then

       ii = min( max( int(xp+0.5), 1 ), nm-1 )

       kk = min( max( int(yp+0.5), 1 ), nn-1 )

       dx = xp - (float(ii)-0.5)

       dy = yp - (float(kk)-0.5)

    elseif (index=='u') then

       ii = min(max(int(xp),1),nm-1)

       kk = min(max(int(yp+0.5),1),nn-1)

       dx = xp-float(ii)

       dy = yp-(float(kk)-.5)

    elseif (index=='v') then

       ii = min(max(int(xp+0.5),1),nm-1)

       kk = min(max(int(yp),1),nn-1)

       dx = xp-(float(ii)-.5)

       dy = yp-float(kk)

    else

       stop 'Index problems when calling value'

    endif

    update = 0.0

    d_z    = 0.0

    j1     = 1

    j2     = 1

    ! Find local vertical coordinate
    i = min(max(int(xp)+1,1),nm)

    k = min(max(int(yp)+1,1),nn)

    if(nlayer(i,k)==0) then ! below the bottom or value = 0.
       return

    elseif(zp>mid_depth(nlayer(i,k))) then

       return

    endif

    if(index=='w') then ! Fields at the top of each layer(ex: w)

       if(zp<=0.0) then ! If above the surface

          d_z = 0.

          j1  = 1

          j2  = 1

       else

          loop1 : do j = 1,nlayer(i,k) ! Find the bracketing layers

             if(zp>up_depth(j).and.zp<=mid_depth(j)) then

                d_z = (zp-up_depth(j))/(mid_depth(j)-up_depth(j)) ! Thickness fraction

                j1  = j 

                j2  = j+1

                exit loop1  ! Exit the loop to save time

             endif

          enddo loop1

          if(j1==nlayer(i,k)) then

             j2  = j1 ! Constante value in the bottom layer

             d_z = 0.

          endif

       endif

       ! Fields at the center of each layer (ex:u)

    elseif(index=='u'.or.index=='v'.or.index=='s') then

       ! Constant value if above middle of top layer

       if(zp<=mid_depth(1)) then

          d_z = 0.

          j1  = 1

          j2  = 1

          ! Constant value if if below the middle of bottom layer

       elseif (nlayer(i,k)>0.and.zp>mid_depth(nlayer(i,k))) then 

          d_z =  0.

          j1  = nlayer(i,k)

          j2  = nlayer(i,k)

       else 

          loop2 : do j = 1,nlayer(i,k)-1 ! Find the bracketing layers

             if(zp>mid_depth(j).and.zp<=mid_depth(j+1)) then

                thick = (mid_depth(j)-up_depth(j))*0.5+(mid_depth(j+1)-up_depth(j+1))*0.5

                d_z   = (zp-mid_depth(j))/thick ! Thickness fraction

                j1    = j 

                j2    = j+1

                exit loop2

             endif

          enddo loop2

       endif

    endif

    !--- Horizontal interpolation (bilinear)

    !    First horizontal plan

    val1 = dx *     dy * A3d(ii+1,kk+1,j1) + (1.-dx) *     dy * A3d(ii,kk+1,j1) + &
         dx * (1.-dy)* A3d(ii+1,kk,j1)   + (1.-dx) * (1.-dy)* A3d(ii,kk,j1)

    !    Second horizontal plan

    val2 = dx *    dy * A3d(ii+1,kk+1,j2) + (1.-dx) *     dy * A3d(ii,kk+1,j2) + &
         dx *(1.-dy)* A3d(ii+1,kk,j2)   + (1.-dx) * (1.-dy)* A3d(ii,kk,j2)

    !    Vertical interpolation (linear)

    update =  d_z*val2 + (1.-d_z)*val1

    return

  end function update

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!72
  
  real function update2d(A2d,nm,nn,index,xp,yp)

    !------------------------------------------------------------------72
    !
    ! This routine interpolates linearly 2D arraye at a given point
    ! where (x,y) are in grid units (i,k)
    !
    !------------------------------------------------------------------72
    !
    ! index can be u,v,w or s(sigmat). One can use s for temperature and 
    ! salinity. Use w for eta.
    !
    !------------------------------------------------------------------72

    integer,intent(in) :: nm,nn
    integer            :: j,j1,j2,ii,kk

    real,   intent(in) :: xp,yp,A2d(nm,nn)
    real               :: dx,dy,thick,val1,val2

    character(len=1),intent(in) :: index


    ! Find local coordinates (dx,dy)
    if(index=='w'.or.index=='s') then

       ii = min(max(int(xp+0.5),1),nm-1)

       kk = min(max(int(yp+0.5),1),nn-1)

       dx = xp-(float(ii)-.5)

       dy = yp-(float(kk)-.5)

    elseif (index=='u') then

       ii = min(max(int(xp),1),nm-1)

       kk = min(max(int(yp+0.5),1),nn-1)

       dx = xp-float(ii)

       dy = yp-(float(kk)-.5)

    elseif (index=='v') then

       ii = min(max(int(xp+0.5),1),nm-1)

       kk = min(max(int(yp),1),nn-1)

       dx = xp-(float(ii)-.5)

       dy = yp-float(kk)

    endif

    update2d = 0.0

    !--- Horizontal interpolation (bilinear)

    val1 = dx *    dy * A2d(ii+1,kk+1) + (1.-dx) *     dy * A2d(ii,kk+1) + &
         dx *(1.-dy)* A2d(ii+1,kk)   + (1.-dx) * (1.-dy)* A2d(ii,kk)

    update2d = val1

    return

  end function update2d
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!72
  
  subroutine trajectory(mylist, part_spec)

    class(partList)   :: mylist
    class(*), pointer :: curr
    integer, intent(in)  :: part_spec

    integer :: it, i, j, k, ii, kk, dt2

    integer :: kstep_diff, istart, isteps, mode

    real    :: up, vp, wp
    real    :: xp_temp, yp_temp, zp_temp

    real    :: kdiff = 20 !m2/s; ramdom walk

    real    :: rn, xdiffus, ydiffus

    real    :: dayhour

    real    :: xpo, ypo, zpo,      &
               xk1, xk2, xk3, xk4, &
               yk1, yk2, yk3, yk4, &
               zk1, zk2, zk3, zk4, &
               zz,  zzz, sigma_local, zday, znight, bdep 

    logical :: light_now, dvm, cell_out
    integer :: part_status

    ! Pointer array declaration for parallelization
    type parray
       class(Zoo), pointer :: ptr
    end type parray
    
    type(parray), allocatable :: ptrarray(:)

    integer :: nptr, thrid, thrs

    !------------------------------------------------------------------72

    ! Compute hour of day (hh.hhh)
    dayhour = time(3)+time(2)/60

    ! Identify daytime | nighttime        
    if (dayhour>=time_sunrise.and.dayhour<=time_sunset) then
    !if (dayhour>=time_go_down.and.dayhour<=time_go_up) then ! KB

       light_now = .true. ! day

    else

       light_now = .false. ! night

    endif

    ! Fancy output
    if (light_now.eqv.light_day) then
       write(*,'($a)') '-'
    else
       if (light_now) then
#ifdef BACK
          write(*,'($a)') '°'
#else
          write(*,'($a)') '*'
#endif
       else
#ifdef BACK
          write(*,'($a)') '*'
#else
          write(*,'($a)') '°'
#endif
       endif
    endif

    !------------------------------------------------------------------72

    dt2 = 2.*dt

    ! Create pointer array for OMP parallelization

    nptr = mylist%n_part()
    allocate(ptrarray(nptr))

    call mylist%reset()

    i = 1
    do while(mylist%moreValues())

       ! Select parent Zoo class for advection
       select type(curr => mylist%currentValue())
       class is (Zoo)
          ptrarray(i)%ptr => curr
       class is (Krill)
          ptrarray(i)%ptr => curr
       class is (Calanus)
          ptrarray(i)%ptr => curr
       end select

       i = i + 1
       call mylist%next()

    enddo

    call mylist%reset()

#ifdef PAR
#ifdef DEBUG
    !$OMP PARALLEL DO &
    !$OMP PRIVATE(xpo,ypo,zpo,xk1,xk2,xk3,xk4,yk1,yk2,yk3,yk4,zk1,zk2,zk3,zk4) &
    !$OMP PRIVATE(up,vp,wp,xp_temp,yp_temp,zp_temp,xdiffus,ydiffus) &
    !$OMP PRIVATE(dvm,zz,istart,isteps,cell_out,i,k,kstep_diff, part_status) &
    !$OMP PRIVATE(j,thrid,thrs) &
    !$OMP SHARED(nptr,ptrarray)
#else
    !$OMP PARALLEL DO &
    !$OMP PRIVATE(xpo,ypo,zpo,xk1,xk2,xk3,xk4,yk1,yk2,yk3,yk4,zk1,zk2,zk3,zk4) &
    !$OMP PRIVATE(up,vp,wp,xp_temp,yp_temp,zp_temp,xdiffus,ydiffus) &
    !$OMP PRIVATE(dvm,zz,istart,isteps,cell_out,i,k,kstep_diff) &
    !$OMP PRIVATE(j) &
    !$OMP SHARED(nptr,ptrarray)
#endif
#endif
    do j = 1, nptr

#ifdef PAR
#ifdef DEBUG
       thrs  = omp_get_num_threads()
       thrid = omp_get_thread_num()
       print*, ' Number of threads = ', thrs, &
               ';   Thread ID      = ', thrid, &
               ';   DO loop index  = ', j
#endif
#endif
       
       if (part_spec==3)then
          part_status=ptrarray(j)%ptr%get_part_status()
       else
          part_status=1
       endif
       if (part_status==1) then
          xpo = ptrarray(j)%ptr%get_xpo()
          ypo = ptrarray(j)%ptr%get_ypo()
          xk1 = ptrarray(j)%ptr%get_xk1()
          xk2 = ptrarray(j)%ptr%get_xk2()
          xk3 = ptrarray(j)%ptr%get_xk3()
          xk4 = ptrarray(j)%ptr%get_xk4()
          yk1 = ptrarray(j)%ptr%get_yk1()
          yk2 = ptrarray(j)%ptr%get_yk2()
          yk3 = ptrarray(j)%ptr%get_yk3()
          yk4 = ptrarray(j)%ptr%get_yk4()       
          zpo = ptrarray(j)%ptr%get_zpo()
#ifdef VERT
          zk1 = ptrarray(j)%ptr%get_zk1()
          zk2 = ptrarray(j)%ptr%get_zk2()
          zk3 = ptrarray(j)%ptr%get_zk3()
          zk4 = ptrarray(j)%ptr%get_zk4()
#endif

          dvm = ptrarray(j)%ptr%get_dvm()
          zz  = ptrarray(j)%ptr%get_zday()
          zzz = ptrarray(j)%ptr%get_znight()

          istart = ptrarray(j)%ptr%get_istart()
          isteps = ptrarray(j)%ptr%get_isteps()



          i = int(xpo) + 1
          if (i > m-1 .or. i < 1) then 
             cell_out = .true.
             call ptrarray(j)%ptr%set_cell_out(cell_out)
          end if

          k = int(ypo) + 1
          if (k > n-1 .or. k < 1) then 
             cell_out = .true.
             call ptrarray(j)%ptr%set_cell_out(cell_out)
          end if

          if (ptrarray(j)%ptr%get_cell_out() .eqv. .false.) then 

             i = min(m, max(1,i))
             k = min(n, max(1,k))

             kstep_diff = kstep - istart
!             print*, ' isteps ', isteps  

             !------------------------------------------------------------72
             ! Vertical migration updated at each timestep
              
             if (kstep_diff == 0) then 

             ! Initialisation of k1s
#ifdef ADV
                up =  update(u, m, n, ilo, 'u', xpo, ypo, zpo)
                vp =  update(v, m, n, ilo, 'v', xpo, ypo, zpo)

#endif

#ifdef VERT
                wp = -update(w, m, n, ilo, 'w', xpo, ypo, zpo)
#endif

#ifdef ADV
                xk1 = dt2 * up / dlx (i,k) 
                call ptrarray(j)%ptr%set_xk1(xk1)

                yk1 = dt2 * vp / dly (i,k)
                call ptrarray(j)%ptr%set_yk1(yk1)
#endif

#ifdef VERT
                zk1 = dt2 * wp
                call ptrarray(j)%ptr%set_zk1(zk1)
#endif

              elseif (kstep_diff>0.and.kstep_diff<=isteps) then

                if (mod(kstep_diff,2)==0) then 
                ! Avanced time : computation of k4s and  the new positions, then k1s 

#ifdef ADV
                    up =  update(u, m, n, ilo, 'u', xpo + xk3, ypo + yk3, zpo + zk3)
                    vp =  update(v, m, n, ilo, 'v', xpo + xk3, ypo + yk3, zpo + zk3)
#endif

#ifdef VERT
                    wp = -update(w, m, n, ilo, 'w', xpo + xk3, ypo + yk3, zpo + zk3)
#endif

#ifdef ADV
                    xk4 = dt2 * up / dlx (i,k) 
                    call ptrarray(j)%ptr%set_xk4(xk4)

                    yk4 = dt2 * vp / dly (i,k)
                    call ptrarray(j)%ptr%set_yk4(yk4)
#endif

#ifdef VERT
                    zk4 = dt2 * wp
                    call ptrarray(j)%ptr%set_zk4(zk4)
#endif

#ifdef ADV
                    xp_temp = xpo + (xk1 + 2.*xk2 + 2*xk3 + xk4)/6.
                   ! Random walk
                    call rng_norm(rn)
                    xdiffus = rn / dlx(i,k) * ( 2. * kdiff * dt2 ) ** 0.5
                    xp_temp = xp_temp + xdiffus

                    yp_temp = ypo + (yk1 + 2.*yk2 + 2.*yk3 + yk4)/6.
                   ! Random walk
                    call rng_norm(rn)
                    ydiffus = rn / dly(i,k) * ( 2. * kdiff * dt2 ) ** 0.5
                    yp_temp = yp_temp + ydiffus
#endif

                    zp_temp = zpo
#ifdef VERT
                    zp_temp = zpo + (zk1 + 2.*zk2 + 2.*zk3 + zk4)/6.
#endif

                   ! Move wet particles

                    ii = int(xp_temp) + 1 ! Destination cell
                    kk = int(yp_temp) + 1

                    ii = min(m, max(1,ii)) ! Necessary with diffusion to avoid out of bounds
                    kk = min(n, max(1,kk)) 


                    if(nlayer(ii,kk)>0) then 
                       if(zp_temp<mid_depth(nlayer(ii,kk))) then ! NL
#ifdef ADV
                          xpo = xp_temp
                          call ptrarray(j)%ptr%set_xpo(xpo)
                          ypo = yp_temp
                          call ptrarray(j)%ptr%set_ypo(ypo)
#endif

#ifdef VERT
                          ! vertical advection
                          zpo = max(0., zp_temp)
#endif

                        endif
                     endif
#ifdef ADV
                     up  =  update(u, m, n, ilo, 'u', xpo, ypo, zpo)
                     vp  =  update(v, m, n, ilo, 'v', xpo, ypo, zpo)

#endif

#ifdef VERT
                     wp  = -update(w, m, n, ilo, 'w', xpo, ypo, zpo)
#endif

#ifdef ADV
                     xk1 = dt2 * up / dlx (i,k)
                     call ptrarray(j)%ptr%set_xk1(xk1)

                     yk1 = dt2 * vp / dly (i,k)
                     call ptrarray(j)%ptr%set_yk1(yk1)
#endif

#ifdef VERT
                     zk1 = dt2 * wp
                     call ptrarray(j)%ptr%set_zk1(zk1)
#endif

                else 
                  ! Half time : computation of the k2s and k3s
#ifdef ADV
                     up  =  update(u, m, n, ilo, 'u', xpo+0.5*xk1, ypo+0.5*yk1, zpo+0.5*zk1)
                     vp  =  update(v, m, n, ilo, 'v', xpo+0.5*xk1, ypo+0.5*yk1, zpo+0.5*zk1)
#endif

#ifdef VERT
                     ! WARNING : wp = 0
                     wp  = -update(w, m, n, ilo, 'w',xpo+0.5*xk1, ypo+0.5*yk1, zpo+0.5*zk1)
#endif

#ifdef ADV
                     xk2 = dt2 * up / dlx (i,k)
                     call ptrarray(j)%ptr%set_xk2(xk2)

                     yk2 = dt2 * vp / dly (i,k)
                     call ptrarray(j)%ptr%set_yk2(yk2)
#endif

#ifdef VERT
                     zk2 = dt2 * wp
                     call ptrarray(j)%ptr%set_zk2(zk2)
#endif


#ifdef ADV
                     up  =  update(u, m, n, ilo, 'u', xpo+0.5*xk2, ypo+0.5*yk2, zpo+0.5*zk2)
                     vp  =  update(v, m, n, ilo, 'v', xpo+0.5*xk2, ypo+0.5*yk2, zpo+0.5*zk2)
#endif

#ifdef VERT
                     ! WARNING : wp = 0
                     wp  = -update(w, m, n, ilo, 'w',xpo+0.5*xk2, ypo+0.5*yk2, zpo+0.5*zk2)
#endif

#ifdef ADV
                     xk3 = dt2 * up / dlx (i,k)
                     call ptrarray(j)%ptr%set_xk3(xk3)

                     yk3 = dt2 * vp / dly (i,k)
                     call ptrarray(j)%ptr%set_yk3(yk3)
#endif

#ifdef VERT
                     zk3 = dt2 * wp
                     call ptrarray(j)%ptr%set_zk3(zk3)
#endif

                endif

              endif  

              ! vertical migration



              if (part_spec > 4) then
                 mode         = ptrarray(j)%ptr%get_mode()
                 sigma_local  = ptrarray(j)%ptr%get_sigma_local()
              endif


              if (nlayer(i,k) == 0) then

                 if (part_spec > 4) then
                   zpo = vert_mig(zpo, zz, zzz, 0., &
                               theta, sigma_local, dvm)
                
                 else
                  zpo = vert_mig(zpo, zz, zzz, 0., &
                              theta, sigma, dvm)

                 endif

              else
                 if (part_spec >4) then
                   zpo = vert_mig(zpo, zz, zzz, mid_depth(nlayer(i,k)), &
                               theta, sigma_local, dvm)
                                 
                 else
                   zpo = vert_mig(zpo, zz, zzz, mid_depth(nlayer(i,k)), &
                               theta, sigma, dvm)

                 endif

              endif

              call ptrarray(j)%ptr%set_zpo(zpo)

              call ptrarray(i)%ptr%set_light_day(light_day)

          endif
       endif

    enddo

#ifdef PAR
    !$OMP END PARALLEL DO
#endif

    deallocate(ptrarray)

    ! Update light_day
    light_day = light_now

  end subroutine trajectory

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!72

end module lagrangian_motion_mod
