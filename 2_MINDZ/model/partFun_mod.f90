
module partFun_mod

  use netcdf      ! netCDF file format module
  use rng_par_zig ! random generator module

  use partClass_mod
  use partSubClass_mod

  use lagrangian_motion_mod

#ifdef PAR
  use omp_lib
#endif

  implicit none 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!72

contains

  subroutine evolveKrill(kriList)
    class(partList)           :: kriList
    class(*), pointer         :: curr

    ! Pointer array declaration for parallelization (critical) !
    type parray
       class(Krill), pointer  :: ptr
    end type parray

    type(parray), allocatable :: ptrarray(:)

    integer                   :: i, j, ii, kk
    integer                   :: nptr, thrid, thrs

    real                      :: xpo, ypo, zpo, tp, zday, znight

    ! Create pointer array for OMP parallelization
    nptr = kriList%n_part()
    allocate(ptrarray(nptr))

    call kriList%reset()
    i = 0
    do while(kriList%moreValues())
       select type(curr => kriList%currentValue())
       class is (Krill)
          i = i + 1
          ptrarray(i)%ptr => curr
       end select       
       call kriList%next()
    enddo
    call kriList%reset()

#ifdef PAR
    !$OMP PARALLEL DO &
    !$OMP PRIVATE(i,j,xpo,ypo,zpo,tp) &
    !$OMP SHARED(nptr,ptrarray)
#endif
    do j = 1, i

       if( ptrarray(j)%ptr%get_cell_out() .eqv. .false.) then  

          !------------------------------------------------------------72
          ! Spatial interpolation of forcing fields
          
          xpo = ptrarray(j)%ptr%get_xpo()
          ii  = int(xpo) + 1
          
          ypo = ptrarray(j)%ptr%get_ypo()
          kk  = int(ypo) + 1
          
          zpo = ptrarray(j)%ptr%get_zpo()
          
          tp  = update(temp, m, n, ilo, 's', xpo, ypo, zpo)
          call ptrarray(j)%ptr%set_tp(tp)
          
          call ptrarray(j)%ptr%grow()
          
          call ptrarray(j)%ptr%molt()

          !------------------------------------------------------------72
          ! Diel vertical migrations
          
          zday = ptrarray(j)%ptr%zday_fun(salt2(ii,kk))
          call ptrarray(j)%ptr%set_zday(zday)

          znight = ptrarray(j)%ptr%znight_fun()
          call ptrarray(j)%ptr%set_znight(znight)

       end if

    enddo
#ifdef PAR
    !$OMP END PARALLEL DO
#endif

    deallocate(ptrarray)

  end subroutine evolveKrill

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!72

  subroutine evolveCalanus(calList, runList)
    class(partList)            :: calList
    class(*), pointer          :: curr
    character(*), intent(in)   :: runList

    
    ! Pointer array declaration for parallelization (critical) !
    type parray
       class(Calanus), pointer :: ptr
    end type parray

    type(parray), allocatable  :: ptrarray(:)

    integer                    :: i, j, ii, kk, burn_num, jj, cc
    integer                    :: nptr, thrid, thrs, sex, ioerr
    real                       :: rds,rdvm !random number for calanus vertical behaviour
    real                       :: xpo, ypo, zpo, tp, zday, znight, fd, ph, new_egg
    real,dimension(1:9)        :: Fs, Ts
    real                       :: mass, abund, age !, massp
    integer                    :: serial, mater, birthdate
    real                       :: we
    integer                    :: new_serial  !, new_birthdate 
    logical                    :: light_day, cell_out 
    real                       :: xk1,xk2,xk3,xk4,yk1,yk2,yk3,yk4,zk1,zk2,zk3,zk4
    integer                    :: activity

    open(10, file = runList, status = 'old', action = 'read', iostat = ioerr)
    read(10, nml = calanus_nml)
    close(10)
    we=paramosome(13)
    ! Create pointer array for OMP parallelization
    nptr = calList%n_part()
    allocate(ptrarray(nptr))

    call calList%reset()
    i = 0
    do while(calList%moreValues())
       select type(curr => calList%currentValue())
       class is (Calanus)
          i = i + 1
          ptrarray(i)%ptr => curr
       end select
       call calList%next()
    enddo
    call calList%reset()

#ifdef PAR
    !$OMP PARALLEL DO &
    !$OMP PRIVATE(i,j,xpo,ypo,zpo,tp,fd,jj, xk1, xk2, xk3, xk4,yk1,yk2,yk3,yk4,zk1,zk2,zk3,zk4,light_day,ph) &
    !$OMP SHARED(nptr,ptrarray)
#endif
    do j = 1, i
       xpo = ptrarray(j)%ptr%get_xpo()
       ii  = int(xpo) + 1
       !if (ii<1) then
       !    write(*,'(/a,i4)') 'ii ', ii
       !endif
       ypo = ptrarray(j)%ptr%get_ypo()
       kk  = int(ypo) + 1

       if (ii> m-1 .or. ii < 1)then
           cell_out=.true.
           call ptrarray(j)%ptr%set_cell_out(cell_out)
       endif
       if (kk> n-1 .or. kk < 1)then
           cell_out=.true.
           call ptrarray(j)%ptr%set_cell_out(cell_out)
       endif
       if (ptrarray(j)%ptr%get_cell_out() .eqv. .true.) then
           call ptrarray(j)%ptr%set_part_status(0)
       endif
               
       ! only calculate those particles has status =1 to save computation time
       if( (ptrarray(j)%ptr%get_cell_out() .eqv. .false.) .and. (ptrarray(j)%ptr%get_part_status()==1)) then  

          !------------------------------------------------------------72
          ! Spatial interpolation of forcing fields
          
          
          zpo = ptrarray(j)%ptr%get_zpo()
  
          xk1 = ptrarray(j)%ptr%get_xk1()
          xk2 = ptrarray(j)%ptr%get_xk2()
          xk3 = ptrarray(j)%ptr%get_xk3()
          xk4 = ptrarray(j)%ptr%get_xk4()
          yk1 = ptrarray(j)%ptr%get_yk1()
          yk2 = ptrarray(j)%ptr%get_yk2()
          yk3 = ptrarray(j)%ptr%get_yk3()
          yk4 = ptrarray(j)%ptr%get_yk4()       
          zk1 = ptrarray(j)%ptr%get_zk1()
          zk2 = ptrarray(j)%ptr%get_zk2()
          zk3 = ptrarray(j)%ptr%get_zk3()
          zk4 = ptrarray(j)%ptr%get_zk4()
          light_day = ptrarray(j)%ptr%get_light_day()
          activity  = ptrarray(j)%ptr%get_activity()
                
          tp  = update(temp, m, n, ilo, 's', xpo, ypo, zpo)
  
#ifdef YUMY  
          fd  = update(food, m, n, ilo, 's', xpo, ypo, zpo)
  
          ph  = update(phyto3d, m, n, ilo, 's', xpo, ypo, zpo)    

#endif
          if( nlayer(ii,kk) .ne. 0) then
#ifdef YUMY
              Fs  = food(ii,kk,1:9)
#endif
              Ts  = temp(ii,kk,1:9)   
          endif
  
          call ptrarray(j)%ptr%set_pre_activity(activity)

          call ptrarray(j)%ptr%set_tp(tp)
  
          call ptrarray(j)%ptr%set_fd(fd)

          call ptrarray(j)%ptr%set_part_status(1)
          
          !birthdate = time(4) + time(5)*100 + time(6)*10000  !birthdate = d_Astart + m_Astart*100 + y_Astart*10000
          !new_serial = j   
#ifdef DEBUG
          if (j .eq. 2) then
              write(*,'(/a,i8)') 'Check: d_Astart + m_Astart*100 + y_Astart*10000 = ', birthdate
              write(*,'(/a,i8)') 'new_serial : ', new_serial
          endif
#endif
          !call ptrarray(j)%ptr%set_serial(new_serial)
          
          !call ptrarray(j)%ptr%set_birthdate(birthdate)
          
          call ptrarray(j)%ptr%grow()
          
          call ptrarray(j)%ptr%molt()
  
          call ptrarray(j)%ptr%set_ph(ph)
          
          call ptrarray(j)%ptr%diapause_entry
          
          call ptrarray(j)%ptr%diapause_exit
          
          call ptrarray(j)%ptr%egg_pro
          !------------------------------------------------------------72
          ! Diel vertical migrations for calanus
          
          rdvm= rand()
          !if (nlayer(ii,kk) .ne. 0) then
              zday = ptrarray(j)%ptr%zday_fun(Fs,Ts,mid_depth(nlayer(ii,kk)),rdvm)
          !endif
          call ptrarray(j)%ptr%set_zday(zday)

          !if (nlayer(ii,kk) .ne. 0) then
              znight = ptrarray(j)%ptr%znight_fun(Fs,Ts,mid_depth(nlayer(ii,kk)))
          !endif
          call ptrarray(j)%ptr%set_znight(znight)
  
           !!!!!!!!!!!!!!!!!!!!!!add new particles!!!!!!!!!!!!!!!!!!!!!!!!!
          abund = ptrarray(j)%ptr%get_abund()
          !new_egg=ptrarray(j)%ptr%get_EP()* abund
          new_egg=ptrarray(j)%ptr%get_EP_CLUTCH()*abund
          serial= ptrarray(j)%ptr%get_serial()
          mater = ptrarray(j)%ptr%get_mater()
 
          ! set the values of new created particles the same as the particle produced them
          if (new_egg>0)then
             jj=next_available_part
             call ptrarray(jj)%ptr%set_xpo(xpo)
             call ptrarray(jj)%ptr%set_ypo(ypo)
             call ptrarray(jj)%ptr%set_zpo(zpo)
             call ptrarray(jj)%ptr%set_xk1(xk1)
             call ptrarray(jj)%ptr%set_xk2(xk2)
             call ptrarray(jj)%ptr%set_xk3(xk3)
             call ptrarray(jj)%ptr%set_xk4(xk4)
             call ptrarray(jj)%ptr%set_yk1(yk1)
             call ptrarray(jj)%ptr%set_yk2(yk2)
             call ptrarray(jj)%ptr%set_yk3(yk3)
             call ptrarray(jj)%ptr%set_yk4(yk4)
             call ptrarray(jj)%ptr%set_zk1(zk1)
             call ptrarray(jj)%ptr%set_zk2(zk2)
             call ptrarray(jj)%ptr%set_zk3(zk3)
             call ptrarray(jj)%ptr%set_zk4(zk4)
             call ptrarray(jj)%ptr%set_zday(zday)
             call ptrarray(jj)%ptr%set_znight(znight)

            ! CALANUS attributes
             call ptrarray(jj)%ptr%set_tp(tp)
             call ptrarray(jj)%ptr%set_fd(fd)
             call ptrarray(jj)%ptr%set_ph(ph)
             call ptrarray(jj)%ptr%set_sex(0)
             call ptrarray(jj)%ptr%set_mass(we)
             call ptrarray(jj)%ptr%set_lipid(0.0)
             call ptrarray(jj)%ptr%set_stage(0.0)
             call ptrarray(jj)%ptr%set_activity(0)
            ! newborn egg's lipid, stage, activity should be 0
             call ptrarray(jj)%ptr%set_abund(new_egg)
             call ptrarray(jj)%ptr%set_part_status(1)  
             call ptrarray(jj)%ptr%set_lipid_floor(0.0)
             call ptrarray(jj)%ptr%set_ep_daily(0.0)
             call ptrarray(jj)%ptr%set_EP(0)
             call ptrarray(jj)%ptr%set_egg_clutch(0)
             call ptrarray(jj)%ptr%set_EP_CLUTCH(0)
             call ptrarray(jj)%ptr%set_start_ep(0.0)
             call ptrarray(jj)%ptr%set_reproduce_info([0.0,0.0,0.0,0.0,0.0,0.0])
             call ptrarray(jj)%ptr%set_duration_ep(0.0)
             call ptrarray(jj)%ptr%set_counter_eggs(0.0)
             call ptrarray(jj)%ptr%set_cum_eggs(0)
             call ptrarray(jj)%ptr%set_max_diap_lipid(0.0)
             call ptrarray(jj)%ptr%set_lipid_exit_diap(0.0)
             call ptrarray(jj)%ptr%set_moult_mass_limit(0.0)
             call ptrarray(jj)%ptr%set_light_day(light_day)
             call ptrarray(jj)%ptr%set_a_molt(a_devlp(0.0))
             call ptrarray(jj)%ptr%set_b_molt(t_devlp)  
             call ptrarray(jj)%ptr%set_age(0.0)  
             call ptrarray(jj)%ptr%set_massp(0.0)

             !new_birthdate = time(4) + time(5)*100 + time(6)*10000
             birthdate = time(4) + time(5)*100 + time(6)*10000
             new_serial = jj       
#ifdef DEBUG
            if (jj .eq. 51) then
                write(*,'(/a,i8)') 'Check: birthdate(jj is 51) = ', birthdate
                write(*,'(/a,i8)') 'new_serial : ', new_serial
                write(*,'(/a,i2,a,i4,a,i8)') 'time(4), time(5)*100, time(6)*10000: ',time(4),'-',time(5)*100,'-',time(6)*10000
            endif
#endif
             call ptrarray(jj)%ptr%set_serial(new_serial)       ! next available part #
             call ptrarray(jj)%ptr%set_birthdate(birthdate) ! yyyymmdd
             call ptrarray(jj)%ptr%set_mater(serial)            ! serial # of mother
             next_available_part=next_available_part+1
          endif
  
          !!! calculate mortality
  
          call ptrarray(j)%ptr%starve
          !call ptrarray(j)%ptr%predation(salt2(ii,kk), zpo, parsurf(ii,kk), temp2(ii,kk), phyto(ii,kk))
          !call ptrarray(j)%ptr%mort_maps2012
          call ptrarray(j)%ptr%mort_maps2012(zpo, parsurf(ii,kk))
          !!!!!!!!!!!!!!!!!!!!remove 0 abundance particles!!!!!!!!!!!!!!!!!!!!!!!!
          abund = ptrarray(j)%ptr%get_abund()  
          if (abund<1.0) then
             call ptrarray(j)%ptr%set_abund(0.0)
             call ptrarray(j)%ptr%set_part_status(0)
          endif
       end if
    enddo
#ifdef PAR
    !$OMP END PARALLEL DO
#endif

    deallocate(ptrarray)

  end subroutine evolveCalanus

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!72




  subroutine evolveLate_hyperboreus(LHList, in_diapause)
    class(partList)            :: LHList
    class(*), pointer          :: curr
    
    ! Pointer array declaration for parallelization (critical) !
    type parray
       class(Late_hyperboreus), pointer :: ptr
    end type parray

    type(parray), allocatable  :: ptrarray(:)

    integer, intent(in)        :: in_diapause
    integer                    :: nptr
    integer                    :: i, j, ii, kk, jj
    real                       :: xpo, ypo, zpo, bdep, zday, znight, sigma_local
    logical                    :: cell_out 
    real                       :: xk1,xk2,xk3,xk4,yk1,yk2,yk3,yk4,zk1,zk2,zk3,zk4
    integer                    :: mode



    nptr = LHList%n_part()
    allocate(ptrarray(nptr))

    call LHList%reset()
    i = 0
    do while(LHList%moreValues())
       select type(curr => LHList%currentValue())
       class is (Late_hyperboreus)
          i = i + 1
          ptrarray(i)%ptr => curr
       end select
       call LHList%next()
    enddo
    call LHList%reset()


    do j = 1, i
       xpo = ptrarray(j)%ptr%get_xpo()
       ii  = int(xpo) + 1
       !if (ii<1) then
       !    write(*,'(/a,i4)') 'ii ', ii
       !endif
       ypo = ptrarray(j)%ptr%get_ypo()
       kk  = int(ypo) + 1

       if (ii> m-1 .or. ii < 1)then
           cell_out=.true.
           call ptrarray(j)%ptr%set_cell_out(cell_out)
       endif
       if (kk> n-1 .or. kk < 1)then
           cell_out=.true.
           call ptrarray(j)%ptr%set_cell_out(cell_out)
       endif
       if (ptrarray(j)%ptr%get_cell_out() .eqv. .true.) then
           call ptrarray(j)%ptr%set_part_status(0)
       endif
       


       if( (ptrarray(j)%ptr%get_cell_out() .eqv. .false.) ) then  

          !------------------------------------------------------------72
          ! Diel vertical migrations for calanus
        
          zpo = ptrarray(j)%ptr%get_zpo()

! little part Not ready for now, aims to actualise the mode associated to a particle.
! It is important if the mode can change, depending on the time or the bottom depth for exemple.
! It is not the case here for now.

!            mode_prop(i) = LHPart%mode_prop_fun( LHPart%mode_prop, mid_depth(nlayer(xgrid(i),ygrid(i))), in_diapause )

            !here we determine the difference of mode proportion between t-1 and t
!            dif_prop = mode_prop(i) - LHPart%mode_prop

!            if (dif_prop < 0 .and. LHPart%mode == 1) then !if prop decreases and the particle mode is top
!               call rng_uni(r3) !random number between 0 and 1
!               if ( r3 <= abs(dif_prop) ) then
!                  mode(i) = 2
!               endif

!            elseif (dif_prop > 0 .and. LHPart%mode == 2) then !if prop increases and the particle mode is down
!               call rng_uni(r3) !random number between 0 and 1
!               if ( r3 <= dif_prop ) then
!                  mode(i) = 1
!               endif

!            else
!               mode(i) = LHPart%mode

!            endif

          mode  = ptrarray(j)%ptr%get_mode()
          if (nlayer(ii,kk) .ne. 0) then
              zday = ptrarray(j)%ptr%zday_fun(mid_depth(nlayer(ii,kk)), mode, in_diapause)
          endif
          call ptrarray(j)%ptr%set_zday(zday)

          if (nlayer(ii,kk) .ne. 0) then
              znight = ptrarray(j)%ptr%znight_fun(mid_depth(nlayer(ii,kk)), mode, in_diapause)
          endif
          call ptrarray(j)%ptr%set_znight(znight)
 
          sigma_local = ptrarray(j)%ptr%sigma_fun(mid_depth(nlayer(ii, kk)), mode, in_diapause)
          call ptrarray(j)%ptr%set_sigma_local(sigma_local)

          bdep = mid_depth(nlayer(ii,kk))
          call ptrarray(j)%ptr%set_bdep(bdep)


       endif
    enddo
#ifdef PAR
    !$OMP END PARALLEL DO
#endif

    deallocate(ptrarray)

  end subroutine evolveLate_hyperboreus






  subroutine evolveYoung_hyperboreus(YHList, in_diapause)
    class(partList)            :: YHList
    class(*), pointer          :: curr
    
    ! Pointer array declaration for parallelization (critical) !
    type parray
       class(Young_hyperboreus), pointer :: ptr
    end type parray

    type(parray), allocatable  :: ptrarray(:)

    integer, intent(in)        :: in_diapause
    integer                    :: nptr
    integer                    :: i, j, ii, kk, jj
    real                       :: xpo, ypo, zpo, bdep, zday, znight, sigma_local
    logical                    :: cell_out 
    real                       :: xk1,xk2,xk3,xk4,yk1,yk2,yk3,yk4,zk1,zk2,zk3,zk4
    integer                    :: mode



    nptr = YHList%n_part()
    allocate(ptrarray(nptr))

    call YHList%reset()
    i = 0
    do while(YHList%moreValues())
       select type(curr => YHList%currentValue())
       class is (Young_hyperboreus)
          i = i + 1
          ptrarray(i)%ptr => curr
       end select
       call YHList%next()
    enddo
    call YHList%reset()


    do j = 1, i
       xpo = ptrarray(j)%ptr%get_xpo()
       ii  = int(xpo) + 1
       !if (ii<1) then
       !    write(*,'(/a,i4)') 'ii ', ii
       !endif
       ypo = ptrarray(j)%ptr%get_ypo()
       kk  = int(ypo) + 1

       if (ii> m-1 .or. ii < 1)then
           cell_out=.true.
           call ptrarray(j)%ptr%set_cell_out(cell_out)
       endif
       if (kk> n-1 .or. kk < 1)then
           cell_out=.true.
           call ptrarray(j)%ptr%set_cell_out(cell_out)
       endif
       if (ptrarray(j)%ptr%get_cell_out() .eqv. .true.) then
           call ptrarray(j)%ptr%set_part_status(0)
       endif





       if( (ptrarray(j)%ptr%get_cell_out() .eqv. .false.) ) then  
  
          !------------------------------------------------------------72
          ! Diel vertical migrations for calanus
        
          zpo = ptrarray(j)%ptr%get_zpo()

! little part Not ready for now, aims to actualise the mode associated to a particle.
! It is important if the mode can change, depending on the time or the bottom depth for exemple.
! It is not the case here for now.

!            mode_prop(i) = LHPart%mode_prop_fun( LHPart%mode_prop, mid_depth(nlayer(xgrid(i),ygrid(i))), in_diapause )

            !here we determine the difference of mode proportion between t-1 and t
!            dif_prop = mode_prop(i) - LHPart%mode_prop

!            if (dif_prop < 0 .and. LHPart%mode == 1) then !if prop decreases and the particle mode is top
!               call rng_uni(r3) !random number between 0 and 1
!               if ( r3 <= abs(dif_prop) ) then
!                  mode(i) = 2
!               endif

!            elseif (dif_prop > 0 .and. LHPart%mode == 2) then !if prop increases and the particle mode is down
!               call rng_uni(r3) !random number between 0 and 1
!               if ( r3 <= dif_prop ) then
!                  mode(i) = 1
!               endif

!            else
!               mode(i) = LHPart%mode

!            endif

          mode  = ptrarray(j)%ptr%get_mode()
          if (nlayer(ii,kk) .ne. 0) then
              zday = ptrarray(j)%ptr%zday_fun(mid_depth(nlayer(ii,kk)), in_diapause)
          endif
          call ptrarray(j)%ptr%set_zday(zday)

          if (nlayer(ii,kk) .ne. 0) then
              znight = ptrarray(j)%ptr%znight_fun(mid_depth(nlayer(ii,kk)), in_diapause)
          endif
          call ptrarray(j)%ptr%set_znight(znight)
 
          sigma_local = ptrarray(j)%ptr%sigma_fun(mid_depth(nlayer(ii, kk)), in_diapause)
          call ptrarray(j)%ptr%set_sigma_local(sigma_local)

          bdep = mid_depth(nlayer(ii,kk))
          call ptrarray(j)%ptr%set_bdep(bdep)


       endif
    enddo
#ifdef PAR
    !$OMP END PARALLEL DO
#endif

    deallocate(ptrarray)

  end subroutine evolveYoung_hyperboreus





  subroutine evolveLate_finmarchicus(LFList, in_diapause)
    class(partList)            :: LFList
    class(*), pointer          :: curr
    
    ! Pointer array declaration for parallelization (critical) !
    type parray
       class(Late_finmarchicus), pointer :: ptr
    end type parray

    type(parray), allocatable  :: ptrarray(:)

    integer, intent(in)        :: in_diapause
    integer                    :: nptr
    integer                    :: i, j, ii, kk, jj
    real                       :: xpo, ypo, zpo, bdep, zday, znight, sigma_local
    logical                    :: cell_out 
    real                       :: xk1,xk2,xk3,xk4,yk1,yk2,yk3,yk4,zk1,zk2,zk3,zk4
    integer                    :: mode



    nptr = LFList%n_part()
    allocate(ptrarray(nptr))

    call LFList%reset()
    i = 0
    do while(LFList%moreValues())
       select type(curr => LFList%currentValue())
       class is (Late_finmarchicus)
          i = i + 1
          ptrarray(i)%ptr => curr
       end select
       call LFList%next()
    enddo
    call LFList%reset()


    do j = 1, i
       xpo = ptrarray(j)%ptr%get_xpo()
       ii  = int(xpo) + 1
       !if (ii<1) then
       !    write(*,'(/a,i4)') 'ii ', ii
       !endif
       ypo = ptrarray(j)%ptr%get_ypo()
       kk  = int(ypo) + 1

       if (ii> m-1 .or. ii < 1)then
           cell_out=.true.
           call ptrarray(j)%ptr%set_cell_out(cell_out)
       endif
       if (kk> n-1 .or. kk < 1)then
           cell_out=.true.
           call ptrarray(j)%ptr%set_cell_out(cell_out)
       endif
       if (ptrarray(j)%ptr%get_cell_out() .eqv. .true.) then
           call ptrarray(j)%ptr%set_part_status(0)
       endif
       




       if( (ptrarray(j)%ptr%get_cell_out() .eqv. .false.) ) then  

          !------------------------------------------------------------72
          ! Diel vertical migrations for calanus
        
          zpo = ptrarray(j)%ptr%get_zpo()

! little part Not ready for now, aims to actualise the mode associated to a particle.
! It is important if the mode can change, depending on the time or the bottom depth for exemple.
! It is not the case here for now.

!            mode_prop(i) = LHPart%mode_prop_fun( LHPart%mode_prop, mid_depth(nlayer(xgrid(i),ygrid(i))), in_diapause )

            !here we determine the difference of mode proportion between t-1 and t
!            dif_prop = mode_prop(i) - LHPart%mode_prop

!            if (dif_prop < 0 .and. LHPart%mode == 1) then !if prop decreases and the particle mode is top
!               call rng_uni(r3) !random number between 0 and 1
!               if ( r3 <= abs(dif_prop) ) then
!                  mode(i) = 2
!               endif

!            elseif (dif_prop > 0 .and. LHPart%mode == 2) then !if prop increases and the particle mode is down
!               call rng_uni(r3) !random number between 0 and 1
!               if ( r3 <= dif_prop ) then
!                  mode(i) = 1
!               endif

!            else
!               mode(i) = LHPart%mode

!            endif

          mode  = ptrarray(j)%ptr%get_mode()
          if (nlayer(ii,kk) .ne. 0) then
              zday = ptrarray(j)%ptr%zday_fun(mid_depth(nlayer(ii,kk)), mode, in_diapause)
          endif
          call ptrarray(j)%ptr%set_zday(zday)

          if (nlayer(ii,kk) .ne. 0) then
              znight = ptrarray(j)%ptr%znight_fun(mid_depth(nlayer(ii,kk)), mode, in_diapause)
          endif
          call ptrarray(j)%ptr%set_znight(znight)
 
          sigma_local = ptrarray(j)%ptr%sigma_fun(mid_depth(nlayer(ii, kk)), in_diapause)
          call ptrarray(j)%ptr%set_sigma_local(sigma_local)

          bdep = mid_depth(nlayer(ii,kk))
          call ptrarray(j)%ptr%set_bdep(bdep)


       endif
    enddo
#ifdef PAR
    !$OMP END PARALLEL DO
#endif

    deallocate(ptrarray)

  end subroutine evolveLate_finmarchicus







  
  subroutine evolveYoung_finmarchicus(YFList, in_diapause)
    class(partList)            :: YFList
    class(*), pointer          :: curr
    
    ! Pointer array declaration for parallelization (critical) !
    type parray
       class(Young_finmarchicus), pointer :: ptr
    end type parray

    type(parray), allocatable  :: ptrarray(:)

    integer, intent(in)        :: in_diapause
    integer                    :: nptr
    integer                    :: i, j, ii, kk, jj
    real                       :: xpo, ypo, zpo, bdep, zday, znight, sigma_local
    logical                    :: cell_out 
    real                       :: xk1,xk2,xk3,xk4,yk1,yk2,yk3,yk4,zk1,zk2,zk3,zk4
    integer                    :: mode



    nptr = YFList%n_part()
    allocate(ptrarray(nptr))

    call YFList%reset()
    i = 0
    do while(YFList%moreValues())
       select type(curr => YFList%currentValue())
       class is (Young_finmarchicus)
          i = i + 1
          ptrarray(i)%ptr => curr
       end select
       call YFList%next()
    enddo
    call YFList%reset()


    do j = 1, i
       xpo = ptrarray(j)%ptr%get_xpo()
       ii  = int(xpo) + 1
       !if (ii<1) then
       !    write(*,'(/a,i4)') 'ii ', ii
       !endif
       ypo = ptrarray(j)%ptr%get_ypo()
       kk  = int(ypo) + 1

       if (ii> m-1 .or. ii < 1)then
           cell_out=.true.
           call ptrarray(j)%ptr%set_cell_out(cell_out)
       endif
       if (kk> n-1 .or. kk < 1)then
           cell_out=.true.
           call ptrarray(j)%ptr%set_cell_out(cell_out)
       endif
       if (ptrarray(j)%ptr%get_cell_out() .eqv. .true.) then
           call ptrarray(j)%ptr%set_part_status(0)
       endif
          





       if( (ptrarray(j)%ptr%get_cell_out() .eqv. .false.) ) then  

          !------------------------------------------------------------72
          ! Diel vertical migrations for calanus
        
          zpo = ptrarray(j)%ptr%get_zpo()

! little part Not ready for now, aims to actualise the mode associated to a particle.
! It is important if the mode can change, depending on the time or the bottom depth for exemple.
! It is not the case here for now.

!            mode_prop(i) = LHPart%mode_prop_fun( LHPart%mode_prop, mid_depth(nlayer(xgrid(i),ygrid(i))), in_diapause )

            !here we determine the difference of mode proportion between t-1 and t
!            dif_prop = mode_prop(i) - LHPart%mode_prop

!            if (dif_prop < 0 .and. LHPart%mode == 1) then !if prop decreases and the particle mode is top
!               call rng_uni(r3) !random number between 0 and 1
!               if ( r3 <= abs(dif_prop) ) then
!                  mode(i) = 2
!               endif

!            elseif (dif_prop > 0 .and. LHPart%mode == 2) then !if prop increases and the particle mode is down
!               call rng_uni(r3) !random number between 0 and 1
!               if ( r3 <= dif_prop ) then
!                  mode(i) = 1
!               endif

!            else
!               mode(i) = LHPart%mode

!            endif

          if (nlayer(ii,kk) .ne. 0) then
              zday = ptrarray(j)%ptr%zday_fun(mid_depth(nlayer(ii,kk)), in_diapause)
          endif
          call ptrarray(j)%ptr%set_zday(zday)

          if (nlayer(ii,kk) .ne. 0) then
              znight = ptrarray(j)%ptr%znight_fun(mid_depth(nlayer(ii,kk)), in_diapause)
          endif
          call ptrarray(j)%ptr%set_znight(znight)
 
          sigma_local = ptrarray(j)%ptr%sigma_fun(mid_depth(nlayer(ii, kk)), in_diapause)
          call ptrarray(j)%ptr%set_sigma_local(sigma_local)

          bdep = mid_depth(nlayer(ii,kk))
          call ptrarray(j)%ptr%set_bdep(bdep)


       endif
    enddo
#ifdef PAR
    !$OMP END PARALLEL DO
#endif

    deallocate(ptrarray)

  end subroutine evolveYoung_finmarchicus









!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!72




  subroutine generateZoo(zooList, runList)
    type(Zoo)                    :: zooPart
    class(partList)              :: zooList


    character(*), intent(in)     :: runList

    integer                      :: i, j, ii, k, l, ioerr, t, IOstatus ! counters or iterators

    integer                      :: kd_start, kd_end, kd

    integer, dimension(6)        :: part_start, part_end

    integer, dimension(max_part) :: xgrid, ygrid, zgrid

    real                         :: r1, r2

    real,    dimension(max_part) :: xpo, ypo, zpo, zz, zzz, tp


    !------------------------------------------------------------------72

    open(10,file = runList, status = 'old', action = 'read', err = 1, iostat = ioerr)
1   if (ioerr /= 0) STOP 'problem opening file'
    read(10, nml = part_nml)
    close(10)

    !------------------------------------------------------------------72
    ! Initial particles spatial distribution
    !
    ! Generic case creates one particle per wet cell ito provide 
    ! the # of particles for Lyapunov exponents simulations

    ! nxy marks the cells with particles
                       nxy = 1
    where( nlayer==0 ) nxy = 0

    ! npart0 is the number of cells with particles
    npart0 = min( sum(nxy), max_part )

    ! allocate the vectors used for the cells with particles
    allocate(ipart(npart0))
    allocate(jpart(npart0))

    ! record the coordinates of the cells with paticles
    ii = 0
    outer: do j = 1, n    ! n corresponds to the y horizontal dimension of the grid
       inner :do i = 1, m ! m corresponds to the x horizontal dimension of the grid

          if( nlayer(i,j) .ne. 0 ) then
             ii = min( ii+1, npart0 )
             if( ii .gt. npart0 ) exit outer
             ipart(ii) = i
             jpart(ii) = j
          endif

       end do inner
    end do outer

    write(*,'(/a,i9,a)') ' Cells with particles = ', npart0, '; now generating all the particles variables'


    !------------------------------------------------------------------72
    ! Initialization of individual particles

    ! Initializing time (time would be tzero + one dt)
    call update2(initial)

    ! Initialize time loop for periodic release of particles:
    ! start of particles release -> end of release (by frequency of release)

    part_start = (/ 0, 0, 0, d_Astart, m_Astart, y_Astart /) ! Start of particles release
    call cday2(part_start(4),part_start(5),part_start(6),kd_start)

    part_end   = (/ 0, 0, 0, d_Aend, m_Aend, y_Aend /)       ! End of particles release
    call cday2(part_end(4),part_end(5),part_end(6),kd_end)

#ifdef DEBUG
    write(*,'(/a,i2,1x,i2,1x,i4,a,i6)') ' part start              : ', part_start(4), part_start(5), part_start(6), &
                                        ';   start timestep : ', kd_start
    write(*,'(a,i2,1x,i2,1x,i4,a,i6)')  ' part end                : ', part_end(4),   part_end(5),   part_end(6),   &
                                        ';   end timestep   : ', kd_end
#endif

    ! Duration of particles advection; same for all particles
    isteps_part = part_duration * 86400/dt

#ifdef DEBUG                 
    write(*,'(a,i10,a)') ' part advection duration : ', isteps_part, ' (inner time steps)'
#endif

    ! Periodic initialization
#ifdef BACK
    tpart = ceiling( real(kd_start - kd_end+1)/part_freq)
#else
    tpart = ceiling( real(kd_end - kd_start+1)/part_freq )
#endif    

#ifdef DEBUG                 
    write(*,'(a,i10)') ' part seeding events #   : ', tpart
#endif

    ! NOW we initialize each particle's values

    ! FIRST open sequence of repeated initialization...
    do t = 1, tpart

       i = 0

       ! ... then get random positions on the grid...
       do k = 1, npart0      ! Cell position
          do l = 1, part_fac ! Individual particle within cell

             ! ... and duplicate particles for each sequence of initialization
             i = i + 1

             if (t==1) then

                !initialize position
3               call rng_uni(r1)
                call rng_uni(r2)

                xpo(i) = real(ipart(k)) - 1. + r1
                ypo(i) = real(jpart(k)) - 1. + r2

                xgrid(i) = max(1, min(m, ceiling(xpo(i))))
                ygrid(i) = max(1, min(n, ceiling(ypo(i))))

                do while (nlayer(xgrid(i), ygrid(i)) == 0)
                   goto 3 ! check that the particle is in a wet cell !
                end do

                zz(i)  = zooPart%zday_funZ()
                zzz(i) = zooPart%znight_funZ()
                zpo(i) = vert_mig(0., zz(i), zzz(i), mid_depth(nlayer(xgrid(i), ygrid(i))), &
                                  theta, sigma, dvm)


                ! Loop to find zgrid(i), the depth layer corresponding to the depth zpo(i) :
                zloop: do j = 1, min( ilo, nlayer(xgrid(i),ygrid(i)) )
                   zgrid(i) = j
                   if (zpo(i)<=mid_depth(j)) exit zloop
                enddo zloop

#ifdef BACK
                ! Days from the run start
                kd = kd_ini + istart*dtnc/86400 - (kd_start+(t-1)*part_freq)

                ! Initial particle timestep (kstep); use dt in sec
                istart_part = int( (kd*86400 - part_start(3)*3600 - part_start(2)*60 - part_start(1))/dt )
#else
                ! Days from the run start
                kd = kd_start - kd_ini - istart*dtnc/86400

                ! Initial particle timestep (kstep); use dt in sec
                istart_part = int( (kd*86400 + part_start(3)*3600 + part_start(2)*60 + part_start(1))/dt )
#endif
                tp(i) = temp(xgrid(i),ygrid(i),zgrid(i))
                
             endif

             zooPart%part_spec = part_spec

             zooPart%istart    = istart_part
             zooPart%isteps    = isteps_part

             zooPart%dvm       = dvm
             zooPart%light_day = light_day
             zooPart%zday      = zz(i)
             zooPart%znight    = zzz(i)

             zooPart%cell_out  = cell_out

             zooPart%xpo       = xpo(i)
             zooPart%ypo       = ypo(i)
             zooPart%zpo       = zpo(i)

             zooPart%xk1       = xk1
             zooPart%xk2       = xk2
             zooPart%xk3       = xk3
             zooPart%xk4       = xk4

             zooPart%yk1       = yk1
             zooPart%yk2       = yk2
             zooPart%yk3       = yk3
             zooPart%yk4       = yk4

             zooPart%zk1       = zk1
             zooPart%zk2       = zk2
             zooPart%zk3       = zk3
             zooPart%zk4       = zk4                

             call zooList%addPart(zooPart)

          end do
       end do
    end do

    npart = npart0 * part_fac

    deallocate(ipart,jpart)

    write(*,'(/a,i6,a)') ' Generated ', npart, ' particles'

  end subroutine generateZoo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!72

  subroutine generateKrill(kriList, runList)
    type(Krill)                  :: kriPart
    class(partList)              :: kriList


    character(*), intent(in)     :: runList

    integer                      :: i, j, ii, k, l, ioerr, t, IOstatus

    integer                      :: kd_start, kd_end, kd

    integer, dimension(6)        :: part_start, part_end

    integer, dimension(max_part) :: xgrid, ygrid, zgrid

    real                         :: r1, r2

    real,    dimension(max_part) :: xpo, ypo, zpo, zz, zzz, tp

    real,    dimension(max_part) :: mass0

    real,    dimension(14, 6576) :: pars ! very specific...

    !------------------------------------------------------------------72

    open(10,file = runList, status = 'old', action = 'read', err = 4, iostat = ioerr)
4   if (ioerr /= 0) STOP 'problem opening file'
    read(10, nml = part_nml)
    close(10)

    !------------------------------------------------------------------72
    ! Initial particles distribution (spatial)
    !
    ! Works with input_particles (below) to provide the # of particles
    ! added at specific time steps in case of periodic seeding of particles

    open(10,file = partfile, status = 'old', action = 'read', err = 5, iostat = ioerr)
5   if (ioerr /= 0) STOP 'problem opening file, particle input'

    do i = 1, m
       read(10, *) (part_grid_init(i, j), j = 1, n)
    end do

    npart0  = sum(part_grid_init) !npart0 is the number of cells with particles

    ii = 0

    do j = 1, n    !n is a horizontal dimension of the grid
       do i = 1, m !m is the other horizontal dimension of the grid
          if(part_grid_init(i, j) == 1) then
             ii = min(ii + 1, npart0)
             ipart(ii) = i
             jpart(ii) = j
          endif
       end do
    end do

    close(10)

    ! Number of individuals generated
    write(*,'(/a,i9,a,/)') ' Cells with particles = ', npart0, '; now generating all the particles variables'

    !------------------------------------------------------------------72
    ! Initialization of individual particles

    ! Open namelist that contains physiological parameters  
    open(10, file = runList, status = 'old', action = 'read', err = 4, iostat = ioerr)
    read(10, nml = krill_nml)

    close(10)

    ! Open matrix of individual parameters
    if (INTRA_SPE .eqv. .true.) then

       open(10, file=parsfile, status = 'old', action = 'read', err = 6, iostat = ioerr)
6      if (ioerr/=0) stop '!!! Problem opening parameters matrix !!!'

       do i = 1, nb_pars 
          read(10,*) (pars(i,j), j = 1, npart0)
       end do

       close(10)

    endif

    ! Initializing time (time would be tzero + one dt)
    call update2(initial)

    ! Initialize time loop for periodic release of particles:
    ! start of particles release -> end of release (by frequency of release)

    part_start = (/ 0, 0, 0, d_Astart, m_Astart, y_Astart /) ! Start of particles release
    call cday2(part_start(4),part_start(5),part_start(6),kd_start)

    part_end   = (/ 0, 0, 0, d_Aend, m_Aend, y_Aend /)       ! End of particles release
    call cday2(part_end(4),part_end(5),part_end(6),kd_end)

#ifdef DEBUG
    write(*,'(/a,i2,1x,i2,1x,i4,a,i6)') ' part start              : ', part_start(4), part_start(5), part_start(6), &
         '; start timestep : ', kd_start
    write(*,'(a,i2,1x,i2,1x,i4,a,i6)') ' part end                : ', part_end(4),   part_end(5),   part_end(6),   &
         '; end timestep   : ', kd_end
#endif

    ! Duration of particles advection; same for all particles
    isteps_part = part_duration * 86400/dt

#ifdef DEBUG                 
    write(*,'(a,i10,a)') ' part advection duration : ', isteps_part, ' (inner time steps)'
#endif

    ! Periodic initialization
#ifdef BACK
    tpart = ceiling( real(kd_start - kd_end+1)/part_freq)
#else
    tpart = ceiling( real(kd_end-kd_start+1)/part_freq )
#endif    


#ifdef DEBUG                 
    write(*,'(a,i10)') ' part seeding events #   : ', tpart
#endif

    ! NOW we initialize our krill particles

    ! FIRST open sequence of repeated initialization...
    do t = 1, tpart

       i = 0

       ! ... then get random positions on the grid...
       do k = 1, npart0         ! Cell position
          do l = 1, part_fac ! Individual particle within cell

             i = i + 1

             if (t==1) then
                !initialize position
7               call rng_uni(r1)
                call rng_uni(r2)

                xpo(i) = real(ipart(k)) - 1. + r1
                ypo(i) = real(jpart(k)) - 1. + r2

                xgrid(i) = max(1, min(m, ceiling(xpo(i))))
                ygrid(i) = max(1, min(n, ceiling(ypo(i))))

                do while (nlayer(xgrid(i), ygrid(i)) == 0)
                   goto 7 ! check that the particle is in a wet cell !
                end do

                zz(i)  = kriPart%zday_fun(salt2(xgrid(i),ygrid(i)))
                zzz(i) = kriPart%znight_fun()
                zpo(i) = vert_mig(0., zz(i), zzz(i), mid_depth(nlayer(xgrid(i), ygrid(i))), &
                                  theta, sigma, dvm)

                zloop: do j = 1, min( 20, nlayer(xgrid(i),ygrid(i)) )
                   zgrid(i) = j
                   if (zpo(i)<=mid_depth(j)) exit zloop
                enddo zloop

#ifdef BACK
                ! Days from the run start
                kd = kd_ini + istart*dtnc/86400 - (kd_start+(t-1)*part_freq)

                ! Initial particle timestep (kstep); use dt in sec
                istart_part = int( (kd*86400 - part_start(3)*3600 - part_start(2)*60 - part_start(1))/dt )
#else
                ! Days from the run start
                kd = kd_start - kd_ini - istart*dtnc/86400

                ! Initial particle timestep (kstep); use dt in sec
                istart_part = int( (kd*86400 + part_start(3)*3600 + part_start(2)*60 + part_start(1))/dt )
#endif
                tp(i)   = temp(xgrid(i),ygrid(i),zgrid(i))

                if (clone .eqv. .true.) then 
                   mass0(i) = mean_mass
                else 
                   call rng_uni(r1)
                   mass0(i) = min_mass + r1 * (max_mass - min_mass)
                end if
                
                kriPart%part_spec = part_spec

                kriPart%istart    = istart_part
                kriPart%isteps    = isteps_part

                kriPart%dvm       = dvm
                kriPart%light_day = light_day
                kriPart%zday      = zz(i)
                kriPart%znight    = zzz(i)

                kriPart%cell_out  = cell_out

                kriPart%xpo       = xpo(i)
                kriPart%ypo       = ypo(i)
                kriPart%zpo       = zpo(i)

                kriPart%xk1       = xk1
                kriPart%xk2       = xk2
                kriPart%xk3       = xk3
                kriPart%xk4       = xk4

                kriPart%yk1       = yk1
                kriPart%yk2       = yk2
                kriPart%yk3       = yk3
                kriPart%yk4       = yk4

                kriPart%zk1       = zk1
                kriPart%zk2       = zk2
                kriPart%zk3       = zk3
                kriPart%zk4       = zk4

                ! KRILL attributes

                kriPart%tp        = tp(i)

                kriPart%sex       = sex
                kriPart%mass      = mass0(i) * 1e3
                kriPart%dev_frac  = 0.0

                kriPart%a_molt    = a_molt
                kriPart%b_molt    = b_molt

                kriPart%er        = er
                kriPart%t_lim     = t_lim

                kriPart%r0        = r0
                kriPart%rb        = rb

                ! ... and finally, duplicate particles for each sequence of initialization
             else
                
                kriPart%part_spec = part_spec

                kriPart%istart    = istart_part
                kriPart%isteps    = isteps_part

                kriPart%dvm       = dvm
                kriPart%light_day = light_day
                kriPart%zday      = zz(i)
                kriPart%znight    = zzz(i)

                kriPart%cell_out  = cell_out

                kriPart%xpo       = xpo(i)
                kriPart%ypo       = ypo(i)
                kriPart%zpo       = zpo(i)

                kriPart%xk1       = xk1
                kriPart%xk2       = xk2
                kriPart%xk3       = xk3
                kriPart%xk4       = xk4

                kriPart%yk1       = yk1
                kriPart%yk2       = yk2
                kriPart%yk3       = yk3
                kriPart%yk4       = yk4

                kriPart%zk1       = zk1
                kriPart%zk2       = zk2
                kriPart%zk3       = zk3
                kriPart%zk4       = zk4

                ! KRILL attributes

                kriPart%tp        = tp(i)

                kriPart%sex       = sex
                kriPart%mass      = mass0(i) * 1e3
                kriPart%dev_frac  = 0.0

                kriPart%a_molt    = a_molt
                kriPart%b_molt    = b_molt

                kriPart%er        = er
                kriPart%t_lim     = t_lim

                kriPart%r0        = r0
                kriPart%rb        = rb

             endif

             call kriList%addPart(kriPart)

          end do
       end do
    end do

    npart = npart0 * part_fac

    write(*,'(/a,i6,a)') ' Generated ', npart, ' particles'

  end subroutine generateKrill

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!72

  subroutine generateCalanus(calList, runList)

    type(Calanus)                   :: calPart
  
    class(partList)                 :: calList


    character(*), intent(in)        :: runList

    integer                         :: i, j, ii, k, l, ioerr, t, IOstatus ! counters or iterators

    integer                         :: kd_start, kd_end, kd

    integer, dimension(6)           :: part_start, part_end

    integer, dimension(max_part)    :: xgrid, ygrid, zgrid

    real                            :: r1, r2, re1, re2

    real,    dimension(max_part)    :: xpo, ypo, zpo, zz, zzz, tp, fd, ph

    real,    dimension(max_part)    :: mass0, stage0, abund0, age0

    integer, dimension(max_part)    :: part_serial, part_mater, part_birthdate
    
    integer, dimension(max_part)    :: activity0

    real,    dimension(14, 6576)    :: pars ! very specific...

    real ,dimension(1:9)            :: Ts, Fs
    
    ! hard coding epsilon range for diapause, since it need to be round to second decimal ! prev dimension(1:5) 
    real, dimension(1:8),parameter  :: ep_range1=[0.33,0.34,0.35,0.36,0.37,0.38,0.39,0.40] !
    
    real, dimension(1:21),parameter :: ep_range2=[0.4 ,0.41,0.42,0.43,0.44,0.45,0.46,0.47,0.48,0.49,0.5,0.51,0.52, &
                                                      0.53,0.54,0.55,0.56,0.57,0.58,0.59,0.6]
    
    integer                         :: re1_index, re2_index
    !------------------------------------------------------------------72
 
    open(10,file = runList, status = 'old', action = 'read', err = 8, iostat = ioerr)
8   if (ioerr /= 0) STOP 'problem opening namelist file'
    read(10, nml = part_nml)
    close(10)

    !------------------------------------------------------------------72
    ! Initial particles distribution (spatial)
    !
    ! Works with input_particles (below) to provide the # of particles
    ! added at specific time steps in case of periodic seeding of particles

    open(10,file = partfile, status = 'old', action = 'read', err = 9, iostat = ioerr)
9   if (ioerr /= 0) STOP 'problem opening file, particle input'

    do i = 1, m
       read(10, *) (part_grid_init(i, j), j = 1, n)
    end do

    npart0  = sum(part_grid_init) !npart0 is the number of cells with particles

    ii = 0

    do j = 1, n    !n is a horizontal dimension of the grid
       do i = 1, m !m is the other horizontal dimension of the grid
          if(part_grid_init(i, j) == 1) then
             ii = min(ii + 1, npart0)
             ipart(ii) = i
             jpart(ii) = j
          endif
       end do
    end do

    close(10)

    ! Number of individuals generated
    write(*,'(/a,i9,a,/)') ' Cells with particles = ', npart0, '; now generating all the particles variables'

    !------------------------------------------------------------------72
    ! Initialization of individual particles

    ! Open namelist that contains physiological parameters  
    open(10, file = runList, status = 'old', action = 'read', err = 8, iostat = ioerr)
    read(10, nml = calanus_nml)

    close(10)


    ! Open matrix of individual parameters
    if (INTRA_SPE .eqv. .true.) then

       open(10, file=parsfile, status = 'old', action = 'read', err = 10, iostat = ioerr)
10     if (ioerr/=0) stop '!!! Problem opening parameters matrix !!!'

       do i = 1, nb_pars 
          read(10,*) (pars(i,j), j = 1, npart0)
       end do

       close(10)

    endif

    ! Initializing time (time would be tzero + one dt)
    call update2(initial)

    ! Initialize time loop for periodic release of particles:
    ! start of particles release -> end of release (by frequency of release)

    part_start = (/ 0, 0, 0, d_Astart, m_Astart, y_Astart /) ! Start of particles release
    call cday2(part_start(4),part_start(5),part_start(6),kd_start)

    part_end   = (/ 0, 0, 0, d_Aend, m_Aend, y_Aend /)       ! End of particles release
    call cday2(part_end(4),part_end(5),part_end(6),kd_end)

#ifdef DEBUG
    write(*,'(/a,i2,1x,i2,1x,i4,a,i6)') ' part start              : ', part_start(4), part_start(5), part_start(6), &
         '; start timestep : ', kd_start
    write(*,'(a,i2,1x,i2,1x,i4,a,i6)') ' part end                : ', part_end(4),   part_end(5),   part_end(6),   &
         '; end timestep   : ', kd_end
#endif

    ! Duration of particles advection; same for all particles
    isteps_part = part_duration * 86400/dt

#ifdef DEBUG                 
    write(*,'(a,i10,a)') ' part advection duration : ', isteps_part, ' (inner time steps)'
#endif

    ! Periodic initialization
#ifdef BACK
    tpart = ceiling( real(kd_start - kd_end+1)/part_freq)
#else
    tpart = ceiling( real(kd_end-kd_start+1)/part_freq )
#endif    


#ifdef DEBUG                 
    write(*,'(a,i10)') ' part seeding events #   : ', tpart
#endif

    ! NOW we initialize our Calanus copepods

    ! FIRST open sequence of repeated initialization...
    do t = 1, tpart

       i = 0
     
#ifdef ADV
       npart0=npart0
#else
       npart0=1
#endif
 
       ! ... then get random positions on the grid...
       do k = 1, npart0         ! Cell position
          do l = 1, part_fac ! Individual particle within cell

             i = i + 1

             if (t==1) then
                !initialize position
#ifdef FTLE
! in case of ftle, there is no randomness added to the particles positions :

                xpo(i) = real(ipart(k))
                ypo(i) = real(jpart(k))

                xgrid(i) = max(1, min(m, ceiling(xpo(i))))
                ygrid(i) = max(1, min(n, ceiling(ypo(i))))

                if (nlayer(xgrid(i), ygrid(i)) == 0) stop '!!! STOP: one partile in the initial position is in a dry cell !!!'

#else
11              call rng_uni(r1)
                call rng_uni(r2)
               
#ifdef ADV
                xpo(i) = real(ipart(k)) - 1. + r1
                ypo(i) = real(jpart(k)) - 1. + r2
#else
                ! when not define ADV initilize x,y position from user defined init_xpo and init_ypo from run.list 
                xpo(i)=init_xpo
                ypo(i)=init_ypo
#endif

                xgrid(i) = max(1, min(m, ceiling(xpo(i))))
                ygrid(i) = max(1, min(n, ceiling(ypo(i))))

                do while (nlayer(xgrid(i), ygrid(i)) == 0)
                   goto 11 ! check that the particle is in a wet cell !
                end do
#endif

#ifdef YUMMY
                Fs  = food(xgrid(i),ygrid(i),1:9)
#else
                Fs  = [1, 1, 1, 1, 1, 1, 1, 1, 1]
#endif
                Ts  = temp(xgrid(i),ygrid(i),1:9)

                zz(i)  = calPart%zday_fun(Fs,Ts,mid_depth(nlayer(xgrid(i),ygrid(i))),0.0)
                zzz(i) = calPart%znight_fun(Fs,Ts,mid_depth(nlayer(xgrid(i),ygrid(i))))
                zpo(i) = vert_mig(0., zz(i), zzz(i), mid_depth(nlayer(xgrid(i), ygrid(i))), &
                                  theta, sigma, dvm)

                zloop: do j = 1, min( 20, nlayer(xgrid(i),ygrid(i)) )
                   zgrid(i) = j
                   if (zpo(i)<=mid_depth(j)) exit zloop
                enddo zloop

#ifdef BACK
                ! Days from the run start
                kd = kd_ini + istart*dtnc/86400 - (kd_start+(t-1)*part_freq)

                ! Initial particle timestep (kstep); use dt in sec
                istart_part = int( (kd*86400 - part_start(3)*3600 - part_start(2)*60 - part_start(1))/dt )
#else
                ! Days from the run start
                kd = kd_start - kd_ini - istart*dtnc/86400

                ! Initial particle timestep (kstep); use dt in sec
                istart_part = int( (kd*86400 + part_start(3)*3600 + part_start(2)*60 + part_start(1))/dt )
#endif
                tp(i)   = temp(xgrid(i),ygrid(i),zgrid(i))

#ifdef YUMMY
                fd(i)   = food(xgrid(i),ygrid(i),zgrid(i))

                ph(i)   = phyto3d(xgrid(i),ygrid(i),zgrid(i))
#endif

                stage0(i) = init_stage
                age0(i) = init_age
                abund0(i) = init_abund
                activity0(i) = init_activity
                
                part_birthdate(i) = d_Astart + m_Astart*100 + y_Astart*10000
                part_serial(i) = i
                part_mater(i) = 0
                

#ifdef DEBUG
                if (i .eq. 1) then
                    write(*,'(/a,i8)') 'part_birthdate(i) : ', part_birthdate(i)
                    write(*,'(/a,i6)') 'part_serial(i) : ', part_serial(i)
                endif
#endif
                if (clone .eqv. .true.) then 
                   mass0(i) = mean_mass
                else 
                   call rng_uni(r1)
                   mass0(i) = min_mass + r1 * (max_mass - min_mass)
                end if
                
                calPart%part_spec = part_spec

                calPart%istart    = istart_part
                calPart%isteps    = isteps_part

                calPart%dvm       = dvm
                calPart%light_day = light_day
                calPart%zday      = zz(i)
                calPart%znight    = zzz(i)

                calPart%cell_out  = cell_out

                calPart%xpo       = xpo(i)
                calPart%ypo       = ypo(i)
                calPart%zpo       = zpo(i)

                calPart%xk1       = xk1
                calPart%xk2       = xk2
                calPart%xk3       = xk3
                calPart%xk4       = xk4

                calPart%yk1       = yk1
                calPart%yk2       = yk2
                calPart%yk3       = yk3
                calPart%yk4       = yk4

                calPart%zk1       = zk1
                calPart%zk2       = zk2
                calPart%zk3       = zk3
                calPart%zk4       = zk4

                ! CALANUS attributes

                calPart%tp        = tp(i)
#ifdef YUMMY
                calPart%fd        = fd(i)
                calPart%ph        = ph(i)
#endif
                calPart%sex       = sex
                calPart%mass      = mass0(i)
                calPart%massp     = 0.
                calPart%lipid     = 0.
                calPart%lipid_floor=0.
                calPart%counter_eggs=0.
                calPart%cum_eggs  =0
                calPart%start_ep  = 0.
                calPart%reproduce_info  = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
                calPart%duration_ep=0.
                calPart%ep_daily  = 0.
                calPart%EP        = 0
                calPart%egg_clutch= 0
                calPart%EP_CLUTCH = 0
                calPart%max_diap_lipid=0.
                calPart%lipid_exit_diap=0.
                calPart%mass_limit= 0.
                calPart%stage     = stage0(i)
                calPart%age       = age0(i)
                calPart%abund     = abund0(i)
                calPart%activity  = activity0(i)
                calPart%part_status  = 1
                calPart%a_molt    = a_devlp(int(stage0(i)))
                calPart%b_molt    = t_devlp

                calPart%er        = er
                calPart%t_lim     = t_lim

                calPart%r0        = r0
                calPart%rb        = rb
                !initialize epsilon for diapause entry and exit
                re1=rand()
                re2=rand()
                re1_index         = int(re1*5)+1
                re2_index         = int(re2*21)+1
                calPart%ep1       = ep_range1(re1_index)
                calPart%ep2       = ep_range2(re2_index)
                calPart%counter_c6_active  = 0.0
                calPart%flag_enter_diap    = .false.
                calPart%pre_activity       = 0
                calPart%serial    = part_serial(i)
                calPart%birthdate = part_birthdate(i)
                calPart%mater     = 0
                ! ... and finally, duplicate particles for each sequence of initialization
             else
                             
                calPart%part_spec = part_spec

                calPart%istart    = istart_part
                calPart%isteps    = isteps_part

                calPart%dvm       = dvm
                calPart%light_day = light_day
                calPart%zday      = zz(i)
                calPart%znight    = zzz(i)

                calPart%cell_out  = cell_out

                calPart%xpo       = xpo(i)
                calPart%ypo       = ypo(i)
                calPart%zpo       = zpo(i)

                calPart%xk1       = xk1
                calPart%xk2       = xk2
                calPart%xk3       = xk3
                calPart%xk4       = xk4

                calPart%yk1       = yk1
                calPart%yk2       = yk2
                calPart%yk3       = yk3
                calPart%yk4       = yk4

                calPart%zk1       = zk1
                calPart%zk2       = zk2
                calPart%zk3       = zk3
                calPart%zk4       = zk4
                calPart%part_status  = 1
                ! CALANUS attributes

                calPart%tp        = tp(i)
#ifdef YUMMY
                calPart%fd        = fd(i)
                calPart%ph        = ph(i)
#endif
                calPart%sex       = sex
                calPart%mass      = mass0(i)
                calPart%massp     = 0.
                calPart%lipid     = 0.
                calPart%lipid_floor=0.
                calPart%counter_eggs=0.
                calPart%cum_eggs  = 0
                calPart%start_ep  = 0.
                calPart%duration_ep=0.
                calPart%ep_daily  = 0.
                calPart%egg_clutch= 0
                calPart%reproduce_info  = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
                calPart%EP        = 0
                calPart%EP_CLUTCH = 0
                calPart%max_diap_lipid=0.
                calPart%lipid_exit_diap=0.
                calPart%mass_limit= 0.
                calPart%stage     = stage0(i)
                calPart%age       = age0(i)
                calPart%abund     = abund0(i)
                calPart%activity  = activity0(i)
                calPart%a_molt    = a_devlp(int(stage0(i)))
                calPart%b_molt    = t_devlp

                calPart%er        = er
                calPart%t_lim     = t_lim

                calPart%r0        = r0
                calPart%rb        = rb
                re1=rand()
                re2=rand()
                re1_index         = int(re1*5)+1
                re2_index         = int(re2*21)+1
                calPart%ep1       = ep_range1(re1_index)
                calPart%ep2       = ep_range2(re2_index)
                calPart%counter_c6_active  = 0.0
                calPart%flag_enter_diap    = .false.
                calPart%pre_activity       = 0
                calPart%serial    = part_serial(i)
                calPart%birthdate = part_birthdate(i)
                calPart%mater     = 0

             endif

             call calList%addPart(calPart)

          end do
       end do
   
           ! generate additional status=0 particles for preparation of adding particles for egg production
       do k=1,init_part*part_fac*npart0
           calPart%part_spec = part_spec

           calPart%istart    = istart_part
           calPart%isteps    = isteps_part

           calPart%dvm       = dvm
           calPart%light_day = 0
           calPart%zday      = 0
           calPart%znight    = 0

           calPart%cell_out  = .false.

           calPart%xpo       = 0
           calPart%ypo       = 0
           calPart%zpo       = 0

           calPart%xk1       = 0
           calPart%xk2       = 0
           calPart%xk3       = 0
           calPart%xk4       = 0

           calPart%yk1       = 0
           calPart%yk2       = 0
           calPart%yk3       = 0
           calPart%yk4       = 0

           calPart%zk1       = 0
           calPart%zk2       = 0
           calPart%zk3       = 0
           calPart%zk4       = 0
           ! give these particles 0 status
           calPart%part_status  = 0
           ! CALANUS attributes

           calPart%tp        = 0
#ifdef YUMMY
           calPart%fd        = 0
           calPart%ph        = 0
#endif
           calPart%sex       = 0
           calPart%mass      = 0
           calPart%massp     = 0.
           calPart%lipid     = 0
           calPart%stage     = 0
           calPart%age       = 0
           calPart%abund     = 0
           calPart%activity  = 0
           calPart%a_molt    = 0
           calPart%b_molt    = 0
           calPart%lipid_floor=0.
           calPart%counter_eggs=0.
           calPart%cum_eggs  =0
           calPart%start_ep  = 0.
           calPart%reproduce_info  = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
           calPart%duration_ep=0.
           calPart%ep_daily  = 0.
           calPart%EP        = 0
           calPart%egg_clutch= 0
           calPart%EP_CLUTCH = 0
           calPart%max_diap_lipid=0.
           calPart%lipid_exit_diap=0.
           calPart%mass_limit= 0.
           calPart%er        = 0
           calPart%t_lim     = 0

           calPart%r0        = 0
           calPart%rb        = 0
   
           re1=rand()
           re2=rand()

           re1_index         = int(re1*5)+1
           re2_index         = int(re2*21)+1
           calPart%ep1       = ep_range1(re1_index)
           calPart%ep2       = ep_range2(re2_index)

           calPart%counter_c6_active  = 0.0
           calPart%flag_enter_diap    = .false.
           calPart%pre_activity       = 0
           calPart%serial    = 0
           calPart%birthdate = 0
           calPart%mater     = 0
   
           call calList%addPart(calPart)
       end do   
    end do
    
    npart = npart0 * part_fac*(init_part+1)
    next_available_part=npart0*part_fac+1

    write(*,'(/a,i6,a)') ' Generated ', npart, ' particles'

  end subroutine generateCalanus


  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!72




subroutine generateLate_hyperboreus(LHList, runList, in_diapause)
  type(Late_hyperboreus)                    :: LHPart
  class(partList)              :: LHList


  character(*), intent(in)     :: runList

  integer, intent(in)          :: in_diapause

  integer                      :: i, j, ii, k, l, ioerr, t, IOstatus ! counters or iterators

  integer                      :: kd_start, kd_end, kd

  integer, dimension(6)        :: part_start, part_end

  integer, dimension(max_part) :: xgrid, ygrid, zgrid

  integer                      :: buffer_lim

  real                         :: r1, r2, r3, r4, r5

  real,    dimension(max_part) :: xpo, ypo, zpo, zz, zzz, tp, sigma_local
  ! sigma local is the sigma variable that can change depending on the local environmental conditions

  integer                      :: ipo, jpo

  integer, dimension(max_part) :: mode

  real,    dimension(max_part) :: mode_prop


  !------------------------------------------------------------------72

  open(10,file = runList, status = 'old', action = 'read', err = 1, iostat = ioerr)
1   if (ioerr /= 0) STOP 'problem opening file'
  read(10, nml = part_nml)
  close(10)

!------------------------------------------------------------------72
! Initial particles distribution (spatial)
!
! Works with input_particles (below) to provide the # of particles
! added at specific time steps in case of periodic seeding of particles

  open(10,file = partfile, status = 'old', action = 'read', err = 5, iostat = ioerr)
  5   if (ioerr /= 0) STOP 'problem opening file, particle input'

  do i = 1, m
     read(10, *) (part_grid_init(i, j), j = 1, n)
  end do

  npart0  = sum(part_grid_init) !npart0 is the number of cells with particles

  ! allocate the vectors used for the cells with particles
  allocate(ipart(npart0))
  allocate(jpart(npart0))

  ! record the coordinates of the cells with particles
  ii = 0
  outer : do j = 1, n    !n is a horizontal dimension of the grid
     inner : do i = 1, m !m is the other horizontal dimension of the grid
        if(part_grid_init(i, j) >= 1) then
           ii = min(ii + 1, npart0)
           ipart(ii) = i
           jpart(ii) = j
        endif
     end do inner
  end do outer

  close(10)

  ! Number of individuals generated
  write(*,'(/a,i9,a,/)') ' Cells with particles = ', npart0, '; now generating all the particles variables'

  !------------------------------------------------------------------72
  ! Initialization of individual particles

  ! Initializing time (time would be tzero + one dt)
  call update2(initial)

  ! Initialize time loop for periodic release of particles:
  ! start of particles release -> end of release (by frequency of release)

  part_start = (/ 0, 0, 0, d_Astart, m_Astart, y_Astart /) ! Start of particles release
  call cday2(part_start(4),part_start(5),part_start(6),kd_start)

  part_end   = (/ 0, 0, 0, d_Aend, m_Aend, y_Aend /)       ! End of particles release
  call cday2(part_end(4),part_end(5),part_end(6),kd_end)

#ifdef DEBUG
  write(*,'(/a,i2,1x,i2,1x,i4,a,i6)') ' part start              : ', part_start(4), part_start(5), part_start(6), &
                                      ';   start timestep : ', kd_start
  write(*,'(a,i2,1x,i2,1x,i4,a,i6)')  ' part end                : ', part_end(4),   part_end(5),   part_end(6),   &
                                      ';   end timestep   : ', kd_end
#endif

  ! Duration of particles advection; same for all particles
  isteps_part = part_duration * 86400/dt

#ifdef iEBUG
  write(*,'(a,i10,a)') ' part advection duration : ', isteps_part, ' (inner time steps)'
#endif

  ! Periodic initialization
#ifdef BACK
  tpart = ceiling( real(kd_start - kd_end+1)/part_freq)
#else
  tpart = ceiling( real(kd_end - kd_start+1)/part_freq )
#endif

#ifdef DEBUG
  write(*,'(a,i10)') ' part seeding events #   : ', tpart
#endif

  ! NOW we initialize each particle's values

  ! FIRST open sequence of repeated initialization...
  do t = 1, tpart

     i = 0

     ! ... then get random positions on the grid...
     do k = 1, npart0      ! Cell position
        ! write(*,'(a,i10)') ' k  : ', k , ' / ', npart0 !just to debug
        do l = 1, part_fac ! Individual particle within cell

           ! ... and duplicate particles for each sequence of initialization
           i = i + 1

           if (t==1) then

              !initialize position
#ifdef FTLE
! in case of ftle, there is no randomness added to the particles positions :

              xpo(i) = real(ipart(k))
              ypo(i) = real(jpart(k))

              xgrid(i) = max(1, min(m, ceiling(xpo(i))))
              ygrid(i) = max(1, min(n, ceiling(ypo(i))))

              if (nlayer(xgrid(i), ygrid(i)) == 0) stop '!!! STOP: one partile in the initial position is in a dry cell !!!'

#else
3             call rng_uni(r1)
              call rng_uni(r2)

              xpo(i) = real(ipart(k)) - 1. + r1
              ypo(i) = real(jpart(k)) - 1. + r2

              xgrid(i) = max(1, min(m, ceiling(xpo(i))))
              ygrid(i) = max(1, min(n, ceiling(ypo(i))))

              !write(*,'(a,i10)') ' wet ? '  !just to debug
              
              do while (nlayer(xgrid(i), ygrid(i)) == 0)
                 goto 3 ! check that the particle is in a wet cell !
              end do
#endif
              
              !write(*,'(a,i10)') ' wet ! '  !just to debug

              ! Find the mode proportion of each particle place.  The mode proportion btw 0 and 1 is the proportion of particles
              ! that should be associated to the surface mode. This mode proportion can depend or not on
              ! The the local bottom depth. If not, this proportion just depends on the type of particle of the model,
              ! and will therefore be the same for each particle.
              mode_prop(i) = LHPart%mode_prop_fun( LHPart%mode_prop, mid_depth(nlayer(xgrid(i),ygrid(i))), in_diapause )

              call rng_uni(r3) !random number between 0 and 1

              ! if random number is < mode proportion, -> the particle is in the the upper mode (1).
              ! Else, is is in the lower mode (2)
              if (r3 < mode_prop(i)) then
                 mode(i) = 1
              else
                 mode(i) = 2
              endif


              zz(i)  = LHPart%zday_fun  (mid_depth(nlayer(xgrid(i), ygrid(i))), mode(i), in_diapause)
              zzz(i) = LHPart%znight_fun(mid_depth(nlayer(xgrid(i), ygrid(i))), mode(i), in_diapause)
              sigma_local(i) = LHPart%sigma_fun(mid_depth(nlayer(xgrid(i), ygrid(i))), mode(i), in_diapause)
              zpo(i) = vert_mig(0., zz(i), zzz(i), mid_depth(nlayer(xgrid(i), ygrid(i))), &
                                theta, sigma_local(i), dvm)

              ! If the particle's depth is under bottom depth, put it randomly in a buffer 50m under bottom depth :
              ! Or in a smaller buffer is the bdeph is shallower than 50m
              if ( zpo(i) > mid_depth(nlayer(xgrid(i), ygrid(i))) ) then
                 call rng_uni(r4) !random number between 0 and 1
                 buffer_lim = min( mid_depth(nlayer(xgrid(i), ygrid(i))), 10.0 )
                 zpo(i) = (r4*buffer_lim) + ( mid_depth(nlayer(xgrid(i), ygrid(i)))-buffer_lim )
              endif

              ! If the particle's depth is over the surface, put it randomly in a buffer 50m under surface :
              ! Or in a smaller buffer is the bdeph is shallower than 50m
              if ( zpo(i) < 0 ) then
                 call rng_uni(r5) !random number between 0 and 1
                 buffer_lim = min( mid_depth(nlayer(xgrid(i), ygrid(i))), 10.0 )
                 zpo(i) = r5*buffer_lim
              endif

              ! Loop to find zgrid(i), the depth layer corresponding to the depth zpo(i) :
              zloop: do j = 1, min( ilo, nlayer(xgrid(i),ygrid(i)) )
                 zgrid(i) = j
                 if (zpo(i)<=mid_depth(j)) exit zloop
              enddo zloop

#ifdef BACK
              ! Days from the run start
              kd = kd_ini + istart*dtnc/86400 - (kd_start+(t-1)*part_freq)

              ! Initial particle timestep (kstep); use dt in sec
              istart_part = int( (kd*86400 - part_start(3)*3600 - part_start(2)*60 - part_start(1))/dt )
#else
              ! Days from the run start
              kd = kd_start - kd_ini - istart*dtnc/86400

              ! Initial particle timestep (kstep); use dt in sec
              istart_part = int( (kd*86400 + part_start(3)*3600 + part_start(2)*60 + part_start(1))/dt )
#endif
              tp(i) = temp(xgrid(i),ygrid(i),zgrid(i))
              
           endif

           LHPart%part_spec = part_spec

           LHPart%istart    = istart_part
           LHPart%isteps    = isteps_part

           LHPart%dvm       = dvm
           LHPart%light_day = light_day
           LHPart%zday      = zz(i)
           LHPart%znight    = zzz(i)

           LHPart%cell_out  = cell_out

           LHPart%xpo       = xpo(i)
           LHPart%ypo       = ypo(i)
           LHPart%zpo       = zpo(i)
           LHPart%sigma_local = 3 !sigma_local(i)
           LHPart%bdep      = mid_depth(nlayer(xgrid(i), ygrid(i)))

#ifdef FTLE
           LHPart%ipo       = ipart(i)
           LHPart%jpo       = jpart(i)
#endif
           LHPart%mode      = mode(i)
           LHPart%mode_prop = mode_prop(i)

           LHPart%xk1       = xk1
           LHPart%xk2       = xk2
           LHPart%xk3       = xk3
           LHPart%xk4       = xk4

           LHPart%yk1       = yk1
           LHPart%yk2       = yk2
           LHPart%yk3       = yk3
           LHPart%yk4       = yk4

           LHPart%zk1       = zk1
           LHPart%zk2       = zk2
           LHPart%zk3       = zk3
           LHPart%zk4       = zk4

           call LHList%addPart(LHPart)

        end do
     end do
  end do

  npart = npart0 * part_fac

  deallocate(ipart,jpart)

  write(*,'(/a,i6,a)') ' Generated ', npart, ' particles'

end subroutine generateLate_hyperboreus


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!72




subroutine generateYoung_hyperboreus(YHList, runList)
  type(Young_hyperboreus)                    :: YHPart
  class(partList)              :: YHList


  character(*), intent(in)     :: runList

  integer                      :: i, j, ii, k, l, ioerr, t, IOstatus ! counters or iterators

  integer                      :: kd_start, kd_end, kd

  integer, dimension(6)        :: part_start, part_end

  integer, dimension(max_part) :: xgrid, ygrid, zgrid

  integer                     :: buffer_lim

  real                         :: r1, r2, r3, r4, r5

  real,      dimension(max_part) :: xpo, ypo, zpo, zz, zzz, tp, sigma_local
  integer,    dimension(max_part) :: mode
  real,       dimension(max_part) :: mode_prop


  !------------------------------------------------------------------72

  open(10,file = runList, status = 'old', action = 'read', err = 1, iostat = ioerr)
1   if (ioerr /= 0) STOP 'problem opening file'
  read(10, nml = part_nml)
  close(10)

!------------------------------------------------------------------72
! Initial particles distribution (spatial)
!
! Works with input_particles (below) to provide the # of particles
! added at specific time steps in case of periodic seeding of particles

  open(10,file = partfile, status = 'old', action = 'read', err = 5, iostat = ioerr)
  5   if (ioerr /= 0) STOP 'problem opening file, particle input'

  do i = 1, m
     read(10, *) (part_grid_init(i, j), j = 1, n)
  end do

  npart0  = sum(part_grid_init) !npart0 is the number of cells with particles

  ! allocate the vectors used for the cells with particles
  allocate(ipart(npart0))
  allocate(jpart(npart0))

  ! record the coordinates of the cells with particles
  ii = 0
  outer : do j = 1, n    !n is a horizontal dimension of the grid
     inner : do i = 1, m !m is the other horizontal dimension of the grid
        if(part_grid_init(i, j) == 1) then
           ii = min(ii + 1, npart0)
           ipart(ii) = i
           jpart(ii) = j
        endif
     end do inner
  end do outer

  close(10)

  ! Number of individuals generated
  write(*,'(/a,i9,a,/)') ' Cells with particles = ', npart0, '; now generating all the particles variables'

  !------------------------------------------------------------------72
  ! Initialization of individual particles

  ! Initializing time (time would be tzero + one dt)
  call update2(initial)

  ! Initialize time loop for periodic release of particles:
  ! start of particles release -> end of release (by frequency of release)

  part_start = (/ 0, 0, 0, d_Astart, m_Astart, y_Astart /) ! Start of particles release
  call cday2(part_start(4),part_start(5),part_start(6),kd_start)

  part_end   = (/ 0, 0, 0, d_Aend, m_Aend, y_Aend /)       ! End of particles release
  call cday2(part_end(4),part_end(5),part_end(6),kd_end)

#ifdef DEBUG
  write(*,'(/a,i2,1x,i2,1x,i4,a,i6)') ' part start              : ', part_start(4), part_start(5), part_start(6), &
                                      ';   start timestep : ', kd_start
  write(*,'(a,i2,1x,i2,1x,i4,a,i6)')  ' part end                : ', part_end(4),   part_end(5),   part_end(6),   &
                                      ';   end timestep   : ', kd_end
#endif

  ! Duration of particles advection; same for all particles
  isteps_part = part_duration * 86400/dt

#ifdef DEBUG
  write(*,'(a,i10,a)') ' part advection duration : ', isteps_part, ' (inner time steps)'
#endif

  ! Periodic initialization
#ifdef BACK
  tpart = ceiling( real(kd_start - kd_end+1)/part_freq)
#else
  tpart = ceiling( real(kd_end - kd_start+1)/part_freq )
#endif

#ifdef DEBUG
  write(*,'(a,i10)') ' part seeding events #   : ', tpart
#endif

  ! NOW we initialize each particle's values

  ! FIRST open sequence of repeated initialization...
  do t = 1, tpart

     i = 0

     ! ... then get random positions on the grid...
     do k = 1, npart0      ! Cell position
        !write(*,'(a,i10)') ' k  : ', k !just to debug
        do l = 1, part_fac ! Individual particle within cell

           ! ... and duplicate particles for each sequence of initialization
           i = i + 1

           if (t==1) then

              !initialize position
#ifdef FTLE
! in case of ftle, there is no randomness added to the particles positions :

              xpo(i) = real(ipart(k))
              ypo(i) = real(jpart(k))

              xgrid(i) = max(1, min(m, ceiling(xpo(i))))
              ygrid(i) = max(1, min(n, ceiling(ypo(i))))

              if (nlayer(xgrid(i), ygrid(i)) == 0) stop '!!! STOP: one partile in the initial position is in a dry cell !!!'

#else
3             call rng_uni(r1)
              call rng_uni(r2)

              xpo(i) = real(ipart(k)) - 1. + r1
              ypo(i) = real(jpart(k)) - 1. + r2

              xgrid(i) = max(1, min(m, ceiling(xpo(i))))
              ygrid(i) = max(1, min(n, ceiling(ypo(i))))

              do while (nlayer(xgrid(i), ygrid(i)) == 0)
                 goto 3 ! check that the particle is in a wet cell !
              end do
#endif
              !write(*,'(a,i10)') ' wet '  !just to debug

              ! Find the mode proportion of each particle place.  The mode proportion btw 0 and 1 is the proportion of particles
              ! that should be associated to the surface mode. This mode proportion can depend or not on
              ! The the local bottom depth. If not, this proportion just depends on the type of particle of the model,
              ! and will therefore be the same for each particle.
              mode_prop(i) = YHPart%mode_prop_fun( YHPart%mode_prop, mid_depth(nlayer(xgrid(i),ygrid(i))), in_diapause )

              call rng_uni(r3) !random number between 0 and 1

              ! if random number is < mode proportion, -> the particle is in the the upper mode (1).
              ! Else, is is in the lower mode (2)
              if (r3 < mode_prop(i)) then
                 mode(i) = 1
              else
                 mode(i) = 2
              endif


              zz(i)  = YHPart%zday_fun  (mid_depth(nlayer(xgrid(i), ygrid(i))), mode(i))
              zzz(i) = YHPart%znight_fun(mid_depth(nlayer(xgrid(i), ygrid(i))), mode(i))
              sigma_local(i) = YHPart%sigma_fun(mid_depth(nlayer(xgrid(i), ygrid(i))), mode(i))
              zpo(i) = vert_mig(0., zz(i), zzz(i), mid_depth(nlayer(xgrid(i), ygrid(i))), &
                                theta, sigma_local(i), dvm)

              ! If the particle's depth is under bottom depth, put it randomly in a buffer 50m under bottom depth :
              ! Or in a smaller buffer is the bdeph is shallower than 50m
              if ( zpo(i) > mid_depth(nlayer(xgrid(i), ygrid(i))) ) then
                 call rng_uni(r4) !random number between 0 and 1
                 buffer_lim = min( mid_depth(nlayer(xgrid(i), ygrid(i))), 10.0 )
                 zpo(i) = (r4*buffer_lim) + ( mid_depth(nlayer(xgrid(i), ygrid(i)))-buffer_lim )
              endif

              ! If the particle's depth is over the surface, put it randomly in a buffer 50m under surface :
              ! Or in a smaller buffer is the bdeph is shallower than 50m
              if ( zpo(i) < 0 ) then
                 call rng_uni(r5) !random number between 0 and 1
                 buffer_lim = min( mid_depth(nlayer(xgrid(i), ygrid(i))), 10.0 )
                 zpo(i) = r5*buffer_lim
              endif

              ! Loop to find zgrid(i), the depth layer corresponding to the depth zpo(i) :
              zloop: do j = 1, min( ilo, nlayer(xgrid(i),ygrid(i)) )
                 zgrid(i) = j
                 if (zpo(i)<=mid_depth(j)) exit zloop
              enddo zloop

#ifdef BACK
              ! Days from the run start
              kd = kd_ini + istart*dtnc/86400 - (kd_start+(t-1)*part_freq)

              ! Initial particle timestep (kstep); use dt in sec
              istart_part = int( (kd*86400 - part_start(3)*3600 - part_start(2)*60 - part_start(1))/dt )
#else
              ! Days from the run start
              kd = kd_start - kd_ini - istart*dtnc/86400

              ! Initial particle timestep (kstep); use dt in sec
              istart_part = int( (kd*86400 + part_start(3)*3600 + part_start(2)*60 + part_start(1))/dt )
#endif
              tp(i) = temp(xgrid(i),ygrid(i),zgrid(i))
              
           endif

           YHPart%part_spec = part_spec

           YHPart%istart    = istart_part
           YHPart%isteps    = isteps_part

           YHPart%dvm       = dvm
           YHPart%light_day = light_day
           YHPart%zday      = zz(i)
           YHPart%znight    = zzz(i)

           YHPart%cell_out  = cell_out

           YHPart%xpo       = xpo(i)
           YHPart%ypo       = ypo(i)
           YHPart%zpo       = zpo(i)
           YHPart%bdep      = mid_depth(nlayer(xgrid(i), ygrid(i)))

#ifdef FTLE
           YHPart%ipo       = ipart(i)
           YHPart%jpo       = jpart(i)
#endif

           YHPart%mode      = mode(i)
           YHPart%mode_prop = mode_prop(i)

           YHPart%xk1       = xk1
           YHPart%xk2       = xk2
           YHPart%xk3       = xk3
           YHPart%xk4       = xk4

           YHPart%yk1       = yk1
           YHPart%yk2       = yk2
           YHPart%yk3       = yk3
           YHPart%yk4       = yk4

           YHPart%zk1       = zk1
           YHPart%zk2       = zk2
           YHPart%zk3       = zk3
           YHPart%zk4       = zk4

           call YHList%addPart(YHPart)

        end do
     end do
  end do

  npart = npart0 * part_fac

  deallocate(ipart,jpart)

  write(*,'(/a,i6,a)') ' Generated ', npart, ' particles'

end subroutine generateYoung_hyperboreus




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!72




subroutine generateLate_finmarchicus(LFList, runList, in_diapause)
  type(Late_finmarchicus)                    :: LFPart
  class(partList)              :: LFList


  character(*), intent(in)     :: runList

  integer, intent(in)          :: in_diapause

  integer                      :: i, j, ii, k, l, ioerr, t, IOstatus ! counters or iterators

  integer                      :: kd_start, kd_end, kd

  integer, dimension(6)        :: part_start, part_end

  integer, dimension(max_part) :: xgrid, ygrid, zgrid

  integer                     :: buffer_lim

  real                         :: r1, r2, r3, r4, r5

  real,    dimension(max_part) :: xpo, ypo, zpo, zz, zzz, tp, sigma_local
  integer,    dimension(max_part) :: mode
  real,       dimension(max_part) :: mode_prop


  !------------------------------------------------------------------72

  open(10,file = runList, status = 'old', action = 'read', err = 1, iostat = ioerr)
1   if (ioerr /= 0) STOP 'problem opening file'
  read(10, nml = part_nml)
  close(10)

!------------------------------------------------------------------72
! Initial particles distribution (spatial)
!
! Works with input_particles (below) to provide the # of particles
! added at specific time steps in case of periodic seeding of particles

  open(10,file = partfile, status = 'old', action = 'read', err = 5, iostat = ioerr)
  5   if (ioerr /= 0) STOP 'problem opening file, particle input'

  do i = 1, m
     read(10, *) (part_grid_init(i, j), j = 1, n)
  end do

  npart0  = sum(part_grid_init) !npart0 is the number of cells with particles

  ! allocate the vectors used for the cells with particles
  allocate(ipart(npart0))
  allocate(jpart(npart0))

  ! record the coordinates of the cells with particles
  ii = 0
  outer : do j = 1, n    !n is a horizontal dimension of the grid
     inner : do i = 1, m !m is the other horizontal dimension of the grid
        if(part_grid_init(i, j) == 1) then
           ii = min(ii + 1, npart0)
           ipart(ii) = i
           jpart(ii) = j
        endif
     end do inner
  end do outer

  close(10)

  ! Number of individuals generated
  write(*,'(/a,i9,a,/)') ' Cells with particles = ', npart0, '; now generating all the particles variables'

  !------------------------------------------------------------------72
  ! Initialization of individual particles

  ! Initializing time (time would be tzero + one dt)
  call update2(initial)

  ! Initialize time loop for periodic release of particles:
  ! start of particles release -> end of release (by frequency of release)

  part_start = (/ 0, 0, 0, d_Astart, m_Astart, y_Astart /) ! Start of particles release
  call cday2(part_start(4),part_start(5),part_start(6),kd_start)

  part_end   = (/ 0, 0, 0, d_Aend, m_Aend, y_Aend /)       ! End of particles release
  call cday2(part_end(4),part_end(5),part_end(6),kd_end)

#ifdef DEBUG
  write(*,'(/a,i2,1x,i2,1x,i4,a,i6)') ' part start              : ', part_start(4), part_start(5), part_start(6), &
                                      ';   start timestep : ', kd_start
  write(*,'(a,i2,1x,i2,1x,i4,a,i6)')  ' part end                : ', part_end(4),   part_end(5),   part_end(6),   &
                                      ';   end timestep   : ', kd_end
#endif

  ! Duration of particles advection; same for all particles
  isteps_part = part_duration * 86400/dt

#ifdef DEBUG
  write(*,'(a,i10,a)') ' part advection duration : ', isteps_part, ' (inner time steps)'
#endif

  ! Periodic initialization
#ifdef BACK
  tpart = ceiling( real(kd_start - kd_end+1)/part_freq)
#else
  tpart = ceiling( real(kd_end - kd_start+1)/part_freq )
#endif

#ifdef DEBUG
  write(*,'(a,i10)') ' part seeding events #   : ', tpart
#endif

  ! NOW we initialize each particle's values

  ! FIRST open sequence of repeated initialization...
  do t = 1, tpart

     i = 0

     ! ... then get random positions on the grid...
     do k = 1, npart0      ! Cell position
        !write(*,'(a,i10)') ' k  : ', k !just to debug
        do l = 1, part_fac ! Individual particle within cell

           ! ... and duplicate particles for each sequence of initialization
           i = i + 1

           if (t==1) then

              !initialize position
#ifdef FTLE
! in case of ftle, there is no randomness added to the particles positions :

              xpo(i) = real(ipart(k))
              ypo(i) = real(jpart(k))

              xgrid(i) = max(1, min(m, ceiling(xpo(i))))
              ygrid(i) = max(1, min(n, ceiling(ypo(i))))

              if (nlayer(xgrid(i), ygrid(i)) == 0) stop '!!! STOP: one partile in the initial position is in a dry cell !!!'

#else
3             call rng_uni(r1)
              call rng_uni(r2)

              xpo(i) = real(ipart(k)) - 1. + r1
              ypo(i) = real(jpart(k)) - 1. + r2

              xgrid(i) = max(1, min(m, ceiling(xpo(i))))
              ygrid(i) = max(1, min(n, ceiling(ypo(i))))

              do while (nlayer(xgrid(i), ygrid(i)) == 0)
                 goto 3 ! check that the particle is in a wet cell !
              end do
#endif
              !write(*,'(a,i10)') ' wet '  !just to debug

              ! Find the mode proportion of each particle place.  The mode proportion btw 0 and 1 is the proportion of particles
              ! that should be associated to the surface mode. This mode proportion can depend or not on
              ! The the local bottom depth. If not, this proportion just depends on the type of particle of the model,
              ! and will therefore be the same for each particle.
              mode_prop(i) = LFPart%mode_prop_fun( LFPart%mode_prop, mid_depth(nlayer(xgrid(i),ygrid(i))), in_diapause )

              call rng_uni(r3) !random number between 0 and 1

              ! if random number is < mode proportion, -> the particle is in the the upper mode (1).
              ! Else, is is in the lower mode (2)
              if (r3 < mode_prop(i)) then
                 mode(i) = 1
              else
                 mode(i) = 2
              endif


              zz(i)  = LFPart%zday_fun  (mid_depth(nlayer(xgrid(i), ygrid(i))), mode(i), in_diapause)
              zzz(i) = LFPart%znight_fun(mid_depth(nlayer(xgrid(i), ygrid(i))), mode(i), in_diapause)
              sigma_local(i) = LFPart%sigma_fun(mid_depth(nlayer(xgrid(i), ygrid(i))), mode(i))
              zpo(i) = vert_mig(0., zz(i), zzz(i), mid_depth(nlayer(xgrid(i), ygrid(i))), &
                                theta, sigma_local(i), dvm)

              ! If the particle's depth is under bottom depth, put it randomly in a buffer 50m under bottom depth :
              ! Or in a smaller buffer is the bdeph is shallower than 50m
              if ( zpo(i) > mid_depth(nlayer(xgrid(i), ygrid(i))) ) then
                 call rng_uni(r4) !random number between 0 and 1
                 buffer_lim = min( mid_depth(nlayer(xgrid(i), ygrid(i))), 10.0 )
                 zpo(i) = (r4*buffer_lim) + ( mid_depth(nlayer(xgrid(i), ygrid(i)))-buffer_lim )
              endif

              ! If the particle's depth is over the surface, put it randomly in a buffer 50m under surface :
              ! Or in a smaller buffer is the bdeph is shallower than 50m
              if ( zpo(i) < 0 ) then
                 call rng_uni(r5) !random number between 0 and 1
                 buffer_lim = min( mid_depth(nlayer(xgrid(i), ygrid(i))), 10.0 )
                 zpo(i) = r5*buffer_lim
              endif

              ! Loop to find zgrid(i), the depth layer corresponding to the depth zpo(i) :
              zloop: do j = 1, min( ilo, nlayer(xgrid(i),ygrid(i)) )
                 zgrid(i) = j
                 if (zpo(i)<=mid_depth(j)) exit zloop
              enddo zloop

#ifdef BACK
              ! Days from the run start
              kd = kd_ini + istart*dtnc/86400 - (kd_start+(t-1)*part_freq)

              ! Initial particle timestep (kstep); use dt in sec
              istart_part = int( (kd*86400 - part_start(3)*3600 - part_start(2)*60 - part_start(1))/dt )
#else
              ! Days from the run start
              kd = kd_start - kd_ini - istart*dtnc/86400

              ! Initial particle timestep (kstep); use dt in sec
              istart_part = int( (kd*86400 + part_start(3)*3600 + part_start(2)*60 + part_start(1))/dt )
#endif
              tp(i) = temp(xgrid(i),ygrid(i),zgrid(i))
              
           endif

           LFPart%part_spec = part_spec

           LFPart%istart    = istart_part
           LFPart%isteps    = isteps_part

           LFPart%dvm       = dvm
           LFPart%light_day = light_day
           LFPart%zday      = zz(i)
           LFPart%znight    = zzz(i)

           LFPart%cell_out  = cell_out

           LFPart%xpo       = xpo(i)
           LFPart%ypo       = ypo(i)
           LFPart%zpo       = zpo(i)
           LFPart%bdep      = mid_depth(nlayer(xgrid(i), ygrid(i)))

#ifdef FTLE
           LFPart%ipo       = ipart(i)
           LFPart%jpo       = jpart(i)
#endif

           LFPart%mode      = mode(i)
           LFPart%mode_prop = mode_prop(i)

           LFPart%xk1       = xk1
           LFPart%xk2       = xk2
           LFPart%xk3       = xk3
           LFPart%xk4       = xk4

           LFPart%yk1       = yk1
           LFPart%yk2       = yk2
           LFPart%yk3       = yk3
           LFPart%yk4       = yk4

           LFPart%zk1       = zk1
           LFPart%zk2       = zk2
           LFPart%zk3       = zk3
           LFPart%zk4       = zk4

           call LFList%addPart(LFPart)

        end do
     end do
  end do

  npart = npart0 * part_fac

  deallocate(ipart,jpart)

  write(*,'(/a,i6,a)') ' Generated ', npart, ' particles'

end subroutine generateLate_finmarchicus




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!72




subroutine generateYoung_finmarchicus(YFList, runList)
  type(Young_finmarchicus)                    :: YFPart
  class(partList)              :: YFList


  character(*), intent(in)     :: runList

  integer                      :: i, j, ii, k, l, ioerr, t, IOstatus ! counters or iterators

  integer                      :: kd_start, kd_end, kd

  integer, dimension(6)        :: part_start, part_end

  integer, dimension(max_part) :: xgrid, ygrid, zgrid

  integer                     :: buffer_lim

  real                         :: r1, r2, r3, r4, r5

  real,    dimension(max_part) :: xpo, ypo, zpo, zz, zzz, tp, sigma_local
  integer,    dimension(max_part) :: mode
  real,       dimension(max_part) :: mode_prop


  !------------------------------------------------------------------72

  open(10,file = runList, status = 'old', action = 'read', err = 1, iostat = ioerr)
1   if (ioerr /= 0) STOP 'problem opening file'
  read(10, nml = part_nml)
  close(10)

!------------------------------------------------------------------72
! Initial particles distribution (spatial)
!
! Works with input_particles (below) to provide the # of particles
! added at specific time steps in case of periodic seeding of particles

  open(10,file = partfile, status = 'old', action = 'read', err = 5, iostat = ioerr)
  5   if (ioerr /= 0) STOP 'problem opening file, particle input'

  do i = 1, m
     read(10, *) (part_grid_init(i, j), j = 1, n)
  end do

  npart0  = sum(part_grid_init) !npart0 is the number of cells with particles

  ! allocate the vectors used for the cells with particles
  allocate(ipart(npart0))
  allocate(jpart(npart0))

  ! record the coordinates of the cells with particles
  ii = 0
  outer : do j = 1, n    !n is a horizontal dimension of the grid
     inner : do i = 1, m !m is the other horizontal dimension of the grid
        if(part_grid_init(i, j) == 1) then
           ii = min(ii + 1, npart0)
           ipart(ii) = i
           jpart(ii) = j
        endif
     end do inner
  end do outer

  close(10)

  ! Number of individuals generated
  write(*,'(/a,i9,a,/)') ' Cells with particles = ', npart0, '; now generating all the particles variables'

  !------------------------------------------------------------------72
  ! Initialization of individual particles

  ! Initializing time (time would be tzero + one dt)
  call update2(initial)

  ! Initialize time loop for periodic release of particles:
  ! start of particles release -> end of release (by frequency of release)

  part_start = (/ 0, 0, 0, d_Astart, m_Astart, y_Astart /) ! Start of particles release
  call cday2(part_start(4),part_start(5),part_start(6),kd_start)

  part_end   = (/ 0, 0, 0, d_Aend, m_Aend, y_Aend /)       ! End of particles release
  call cday2(part_end(4),part_end(5),part_end(6),kd_end)

#ifdef DEBUG
  write(*,'(/a,i2,1x,i2,1x,i4,a,i6)') ' part start              : ', part_start(4), part_start(5), part_start(6), &
                                      ';   start timestep : ', kd_start
  write(*,'(a,i2,1x,i2,1x,i4,a,i6)')  ' part end                : ', part_end(4),   part_end(5),   part_end(6),   &
                                      ';   end timestep   : ', kd_end
#endif

  ! Duration of particles advection; same for all particles
  isteps_part = part_duration * 86400/dt

#ifdef DEBUG
  write(*,'(a,i10,a)') ' part advection duration : ', isteps_part, ' (inner time steps)'
#endif

  ! Periodic initialization
#ifdef BACK
  tpart = ceiling( real(kd_start - kd_end+1)/part_freq)
#else
  tpart = ceiling( real(kd_end - kd_start+1)/part_freq )
#endif

#ifdef DEBUG
  write(*,'(a,i10)') ' part seeding events #   : ', tpart
#endif

  ! NOW we initialize each particle's values

  ! FIRST open sequence of repeated initialization...
  do t = 1, tpart

     i = 0

     ! ... then get random positions on the grid...
     do k = 1, npart0      ! Cell position
        !write(*,'(a,i10)') ' k  : ', k !just to debug
        do l = 1, part_fac ! Individual particle within cell

           ! ... and duplicate particles for each sequence of initialization
           i = i + 1

           if (t==1) then

              !initialize position
#ifdef FTLE
! in case of ftle, there is no randomness added to the particles positions :

              xpo(i) = real(ipart(k))
              ypo(i) = real(jpart(k))

              xgrid(i) = max(1, min(m, ceiling(xpo(i))))
              ygrid(i) = max(1, min(n, ceiling(ypo(i))))

              if (nlayer(xgrid(i), ygrid(i)) == 0) stop '!!! STOP: one partile in the initial position is in a dry cell !!!'

#else
3             call rng_uni(r1)
              call rng_uni(r2)

              xpo(i) = real(ipart(k)) - 1. + r1
              ypo(i) = real(jpart(k)) - 1. + r2

              xgrid(i) = max(1, min(m, ceiling(xpo(i))))
              ygrid(i) = max(1, min(n, ceiling(ypo(i))))

              do while (nlayer(xgrid(i), ygrid(i)) == 0)
                 goto 3 ! check that the particle is in a wet cell !
              end do
#endif
              !write(*,'(a,i10)') ' wet '  !just to debug

              ! Find the mode proportion of each particle place.  The mode proportion btw 0 and 1 is the proportion of particles
              ! that should be associated to the surface mode. This mode proportion can depend or not on
              ! The the local bottom depth. If not, this proportion just depends on the type of particle of the model,
              ! and will therefore be the same for each particle.
              mode_prop(i) = YFPart%mode_prop_fun( YFPart%mode_prop, mid_depth(nlayer(xgrid(i),ygrid(i))) )

              call rng_uni(r3) !random number between 0 and 1

              ! if random number is < mode proportion, -> the particle is in the the upper mode (1).
              ! Else, is is in the lower mode (2)
              if (r3 < mode_prop(i)) then
                 mode(i) = 1
              else
                 mode(i) = 2
              endif


              zz(i)  = YFPart%zday_fun  (mid_depth(nlayer(xgrid(i), ygrid(i))), mode(i))
              zzz(i) = YFPart%znight_fun(mid_depth(nlayer(xgrid(i), ygrid(i))), mode(i))
              sigma_local(i) = YFPart%sigma_fun(mid_depth(nlayer(xgrid(i), ygrid(i))), mode(i))
              zpo(i) = vert_mig(0., zz(i), zzz(i), mid_depth(nlayer(xgrid(i), ygrid(i))), &
                                theta, sigma_local(i), dvm)

              ! If the particle's depth is under bottom depth, put it randomly in a buffer 50m over bottom depth :
              if ( zpo(i) > mid_depth(nlayer(xgrid(i), ygrid(i))) ) then
                 call rng_uni(r4) !random number between 0 and 1
                 buffer_lim = min( mid_depth(nlayer(xgrid(i), ygrid(i))), 10.0 )
                 zpo(i) = (r4*buffer_lim) + ( mid_depth(nlayer(xgrid(i), ygrid(i)))-buffer_lim )
              endif

              ! If the particle's depth is over the surface, put it randomly in a buffer 50m under surface :
              if ( zpo(i) < 0 ) then
                 call rng_uni(r5) !random number between 0 and 1
                 buffer_lim = min( mid_depth(nlayer(xgrid(i), ygrid(i))), 10.0 )
                 zpo(i) = r5*buffer_lim
              endif


              ! Loop to find zgrid(i), the depth layer corresponding to the depth zpo(i) :
              zloop: do j = 1, min( ilo, nlayer(xgrid(i),ygrid(i)) )
                 zgrid(i) = j
                 if (zpo(i)<=mid_depth(j)) exit zloop
              enddo zloop

#ifdef BACK
              ! Days from the run start
              kd = kd_ini + istart*dtnc/86400 - (kd_start+(t-1)*part_freq)

              ! Initial particle timestep (kstep); use dt in sec
              istart_part = int( (kd*86400 - part_start(3)*3600 - part_start(2)*60 - part_start(1))/dt )
#else
              ! Days from the run start
              kd = kd_start - kd_ini - istart*dtnc/86400

              ! Initial particle timestep (kstep); use dt in sec
              istart_part = int( (kd*86400 + part_start(3)*3600 + part_start(2)*60 + part_start(1))/dt )
#endif
              tp(i) = temp(xgrid(i),ygrid(i),zgrid(i))
              
           endif

           YFPart%part_spec = part_spec

           YFPart%istart    = istart_part
           YFPart%isteps    = isteps_part

           YFPart%dvm       = dvm
           YFPart%light_day = light_day
           YFPart%zday      = zz(i)
           YFPart%znight    = zzz(i)

           YFPart%cell_out  = cell_out

           YFPart%xpo       = xpo(i)
           YFPart%ypo       = ypo(i)
           YFPart%zpo       = zpo(i)
           YFPart%bdep      = mid_depth(nlayer(xgrid(i), ygrid(i)))

#ifdef FTLE
           YFPart%ipo       = ipart(i)
           YFPart%jpo       = jpart(i)
#endif
           YFPart%mode      = mode(i)
           YFPart%mode_prop = mode_prop(i)

           YFPart%xk1       = xk1
           YFPart%xk2       = xk2
           YFPart%xk3       = xk3
           YFPart%xk4       = xk4

           YFPart%yk1       = yk1
           YFPart%yk2       = yk2
           YFPart%yk3       = yk3
           YFPart%yk4       = yk4

           YFPart%zk1       = zk1
           YFPart%zk2       = zk2
           YFPart%zk3       = zk3
           YFPart%zk4       = zk4

           call YFList%addPart(YFPart)

        end do
     end do
  end do

  npart = npart0 * part_fac

  deallocate(ipart,jpart)

  write(*,'(/a,i6,a)') ' Generated ', npart, ' particles'

end subroutine generateYoung_finmarchicus





!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine output_ftle_netcdf(list, part_spec, action)

  ! Output particles trajectories in netcdf format, according to the
  ! requiremetns of Finite Time Lyapunov Exponents computation, i.e.
  ! 1 particle at the center of each grid cell.
  !   Argument :
  !   action = 0 => initialization
  !   action = 1 => Record data
  !   action = 2 => Finalization



  use netcdf

  implicit none

  class(partList)     :: list
  class(*), pointer   :: curr

  integer, intent(in) :: part_spec, action

  integer :: ioerr, xdimid, ydimid, zdimid, timedimid, &
             lonid, latid, xid, yid, zid, nb

  !integer :: ncid, timeid, xtid, ytid, ncinc, nctlen

  integer, save :: ncid, timeid, xtid, ytid, ncinc, nctlen

  integer  :: j

  real     :: xx, yy
  integer  :: ii, jj

  real, dimension(m,n) :: xt, yt

  real, save :: ncfreq, nctime

  character(len=100) :: date_print, filename_nc

  character(len=3), dimension(12) :: month_nc = (/'JAN','FEB','MAR','APR','MAY','JUN','JUL','AUG','SEP','OCT','NOV','DEC'/)

  select case (action)

  case (0)

       ! First check that the # of particle in the list is correct
       nb = list%n_part()
       if (nb/=npart*tpart) then
          write(*,'(/a)') '!!! ERROR: particles number not correct !!!'
          write(*,'(3(a,i6))') 'nb = ', nb, ' npart = ', npart, ' tpart = ', tpart
          stop
       endif

     ! Create file

     filename_nc = trim(outfile)//"_traj.nc"
#ifdef DEBUG
     write(*,'(/a,a)'), ' Creation of netcdf output file : ', trim(filename_nc)
#endif

     ioerr = nf90_create(trim(filename_nc),nf90_clobber,ncid)
     if (ioerr/=nf90_noerr) print*,'Error opening netcdf file ',trim(nf90_strerror(ioerr))


     ! Dimensions

     ioerr = nf90_def_dim(ncid,"x",m,xdimid)
     if (ioerr/=nf90_noerr) print*,'Error defining X dim in netcdf file ',trim(nf90_strerror(ioerr))

     ioerr = nf90_def_dim(ncid,"y",n,ydimid)
     if (ioerr/=nf90_noerr) print*,'Error defining Y dim in netcdf file ',trim(nf90_strerror(ioerr))

     ioerr = nf90_def_dim(ncid,"time",nf90_unlimited,timedimid)
     if (ioerr/=nf90_noerr) print*,'Error defining T dim in netcdf file ',trim(nf90_strerror(ioerr))



     ! Variables

     ioerr = nf90_def_var(ncid,"time",nf90_double,(/timedimid/),timeid)
     if (ioerr/=nf90_noerr) print*,'Error defining T var ',trim(nf90_strerror(ioerr))

     ioerr = nf90_def_var(ncid,"lon",nf90_real,(/xdimid,ydimid/),lonid)
     if (ioerr/=nf90_noerr) print*,'Error defining Lon var ',trim(nf90_strerror(ioerr))

     ioerr = nf90_def_var(ncid,"lat",nf90_real,(/xdimid,ydimid/),latid)
     if (ioerr/=nf90_noerr) print*,'Error defining Lat var ',trim(nf90_strerror(ioerr))

     ioerr = nf90_def_var(ncid,"dx",nf90_real,(/xdimid,ydimid/),xid)
     if (ioerr/=nf90_noerr) print*,'Error defining dX var ',trim(nf90_strerror(ioerr))

     ioerr = nf90_def_var(ncid,"dy",nf90_real,(/xdimid,ydimid/),yid)
     if (ioerr/=nf90_noerr) print*,'Error defining dY var ',trim(nf90_strerror(ioerr))

     ioerr = nf90_def_var(ncid,"xt",nf90_real,(/xdimid,ydimid,timedimid/),xtid)
     if (ioerr/=nf90_noerr) print*,'Error defining XT var ',trim(nf90_strerror(ioerr))

     ioerr = nf90_def_var(ncid,"yt",nf90_real,(/xdimid,ydimid,timedimid/),ytid)
     if (ioerr/=nf90_noerr) print*,'Error defining YT var ',trim(nf90_strerror(ioerr))



     ! Attributes

     ioerr = nf90_put_att(ncid,NF90_GLOBAL,"File_name",filename_nc)
     if (ioerr/=nf90_noerr) print*,'Error defining FileName att ',trim(nf90_strerror(ioerr))

     ioerr = nf90_put_att(ncid,NF90_GLOBAL,"Description",'Part lagrangian experiment')
     if (ioerr/=nf90_noerr) print*,'Error defining Description att ',trim(nf90_strerror(ioerr))

     ioerr = nf90_put_att(ncid,NF90_GLOBAL,"Version",'1.0')
     if (ioerr/=nf90_noerr) print*,'Error defining Version att ',trim(nf90_strerror(ioerr))

     write(date_print,'(i3,i3,i5,a,3(i2,a1))'),today(2),today(1),today(3),' - ',now(1),'h',now(2),'m',now(3),'s' ! JJ MM AAAA -- HHhMMmSSs

     ioerr = nf90_put_att(ncid,NF90_GLOBAL,"Date",trim(date_print))
     if (ioerr/=nf90_noerr) print*,'Error defining Date att ',trim(nf90_strerror(ioerr))


     ioerr = nf90_put_att(ncid,timedimid,"point_spacing",'even')
     if (ioerr/=nf90_noerr) print*,'Error defining point_spacing att for T var ',trim(nf90_strerror(ioerr))

     ioerr = nf90_put_att(ncid,timeid,"Long_name",'Days')
     if (ioerr/=nf90_noerr) print*,'Error defining Long_name att for T var ',trim(nf90_strerror(ioerr))

     ioerr = nf90_put_att(ncid,timeid,"Units",trim(tnc_units))
     if (ioerr/=nf90_noerr) print*,'Error defining Units att for T var ',trim(nf90_strerror(ioerr))

     ioerr = nf90_put_att(ncid,timeid,"time_origin",trim(tnc_origin_out))
     if (ioerr/=nf90_noerr) print*,'Error defining time_origin att for T var ',trim(nf90_strerror(ioerr))


     !ioerr = nf90_put_att(ncid,zdimid,"Positive",'down')
     !if (ioerr/=nf90_noerr) print*,'Error defining Positive att for Z var ',trim(nf90_strerror(ioerr))


     ioerr = nf90_put_att(ncid,lonid,"Long_name",'Longitude')
     if (ioerr/=nf90_noerr) print*,'Error defining Long_name att for Lon var ',trim(nf90_strerror(ioerr))

     ioerr = nf90_put_att(ncid,lonid,"Units",'degrees_east')
     if (ioerr/=nf90_noerr) print*,'Error defining Units att for Lon var ',trim(nf90_strerror(ioerr))

     !ioerr = nf90_put_att(ncid,lonid,"valid_min",-7.138476e+01)
     !if (ioerr/=nf90_noerr) print*,'Error defining Units att for Lon var ',trim(nf90_strerror(ioerr))

     !ioerr = nf90_put_att(ncid,lonid,"valid_max",-6.299908e+01)
     !if (ioerr/=nf90_noerr) print*,'Error defining valid_max att for Lon var ',trim(nf90_strerror(ioerr))


     ioerr = nf90_put_att(ncid,latid,"Long_name",'Latitude')
     if (ioerr/=nf90_noerr) print*,'Error defining Long_name att for Lat var ',trim(nf90_strerror(ioerr))

     ioerr = nf90_put_att(ncid,latid,"Units",'degrees_north')
     if (ioerr/=nf90_noerr) print*,'Error defining Units att for Lat var ',trim(nf90_strerror(ioerr))

     !ioerr = nf90_put_att(ncid,latid,"valid_min",3.852767e+01)
     !if (ioerr/=nf90_noerr) print*,'Error defining Units att for Lat var ',trim(nf90_strerror(ioerr))

     !ioerr = nf90_put_att(ncid,latid,"valid_max",4.223735e+01)
     !if (ioerr/=nf90_noerr) print*,'Error defining valid_max att for Lat var ',trim(nf90_strerror(ioerr))


     ioerr = nf90_put_att(ncid,xid,"Long_name",'Zonal cell width')
     if (ioerr/=nf90_noerr) print*,'Error defining Long_name att for dX var ',trim(nf90_strerror(ioerr))

     ioerr = nf90_put_att(ncid,xid,"Units",'m')
     if (ioerr/=nf90_noerr) print*,'Error defining Units att for dX var ',trim(nf90_strerror(ioerr))


     ioerr = nf90_put_att(ncid,yid,"Long_name",'Meridional cell width')
     if (ioerr/=nf90_noerr) print*,'Error defining Long_name att for dY var ',trim(nf90_strerror(ioerr))

     ioerr = nf90_put_att(ncid,yid,"Units",'m')
     if (ioerr/=nf90_noerr) print*,'Error defining Units att for dY var ',trim(nf90_strerror(ioerr))


     ioerr = nf90_put_att(ncid,xtid,"Long_name",'Zonal position')
     if (ioerr/=nf90_noerr) print*,'Error defining Long_name att for XT var ',trim(nf90_strerror(ioerr))

     ioerr = nf90_put_att(ncid,xtid,"Units",'Grid coordinates')
     if (ioerr/=nf90_noerr) print*,'Error defining Units att for XT var ',trim(nf90_strerror(ioerr))

     ioerr = nf90_put_att(ncid,xtid,"_FillValue",nf90_fill_real)
     if (ioerr/=nf90_noerr) print*,'Error defining fillvalue att for XT var ',trim(nf90_strerror(ioerr))


     ioerr = nf90_put_att(ncid,ytid,"Long_name",'Meridional position')
     if (ioerr/=nf90_noerr) print*,'Error defining Long_name att for YT var ',trim(nf90_strerror(ioerr))

     ioerr = nf90_put_att(ncid,ytid,"Units",'Grid coordinates')
     if (ioerr/=nf90_noerr) print*,'Error defining Units att for YT var ',trim(nf90_strerror(ioerr))

     ioerr = nf90_put_att(ncid,ytid,"_FillValue",nf90_fill_real)
     if (ioerr/=nf90_noerr) print*,'Error defining fillvalue att for YT var ',trim(nf90_strerror(ioerr))


     ! End definition

     ioerr = nf90_enddef(ncid)
     if (ioerr/=nf90_noerr) print*,'Error exiting define mode in netcdf file ',trim(nf90_strerror(ioerr))



     ! Put first timestep
      ncinc = 1

      ncfreq = output_freq2*dt
#ifdef DEBUG
      write(*,'(/a,f6.0,/)') ' Output frequency (in sec) : ', ncfreq
#endif

      nctime = dble(istart) * dble(dtnc)

      ioerr = nf90_put_var(ncid, timeid, nctime, (/ncinc/))
      if (ioerr/=nf90_noerr) print*, ' Error writing T var ', trim(nf90_strerror(ioerr))


     ! Put lon/lat

     ioerr = nf90_put_var(ncid,lonid,lon,(/ 1, 1 /),(/ m, n /))
     if (ioerr/=nf90_noerr) print*,'Error writing Lon var ',trim(nf90_strerror(ioerr))

     ioerr = nf90_put_var(ncid,latid,lat,(/ 1, 1 /),(/ m, n /))
     if (ioerr/=nf90_noerr) print*,'Error writing Lat var ',trim(nf90_strerror(ioerr))

     ! Put cells width

     ioerr = nf90_put_var(ncid,xid,dlx,(/ 1, 1 /),(/ m, n /))
     if (ioerr/=nf90_noerr) print*,'Error writing dX var ',trim(nf90_strerror(ioerr))

     ioerr = nf90_put_var(ncid,yid,dly,(/ 1, 1 /),(/ m, n /))
     if (ioerr/=nf90_noerr) print*,'Error writing dY var ',trim(nf90_strerror(ioerr))




     ! Put initial particles positions
     
     xt = nf90_fill_real
     
     yt = nf90_fill_real



     call list%reset()

     do while(list%moreValues())
        curr => list%currentValue()

        select type(curr)
        
        type is (Late_hyperboreus)

           xx = curr%get_xpo()
           yy = curr%get_ypo()
           ii = curr%get_ipo()
           jj = curr%get_jpo()


        type is (Late_finmarchicus)

           xx = curr%get_xpo()
           yy = curr%get_ypo()
           ii = curr%get_ipo()
           jj = curr%get_jpo()


        type is (Young_finmarchicus)

           xx = curr%get_xpo()
           yy = curr%get_ypo()
           ii = curr%get_ipo()
           jj = curr%get_jpo()

        end select

           xt(ii,jj) = xx
           yt(ii,jj) = yy

        call list%next()
     end do


     ioerr = nf90_put_var(ncid,xtid,xt,(/ 1, 1, ncinc /),(/ m, n, 1 /))
     if (ioerr/=nf90_noerr) print*,'Error writing XT var ',trim(nf90_strerror(ioerr))

     ioerr = nf90_put_var(ncid,ytid,yt,(/ 1, 1, ncinc /),(/ m, n, 1 /))
     if (ioerr/=nf90_noerr) print*,'Error writing YT var ',trim(nf90_strerror(ioerr))





!________________________________________________________________________________________________




  case (1)

     !     print*,'!V! arrived in case 1'

     ncinc = ncinc+1

     ! Put timestep

#ifdef BACK
       nctime = nctime - ncfreq
#else
       nctime = nctime + ncfreq
#endif
     !print*,'!V! actualised nctime'

     ioerr = nf90_put_var(ncid,timeid,nctime,(/ncinc/))
     if (ioerr/=nf90_noerr) print*,'Error writing T var ',trim(nf90_strerror(ioerr))

     !print*,'!V! saved time'

     ! Put particles positions

     xt = nf90_fill_real
     
     yt = nf90_fill_real


     call list%reset()

     !print*,'!V! before output loop'

     do while(list%moreValues())
        curr => list%currentValue()

        select type(curr)
        
        type is (Late_hyperboreus)

           xx = curr%get_xpo()
           yy = curr%get_ypo()
           ii = curr%get_ipo()
           jj = curr%get_jpo()


        type is (Late_finmarchicus)

           xx = curr%get_xpo()
           yy = curr%get_ypo()
           ii = curr%get_ipo()
           jj = curr%get_jpo()


        type is (Young_finmarchicus)

           xx = curr%get_xpo()
           yy = curr%get_ypo()
           ii = curr%get_ipo()
           jj = curr%get_jpo()

        end select

           j = j + 1

           xt(ii,jj) = xx
           yt(ii,jj) = yy

        call list%next()
     end do

     !print*,'!V! after output loop'

     ioerr = nf90_put_var(ncid,xtid,xt,(/ 1, 1, ncinc /),(/ m, n, 1 /))
     if (ioerr/=nf90_noerr) print*,'Error writing XT var ',trim(nf90_strerror(ioerr))

     ioerr = nf90_put_var(ncid,ytid,yt,(/ 1, 1, ncinc /),(/ m, n, 1 /))
     if (ioerr/=nf90_noerr) print*,'Error writing YT var ',trim(nf90_strerror(ioerr))

     
     !print*,'!V! after save'



!________________________________________________________________________________________________




  case (2) ! Close

     ioerr = nf90_close(ncid)

  end select


  return

end subroutine output_ftle_netcdf



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!72




  subroutine netcdfPart(list, part_spec, action)

    ! Output particles trajectories in netcdf format
    !   Argument :
    !
    !   part_spec  => class of particle (species)
    !
    !   action = 0 => initialization
    !   action = 1 => Record data
    !   action = 2 => Finalization

    use netcdf

    implicit none

    class(partList)     :: list
    class(*), pointer   :: curr

    integer, intent(in) :: part_spec, action

    ! variables for the netcdf procedures
    character(100)      :: filename_nc, date_print

    integer             :: ierr, i, j, nb

    integer, save       :: xdimid, ydimid, zdimid, timedimid, reprodimid

    integer, save       :: ncinc, nctlen 

    !------------------------------------------------------------------72
    ! Particle-specific variables IDs
    integer, save       :: file_id, attributesPartid, timeid, &
                           xpoid, ypoid, zpoid, bdepid, lonid, latid, &
                           znightid, speciesid,               &
    !<><><><><>
    ! USERS MAY CHANGE THE VARIABLES TO SAVE HERE 
                           modeid, mode_propid, tpid, sexid, massid, lipidid,      &
                           stageid, foodid, part_statusid,    &
                           abundid, activityid, ageid,        &
                           reproid, cum_eggsid, phyto3did,    &
                           epid, serialid, materid, birthdateid
    !<><><><><>
    !------------------------------------------------------------------72

    character(len=3), dimension(12)    :: month_nc = (/'JAN','FEB','MAR','APR','MAY','JUN','JUL','AUG','SEP','OCT','NOV','DEC'/)

    ! buffer arrays to allocate values to netcdf variables
    integer, dimension(:), allocatable :: all_species
    real,    save                      :: ncfreq 
    real*8,  save                      :: nctime 
    real,    dimension(:), allocatable :: all_xpo
    real,    dimension(:), allocatable :: all_lon
    real,    dimension(:), allocatable :: all_ypo
    real,    dimension(:), allocatable :: all_lat
    real,    dimension(:), allocatable :: all_zpo
    real,    dimension(:), allocatable :: all_bdep
    real,    dimension(:), allocatable :: all_znight

    !<><><><><>
    !!! DEFINE BELOW CLASS-SPECIFIC VARIABLES TO SAVE !!!
    integer, dimension(:), allocatable :: all_mode
    integer, dimension(:), allocatable :: all_mode_prop
    integer, dimension(:), allocatable :: all_sex
    integer, dimension(:), allocatable :: all_part_status
    real,    dimension(:), allocatable :: all_tp
    real,    dimension(:), allocatable :: all_mass
    real,    dimension(:), allocatable :: all_lipid
    real,    dimension(:), allocatable :: all_stage
    real,    dimension(:), allocatable :: all_abund
    integer, dimension(:), allocatable :: all_activity
#ifdef YUMMY
    real,    dimension(:), allocatable :: all_fd
    real,    dimension(:), allocatable :: all_ph
#endif
    real,    dimension(:), allocatable :: all_age
    real,  dimension(:,:), allocatable :: all_repro
    real,    dimension(:), allocatable :: all_cum_eggs
    real,    dimension(:), allocatable :: all_ep
    integer, dimension(:), allocatable :: all_serial
    integer, dimension(:), allocatable :: all_mater
    integer, dimension(:), allocatable :: all_birthdate
    !!! DEFINE ABOVE CLASS-SPECIFIC VARIABLES TO SAVE !!!
    !<><><><><>

    select case (action)

    case (0)

       ! First check that the # of particle in the list is correct
       nb = list%n_part()
       if (nb/=npart*tpart) then
          write(*,'(/a)') '!!! ERROR: particles number not correct !!!'
          write(*,'(3(a,i6))') 'nb = ', nb, ' npart = ', npart, ' tpart = ', tpart
          stop
       endif

       ! Creation of the netcdf dataset
       filename_nc = trim(outfile)//"_traj.nc"
#ifdef DEBUG
       write(*,'(/a,a)'), ' Creation of netcdf output file : ', trim(filename_nc)
#endif

       ierr = nf90_create(trim(filename_nc), nf90_clobber, file_id)
       if (ierr/=nf90_noerr) print*, 'ERROR NETCDF CREATION ', trim(nf90_strerror(ierr))

       ! Dimensions

       ! X axis here is for each particles
       ierr = nf90_def_dim(file_id, "part", npart, xdimid)
       if (ierr/=nf90_noerr) print*, 'ERROR PARTICLES DIMENSION DEFINITION', trim(nf90_strerror(ierr))

       ! Y axis here is for each sequence of particles (initialized at a given period)
       ierr = nf90_def_dim(file_id, "sequence", tpart, ydimid)
       if (ierr/=nf90_noerr) print*, 'ERROR SEQUENCE DIMENSION DEFINITION', trim(nf90_strerror(ierr))

       ! Z axis here is for the 4 values of "pref" only
       !ierr = nf90_def_dim(file_id, "pref length", 4, zdimid)
       !if (ierr/=nf90_noerr) print*, 'ERROR PREF LENGTH DIMENSION DEFINTION', trim(nf90_strerror(ierr))

       ! output_freq not always = input_freq
       nctlen = final/output_freq2
   
       !dimension for reproduce information
       if (part_spec.eq.3) then
           ierr = nf90_def_dim(file_id, "repro", 6, reprodimid)
           if (ierr/=nf90_noerr) print*, 'ERROR REPRODUCE DIMENSION DEFINITION', trim(nf90_strerror(ierr))
       endif

       ! Time axis (unlimited)
       ierr = nf90_def_dim(file_id, 'time', nf90_unlimited, timedimid)    
       if (ierr/=nf90_noerr) print*, 'ERROR TIME DIMENSION DEFINITION HERE', trim(nf90_strerror(ierr))   


       ! creation of the variables for the files

       ierr = nf90_def_var(file_id, "species",  nf90_int,  (/xdimid, ydimid/), speciesid)
       if (ierr/=nf90_noerr) print*, 'ERROR SPECIES VARIABLE DEFINITION ', nf90_strerror(ierr)

       !<><><><><>
       !!! DEFINE BELOW CLASS-SPECIFIC VARIABLES TO SAVE
       if (part_spec.gt.0 ) then
          ierr = nf90_def_var(file_id, "sex",    nf90_int,   (/xdimid, ydimid, timedimid/), sexid)
          if (ierr/=nf90_noerr) print*, 'ERROR SEX VARIABLE DEFINITION ', nf90_strerror(ierr)

          ierr = nf90_def_var(file_id, "serial", nf90_int,   (/xdimid, ydimid/), serialid)
          if (ierr/=nf90_noerr) print*, 'ERROR SERIAL VARIABLE DEFINITION ', nf90_strerror(ierr)
          
          ierr = nf90_def_var(file_id, "birthdate", nf90_int,   (/xdimid, ydimid/), birthdateid)
          if (ierr/=nf90_noerr) print*, 'ERROR BIRTHDATE VARIABLE DEFINITION ', nf90_strerror(ierr)
          
          ierr = nf90_def_var(file_id, "mater",  nf90_int,   (/xdimid, ydimid/), materid)
          if (ierr/=nf90_noerr) print*, 'ERROR MATER VARIABLE DEFINITION ', nf90_strerror(ierr)
       end if
       !<><><><><>
       
       ! Time-dependent variables

       ierr = nf90_def_var(file_id, "time", nf90_double, (/timedimid/), timeid)
       if (ierr/=nf90_noerr) print*, "ERROR TIME VARIABLE DEFINITION ", nf90_strerror(ierr)

       ierr = nf90_def_var(file_id, "lon", nf90_float, (/xdimid, ydimid, timedimid/), lonid)
       if (ierr/=nf90_noerr) print*, 'ERROR LON VARIABLE DEFINITION ', trim(nf90_strerror(ierr))

       ierr = nf90_def_var(file_id, "lat", nf90_float, (/xdimid, ydimid, timedimid/), latid)
       if (ierr/=nf90_noerr) print*, 'ERROR LAT VARIABLE DEFINITION ', trim(nf90_strerror(ierr))

       ierr = nf90_def_var(file_id, "xpo",  nf90_float, (/xdimid, ydimid, timedimid/), xpoid)
       if (ierr/=nf90_noerr) print*, "ERROR XPO VARIABLE DEFINITION ", nf90_strerror(ierr)

       ierr = nf90_def_var(file_id, "ypo",  nf90_float, (/xdimid, ydimid, timedimid/), ypoid)
       if (ierr/=nf90_noerr) print*, "ERROR YPO VARIABLE DEFINITION ", nf90_strerror(ierr)

       ierr = nf90_def_var(file_id, "zpo",  nf90_float, (/xdimid, ydimid, timedimid/), zpoid)
       if (ierr/=nf90_noerr) print*, "ERROR ZPO VARIABLE DEFINITION ", nf90_strerror(ierr)

       ierr = nf90_def_var(file_id, "bdep",  nf90_float, (/xdimid, ydimid, timedimid/), bdepid)
       if (ierr/=nf90_noerr) print*, "ERROR bdep VARIABLE DEFINITION ", nf90_strerror(ierr)

       ierr = nf90_def_var(file_id, "znight",  nf90_float, (/xdimid, ydimid, timedimid/), znightid)
       if (ierr/=nf90_noerr) print*, "ERROR ZNIGHT VARIABLE DEFINITION ", nf90_strerror(ierr)

       !!! DEFINE BELOW CLASS-SPECIFIC VARIABLES TO SAVE
       if (part_spec.gt.0 ) then
          ierr = nf90_def_var(file_id, "tp",  nf90_float,   (/xdimid, ydimid, timedimid/), tpid)
          if (ierr/=nf90_noerr) print*, 'ERROR TP VARIABLE DEFINITION ',    nf90_strerror(ierr)

          ierr = nf90_def_var(file_id, "mass", nf90_float, (/xdimid, ydimid,timedimid/), massid)
          if (ierr/=nf90_noerr) print*, 'ERROR MASS VARIABLE DEFINITION ', nf90_strerror(ierr)
       end if

       if (part_spec.eq.3) then
          ierr = nf90_def_var(file_id, "lipid", nf90_float, (/xdimid, ydimid,timedimid/), lipidid)
          if (ierr/=nf90_noerr) print*, 'ERROR LIPID VARIABLE DEFINITION ', nf90_strerror(ierr)
   
          ierr = nf90_def_var(file_id, "stage", nf90_float, (/xdimid, ydimid,timedimid/), stageid)
          if (ierr/=nf90_noerr) print*, 'ERROR STAGE VARIABLE DEFINITION ', nf90_strerror(ierr)
  
#ifdef YUMY
          ierr = nf90_def_var(file_id, "food", nf90_float, (/xdimid, ydimid,timedimid/), foodid)
          if (ierr/=nf90_noerr) print*, 'ERROR FOOD VARIABLE DEFINITION ', nf90_strerror(ierr)
#endif
  
          ierr = nf90_def_var(file_id, "part_status", nf90_int,   (/xdimid, ydimid,timedimid/), part_statusid)
          if (ierr/=nf90_noerr) print*, 'ERROR PARTICLE STATUS VARIABLE DEFINITION ', nf90_strerror(ierr)
  
          ierr = nf90_def_var(file_id, "abund", nf90_float, (/xdimid, ydimid,timedimid/), abundid)
          if (ierr/=nf90_noerr) print*, 'ERROR ABUNDANCE VARIABLE DEFINITION ', nf90_strerror(ierr)
  
          ierr = nf90_def_var(file_id, "activity", nf90_int,   (/xdimid, ydimid, timedimid/), activityid)
          if (ierr/=nf90_noerr) print*, 'ERROR ACTIVITY VARIABLE DEFINITION ', nf90_strerror(ierr)

          ierr = nf90_def_var(file_id, "age", nf90_float, (/xdimid, ydimid,timedimid/), ageid)
          if (ierr/=nf90_noerr) print*, 'ERROR AGE VARIABLE DEFINITION ', nf90_strerror(ierr)
  
          ierr = nf90_def_var(file_id, "repro", nf90_float, (/reprodimid, xdimid, ydimid,timedimid/), reproid)
          if (ierr/=nf90_noerr) print*, 'ERROR REPRODUCE VARIABLE DEFINITION ', nf90_strerror(ierr)
  
          ierr = nf90_def_var(file_id, "accum_eggs", nf90_float, (/xdimid, ydimid,timedimid/), cum_eggsid)
          if (ierr/=nf90_noerr) print*, 'ERROR ACCUMULATED EGGS VARIABLE DEFINITION ', nf90_strerror(ierr)
  
          ierr = nf90_def_var(file_id, "daily_eggs", nf90_float, (/xdimid, ydimid,timedimid/), epid)
          if (ierr/=nf90_noerr) print*, 'ERROR Daily EGGS VARIABLE DEFINITION ', nf90_strerror(ierr)
  
#ifdef YUMMY
          ierr = nf90_def_var(file_id, "phyto", nf90_float, (/xdimid, ydimid,timedimid/), phyto3did)
          if (ierr/=nf90_noerr) print*, 'ERROR PHYTOPLANKTON VARIABLE DEFINITION ', nf90_strerror(ierr)
#endif
       end if
       !!! DEFINE ABOVE CLASS-SPECIFIC VARIABLES TO SAVE

       ! Attributes

       ierr = nf90_put_att(file_id, NF90_GLOBAL, "File_name", filename_nc)
       if (ierr/=nf90_noerr) print*, 'Error defining FileName att ', trim(nf90_strerror(ierr))

       ierr = nf90_put_att(file_id, NF90_GLOBAL, "Description", 'DFO Calanus Population Model')
       if (ierr/=nf90_noerr) print*, 'Error defining Description att ', trim(nf90_strerror(ierr))

       ierr = nf90_put_att(file_id, NF90_GLOBAL, "Version", '1.0')
       if (ierr/=nf90_noerr) print*, 'Error defining Version att ', trim(nf90_strerror(ierr))

       ! JJ MM AAAA -- HHhMMmSSs
       write(date_print, '(i3,i3,i5,a,3(i2,a1))') today(2), today(1), today(3), ' - ', now(1), &
            'h', now(2), 'm', now(3), 's'

       ierr = nf90_put_att(file_id, NF90_GLOBAL, "Date", trim(date_print))
       if (ierr/=nf90_noerr) print*, 'Error defining Date att ', trim(nf90_strerror(ierr))

       ierr = nf90_put_att(file_id, timeid, "Long_name", 'Days')
       if (ierr/=nf90_noerr) print*, 'Error defining Long_name att for T var ', trim(nf90_strerror(ierr))

       ierr = nf90_put_att(file_id, timeid, "Units", trim(tnc_units))
       if (ierr/=nf90_noerr) print*,'Error defining Units att for T var ', trim(nf90_strerror(ierr))

       ierr = nf90_put_att(file_id,timeid, "time_origin", trim(tnc_origin_out))
       if (ierr/=nf90_noerr) print*,'Error defining time_origin att for T var ', trim(nf90_strerror(ierr))

       ierr = nf90_put_att(file_id, xpoid, "Long_name", 'Zonal index')
       if (ierr/=nf90_noerr) print*, 'Error defining Long_name att for XPO var ', trim(nf90_strerror(ierr))

       ierr = nf90_put_att(file_id, ypoid, "Long_name", 'Meridional index')
       if (ierr/=nf90_noerr) print*, 'Error defining Long_name att for YPO var ', trim(nf90_strerror(ierr))

       ierr = nf90_put_att(file_id, lonid, "Long_name", 'Longitude')
       if (ierr/=nf90_noerr) print*, 'Error defining Long_name att for LON var ', trim(nf90_strerror(ierr))

       ierr = nf90_put_att(file_id, lonid, "Units", 'degrees_east')
       if (ierr/=nf90_noerr) print*, 'Error defining Units att for LON var ', trim(nf90_strerror(ierr))

       ierr = nf90_put_att(file_id, latid, "Long_name", 'Latitude')
       if (ierr/=nf90_noerr) print*, 'Error defining Long_name att for LAT var ', trim(nf90_strerror(ierr))

       ierr = nf90_put_att(file_id, latid, "Units", 'degrees_north')
       if (ierr/=nf90_noerr) print*, 'Error defining Units att for LAT var ', trim(nf90_strerror(ierr))

       ierr = nf90_put_att(file_id, znightid, "Long_name", 'Surface depth optimal')
       if (ierr/=nf90_noerr) print*,'Error defining Long_name att for ZNIGHT var ',trim(nf90_strerror(ierr))

       ierr = nf90_put_att(file_id, znightid, "Units", 'm')
       if (ierr/=nf90_noerr) print*, 'Error defining Units att for ZNIGHT var ', trim(nf90_strerror(ierr))

       ierr = nf90_put_att(file_id, zpoid, "Long_name", 'Depth')
       if (ierr/=nf90_noerr) print*,'Error defining Long_name att for ZPO var ',trim(nf90_strerror(ierr))

       ierr = nf90_put_att(file_id, bdepid, "Long_name", 'Bottom Depth')
       if (ierr/=nf90_noerr) print*,'Error defining Long_name att for bdep var ',trim(nf90_strerror(ierr))

       ierr = nf90_put_att(file_id, zpoid, "Units", 'm')
       if (ierr/=nf90_noerr) print*, 'Error defining Units att for ZPO var ', trim(nf90_strerror(ierr))

       ierr = nf90_put_att(file_id, bdepid, "Units", 'm')
       if (ierr/=nf90_noerr) print*, 'Error defining Units att for bdep var ', trim(nf90_strerror(ierr))

       ierr = nf90_put_att(file_id, speciesid, "Long_name", 'Species')
       if (ierr/=nf90_noerr) print*,'Error defining Long_name att for SPECIES var ',trim(nf90_strerror(ierr))

       !!! DEFINE BELOW CLASS-SPECIFIC VARIABLES TO SAVE
       if (part_spec.gt.0 ) then
          ierr = nf90_put_att(file_id,tpid,"Long_name",'Temperature')
          if (ierr/=nf90_noerr) print*,'Error defining Long_name att for Temp var',trim(nf90_strerror(ierr))

          ierr = nf90_put_att(file_id,tpid,"Units",'°C')
          if (ierr/=nf90_noerr) print*,'Error defining Units att for Temp var',trim(nf90_strerror(ierr))

          ierr = nf90_put_att(file_id,tpid,"_FillValue",nf90_fill_real)
          if (ierr/=nf90_noerr) print*,'Error defining fillvalue att for Temp var',trim(nf90_strerror(ierr))

          ierr = nf90_put_att(file_id, sexid, "Long_name", 'Sex')
          if (ierr/=nf90_noerr) print*,'Error defining Long_name att for SEX var ',trim(nf90_strerror(ierr))

          ierr = nf90_put_att(file_id, massid, "Long_name", 'Body mass')
          if (ierr/=nf90_noerr) print*,'Error defining Long_name att for MASS var ',trim(nf90_strerror(ierr))

          ierr = nf90_put_att(file_id, massid, "Units", 'ug C')
          if (ierr/=nf90_noerr) print*, 'Error defining Units att for MASS var ', trim(nf90_strerror(ierr))
       end if

       if (part_spec.ge.3  .and. part_spec.lt.5) then
          ierr = nf90_put_att(file_id, lipidid, "Long_name", 'Body lipid')
          if (ierr/=nf90_noerr) print*,'Error defining Long_name att for LIPID var ',trim(nf90_strerror(ierr))

          ierr = nf90_put_att(file_id, lipidid, "Units", 'ug C')
          if (ierr/=nf90_noerr) print*, 'Error defining Units att for LIPID var ', trim(nf90_strerror(ierr))   
   
          ierr = nf90_put_att(file_id, stageid, "Long_name", 'Development stage')
          if (ierr/=nf90_noerr) print*,'Error defining Long_name att for STAGE var ',trim(nf90_strerror(ierr))

          ierr = nf90_put_att(file_id, stageid, "Units", '-')
          if (ierr/=nf90_noerr) print*, 'Error defining Units att for STAGE var ', trim(nf90_strerror(ierr))
  
#ifdef YUMY
          ierr = nf90_put_att(file_id, foodid, "Long_name", 'Environmental food')
          if (ierr/=nf90_noerr) print*,'Error defining Long_name att for FOOD var ',trim(nf90_strerror(ierr))

          ierr = nf90_put_att(file_id, foodid, "Units", 'gCHLa')
          if (ierr/=nf90_noerr) print*, 'Error defining Units att for FOOD var ', trim(nf90_strerror(ierr))
  
          ierr = nf90_put_att(file_id,foodid,"_FillValue",nf90_fill_real)
          if (ierr/=nf90_noerr) print*,'Error defining fillvalue att for FOOD var',trim(nf90_strerror(ierr))
#endif
  
          ierr = nf90_put_att(file_id, part_statusid, "Long_name", 'Particle status')
          if (ierr/=nf90_noerr) print*,'Error defining Long_name att for particle status var ',trim(nf90_strerror(ierr))
  
          ierr = nf90_put_att(file_id, abundid, "Long_name", 'Abundance')
          if (ierr/=nf90_noerr) print*,'Error defining Long_name att for ABUNDANCE var ',trim(nf90_strerror(ierr))

          ierr = nf90_put_att(file_id, abundid, "Units", '-')
          if (ierr/=nf90_noerr) print*, 'Error defining Units att for ABUNDANCE var ', trim(nf90_strerror(ierr))
  
          ierr = nf90_put_att(file_id, activityid, "Long_name", 'Activity')
          if (ierr/=nf90_noerr) print*,'Error defining Long_name att for ACTIVITY var ',trim(nf90_strerror(ierr))
  
          ierr = nf90_put_att(file_id, ageid, "Long_name", 'Copepod age')
          if (ierr/=nf90_noerr) print*,'Error defining Long_name att for AGE var ',trim(nf90_strerror(ierr))

          ierr = nf90_put_att(file_id, ageid, "Units", 'days')
          if (ierr/=nf90_noerr) print*, 'Error defining Units att for AGE var ', trim(nf90_strerror(ierr))
  
          ierr = nf90_put_att(file_id, reproid, "Long_name", 'Reproduce information')
          if (ierr/=nf90_noerr) print*,'Error defining Long_name att for REPRODUCE var ',trim(nf90_strerror(ierr))
  
          ierr = nf90_put_att(file_id, reproid, "Variables", 'start_ep-duration_ep-counter_eggs-start_ep-duration_ep-counter_egg')
          if (ierr/=nf90_noerr) print*,'Error defining Long_name att for REPRODUCE var ',trim(nf90_strerror(ierr))
  
          ierr = nf90_put_att(file_id, cum_eggsid, "Long_name", 'Accumulated eggs')
          if (ierr/=nf90_noerr) print*,'Error defining Long_name att for ACCUMULATED EGGS var ',trim(nf90_strerror(ierr))
  
          ierr = nf90_put_att(file_id, epid, "Long_name", 'Daily eggs')
          if (ierr/=nf90_noerr) print*,'Error defining Long_name att for DAILY EGGS var ',trim(nf90_strerror(ierr))
  
#ifdef YUMY
          ierr = nf90_put_att(file_id, phyto3did, "Long_name", 'Environmental Phytoplankton')
          if (ierr/=nf90_noerr) print*,'Error defining Long_name att for PHYTOPLANKTON var ',trim(nf90_strerror(ierr))

          ierr = nf90_put_att(file_id, phyto3did, "Units", 'gCHLa')
          if (ierr/=nf90_noerr) print*, 'Error defining Units att for PHYTOPLANKTON var ', trim(nf90_strerror(ierr))
#endif
       
          ierr = nf90_put_att(file_id, serialid, "Long_name", 'Serial Number')
          if (ierr/=nf90_noerr) print*, 'Error defining Long_name att for SERIAL var ', trim(nf90_strerror(ierr))
       
          ierr = nf90_put_att(file_id, birthdateid, "Long_name", 'Birth Date (yyyymmdd)')
          if (ierr/=nf90_noerr) print*, 'Error defining Long_name att for BIRTHDATE var ', trim(nf90_strerror(ierr))
       
          ierr = nf90_put_att(file_id, materid, "Long_name", 'Mother Number')
          if (ierr/=nf90_noerr) print*, 'Error defining Long_name att for MATER var ', trim(nf90_strerror(ierr))
       end if
       
       ! end of file definition
       ierr = nf90_enddef(file_id)     
       if (ierr/=nf90_noerr) print*, 'ERROR END NETCDF FILE ', nf90_strerror(ierr)  


       ! Put first timestep
       ncinc = 1

       ncfreq = output_freq2*dt
#ifdef DEBUG
       write(*,'(/a,f6.0,/)') ' Output frequency (in sec) : ', ncfreq
#endif

       nctime = dble(istart) * dble(dtnc)

       ierr = nf90_put_var(file_id, timeid, nctime, (/ncinc/))
       if (ierr/=nf90_noerr) print*, ' Error writing T var ', trim(nf90_strerror(ierr))

       ! Allocate size of list variables
       allocate(all_xpo(npart))
       allocate(all_lon(npart))
       allocate(all_ypo(npart))
       allocate(all_lat(npart))
       allocate(all_zpo(npart))
       allocate(all_bdep(npart))
       allocate(all_species(npart))
       allocate(all_znight(npart))   

       if (part_spec .ge. 5) then
          allocate(all_mode(npart))
          allocate(all_mode_prop(npart))
       end if

       if (part_spec.gt.0 ) then
          allocate(all_tp(npart))
          allocate(all_sex(npart))
          allocate(all_mass(npart))          
       end if
       
       if (part_spec.ge.3  .and. part_spec.lt.5) then
          allocate(all_lipid(npart)) 
          allocate(all_stage(npart))
#ifdef YUMY
          allocate(all_fd(npart))
#endif
          allocate(all_part_status(npart))
          allocate(all_abund(npart)) 
          allocate(all_activity(npart))
          allocate(all_age(npart))  
          allocate(all_repro(6,npart)) 
          allocate(all_cum_eggs(npart))   
          allocate(all_ep(npart))
#ifdef YUMY
          allocate(all_ph(npart))
#endif
          allocate(all_serial(npart))
          allocate(all_birthdate(npart))
          allocate(all_mater(npart))
       end if
   

       call list%reset()

       i = 1
       j = 0

       do while(list%moreValues())
          curr => list%currentValue()

          select type(curr)

          type is (Zoo)
             all_species(i) = curr%get_species()
             all_xpo(i)     = curr%get_xpo()
             all_ypo(i)     = curr%get_ypo()
             all_zpo(i)     = curr%get_zpo()
             all_bdep(i)    = curr%get_bdep()
             all_znight(i)  = curr%get_znight()

          type is (Krill)
             all_species(i) = curr%get_species()
             all_xpo(i)     = curr%get_xpo()
             all_ypo(i)     = curr%get_ypo()
             all_zpo(i)     = curr%get_zpo()
             all_bdep(i)    = curr%get_bdep()
             all_znight(i)  = curr%get_znight()

             all_tp(i)      = curr%get_tp()
             all_sex(i)     = curr%get_sex()
             all_mass(i)    = curr%get_mass()
             
          type is (Calanus)
             all_species(i) = curr%get_species()
             all_xpo(i)     = curr%get_xpo()
             all_ypo(i)     = curr%get_ypo()
             all_zpo(i)     = curr%get_zpo()
             all_bdep(i)    = curr%get_bdep()
             all_znight(i)  = curr%get_znight()

             all_tp(i)          = curr%get_tp()
             all_sex(i)         = curr%get_sex()
             all_mass(i)        = curr%get_mass()
             all_lipid(i)       = curr%get_lipid()
             all_stage(i)       = curr%get_stage()
#ifdef YUMY
             all_fd(i)          = curr%get_fd()
#endif
             all_part_status(i) = curr%get_part_status()
             all_abund(i)       = curr%get_abund()
             all_activity(i)    = curr%get_activity()
             all_age(i)         = curr%get_age()
             all_cum_eggs(i)    = curr%get_cum_eggs()
             all_ep(i)          = curr%get_EP()
             all_repro(1:6,i)   = curr%get_reproduce_info()
#ifdef YUMY
             all_ph(i)          = curr%get_ph()
#endif
             all_serial(i)      = curr%get_serial()
             all_birthdate(i)   = curr%get_birthdate()
             all_mater(i)       = curr%get_mater()

          type is (Late_hyperboreus)
             all_species(i)  = curr%get_species()
             all_xpo(i)      = curr%get_xpo()
             all_ypo(i)      = curr%get_ypo()
             all_zpo(i)      = curr%get_zpo()
             all_bdep(i)     = curr%get_bdep()
             all_mode(i)      = curr%get_mode()
             all_mode_prop(i) = curr%get_mode_prop()
             all_znight(i)   = curr%get_znight()

          type is (Young_hyperboreus)
             all_species(i)  = curr%get_species()
             all_xpo(i)      = curr%get_xpo()
             all_ypo(i)      = curr%get_ypo()
             all_zpo(i)      = curr%get_zpo()
             all_bdep(i)     = curr%get_bdep()
             all_mode(i)      = curr%get_mode()
             all_mode_prop(i) = curr%get_mode_prop()
             all_znight(i)   = curr%get_znight()

          type is (Late_finmarchicus)
             all_species(i)  = curr%get_species()
             all_xpo(i)      = curr%get_xpo()
             all_ypo(i)      = curr%get_ypo()
             all_zpo(i)      = curr%get_zpo()
             all_bdep(i)      = curr%get_bdep()
             all_mode(i)      = curr%get_mode()
             all_mode_prop(i) = curr%get_mode_prop()
             all_znight(i)   = curr%get_znight()

          type is (Young_finmarchicus)
             all_species(i)  = curr%get_species()
             all_xpo(i)      = curr%get_xpo()
             all_ypo(i)      = curr%get_ypo()
             all_zpo(i)      = curr%get_zpo()
             all_bdep(i)     = curr%get_bdep()
             all_mode(i)      = curr%get_mode()
             all_mode_prop(i) = curr%get_mode_prop()
             all_znight(i)   = curr%get_znight()



          end select

          ! Write values to variables
          if ( mod(i, npart)==0 ) then

             j = j + 1

             ierr = nf90_put_var(file_id, speciesid, all_species, (/1, j/), (/npart, 1/))
             if (ierr/=nf90_noerr) print*,'ERROR IN THE SAVING OF SPECIES VARIABLE',nf90_strerror(ierr)

             !!! DEFINE BELOW CLASS-SPECIFIC VARIABLES TO SAVE
             if (part_spec.gt.0 ) then
                ierr = nf90_put_var(file_id, sexid, all_sex, (/1, j, ncinc/), (/npart, 1, 1/))
                if (ierr/=nf90_noerr) print*,'ERROR IN THE SAVING OF SEX VARIABLE ',nf90_strerror(ierr)
             end if            

             ierr = nf90_put_var(file_id, xpoid, all_xpo, (/1, j, ncinc/), (/npart, 1, 1/))
             if (ierr/=nf90_noerr) print*,'ERROR IN THE SAVING OF XPO VARIABLE ',nf90_strerror(ierr)

             ierr = nf90_put_var(file_id, lonid, all_lon, (/1, j, ncinc/), (/npart, 1, 1/))
             if (ierr/=nf90_noerr) print*,'ERROR IN THE SAVING OF LON VARIABLE',nf90_strerror(ierr)

             ierr = nf90_put_var(file_id, ypoid, all_ypo, (/1, j, ncinc/), (/npart, 1, 1/))
             if (ierr/=nf90_noerr) print*,'ERROR IN THE SAVING OF YPO VARIABLE',nf90_strerror(ierr)

             ierr = nf90_put_var(file_id, latid, all_lat, (/1, j, ncinc/), (/npart, 1, 1/))
             if (ierr/=nf90_noerr) print*,'ERROR IN THE SAVING OF LAT VARIABLE',nf90_strerror(ierr)

             ierr = nf90_put_var(file_id, zpoid, all_zpo, (/1, j, ncinc/), (/npart, 1, 1/))
             if (ierr/=nf90_noerr) print*,'ERROR IN THE SAVING OF ZPO VARIABLE',nf90_strerror(ierr)

             ierr = nf90_put_var(file_id, bdepid, all_bdep, (/1, j, ncinc/), (/npart, 1, 1/))
             if (ierr/=nf90_noerr) print*,'ERROR IN THE SAVING OF bdep VARIABLE',nf90_strerror(ierr)

             ierr = nf90_put_var(file_id, znightid, all_znight, (/1, j, ncinc/), (/npart, 1, 1/))
             if (ierr/=nf90_noerr) print*,'ERROR IN THE SAVING OF ZNIGHT VARIABLE',nf90_strerror(ierr)

             !!! DEFINE BELOW CLASS-SPECIFIC VARIABLES TO SAVE
             !if (part_spec .ge. 5) then
             !   ierr = nf90_put_var(file_id, modeid, all_mode, (/1, j, ncinc/), (/npart, 1, 1/))
             !   if (ierr/=nf90_noerr) print*,'ERROR IN THE SAVING OF MODE VARIABLE',nf90_strerror(ierr)
             !end if
             if (part_spec.gt.0 ) then
                ierr = nf90_put_var(file_id, tpid, all_tp, (/1, j,ncinc/), (/npart,1, 1/))
                if (ierr/=nf90_noerr) print*,'ERROR IN THE SAVING OF TP VARIABLE',nf90_strerror(ierr)

                ierr = nf90_put_var(file_id, massid, all_mass, (/1, j, ncinc/), (/npart, 1, 1/))
                if (ierr/=nf90_noerr) print*,'ERROR IN THE SAVING OF MASS VARIABLE',nf90_strerror(ierr)
             end if
             if (part_spec.eq.3) then
                ierr = nf90_put_var(file_id, lipidid, all_lipid, (/1, j, ncinc/), (/npart, 1, 1/))
                if (ierr/=nf90_noerr) print*,'ERROR IN THE SAVING OF LIPID VARIABLE',nf90_strerror(ierr) 

                ierr = nf90_put_var(file_id, stageid, all_stage, (/1, j, ncinc/), (/npart, 1, 1/))
                if (ierr/=nf90_noerr) print*,'ERROR IN THE SAVING OF STAGE VARIABLE',nf90_strerror(ierr)

#ifdef YUMY
                ierr = nf90_put_var(file_id, foodid, all_fd, (/1, j, ncinc/), (/npart, 1, 1/))
                if (ierr/=nf90_noerr) print*,'ERROR IN THE SAVING OF FOOD VARIABLE',nf90_strerror(ierr)
#endif

                ierr = nf90_put_var(file_id, part_statusid, all_part_status, (/1, j, ncinc/), (/npart, 1, 1/))
                if (ierr/=nf90_noerr) print*,'ERROR IN THE SAVING OF PARTICLE STATUS VARIABLE',nf90_strerror(ierr)

                ierr = nf90_put_var(file_id, abundid, all_abund, (/1, j, ncinc/), (/npart, 1, 1/))
                if (ierr/=nf90_noerr) print*,'ERROR IN THE SAVING OF ABUNDANCE VARIABLE',nf90_strerror(ierr)

                ierr = nf90_put_var(file_id, activityid, all_activity, (/1, j, ncinc/), (/npart, 1, 1/))
                if (ierr/=nf90_noerr) print*,'ERROR IN THE SAVING OF ACTIVITY VARIABLE ',nf90_strerror(ierr)

                ierr = nf90_put_var(file_id, ageid, all_age, (/1, j, ncinc/), (/npart, 1, 1/))
                if (ierr/=nf90_noerr) print*,'ERROR IN THE SAVING OF AGE VARIABLE',nf90_strerror(ierr)

                ierr = nf90_put_var(file_id, reproid, all_repro, (/1, 1, j, ncinc/), (/6, npart, 1, 1/))
                if (ierr/=nf90_noerr) print*,'ERROR IN THE SAVING OF REPRODUCE VARIABLE',nf90_strerror(ierr)

                ierr = nf90_put_var(file_id, cum_eggsid, all_cum_eggs, (/1, j, ncinc/), (/npart, 1, 1/))
                if (ierr/=nf90_noerr) print*,'ERROR IN THE SAVING OF ACCUMULATED EGGS VARIABLE',nf90_strerror(ierr)

                ierr = nf90_put_var(file_id, epid, all_ep, (/1, j, ncinc/), (/npart, 1, 1/))
                if (ierr/=nf90_noerr) print*,'ERROR IN THE SAVING OF DAILY EGGS VARIABLE',nf90_strerror(ierr)

#ifdef YUMY
                ierr = nf90_put_var(file_id, phyto3did, all_ph, (/1, j, ncinc/), (/npart, 1, 1/))
                if (ierr/=nf90_noerr) print*,'ERROR IN THE SAVING OF PHYTOPLANKTON VARIABLE',nf90_strerror(ierr)
#endif
             
                ierr = nf90_put_var(file_id, serialid, all_serial, (/1, j/), (/npart, 1/))
                if (ierr/=nf90_noerr) print*,'ERROR IN THE SAVING OF SERIAL VARIABLE',nf90_strerror(ierr)
             
                ierr = nf90_put_var(file_id, birthdateid, all_birthdate, (/1, j/), (/npart, 1/))
                if (ierr/=nf90_noerr) print*,'ERROR IN THE SAVING OF BIRTHDATE VARIABLE',nf90_strerror(ierr)
             
                ierr = nf90_put_var(file_id, materid, all_mater, (/1, j/), (/npart, 1/))
                if (ierr/=nf90_noerr) print*,'ERROR IN THE SAVING OF MATER VARIABLE',nf90_strerror(ierr)
             end if
             
          end if

          i = mod(i+1, npart)
          if (i==0) i = npart

          call list%next()
       end do

       ! Deallocating all allocatable arrays       allocate(all_xpo(npart))
       deallocate(all_xpo)
       deallocate(all_lon)
       deallocate(all_ypo)
       deallocate(all_lat)
       deallocate(all_zpo)
       deallocate(all_bdep)

       deallocate(all_species)

       deallocate(all_znight)

       !!! DEFINE BELOW CLASS-SPECIFIC VARIABLES TO SAVE
       if (part_spec .ge. 5) then
          deallocate(all_mode)
          deallocate(all_mode_prop)
       end if

       if (part_spec.gt.0 ) then
          deallocate(all_sex)
          deallocate(all_tp)
          deallocate(all_mass)
       end if
       
       if (part_spec.ge.3 .and. part_spec.lt.5) then
          deallocate(all_lipid)
          deallocate(all_stage)
#ifdef YUMY
          deallocate(all_fd)
#endif
          deallocate(all_part_status)
          deallocate(all_abund)
          deallocate(all_activity)
          deallocate(all_age)
          deallocate(all_repro)
          deallocate(all_cum_eggs)
          deallocate(all_ep)
#ifdef YUMY
          deallocate(all_ph)
#endif
          deallocate(all_serial)
          deallocate(all_birthdate)
          deallocate(all_mater)
       end if
       
    case (1)

       ncinc = ncinc + 1

       ! Put timestep

#ifdef BACK
       nctime = nctime - ncfreq
#else
       nctime = nctime + ncfreq
#endif

       ierr = nf90_put_var(file_id, timeid, nctime, (/ncinc/))
       if (ierr/=nf90_noerr) print*, 'Error writing T var : ', trim(nf90_strerror(ierr))
       ! Allocate size of list variables
       allocate(all_xpo(npart))
       allocate(all_lon(npart))
       allocate(all_ypo(npart))
       allocate(all_lat(npart))
       allocate(all_zpo(npart))
       allocate(all_bdep(npart))
       allocate(all_mode(npart))
       allocate(all_mode_prop(npart))
       allocate(all_species(npart))
       allocate(all_znight(npart))


       !!! DEFINE BELOW CLASS-SPECIFIC VARIABLES TO SAVE
       if (part_spec.gt.0 ) then
          allocate(all_tp(npart))
          allocate(all_sex(npart))
          allocate(all_mass(npart))
       end if
       
       if (part_spec.ge.3 .and. part_spec.lt.5) then
          allocate(all_lipid(npart))
          allocate(all_stage(npart))
#ifdef YUMY
          allocate(all_fd(npart))
#endif
          allocate(all_part_status(npart))
          allocate(all_abund(npart))
          allocate(all_activity(npart))
          allocate(all_age(npart))
          allocate(all_repro(6,npart))
          allocate(all_cum_eggs(npart))
          allocate(all_ep(npart))
#ifdef YUMY
          allocate(all_ph(npart))
#endif
          allocate(all_serial(npart))
          allocate(all_birthdate(npart))
          allocate(all_mater(npart))
       end if

       call list%reset()

       i = 1
       j = 0

       do while(list%moreValues())
          curr => list%currentValue()

          select type(curr)

          type is (Zoo)
             all_species(i) = curr%get_species()
             all_xpo(i)     = curr%get_xpo()
             all_ypo(i)     = curr%get_ypo()
             all_zpo(i)     = curr%get_zpo()
             all_bdep(i)    = curr%get_bdep()
             all_znight(i)  = curr%get_znight()

          type is (Krill)
             all_species(i) = curr%get_species()
             all_xpo(i)     = curr%get_xpo()
             all_ypo(i)     = curr%get_ypo()
             all_zpo(i)     = curr%get_zpo()
             all_bdep(i)    = curr%get_bdep()
             all_znight(i)  = curr%get_znight()

             all_tp(i)      = curr%get_tp()
             all_sex(i)     = curr%get_sex()
             all_mass(i)    = curr%get_mass()

          type is (Calanus)
             all_species(i) = curr%get_species()
             all_xpo(i)     = curr%get_xpo()
             all_ypo(i)     = curr%get_ypo()
             all_zpo(i)     = curr%get_zpo()
             all_bdep(i)    = curr%get_bdep()
             all_znight(i)  = curr%get_znight()

             all_tp(i)           = curr%get_tp()
             all_sex(i)          = curr%get_sex()
             all_mass(i)         = curr%get_mass()
             all_lipid(i)        = curr%get_lipid()
             all_stage(i)        = curr%get_stage()
#ifdef YUMY
             all_fd(i)           = curr%get_fd()
#endif
             all_part_status(i)  = curr%get_part_status()
             all_abund(i)        = curr%get_abund()
             all_activity(i)     = curr%get_activity()
             all_age(i)          = curr%get_age()
             all_repro(1:6,i)    = curr%get_reproduce_info()
             all_cum_eggs(i)     = curr%get_cum_eggs()
             all_ep(i)           = curr%get_EP()
#ifdef YUMY
             all_ph(i)           = curr%get_ph()
#endif
             all_serial(i)       = curr%get_serial()
             all_birthdate(i)    = curr%get_birthdate()
             all_mater(i)        = curr%get_mater()

          type is (Late_hyperboreus)
             all_species(i)   = curr%get_species()
             all_xpo(i)       = curr%get_xpo()
             all_ypo(i)       = curr%get_ypo()
             all_zpo(i)       = curr%get_zpo()
             all_bdep(i)      = curr%get_bdep()
             all_mode(i)      = curr%get_mode()
             all_mode_prop(i) = curr%get_mode_prop()
            all_znight(i)     = curr%get_znight()

          type is (Young_hyperboreus)
             all_species(i)   = curr%get_species()
             all_xpo(i)       = curr%get_xpo()
             all_ypo(i)       = curr%get_ypo()
             all_zpo(i)       = curr%get_zpo()
             all_bdep(i)      = curr%get_bdep()
             all_mode(i)      = curr%get_mode()
             all_mode_prop(i) = curr%get_mode_prop()
             all_znight(i)     = curr%get_znight()

          type is (Late_finmarchicus)
             all_species(i)   = curr%get_species()
             all_xpo(i)       = curr%get_xpo()
             all_ypo(i)       = curr%get_ypo()
             all_zpo(i)       = curr%get_zpo()
             all_bdep(i)      = curr%get_bdep()
             all_mode(i)      = curr%get_mode()
             all_mode_prop(i) = curr%get_mode_prop()
             all_znight(i)     = curr%get_znight()

          type is (Young_finmarchicus)
             all_species(i)   = curr%get_species()
             all_xpo(i)       = curr%get_xpo()
             all_ypo(i)       = curr%get_ypo()
             all_zpo(i)       = curr%get_zpo()
             all_bdep(i)      = curr%get_bdep()
             all_mode(i)      = curr%get_mode()
             all_mode_prop(i) = curr%get_mode_prop()
             all_znight(i)     = curr%get_znight()


          end select

          ! Compute lat/lon coordinates
          call xy_to_ll_NEMO(all_xpo(i), all_ypo(i), all_lon(i), all_lat(i))

          ! Write values to variables
          if ( mod(i, npart)==0 ) then

             j = j + 1

             ierr = nf90_put_var(file_id, xpoid, all_xpo, (/1, j, ncinc/), (/npart, 1, 1/))
             if (ierr/=nf90_noerr) print*,'ERROR IN THE SAVING OF XPO VARIABLE ',nf90_strerror(ierr)

             ierr = nf90_put_var(file_id, lonid, all_lon, (/1, j, ncinc/), (/npart, 1, 1/))
             if (ierr/=nf90_noerr) print*,'ERROR IN THE SAVING OF LON VARIABLE',nf90_strerror(ierr)

             ierr = nf90_put_var(file_id, ypoid, all_ypo, (/1, j, ncinc/), (/npart, 1, 1/))
             if (ierr/=nf90_noerr) print*,'ERROR IN THE SAVING OF YPO VARIABLE',nf90_strerror(ierr)

             ierr = nf90_put_var(file_id, latid, all_lat, (/1, j, ncinc/), (/npart, 1, 1/))
             if (ierr/=nf90_noerr) print*,'ERROR IN THE SAVING OF LAT VARIABLE',nf90_strerror(ierr)

             ierr = nf90_put_var(file_id, zpoid, all_zpo, (/1, j, ncinc/), (/npart, 1, 1/))
             if (ierr/=nf90_noerr) print*,'ERROR IN THE SAVING OF ZPO VARIABLE',nf90_strerror(ierr)

             ierr = nf90_put_var(file_id, bdepid, all_bdep, (/1, j, ncinc/), (/npart, 1, 1/))
             if (ierr/=nf90_noerr) print*,'ERROR IN THE SAVING OF ZPO VARIABLE',nf90_strerror(ierr)

             ierr = nf90_put_var(file_id, znightid, all_znight, (/1, j, ncinc/), (/npart, 1, 1/))
             if (ierr/=nf90_noerr) print*,'ERROR IN THE SAVING OF ZNIGHT VARIABLE',nf90_strerror(ierr)

             !!! DEFINE BELOW CLASS-SPECIFIC VARIABLES TO SAVE
             if (part_spec.gt.0 ) then
                ierr = nf90_put_var(file_id, tpid, all_tp, (/1, j, ncinc/),(/npart, 1, 1/))
                if (ierr/=nf90_noerr) print*,'ERROR IN THE SAVING OF TP VARIABLE',nf90_strerror(ierr)

                ierr = nf90_put_var(file_id, massid, all_mass, (/1, j, ncinc/), (/npart, 1, 1/))
                if (ierr/=nf90_noerr) print*,'ERROR IN THE SAVING OF MASS VARIABLE',nf90_strerror(ierr)


             end if

             !if (part_spec .ge. 5) then
             !   ierr = nf90_put_var(file_id, modeid, mode_propid, all_mode, all_mode_prop, (/1, j, ncinc/), (/npart, 1, 1/))
             !   if (ierr/=nf90_noerr) print*,'ERROR IN THE SAVING OF MODE VARIABLE',nf90_strerror(ierr)
             !end if

             if (part_spec.ge.3 .and. part_spec.lt.5) then
                ierr = nf90_put_var(file_id, lipidid, all_lipid, (/1, j, ncinc/), (/npart, 1, 1/))
                if (ierr/=nf90_noerr) print*,'ERROR IN THE SAVING OF LIPID VARIABLE',nf90_strerror(ierr)

                ierr = nf90_put_var(file_id, stageid, all_stage, (/1, j, ncinc/), (/npart, 1, 1/))
                if (ierr/=nf90_noerr) print*,'ERROR IN THE SAVING OF STAGE VARIABLE',nf90_strerror(ierr)

#ifdef YUMY
                ierr = nf90_put_var(file_id, foodid, all_fd, (/1, j, ncinc/), (/npart, 1, 1/))
                if (ierr/=nf90_noerr) print*,'ERROR IN THE SAVING OF FOOD VARIABLE',nf90_strerror(ierr)
#endif

                ierr = nf90_put_var(file_id, part_statusid, all_part_status, (/1, j, ncinc/), (/npart, 1, 1/))
                if (ierr/=nf90_noerr) print*,'ERROR IN THE SAVING OF PARTICLE STATUS VARIABLE',nf90_strerror(ierr)

                ierr = nf90_put_var(file_id, abundid, all_abund, (/1, j, ncinc/), (/npart, 1, 1/))
                if (ierr/=nf90_noerr) print*,'ERROR IN THE SAVING OF ABUNDANCE VARIABLE',nf90_strerror(ierr)

                ierr = nf90_put_var(file_id, activityid, all_activity, (/1, j, ncinc/), (/npart, 1, 1/))
                if (ierr/=nf90_noerr) print*,'ERROR IN THE SAVING OF ACTIVITY VARIABLE',nf90_strerror(ierr)

                ierr = nf90_put_var(file_id, ageid, all_age, (/1, j, ncinc/), (/npart, 1, 1/))
                if (ierr/=nf90_noerr) print*,'ERROR IN THE SAVING OF AGE VARIABLE',nf90_strerror(ierr)

                ierr = nf90_put_var(file_id, sexid, all_sex, (/1, j, ncinc/), (/npart, 1, 1/))
                if (ierr/=nf90_noerr) print*,'ERROR IN THE SAVING OF SEX VARIABLE',nf90_strerror(ierr)

                ierr = nf90_put_var(file_id, reproid, all_repro, (/1, 1, j, ncinc/), (/6, npart, 1, 1/))
                if (ierr/=nf90_noerr) print*,'ERROR IN THE SAVING OF REPRODUCE VARIABLE',nf90_strerror(ierr)

                ierr = nf90_put_var(file_id, cum_eggsid, all_cum_eggs, (/1, j, ncinc/), (/npart, 1, 1/))
                if (ierr/=nf90_noerr) print*,'ERROR IN THE SAVING OF ACCUMULATED EGGS VARIABLE',nf90_strerror(ierr)

                ierr = nf90_put_var(file_id, epid, all_ep, (/1, j, ncinc/), (/npart, 1, 1/))
                if (ierr/=nf90_noerr) print*,'ERROR IN THE SAVING OF DAILY EGGS VARIABLE',nf90_strerror(ierr)

#ifdef YUMY
                ierr = nf90_put_var(file_id, phyto3did, all_ph, (/1, j, ncinc/), (/npart, 1, 1/))
                if (ierr/=nf90_noerr) print*,'ERROR IN THE SAVING OF PHYTOPLANKTON VARIABLE',nf90_strerror(ierr)
#endif
                
                ierr = nf90_put_var(file_id, serialid, all_serial, (/1, j/), (/npart, 1/))
                if (ierr/=nf90_noerr) print*,'ERROR IN THE SAVING OF SERIAL VARIABLE',nf90_strerror(ierr)

                ierr = nf90_put_var(file_id, birthdateid, all_birthdate, (/1, j/), (/npart, 1/))
                if (ierr/=nf90_noerr) print*,'ERROR IN THE SAVING OF BIRTHDATE VARIABLE',nf90_strerror(ierr)

                ierr = nf90_put_var(file_id, materid, all_mater, (/1, j/), (/npart, 1/))
                if (ierr/=nf90_noerr) print*,'ERROR IN THE SAVING OF MATER VARIABLE',nf90_strerror(ierr)

             end if
             
          end if

          i = mod(i+1, npart)
          if (i==0) i = npart

          call list%next()
       end do

       ! Deallocating all allocatable arrays       
       deallocate(all_xpo)
       deallocate(all_lon)
       deallocate(all_ypo)
       deallocate(all_lat)
       deallocate(all_zpo)
       deallocate(all_bdep)

       deallocate(all_species)

       deallocate(all_znight)

       !!! DEFINE BELOW CLASS-SPECIFIC VARIABLES TO SAVE
       if (part_spec.gt.0 ) then
          deallocate(all_tp)
          deallocate(all_mass)
          deallocate(all_sex)
       end if
       
       if (part_spec .ge. 5) then
          deallocate(all_mode)
          deallocate(all_mode_prop)
       end if

       if (part_spec.ge.3 .and. part_spec.lt.5) then
          deallocate(all_lipid)
          deallocate(all_stage)
#ifdef YUMY
          deallocate(all_fd)
#endif
          deallocate(all_part_status)
          deallocate(all_abund)
          deallocate(all_activity)
          deallocate(all_age)
          deallocate(all_repro)
          deallocate(all_cum_eggs)
          deallocate(all_ep)
#ifdef YUMY
          deallocate(all_ph)
#endif
          deallocate(all_serial)
          deallocate(all_birthdate)
          deallocate(all_mater)
       end if
       
    case(2) ! Close

       ierr = nf90_close(file_id)

       write(*,'(//a,/)') ' Close the NetCDF output file'

    end select

    return

  end subroutine netcdfPart

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!72

end module partFun_mod
