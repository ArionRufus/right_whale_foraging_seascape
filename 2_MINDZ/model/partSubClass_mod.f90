module partSubClass_mod
  
  use rng_par_zig ! random generator module
  use partClass_mod

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! BELOW IS THE USER-DEFINED SUBCLASS(ES) SECTION        !!
  !! MODIFY AS NEEDED                                      !!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !--------------------------------------------------------------------72
  ! Particle physiology

  integer, parameter :: nb_pars=14

  ! Krill
  integer            :: sex
  
  real               :: mean_mass, min_mass, max_mass, &
                        a_molt, b_molt,                &
                        k0, h0, A,                     &
                        r0, rb,                        &
                        er, t_lim, ep1, ep2
   
  real               :: ratio_ing

  logical            :: mature, clone

  logical            :: FORC_PMZA, INTRA_SPE

#ifdef YUMY
  character(len=150) :: infile_food
#endif

  ! Calanus
  real, dimension(0:12) :: a_devlp
  real, dimension(1:14) :: paramosome
  real, dimension(1:9)  :: crit_moult_max
  real                  :: t_devlp, init_stage, init_abund, &
                           init_age, age
  integer               :: init_activity

  !--------------------------------------------------------------------72
  ! I/O

  namelist /krill_nml/ INTRA_SPE,   FORC_PMZA,           &
#ifdef YUMY
                       infile_food,                      &
#endif                       
                       clone,       mature,    sex,      &
                       mean_mass,   min_mass,  max_mass, &
                       a_molt,      b_molt,              &
                       k0,          h0,        A,        &
                       r0,          rb,                  &
                       er,          t_lim 

  namelist /calanus_nml/ INTRA_SPE,   FORC_PMZA,           &
#ifdef YUMY  
                         infile_food,                      &
#endif                         
                         clone,       mature,    sex,      &
                         mean_mass,   min_mass,  max_mass, &
                         init_stage,  a_devlp,   t_devlp,  &
                         k0,          h0,        A,        &
                         r0,          rb,                  &
                         er,          t_lim,   init_abund, &
                         init_activity, init_age,          &
                         paramosome, crit_moult_max

  !--------------------------------------------------------------------72
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!! Particle I (KRILL) !!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  public :: Krill

  type, extends(Zoo) :: Krill

     ! Environment variables
     real    :: tp

     ! Sex
     integer :: sex 

     ! Size
     real    :: mass

     ! Arrhenius relationship
     real    :: er 
     real    :: t_lim

     ! Respiration
     real    :: r0       
     real    :: rb

     ! Moulting
     real    :: dev_frac
     real    :: a_molt
     real    :: b_molt

     ! Optimal depth
     !real    :: znight  ! KB new
     !real    :: r1      ! KB new

   contains
     
     ! Getters
     procedure, public :: get_tp       => get_tpK

     procedure, public :: get_sex      => get_sexK
     procedure, public :: get_mass     => get_massK
     procedure, public :: get_dev_frac

     ! Setters
     procedure, public :: set_tp       => set_tpK
     
     procedure, public :: set_sex      => set_sexK
     procedure, public :: set_mass     => set_massK
     procedure, public :: set_dev_frac

     ! Methods
     procedure, public :: zday_fun     => zday_funK
     procedure, public :: znight_fun   => znight_funK
     
     procedure         :: devlp        => devlpK    

     procedure         :: arrhenius    => arrheniusK
     procedure         :: adj_arr      => adj_arrK
     
     procedure         :: resp         => respK

     procedure, public :: molt         => moltK
     procedure, public :: grow         => growK

  end type Krill
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!! Particle II (CALANUS) !!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  public :: Calanus

  type, extends(Zoo) :: Calanus  

     ! Environment variables
     real    :: tp
     ! Environmental Food
     real    :: fd 
     ! Environmental phytoplankton
     real    :: ph
     ! Sex
     integer :: sex 
     ! Size
     real    :: mass

     ! Arrhenius relationship
     real    :: er 
     real    :: t_lim
     ! epsilon for entry and exit diapause
     real    :: ep1
     real    :: ep2
     ! Respiration
     real    :: r0       
     real    :: rb

     ! Moulting
     real    :: stage
     real    :: a_molt
     real    :: b_molt
 
     !Lipid
     real    :: lipid
 
     !egg production, starvation, diapause
     real    :: lipid_floor 
     real    :: max_diap_lipid, lipid_exit_diap, mass_limit, massp
     real    :: start_ep, duration_ep, ep_daily,counter_eggs ! counter_egg is combining with other real numbers, so it has to be real
     integer :: EP, cum_eggs, egg_clutch, EP_CLUTCH
     real, dimension(1:6)   ::reproduce_info
 
     !Abundance
     real    :: abund
 
     !Activity
     integer :: activity, pre_activity 
     
     !Age
     real    :: age
     
     ! for diapause
     logical :: flag_enter_diap
     real    :: counter_c6_active
     
     ! Accounting (particle ID=serial, parent ID=mater)
     integer    :: serial, mater, birthdate  ! KB new
     
   contains
     
     ! Getters
     procedure, public :: get_tp       => get_tpC
     procedure, public :: get_sex      => get_sexC
     procedure, public :: get_mass     => get_massC
     procedure, public :: get_a_molt   => get_a_moltC
     procedure, public :: get_b_molt   => get_b_moltC
     procedure, public :: get_lipid
     procedure, public :: get_stage
     procedure, public :: get_fd
     procedure, public :: get_ph
     procedure, public :: get_abund
     procedure, public :: get_activity
     procedure, public :: get_pre_activity
     procedure, public :: get_lipid_floor
     procedure, public :: get_ep_daily
     procedure, public :: get_EP
     procedure, public :: get_egg_clutch
     procedure, public :: get_EP_CLUTCH
     procedure, public :: get_counter_eggs
     procedure, public :: get_cum_eggs
     procedure, public :: get_start_ep
     procedure, public :: get_reproduce_info
     procedure, public :: get_duration_ep
     procedure, public :: get_max_diap_lipid
     procedure, public :: get_lipid_exit_diap
     procedure, public :: get_moult_mass_limit
     procedure, public :: get_massp
     procedure, public :: get_age
     procedure, public :: get_serial
     procedure, public :: get_mater
     procedure, public :: get_birthdate
     ! Setters
     procedure, public :: set_tp       => set_tpC
     procedure, public :: set_sex      => set_sexC
     procedure, public :: set_mass     => set_massC
     procedure, public :: set_a_molt   => set_a_moltC
     procedure, public :: set_b_molt   => set_b_moltC
     procedure, public :: set_lipid
     procedure, public :: set_stage
     procedure, public :: set_fd
     procedure, public :: set_ph
     procedure, public :: set_abund
     procedure, public :: set_activity
     procedure, public :: set_pre_activity
     procedure, public :: set_lipid_floor
     procedure, public :: set_ep_daily
     procedure, public :: set_EP
     procedure, public :: set_egg_clutch
     procedure, public :: set_EP_CLUTCH
     procedure, public :: set_reproduce_info
     procedure, public :: set_counter_eggs
     procedure, public :: set_cum_eggs
     procedure, public :: set_start_ep
     procedure, public :: set_duration_ep
     procedure, public :: set_max_diap_lipid
     procedure, public :: set_lipid_exit_diap
     procedure, public :: set_moult_mass_limit
     procedure, public :: set_massp
     procedure, public :: set_age
     procedure, public :: set_serial
     procedure, public :: set_mater
     procedure, public :: set_birthdate
     ! Methods
     procedure, public :: zday_fun     => zday_funC
     procedure, public :: znight_fun   => znight_funC

     procedure         :: devlp        => devlpC
     procedure, public :: grow_older   => grow_olderC

     procedure         :: arrhenius    => arrheniusC
     procedure         :: adj_arr      => adj_arrC
     
     procedure         :: resp         => respC
     procedure, public :: molt         => moltC
     procedure, public :: grow         => growC
     procedure         :: ingest  
     procedure, public :: max_grow_depth     
     procedure, public :: egg_pro
     procedure, public :: starve
     procedure, public :: mort_maps2012
     procedure, public :: predation
 
     procedure, public :: diapause_entry
     procedure, public :: diapause_exit
  
  end type Calanus



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Particle V (Late_hyperboreus) !!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

public :: Late_hyperboreus

type, extends(Zoo) :: Late_hyperboreus

   ! Environment variables
   real    :: tp
   ! Size
   real    :: mass
   
 contains
   
   ! Methods
   procedure, public :: zday_fun       => zday_funLH
   procedure, public :: znight_fun     => znight_funLH
   procedure, public :: sigma_fun      => sigma_funLH
   procedure, public :: mode_prop_fun  => mode_prop_funLH

end type Late_hyperboreus




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Particle VI (Young_hyperboreus) !!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

public :: Young_hyperboreus

type, extends(Zoo) :: Young_hyperboreus

   ! Environment variables
   real    :: tp
   ! Size
   real    :: mass
   
 contains
   
   ! Methods
   procedure, public :: zday_fun       => zday_funYH
   procedure, public :: znight_fun     => znight_funYH
   procedure, public :: sigma_fun      => sigma_funYH
   procedure, public :: mode_prop_fun  => mode_prop_funYH

end type Young_hyperboreus




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Particle VII (Late_finmarchicus) !!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

public :: Late_finmarchicus

type, extends(Zoo) :: Late_finmarchicus

   ! Environment variables
   real    :: tp
   ! Size
   real    :: mass
   
 contains
   
   ! Methods
   procedure, public :: zday_fun       => zday_funLF
   procedure, public :: znight_fun     => znight_funLF
   procedure, public :: sigma_fun      => sigma_funLF
   procedure, public :: mode_prop_fun  => mode_prop_funLF

end type Late_finmarchicus





!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Particle VII (Young_finmarchicus) !!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

public :: Young_finmarchicus

type, extends(Zoo) :: Young_finmarchicus

   ! Environment variables
   real    :: tp
   ! Size
   real    :: mass
   
 contains
   
   ! Methods
   procedure, public :: zday_fun       => zday_funYF
   procedure, public :: znight_fun     => znight_funYF
   procedure, public :: sigma_fun      => sigma_funYF
   procedure, public :: mode_prop_fun  => mode_prop_funYF

end type Young_finmarchicus





contains

  
  !! Implement your methods here !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!72

    
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!72
  !! KRILL class procedures
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!72

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!72
  !! Getters

  ! Environmental variables
  real function get_tpK(this)
    class(Krill) :: this
    get_tpK = this%tp
  end function get_tpK

  ! Sex
  integer function get_sexK(this)
    class(Krill) :: this
    get_sexK = this%sex
  end function get_sexK

  ! Size
  real function get_massK(this)
    class(Krill) :: this
    get_massK = this%mass
  end function get_massK

  ! Arrhenius relationship
  real function get_erK(this)
    class(Krill) :: this
    get_erK = this%er
  end function get_erK

  real function get_t_limK(this)
    class(Krill) :: this
    get_t_limK = this%t_lim
  end function get_t_limK

  ! Respiration process
  real function get_r0K(this)
    class(Krill) :: this
    get_r0K = this%r0
  end function get_r0K

  real function get_rbK(this)
    class(Krill) :: this
    get_rbK = this%rb
  end function get_rbK
  
  ! Moulting process 
  real function get_dev_frac(this)
    class(Krill) :: this
    get_dev_frac = this%dev_frac
  end function get_dev_frac

  real function get_a_moltK(this)
    class(Krill) :: this
    get_a_moltK = this%a_molt
  end function get_a_moltK

  real function get_b_moltK(this)
    class(Krill) :: this
    get_b_moltK = this%b_molt
  end function get_b_moltK

  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!72
  !! Setters

  ! Environmental variables
  subroutine set_tpK(this, tp)
    class(Krill) :: this
    real         :: tp
    this%tp = tp
  end subroutine set_tpK

  ! Sex
  subroutine set_sexK(this, sex)
    class(Krill) :: this
    real         :: sex
    this%sex = sex
  end subroutine set_sexK

  ! Size
  subroutine set_massK(this, mass)
    class(Krill) :: this
    real         :: mass
    this%mass = mass
  end subroutine set_massK

  ! Arrhenius relationship
  subroutine set_erK(this, er)
    class(Krill) :: this
    real         :: er
    this%er = er
  end subroutine set_erK

  subroutine set_t_limK(this, t_lim)
    class(Krill) :: this
    real         :: t_lim
    this%t_lim = t_lim
  end subroutine set_t_limK

  ! Respiration process
  subroutine set_r0K(this, r0)
    class(Krill) :: this
    real         :: r0
    this%r0 = r0
  end subroutine set_r0K

  subroutine set_rbK(this, rb)
    class(Krill) :: this
    real         :: rb
    this%rb = rb
  end subroutine set_rbK
  
  ! Moulting process 
  subroutine set_dev_frac(this, dev_frac)
    class(Krill) :: this
    real         :: dev_frac
    this%dev_frac = dev_frac
  end subroutine set_dev_frac

  subroutine set_a_moltK(this, a_molt)
    class(Krill) :: this
    real         :: a_molt
    this%a_molt = a_molt
  end subroutine set_a_moltK

  subroutine set_b_moltK(this, b_molt)
    class(Krill) :: this
    real         :: b_molt
    this%b_molt = b_molt
  end subroutine set_b_moltK

  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!72
  !! KRILL class procedures

  ! Daytime preferred depth
  real function zday_funK(this, surf) result(zval)
    class(Krill)   :: this
    real, optional :: surf

    if ( present(surf) ) then
       zval = 11.43 * max(24.,surf) - 200.87 ! T. raschii
    else
       zval = 10.
    endif
     
  end function zday_funK
  
  ! Nightime preferred depth
  real function znight_funK(this) result(zval)
    class(Krill)   :: this

    call rng_uni(r1)
    k = max( 1, ceiling( r1 * 10. ) )
    zval = mid_depth(k)
     
  end function znight_funK
  
  ! Arrhenius function of temperature = Universal Temperature Dependance
  ! (Gillooly et al. 2001)
  ! Returns a temperature dependance coefficient (unitless)
  real function arrheniusK(this, e_act)
    class(Krill)    :: this
    real            :: e_act
    real, parameter :: k  = 8.62e-5 ! The Boltzmann's constant in eV.K-1
    real, parameter :: T0 = 273.15  ! Reference temperature (the frozing point of water) in K

    arrheniusK = exp( ( e_act * this%tp ) / ( k * ( this%tp + T0 ) * T0 ) )

  end function arrheniusK

  ! Arrhenius adjustement function ensures an exponential decrease
  ! AFTER some threshold
  real function adj_arrK(this, alpha_temp)
    class(Krill)    :: this
    real            :: alpha_temp, dtemp

    dtemp = this%t_lim - this%tp

    if( dtemp > 0 ) then 
       adj_arrK = 1
    else 
       adj_arrK = exp( dtemp * alpha_temp )
    endif

  end function adj_arrK

  ! Development time (intermoult period) in days
  real function devlpK(this)
    class(Krill)    :: this

    devlpK = 86400. * ( this%a_molt + this%b_molt * this%tp )

    if ( devlpK < 3. * 86400. ) devlpK = 3. * 86400.

  end function devlpK

  ! Development freq of the individual. 
  ! when the freq > 0.4, computes the molt_factor which fixes
  !                      the next individual's length 
  ! When the freq > 1,   updates the individual's length with the
  !                      previously computed molt_factor.
  subroutine moltK(this)
    class(Krill)    :: this

    if ( this%cell_out .eqv. .true. ) then
       this%dev_frac = 0.0
    else 

       this%dev_frac = this%dev_frac + ( dt / this%devlp() )

       ! NOTE: moulting is irreversibly engaged after 40% of the IMP.
       ! The new individuals length is completely decided at this time
       ! according to the allometric relationship
       !       if ( this%dev_frac > 0.4 .and. this%molt_length < EPS ) then
       !          this%molt_length = ( this%mass / this%aw )** ( 1 / this%bw )
       !       endif

       ! At the end of the development phase, the individual molts
       if ( this%dev_frac > 1. ) then
          this%mass     = this%mass
          this%dev_frac = this%dev_frac - 1.
       endif

    endif

  end subroutine moltK

  ! Respiration = basal metabolism + activity (swimming) metabolism
  ! but NOT specific dynamic action as it is included in assimilation
  ! coeffcient
  ! Data from Angelique
  real function respK(this)
    class(Krill)    :: this
    real, parameter :: alpha = 0.2

    respK =  this%r0                 &
           * this%mass ** this%rb    &
           * this%arrhenius(this%er) &
           * this%adj_arr(alpha)     &
           * dt

  end function respK

  ! Updates the new mass of the individual
  subroutine growK(this)
    class(Krill)    :: this

    if ( this%cell_out .eqv. .true. ) then 
       this%mass = this%mass
    else 
       this%mass = this%mass - this%resp() !+ this%A * this%ingest()
    endif

  end subroutine growK

  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!72
  !! CALANUS class procedures
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!72


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!72
  !! Getters

  ! Environmental variables
  real function get_tpC(this)
    class(Calanus) :: this
    get_tpC = this%tp
  end function get_tpC
  
  real function get_fd(this)
    class(Calanus) :: this
    get_fd = this%fd
  end function get_fd
  
  real function get_ph(this)
    class(Calanus) :: this
    get_ph = this%ph
  end function get_ph
  
  ! Sex
  integer function get_sexC(this)
    class(Calanus) :: this
    get_sexC = this%sex
  end function get_sexC

  ! Size
  real function get_massC(this)
    class(Calanus) :: this
    get_massC = this%mass
  end function get_massC

  ! Arrhenius relationship
  real function get_erC(this)
    class(Calanus) :: this
    get_erC = this%er
  end function get_erC

  real function get_t_limC(this)
    class(Calanus) :: this
    get_t_limC = this%t_lim
  end function get_t_limC

  ! Respiration process
  real function get_r0C(this)
    class(Calanus) :: this
    get_r0C= this%r0
  end function get_r0C

  real function get_rbC(this)
    class(Calanus) :: this
    get_rbC= this%rb
  end function get_rbC
  
  ! Lipid
  real function get_lipid(this)
    class(Calanus) :: this
    get_lipid = this%lipid
  end function get_lipid  
  
  ! Moulting process 
  real function get_stage(this)
    class(Calanus) :: this
    get_stage = this%stage
  end function get_stage

  real function get_a_moltC(this)
    class(Calanus) :: this
    get_a_moltC = this%a_molt
  end function get_a_moltC

  real function get_b_moltC(this)
    class(Calanus) :: this
    get_b_moltC = this%b_molt
  end function get_b_moltC
  ! Abundance
  real function get_abund(this)
    class(Calanus) :: this
    get_abund = this%abund
  end function get_abund
  
  real function get_activity(this)
    class(Calanus) :: this
    get_activity = this%activity
  end function get_activity
  
  real function get_pre_activity(this)
    class(Calanus) :: this
    get_pre_activity = this%pre_activity
  end function get_pre_activity
  
  real function get_lipid_floor(this)
    class(Calanus) :: this
    get_lipid_floor = this%lipid_floor
  end function get_lipid_floor 
  
  real function get_ep_daily(this)
    class(Calanus) :: this
    get_ep_daily = this%ep_daily
  end function get_ep_daily
  
  real function get_egg_clutch(this)
    class(Calanus) :: this
    get_egg_clutch = this%egg_clutch
  end function get_egg_clutch
  
  real function get_counter_eggs(this)
    class(Calanus) :: this
    get_counter_eggs = this%counter_eggs
  end function get_counter_eggs
  
  integer function get_cum_eggs(this)
    class(Calanus) :: this
    get_cum_eggs = this%cum_eggs
  end function get_cum_eggs
  
  integer function get_EP(this)
    class(Calanus) :: this
    get_EP = this%EP
  end function get_EP
  
  integer function get_EP_CLUTCH(this)
    class(Calanus) :: this
    get_EP_CLUTCH = this%EP_CLUTCH
  end function get_EP_CLUTCH
  
  real function get_start_ep(this)
    class(Calanus) :: this
    get_start_ep = this%start_ep
  end function get_start_ep
  
  function get_reproduce_info(this)
    class(Calanus) :: this
	real,dimension(1:6) :: get_reproduce_info
    get_reproduce_info = this%reproduce_info
  end function get_reproduce_info
  
  real function get_duration_ep(this)
    class(Calanus) :: this
    get_duration_ep = this%duration_ep
  end function get_duration_ep
  
  real function get_max_diap_lipid(this)
    class(Calanus) :: this
    get_max_diap_lipid = this%max_diap_lipid
  end function get_max_diap_lipid 
  
  real function get_lipid_exit_diap(this)
    class(Calanus) :: this
    get_lipid_exit_diap = this%lipid_exit_diap
  end function get_lipid_exit_diap
   
  real function get_moult_mass_limit(this)
    class(Calanus) :: this
    get_moult_mass_limit = this%mass_limit
  end function get_moult_mass_limit
  
  real function get_massp(this)
    class(Calanus) :: this
    get_massp = this%massp
  end function get_massp
  
  real function get_age(this)
    class(Calanus) :: this
    get_age = this%age
  end function get_age
  
  integer function get_serial(this)
    class(Calanus) :: this
    get_serial = this%serial
  end function get_serial
  
  integer function get_mater(this)
    class(Calanus) :: this
    get_mater = this%mater
  end function get_mater
   
  integer function get_birthdate(this)
    class(Calanus) :: this
    get_birthdate = this%birthdate
  end function get_birthdate

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!72
  !! Setters

  ! Environmental variables
  subroutine set_tpC(this, tp)
    class(Calanus) :: this
    real         :: tp
    this%tp = tp
  end subroutine set_tpC
  
  ! Food
  subroutine set_fd(this, fd)
    class(Calanus) :: this
    real         :: fd
    this%fd = fd
  end subroutine set_fd
  
  !phytoplankton
  subroutine set_ph(this, ph)
    class(Calanus) :: this
    real         :: ph
    this%ph = ph
  end subroutine set_ph

  ! Sex
  subroutine set_sexC(this, sex)
    class(Calanus) :: this
    integer         :: sex
    this%sex = sex
  end subroutine set_sexC

  ! Size
  subroutine set_massC(this, mass)
    class(Calanus) :: this
    real         :: mass
    this%mass = mass
  end subroutine set_massC

  ! Arrhenius relationship
  subroutine set_erC(this, er)
    class(Calanus) :: this
    real         :: er
    this%er = er
  end subroutine set_erC

  subroutine set_t_limC(this, t_lim)
    class(Calanus) :: this
    real         :: t_lim
    this%t_lim = t_lim
  end subroutine set_t_limC

  ! Respiration process
  subroutine set_r0C(this, r0)
    class(Calanus) :: this
    real         :: r0
    this%r0 = r0
  end subroutine set_r0C

  subroutine set_rbC(this, rb)
    class(Calanus) :: this
    real         :: rb
    this%rb = rb
  end subroutine set_rbC
  
  ! Lipid
  subroutine set_lipid(this, lipid)
    class(Calanus) :: this
    real         :: lipid
    this%lipid = lipid
  end subroutine set_lipid
  
  ! Abundance
  
  subroutine set_abund(this, abund)
    class(Calanus) :: this
    real         :: abund
    this%abund = abund
  end subroutine set_abund
  
  !Activity
  
  subroutine set_activity(this, activity)
    class(Calanus) :: this
    integer        :: activity
    this%activity = activity
  end subroutine set_activity
  
  subroutine set_pre_activity(this, pre_activity)
    class(Calanus) :: this
    integer        :: pre_activity
    this%pre_activity = pre_activity
  end subroutine set_pre_activity
  
  subroutine set_lipid_floor(this, lipid_floor)
    class(Calanus) :: this
    real           :: lipid_floor
    this%lipid_floor = lipid_floor
  end subroutine set_lipid_floor  
  
  subroutine set_ep_daily(this, ep_daily)
    class(Calanus) :: this
    real           :: ep_daily
    this%ep_daily = ep_daily
  end subroutine set_ep_daily  
  
  subroutine set_egg_clutch(this, egg_clutch)
    class(Calanus) :: this
    integer        :: egg_clutch
    this%egg_clutch = egg_clutch
  end subroutine set_egg_clutch 
  
  subroutine set_start_ep(this, start_ep)
    class(Calanus) :: this
    real           :: start_ep
    this%start_ep = start_ep
  end subroutine set_start_ep  
  
  subroutine set_reproduce_info(this, reproduce_info)
    class(Calanus) :: this
    real,dimension(1:6)  :: reproduce_info
    this%reproduce_info = reproduce_info
  end subroutine set_reproduce_info
  
  subroutine set_duration_ep(this, duration_ep)
    class(Calanus) :: this
    real           :: duration_ep
    this%duration_ep = duration_ep
  end subroutine set_duration_ep  
  
  subroutine set_EP(this, EP)
    class(Calanus) :: this
    integer        :: EP
    this%EP = EP
  end subroutine set_EP
  
  subroutine set_EP_CLUTCH(this, EP_CLUTCH)
    class(Calanus) :: this
    integer        :: EP_CLUTCH
    this%EP_CLUTCH = EP_CLUTCH
  end subroutine set_EP_CLUTCH
  
  subroutine set_counter_eggs(this, counter_eggs)
    class(Calanus) :: this
    real           :: counter_eggs
    this%counter_eggs = counter_eggs
  end subroutine set_counter_eggs
  
  subroutine set_cum_eggs(this, cum_eggs)
    class(Calanus) :: this
    integer        :: cum_eggs
    this%cum_eggs = cum_eggs
  end subroutine set_cum_eggs
  
  subroutine set_max_diap_lipid(this, max_diap_lipid)
    class(Calanus) :: this
    real           :: max_diap_lipid
    this%max_diap_lipid = max_diap_lipid
  end subroutine set_max_diap_lipid  
  
  subroutine set_lipid_exit_diap(this, lipid_exit_diap)
    class(Calanus) :: this
    real           :: lipid_exit_diap
    this%lipid_exit_diap = lipid_exit_diap
  end subroutine set_lipid_exit_diap  
  
  subroutine set_moult_mass_limit(this, mass_limit)
    class(Calanus) :: this
    real           :: mass_limit
    this%mass_limit = lipid_exit_diap
  end subroutine set_moult_mass_limit  
  
  subroutine set_massp(this, massp)
    class(Calanus) :: this
    real           :: massp
    this%massp = massp
  end subroutine set_massp  
  
  
  ! Moulting process 
  subroutine set_stage(this, stage)
    class(Calanus) :: this
    real           :: stage
    this%stage = stage
  end subroutine set_stage

  subroutine set_a_moltC(this, a_molt)
    class(Calanus) :: this
    real           :: a_molt
    this%a_molt = a_molt
  end subroutine set_a_moltC

  subroutine set_b_moltC(this, b_molt)
    class(Calanus) :: this
    real           :: b_molt
    this%b_molt = b_molt
  end subroutine set_b_moltC

  subroutine set_age(this, age)
    class(Calanus) :: this
    real           :: age
    this%age = age
  end subroutine set_age

  subroutine set_serial(this, serial)
    class(Calanus) :: this
    integer        :: serial
    this%serial = serial
  end subroutine set_serial

  subroutine set_mater(this, mater)
    class(Calanus) :: this
    integer        :: mater
    this%mater = mater
  end subroutine set_mater
  
  subroutine set_birthdate(this, birthdate)
    class(Calanus) :: this
    integer        :: birthdate
    this%birthdate = birthdate
  end subroutine set_birthdate

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!72
  !! CALANUS class procedures

  ! Daytime preferred depth
  
  real function zday_funC(this,ff,tt,bottom, rdvm) result(zval)
     class(Calanus)     :: this
     real               :: cur_stage,bottom,rdvm
     real,dimension(1:9), intent(in) :: ff,tt
 
     cur_stage=this%stage
     if (cur_stage<=3 .and. this%activity==0) then
        zval = 25.0
     endif
     if (cur_stage>3.0 .and. cur_stage<9.0 .and. this%activity==0) then
        zval=min(mid_depth(this%max_grow_depth(ff,tt)),75.0)
     endif
     if (cur_stage>=9.0 .and. this%activity==0 .and. rdvm>=0.5) then 
        zval=100.0
     else  ! (cur_stage>=9.0 .and. this%activity==0 .and. rdvm<0.5) then 
        zval=min(mid_depth(this%max_grow_depth(ff,tt)),75.0)
     endif
     if (cur_stage>=10.0 .and. this%activity>=1) then
        if (mid_depth(this%max_grow_depth(ff,tt)) > 80.0 ) then
          zval=0.8042*mid_depth(this%max_grow_depth(ff,tt))
        else
          zval=bottom 
        endif
     endif

  end function zday_funC


  ! Nightime preferred depth
 
  real function znight_funC(this,ff,tt,bottom) result(zval)
     class(Calanus)     :: this
     real               :: cur_stage,bottom
     real,dimension(1:9), intent(in) :: ff,tt

     cur_stage=this%stage
     if (cur_stage<=3 .and. this%activity==0 ) then
        zval = 25.0
     endif
     if (cur_stage>3 .and. cur_stage<9 .and. this%activity==0) then
        zval=min(mid_depth(this%max_grow_depth(ff,tt)),75.0)
     endif
     if (cur_stage>=9 .and. this%activity==0) then
        zval=min(mid_depth(this%max_grow_depth(ff,tt)),75.0)
     endif
 
     if (cur_stage>=10 .and. this%activity>=1) then
        if (mid_depth(this%max_grow_depth(ff,tt)) > 80.0 ) then
          zval=0.8042*mid_depth(this%max_grow_depth(ff,tt))
        else
          zval=bottom
        endif
     endif
 
  end function znight_funC   


  ! Arrhenius function of temperature = Universal Temperature Dependance
  ! (Gillooly et al. 2001)
  ! Returns a temperature dependance coefficient (unitless)
  real function arrheniusC(this, e_act)
    class(Calanus)  :: this
    real            :: e_act
    real, parameter :: k  = 8.62e-5 ! The Boltzmann's constant in eV.K-1
    real, parameter :: T0 = 273.15  ! Reference temperature (the frozing point of water) in K

    arrheniusC = exp( ( e_act * this%tp ) / ( k * ( this%tp + T0 ) * T0 ) )

  end function arrheniusC

  ! Arrhenius adjustement function ensures an exponential decrease
  ! AFTER some threshold
  real function adj_arrC(this, alpha_temp)
    class(Calanus) :: this
    real           :: alpha_temp, dtemp

    dtemp = this%t_lim - this%tp

    if( dtemp > 0 ) then 
       adj_arrC = 1
    else 
       adj_arrC = exp( dtemp * alpha_temp )
    endif

  end function adj_arrC

  ! age of individual(s) in days
  real function grow_olderC(this, age) ! real function grow_olderC(this) result(age)
    class(Calanus) :: this
    real           :: dt

    age = this%age + (dt / 86400)  ! age in days, increases by half hour 
    !grow_olderC = age ! this%age + (dt / 86400)  ! age in days, increases by half hour 
    
  end function grow_olderC


  real function devlpC(this)
    class(Calanus)  :: this
    real, parameter :: x_molt = -2.05

    devlpC = max( 86400. * 2., &
                  86400. * ( this%a_molt * (this%tp - this%b_molt) ** x_molt ) )

  end function devlpC

  ! Development freq of the individual. 
  subroutine moltC(this)
    class(Calanus) :: this
    real           :: newstage, a_moult_old
    integer        :: stg, moult


    if ( (this%cell_out .eqv. .false.) .and. (this%activity == 0) .and. (this%stage < 12.0)) then

       stg = int(this%stage)
       if (stg .ge. 4) then
           mass_moult = crit_moult_max(stg-3)  ! crit_moult_max provided as input in run.list
       endif

       newstage = this%stage + ( dt / this%devlp() ) ! theoretical new stage if mass criterion met
       moult = int(newstage) - stg             ! 1 = moult occurs (only if mass criterion met)

       if (stg==11) then 
           moult = 0   ! Never allow C5-->C6 moult while active
       endif

       ! At the end of the development phase, the individual molts
       if (moult .eq. 1) then
          a_molt_old = this%a_molt             ! save old a_molt 
          if (this%mass >= mass_moult) then    ! must adjust development occurring past moult
              this%a_molt = a_devlp(stg+1)     ! overshoot past moult = newstage - int(newstage)
              this%stage = int(newstage) + (newstage - int(newstage))*(this%a_molt/a_molt_old) 
          endif
          if ((this%stage>=12.0) .and. (this%mass_limit .eq. 0)) then
              this%mass_limit = this%mass * 1.15  ! set mass limit after terminal moult
              this%abund = this%abund*0.8         ! kill males that molted while active 
              this%sex = 1                        ! set sex to first year female == 1 
          endif

       elseif (moult .eq. 0) then 
       
          if ((stg .ge. 4) .and. (newstage<12) .and. (this%mass >= mass_moult)) then
             this%stage=newstage
          endif
          if (stg .lt. 4) then
             this%stage = newstage 
          endif
       
       endif

       if (this%stage>=12.0) then
           this%stage = 12.0
       endif

    else

       this%stage=this%stage

    endif


  end subroutine moltC

  ! Respiration = basal metabolism + activity (swimming) metabolism
  ! but NOT specific dynamic action as it is included in assimilation
  ! coefficient; based on Bandara et al. (2019)
  real function respC(this)
    class(Calanus)  :: this
    real            :: f, o, g, p
    real            :: B0, Bb, Ba
    f   =paramosome(3)
    o   =paramosome(4)  ! ooo =paramosome(4)
    g   =paramosome(5)
    p   =paramosome(6)

    if ((this%stage>=10.0) .and. (this%stage<11)) then
       ! if (this%activity==0) then
            f = 0.80*f ! lower CIV resp by 20% (act & diap)
       ! else
       !     f = 0.85*f ! lower CIV resp by 15% during diapause
       ! endif
    endif

    B0=f*((this%mass+this%lipid)**o)*0.5 
    Bb=B0*g*exp(p*this%tp);
    Ba=1.5*Bb
    if (this%activity==1) then
        respC =  0.25*Bb
    endif
    if (this%activity==0) then
        respC =  Ba+Bb
    endif
    if (this%activity==2) then
        !respC = 0.62*Bb ! KB turn off, gonad maturation same as diapause
        respC =  0.25*Bb
    endif
    if (this%activity==3) then
        respC = 0.62*Bb
    endif

  end function respC

  ! Updates the new mass of the individual
  subroutine growC(this)
    class(Calanus) :: this
    real           :: grow_new, lipid_new, mass_new
    real           :: rate, g_sup, dt
    real           :: we  !following Plourde et al. 2003; Huntley & Lopez 1992
    real           :: ingest_co
    we        =paramosome(13)
    rate      =paramosome(8)
    ingest_co =paramosome(7)
    

    if ( this%cell_out .eqv. .true. ) then 
       this%mass = this%mass
       this%lipid=this%lipid
    else
       ! CIV grow more using mean of Tande and Slagstad values
       !if ( (this%stage>=10) .and. (this%stage<11) ) then
       !    ingest_co=0.864
       !endif

       grow_new=ingest_co*this%ingest()-this%resp()-this%EP*we
       
       !if (grow_new>0) then

           g_sup=max(grow_new,0.0) 
           
           if (this%mass<38.0) then
               lipid_new=0.0
           endif
           if (this%mass>=38.0 .and. this%mass<=159.0) then
               lipid_new=g_sup*rate/(1.0+exp((60.0-this%mass)/20.0))
           endif
           if (this%mass >159.0) then
               lipid_new=g_sup*rate
           endif
           if (this%stage<9.0) then
               lipid_new=0.0
           endif
   
       ! KB: new stage requirements: CIII up to 10% lipid, 
       !     CIV & CV partition 'rate'% lipid, 
       !     CVI structural growth ceases, all lipid accum, up to max 70% of DW
       
           if ((this%stage>=9.0) .and. (this%stage<10.0)) then
               if (lipid_new+this%lipid <= 0.1*(this%mass+this%lipid)) then
                   this%lipid = this%lipid + lipid_new
                   mass_new = max(g_sup-lipid_new,0.0)
               else
                   lipid_new = max(0.1*(this%mass+this%lipid)-this%lipid, 0.0)
                   this%lipid = this%lipid + lipid_new
                   mass_new = max(g_sup-lipid_new,0.0)
               endif
           elseif ((this%stage>=10.0) .and. (this%stage<12.0)) then 
               if (lipid_new+this%lipid <= 0.7*(this%mass+this%lipid)) then
                   this%lipid = this%lipid + lipid_new
                   mass_new = max(g_sup-lipid_new,0.0)
               else
                   lipid_new = max(0.7*(this%mass+this%lipid)-this%lipid, 0.0)
                   this%lipid = this%lipid + lipid_new
                   mass_new = max(g_sup-lipid_new,0.0)
               endif
           elseif (this%stage>=12.0) then ! KB added capability for newly moulted CVIF to grow 20% larger 
           ! Next: can allow second year Females to add 10% body mass
               if (this%mass<this%mass_limit) then ! CVIs can add 15% mass due to moult 
                   if (lipid_new+this%lipid <= 0.7*(this%mass+this%lipid)) then
                       this%lipid = this%lipid + lipid_new
                       mass_new = max(g_sup-lipid_new,0.0)
                   else
                       lipid_new = max(0.7*(this%mass+this%lipid)-this%lipid, 0.0)
                       this%lipid = this%lipid + lipid_new
                       mass_new = max(g_sup-lipid_new,0.0)
                   endif
                   mass_new = min(mass_new, this%mass_limit-this%mass) ! ensure not growing beyond mass_limit
               else
                   mass_new = 0.0
                   lipid_new = max(g_sup,0.0)
                   if (lipid_new+this%lipid <= 0.7*(this%mass+this%lipid)) then
                       this%lipid = this%lipid + lipid_new
                   else
                       lipid_new = max(0.7*(this%mass+this%lipid)-this%lipid, 0.0)
                       this%lipid = this%lipid + lipid_new
                   endif
               endif
           elseif (this%stage<9.0) then ! describe growth for N3-C2 (no lipid storage)
               mass_new = max(g_sup,0.0) 
           endif

           this%mass=this%mass+mass_new 
            
        !endif
       
        if (this%abund > 0) then
        ! alive individuals grow older by one timestep
            this%age = this%age + 0.0208 ! (dt / 86400)  ! age in days, increases by half hour 
        endif
  endif

  end subroutine growC
  
  
  ! KB: new option for staged temp-dependent mortality (Maps et al. 2012) 
  subroutine mort_maps2012(this, z_depth, par_surf)  
     class(Calanus)  :: this
     real, intent(in):: z_depth, par_surf !temps, sss, ph
    ! real,parameter  :: factor_n2c=79.5729    
    ! real,parameter  :: ratio_c2chl=0.0180
     real,parameter  :: iomax=360 !,max_body_vis=2500 ! ug C 
     real,parameter  :: min1=0, max1=iomax, min2=1, max2=1.5 !, kw=0.04
     real            :: lin_trans_term1 !, phyt_surf, kcdom, kchla
     real            :: iv, iprime
     real,parameter  :: m1=0.5, m2=0.7, Tref=7.5, kd=0.06
     real            :: sd_num, sd_denom, mort_maps, d2ts !vis

     sd_num = (Tref+13.66)**(-2.05)
     sd_denom = (this%tp+13.66)**(-2.05)
     d2ts = 0.020833       ! daily to 30-min conversion 1/48

  ! calculate daily mortality   
     mort_maps = (0.01+m1*exp(-m2*floor(this%stage)))*(sd_num/sd_denom)
     
     if (mort_maps>0.98) then
         mort_maps=0.98
     elseif (mort_maps<0) then
         mort_maps=0
     endif

  ! get surface chla at (x,y)
     lin_trans_term1=(max2-min2)/(max1-min1)
  !   phyt_surf=ph*factor_n2c*ratio_c2chl
  
  ! calculate light attenuation kd - KB: scale down effect of kcdom, otherwise no light below 10 m
  !   kcdom=-0.0364*sss+1.1942
  !   kchla=0.0518*phyt_surf**(0.428)
  !   kd=kw+0.05*kcdom+kchla  ! kd=kw+kcdom+kchla

  ! calculate light dependent mortality 
     iv = par_surf*exp(-kd*z_depth)        ! based on surface light and particle depth
     iprime=lin_trans_term1*(iv-max1)+max2 ! factor 1 to 1.5
  
  ! include light scaling in mort_maps (increases with light up to 50%), convert to 30-min timestep:
     if (this%stage<1) then
         mort_maps = mort_maps*d2ts
     elseif ((this%stage>=1) .and. (this%activity==0)) then
         mort_maps = iprime*mort_maps*d2ts
  ! reduce mort for copepods diapausing or reproducing at depth to 10% of active   
     elseif ((this%stage>=10) .and. (this%activity>0)) then
         mort_maps=0.1*mort_maps*d2ts  
     endif

  ! calculate the abundance change due to temp-dependent mort
     this%abund=this%abund*(1.0-mort_maps)

  end subroutine mort_maps2012
  

  ! mortality risk by visual and non-visual predation module. (Bandara et al. 2019)
  
  subroutine predation (this, sss, z_depth, par_surf, temps, ph)
     class(Calanus)  :: this
     real, intent(in):: sss,z_depth, par_surf, temps, ph
     real,parameter  :: factor_n2c=79.5729    
     real,parameter  :: ratio_c2chl=0.0180
     real,parameter  :: iomax=360
     real,parameter  :: min1=0, max1=iomax, min2=0.1, max2=0.9, kw=0.04
     real            :: lin_trans_term1,phytot, phyt_surf
     real            :: io, kcdom, kchla, iv, iprime
     real,parameter  :: k_exp=0.002, k_prime=0.15
     real            :: kmort, mort_d, mort_nonvis_daily_max, mort_nonvis_daily_min, mort_asym
     real            :: ino, mort_vis_daily, mort_vis_dt  
     real            :: minn1, minn2, lin_trans_term2, maxx1, maxx2
     real, parameter :: tmin=0.0
     real            :: mort_nonvis, mort_nonvis_dt, tmp
  
  ! visual predation
  
  !get surface chla at (x,y)
     lin_trans_term1=(max2-min2)/(max1-min1)
     phytot=ph*factor_n2c
     phyt_surf=phytot*ratio_c2chl
  
  ! calculate light attenuation kd
     io=par_surf
     kcdom=-0.0364*sss+1.1942
     kchla=0.0518*phyt_surf**(0.428)
     kd=kw+kcdom+kchla
  
  
  ! calculate light dependent mortality
  
     iv = iomax*exp(-kd*z_depth)
     iprime=lin_trans_term1*(iv-max1)+max2
  
  
     kmort=k_prime**4.0
     mort_d=kmort*24.0
     mort_asym=mort_d*2.5
     mort_nonvis_daily_max=0.1*mort_asym
     mort_nonvis_daily_min=0.01*mort_asym
  
     mort_vis_daily = iprime*mort_asym*(1.0-exp(-k_exp*this%mass))*0.33 ! KB CUT MORT TO 30% TEMP
     mort_vis_dt=mort_vis_daily/48.0
 
   ! non_visual_predation	  
  
     minn1=tmin
     maxx1=18.59 ! the number was calculated in matlab, gulf region, 20180101 to 20190101, top 4 layer averaged
     minn2=mort_nonvis_daily_min
     maxx2=mort_nonvis_daily_max
     lin_trans_term2=(maxx2-minn2)/(maxx1-minn1)
   
     tmp=max(0.0, temps)
     mort_nonvis=lin_trans_term2*(tmp-maxx1)+maxx2
     mort_nonvis_dt=mort_nonvis/48.0
  
  ! calculate the abundance change due to visual and non visual predation
     this%abund=this%abund*(1.0-mort_vis_dt-mort_nonvis_dt)
   
  
  end subroutine predation
  
  
  
  ! calculate mortality rate caused by starving
  ! Brandara et al. (2019, Sec. 2.2.2.2.7)
  subroutine starve(this)
     class(Calanus) :: this
     real           :: burn_prop, real_grow ! massp
     real           :: we  !following Plourde et al. 2003; Huntley & Lopez 1992
     real           :: ingest_co
     !massp    =0.0
     we       =paramosome(13)
     ingest_co=paramosome(7)
     
     !if (this%stage>=10 .and. this%stage<11) then
     !    ingest_co=0.864
     !endif
     ! calculate the net growth
     if (this%stage < 3.0) then
         real_grow = 0.0
         !this%massp = 0.0
     else
         real_grow=ingest_co*this%ingest()-this%resp()-we*this%EP
     endif
     
     ! get structure mass
     this%massp=max(this%massp, this%mass) ! max mass achieved before starvation
     burn_prop=0.0
     
     if (real_grow<0) then   ! starvation induced loss of lipid, then mass
        if (this%lipid>0) then
           if (abs(real_grow)<=this%lipid) then
              this%lipid=this%lipid+real_grow
           else
              this%mass=this%mass+(this%lipid+real_grow)
              this%lipid=0.0
           endif
        else
           this%mass=this%mass+real_grow
        endif
  
        burn_prop=(this%massp-this%mass)/this%massp
     endif

     if (this%stage >= 3.0 )then   ! calculate starvation mortality
         if (burn_prop<=0.1)then
             mort_starve=0.0
         endif
  
         if(burn_prop>0.1 .and. burn_prop<0.5) then
             mort_starve=2.0*burn_prop
         endif
  
         if(burn_prop>=0.5)then
             mort_starve=1.0
         endif
     else
         mort_starve=0.0
     endif
 

      ! calculate the abundance change due to starvation
      this%abund=this%abund*(1.0-mort_starve)
      
  end subroutine starve

  ! calculate the depth with the max grow rate based on food and temperature
  ! only used the top 9 layers
  !unit of max_grow rate is the sequence of vertical layers
  real function max_grow_depth(this, ff,tt) 
     class(Calanus)                    :: this
     real, dimension(1:9),  intent(in) :: ff,tt! food and temperature of top 9 layers
     real,dimension(1:9)               :: gd
     real                              :: d

     if ( this%cell_out .eqv. .true. ) then 
         max_grow_depth=9 ! limit the max depth value to 75 m
     else 
     ! calculate the index of grow to simplify the equation (not exactly the grow rate)
     ! equation simplified from the function ingest
         d=0.3*(this%mass)**(-0.138)
         gd=exp(0.0761*tt)*d*ff/(1.0+d*ff)
         max_grow_depth=maxloc(gd,DIM=1)
     endif

  end function max_grow_depth
  
  
  real function ingest(this)
    class(Calanus)  :: this
    real            :: b, c, m, n, kf
    !real, parameter :: kf=25.3 !35.3       ! half saturation constant for Holling Type 3
    real            :: ingest_ref, F, holling3_term   ! KB remove 'd' 
    b = paramosome(9)
    m = paramosome(10)
    c = paramosome(11)
    n = paramosome(12)
    kf = paramosome(14)
    
    if ( this%cell_out .eqv. .true. ) then 
       ingest = 0.
    else
        if ( (this%stage>=10) .and. (this%stage<11)) then
            b = b*1.1
        endif
    ! if stage <3 (was 4 until Jun9, 1st feeding stage likely N3); during EP and diapause, no ingestion
       if (this%stage < 3.0 .or. this%activity > 0) then
          ingest = 0.
       else
          ingest_ref=0.5*b*((this%mass))**m  ! 30-min timestep so multiply by 0.5 here
          F=this%fd*79.57  ! convert mmol N/m3 to ug C/L; no food multiplier
          if (time(5)>=9) then  ! reduce model food by 40% in fall when too high
              F=0.6*F
          endif
          holling3_term = (F**2)/(kf**2 + F**2) ! type 3 holling where kf is half saturation
          !ingest=max(ingest_ref*c*exp(n*(this%tp))*((d*F)/(1.0+d*F))*0.5,0.0)
          ingest=max(ingest_ref*c*exp(n*(this%tp))*holling3_term,0.0)
        endif
    endif

  end function ingest
  
  
  
  subroutine egg_pro(this)
     class(Calanus)     :: this
     real               :: Bb, resp
     real               :: we  !following Plourde et al. 2003; Huntley & Lopez 1992
     real               :: f, oo, p, g
     real,parameter     :: f_resp_ep=0.62, clutch_size=75 
     real               :: T_part
     real               :: ep2
     integer            :: day_hour
     
    ! activity=2: start producing egg, and initialize lipid_floor for egg production
    ! activity=3: producing egg, taking mass from lipid initialized above.

     f =paramosome(3)
     oo=paramosome(4)
     g =paramosome(5)
     p =paramosome(6)
    ! different than the equation in slide. this if should be activity==2, 
    ! since the lipid_floor should only be set once at preparation, or the lipid_floor will change
     this%EP = 0.
     this%EP_CLUTCH = 0. ! KB new
     we=paramosome(13)
     if (this%activity == 2) then
         T_part=this%tp
         if (this%sex == 1) then
             this%lipid_floor = 0.2*this%lipid_exit_diap
         endif
 
         if(this%sex == 2) then
             this%lipid_floor = 0.0
         endif
 
      endif
    ! counter_c6active is in unit of minutes.
      if (this%activity ==2 ) then
          ep2=0.
          this%counter_c6_active = this%counter_c6_active + 30.0
          ! gonad maturation for 3 weeks (no eggs produced)
          if (this%counter_c6_active >= (3.0*7.0*24.0*60.0-0.1)) then
              this%counter_eggs = 0.
              this%counter_c6_active = 0.
              this%start_ep = modeltime
              this%duration_ep = 0.
              this%ep_daily = 0.
              this%activity = 3
              this%egg_clutch = 0.
          endif
       endif
 
       day_hour=time(3) ! KB remove ';'
       
       if ((this%activity == 3)  .and. (this%lipid > this%lipid_floor)) then
           Bb = f_resp_ep*f*(this%mass+this%lipid)**oo ! KB Now M+L
           resp = 0.5*Bb*g*exp(p*T_part)
 
           if (day_hour > 0 ) then ! not midnight
               this%ep_daily = this%ep_daily+2.0*resp
               ep2 = 0.
           else
               ep2 = this%ep_daily+ 2.0*resp
               this%EP = int(ep2/we)
               this%ep_daily = max(0.0, ep2 - this%EP*we)
               this%duration_ep =this%duration_ep+1.       ! add a day
               this%counter_eggs=this%counter_eggs+this%EP
               this%cum_eggs = this%cum_eggs + this%EP
               this%lipid = this%lipid - this%EP*we
               this%egg_clutch=this%egg_clutch + this%EP   ! accumulate eggs in clutch
               if (this%egg_clutch>=clutch_size) then
                   this%EP_CLUTCH = clutch_size
                   this%egg_clutch = this%egg_clutch - clutch_size ! extra eggs to next clutch
               endif
           endif
        endif
 
        if (this%activity == 3 .and. this%lipid <=this%lipid_floor) then
            ep2 = 0.
            this%EP = 0.
            this%egg_clutch = 0.
            this%EP_CLUTCH = 0.
            this%activity = 0
            if (this%sex ==1 ) then
                this%reproduce_info = [this%start_ep, this%duration_ep, this%counter_eggs,0.0,0.0,0.0]
            endif
 
            if(this%sex == 2 )then
               this%reproduce_info(4:6) = [this%start_ep, this%duration_ep, this%counter_eggs]
               this%abund=0.
            endif
         endif
 
  
  end subroutine egg_pro
  
  subroutine diapause_entry(this)
       class(Calanus)     :: this
       real               :: lipidfrac
       real               :: first_diap_stage
  
       first_diap_stage=paramosome(1)
       lipidfrac=this%lipid/(this%lipid+this%mass)
       if ((this%stage >= 11.0) .and. (this%ep1 <= 0.46)) then  ! ep1 = 0.33-0.40
           ! rescale lipid threshold upwards for CV-CVI from 0.33-0.40 to 0.41-0.48
           this%ep1 = ((this%ep1 - 0.33)/0.07)*0.07 + 0.41 ! 1st 0.1=0.40-0.33, 2nd 0.1=0.48-0.41
       endif
       if ((this%stage >= first_diap_stage)  .and. (lipidfrac>=this%ep1) .and. (this%activity==0)) then
           this%flag_enter_diap = .true.
       endif
  
       if((this%flag_enter_diap .eqv. .true.) .and. (this%stage<12.0) .and. &
               (this%stage-int(this%stage))>=0.985 .and. (this%activity==0)) then 
          this%activity = 1
          this%flag_enter_diap = .false.
       endif
  
       if ((this%flag_enter_diap .eqv. .true.) .and. (this%stage>=12.0) .and. (this%activity==0)) then
           this%counter_c6_active = this%counter_c6_active + 30.0
       endif
  
       if ((this%counter_c6_active >= 14.0*24.0*60) .and. (this%activity==0)) then  !(2 weeks counting in minutes) 
           this%activity = 1
           this%counter_c6_active = 0.0
           this%flag_enter_diap = .false.
           if ((this%sex.eq.1) .and. (this%cum_eggs.gt.0)) then
               this%sex = 2
           endif
       endif
       
  end subroutine diapause_entry
  
  
  
  subroutine diapause_exit(this)
     class(Calanus)     :: this
     real               :: lipidfrac_burned, new_stage
     logical            :: flag_exit_diap
     flag_exit_diap = .false.
     lipidfrac_burned = 0.0
     if (this%activity == 1) then   

         if (this%activity-this%pre_activity ==1) then
            this%max_diap_lipid = this%lipid
            flag_exit_diap = .false.
         endif
 
         lipidfrac_burned = (this%max_diap_lipid - this%lipid)/this%max_diap_lipid
         
         if ((this%stage >= 10.0) .and. (this%stage < 11) .and. (this%ep2 <= 0.6)) then
             ! rescale lipid threshold upwards for CIV from 0.4-0.6 to 0.75-0.95
             this%ep2 = ((this%ep2 - 0.4)/0.2)*0.2 + 0.75 ! 1st 0.2=0.6-0.4, 2nd 0.2=0.95-0.75
         endif
         if ((this%stage >= 11) .and. (this%ep2 > 0.6)) then
             ! undo scaling for CV-CVI (return back to 0.4-0.6)
             this%ep2 = ((this%ep2 - 0.75)/0.2)*0.2 + 0.4 ! 1st 0.2=0.95-0.75, 2nd 0.2=0.6-0.4
         endif

         if (lipidfrac_burned >= this%ep2) then
             flag_exit_diap = .true.
         endif
         new_stage = this%stage 
         if ((flag_exit_diap .eqv. .true.) .and. (this%stage <= 12.0)) then
             new_stage = min(real(int(this%stage+1.0)),12.0)
             
             if ((this%stage >= 10.0) .and. (this%stage < 11.0)) then
                 this%activity = 0
                 flag_exit_diap=.false.      ! reset flag to default (F)
             endif
  
             !if ((this%stage >=11.0) .and. (this%stage <= 12.0)) then
             if ((this%stage >=11.0) .and. (this%stage < 12.0)) then
                 this%activity = 2
                 this%abund = 0.8*this%abund  ! 80% female sex ratio; kill males immediately
                 this%sex = 1                 ! set sex to Female (first year) 
                 this%stage = 12.0            ! CVs moult to adult when leaving diapause 
                 flag_exit_diap =.false.      ! reset flag to default (F)
                 this%lipid_exit_diap=this%lipid
                 if (this%mass_limit .eq. 0) then
                     this%mass_limit=this%mass * 1.15 ! set CVI mass limit to 15% above current CV mass
                 endif
              endif
          endif
  
  
          if ((flag_exit_diap .eqv. .true.) .and. (this%stage >= 12.0)) then
              new_stage = this%stage
              !this%sex = 2 ! Now increment sex to 2nd year female when enter diapause after EP>0
              this%activity = 2
              
              flag_exit_diap=.false.
          endif
  
  
          this%stage = new_stage
      endif
  
  end subroutine diapause_exit





!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!72
!! Late_hyperboreus class procedures
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!72


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!72
!! Getters

! Environmental variables
real function get_tpLH(this)
  class(Late_hyperboreus) :: this
  get_tpLH = this%tp
end function get_tpLH

! Size
real function get_massLH(this)
  class(Late_hyperboreus) :: this
  get_massLH = this%mass
end function get_massLH


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!72
!! Setters

! Environmental variables
subroutine set_tpLH(this, tp)
  class(Late_hyperboreus) :: this
  real         :: tp
  this%tp = tp
end subroutine set_tpLH

! Size
subroutine set_massLH(this, mass)
  class(Late_hyperboreus) :: this
  real         :: mass
  this%mass = mass
end subroutine set_massLH

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!72
!! Late_hyperboreus class procedures

! Daytime preferred depth
! few meters above bottom depth. If bottom depth is shallow, on bottom depth.
real function zday_funLH(this,bottom, the_mode, in_diapause) result(zval)
   class(Late_hyperboreus)     :: this
   real               :: bottom, coeff_lo, intercept_lo
   integer            :: the_mode
   integer            :: in_diapause

   if (in_diapause == 0) then
      ! 0 = in diapause
      if (the_mode == 1) then
         ! upper mode
         ! optimal depth is XXm except if the bottom depth is shallower
         intercept_lo = 53.4
         zval = min(bottom, intercept_lo)
      else
         ! lower mode
         intercept_lo = 16.3
         coeff_lo     = 0.79
         if (bottom <= 150) then
            zval = (intercept_lo + 150 * coeff_lo) * bottom / 150
         else
            zval = intercept_lo + coeff_lo * bottom
         endif
      endif
   else
      ! 1 = active
      if (the_mode == 1) then
         ! upper mode
         ! optimal depth is XXm except if the bottom depth is shallower
         intercept_lo = 41
         zval = min(bottom, intercept_lo)
      else
         ! lower mode
         intercept_lo = 1.6
         coeff_lo     = 0.69
         if (bottom <= 150) then
            zval = (intercept_lo + 150 * coeff_lo) * bottom / 150
         else
            zval = intercept_lo + coeff_lo * bottom
         endif
      endif
   endif

end function zday_funLH


! Nightime preferred depth
! few meters above bottom depth. If bottom depth is shallow, on bottom depth.
real function znight_funLH(this,bottom, the_mode, in_diapause) result(zval)
   class(Late_hyperboreus)     :: this
   real               :: bottom, coeff_lo, intercept_lo
   integer            :: the_mode
   integer            :: in_diapause

   if (in_diapause == 0) then
      ! 0 = in diapause
      if (the_mode == 1) then
         ! upper mode
         ! optimal depth is XXm except if the bottom depth is shallower
         intercept_lo = 53.4
         zval = min(bottom, intercept_lo)
      else
         ! lower mode
         intercept_lo = 16.3
         coeff_lo     = 0.79
         if (bottom <= 150) then
            zval = (intercept_lo + 150 * coeff_lo) * bottom / 150
         else
            zval = intercept_lo + coeff_lo * bottom
         endif
      endif
   else
      ! 1 = active
      if (the_mode == 1) then
         ! upper mode
         ! optimal depth is XXm except if the bottom depth is shallower
         intercept_lo = 41
         zval = min(bottom, intercept_lo)
      else
         ! lower mode
         intercept_lo = 1.6
         coeff_lo     = 0.69
         if (bottom <= 150) then
            zval = (intercept_lo + 150 * coeff_lo) * bottom / 150
         else
            zval = intercept_lo + coeff_lo * bottom
         endif
      endif
   endif

end function znight_funLH


! diffusivity around the prefered depth
real function sigma_funLH(this,bottom, the_mode, in_diapause) result(sigma)
   class(Late_hyperboreus)     :: this
   real               :: bottom
   integer            :: the_mode
   integer            :: in_diapause

   ! If in diapause
   if (in_diapause == 0) then
      ! If it is the upper mode
      if (the_mode == 1) then
        sigma = 0.83
      ! if it is the lower mode
      else
        sigma = 0.002 * bottom !0.0054
      endif
   else
      if (the_mode == 1) then
        sigma = 0.8
      else
        sigma = 0.0022 * bottom !0.0061
      endif
   endif
end function sigma_funLH


! proportion of importance of upper mode cmpared to lower mode
real function mode_prop_funLH(this,part_prop, bottom, in_diapause) result(propval)
   class(Late_hyperboreus)     :: this
   real               :: bottom
   real               :: part_prop
   integer            :: in_diapause

   if (in_diapause == 0) then
      propval = 0.057
   else
      propval = 0.44
   endif

end function mode_prop_funLH



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!72
!! Young_hyperboreus class procedures
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!72


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!72
!! Getters

! Environmental variables
real function get_tpYH(this)
  class(Young_hyperboreus) :: this
  get_tpYH = this%tp
end function get_tpYH

! Size
real function get_massYH(this)
  class(Young_hyperboreus) :: this
  get_massYH = this%mass
end function get_massYH


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!72
!! Setters

! Environmental variables
subroutine set_tpYH(this, tp)
  class(Young_hyperboreus) :: this
  real         :: tp
  this%tp = tp
end subroutine set_tpYH

! Size
subroutine set_massYH(this, mass)
  class(Young_hyperboreus) :: this
  real         :: mass
  this%mass = mass
end subroutine set_massYH

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!72
!! Young_hyperboreus class procedures

! Daytime preferred depth
real function zday_funYH(this,bottom, the_mode) result(zval)
   class(Young_hyperboreus)     :: this
   real               :: bottom
   integer             :: the_mode

   if (the_mode == 1) then
      zval = 0.23 * bottom
   else
      zval = 0.77 * bottom
   endif
end function zday_funYH


! Nightime preferred depth
real function znight_funYH(this,bottom, the_mode) result(zval)
   class(Young_hyperboreus)     :: this
   real               :: bottom
   integer             :: the_mode

   if (the_mode == 1) then
      zval = 0.23 * bottom
   else
      zval = 0.77 * bottom
   endif
end function znight_funYH


real function sigma_funYH(this,bottom, the_mode) result(sigma)
   class(Young_hyperboreus)     :: this
   real               :: bottom
   integer            :: the_mode

   if (the_mode == 1) then
     sigma = 0.5707538 + 0.01409882 * bottom
   else
     sigma = 0.6450958 + 0.01423661 * bottom
   endif

end function sigma_funYH


! proportion of importance of upper mode cmpared to lower mode
real function mode_prop_funYH(this,part_prop, bottom, in_diapause) result(propval)
   class(Young_hyperboreus)     :: this
   real               :: bottom
   real               :: part_prop
   integer            :: in_diapause

   if (in_diapause == 0) then
      propval = 0.0784
   else
      propval = 0.546
   endif

end function mode_prop_funYH




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!72
!! Late_finmarchicus class procedures
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!72


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!72
!! Getters

! Environmental variables
real function get_tpLF(this)
  class(Late_finmarchicus) :: this
  get_tpLF = this%tp
end function get_tpLF

! Size
real function get_massLF(this)
  class(Late_finmarchicus) :: this
  get_massLF = this%mass
end function get_massLF


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!72
!! Setters

! Environmental variables
subroutine set_tpLF(this, tp)
  class(Late_finmarchicus) :: this
  real         :: tp
  this%tp = tp
end subroutine set_tpLF

! Size
subroutine set_massLF(this, mass)
  class(Late_finmarchicus) :: this
  real         :: mass
  this%mass = mass
end subroutine set_massLF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!72
!! Late_finmarchicus class procedures

! Daytime preferred depth
! few meters above bottom depth. If bottom depth is shallow, on bottom depth.
real function zday_funLF(this,bottom, the_mode, in_diapause) result(zval)
   class(Late_finmarchicus)     :: this
   real               :: bottom, coeff_lo, intercept_lo
   integer            :: the_mode
   integer            :: in_diapause

   if (in_diapause == 0) then
      ! 0 = in diapause
      if (the_mode == 1) then
         ! upper mode
         ! optimal depth is XXm except if the bottom depth is shallower
         intercept_lo = 37
         zval = min(bottom, intercept_lo)
      else
         ! lower mode
         intercept_lo = 34.4
         coeff_lo     = 0.56
         if (bottom <= 150) then
            zval = (intercept_lo + 150 * coeff_lo) * bottom / 150
         else
            zval = intercept_lo + coeff_lo * bottom
         endif
      endif
   else
      ! 1 = active
      if (the_mode == 1) then
         ! upper mode
         ! optimal depth is XXm except if the bottom depth is shallower
         intercept_lo = 41
         zval = min(bottom, intercept_lo)
      else
         ! lower mode
         intercept_lo = 1.6
         coeff_lo     = 0.69
         if (bottom <= 150) then
            zval = (intercept_lo + 150 * coeff_lo) * bottom / 150
         else
            zval = intercept_lo + coeff_lo * bottom
         endif
      endif
   endif
end function zday_funLF


! Nightime preferred depth
! few meters above bottom depth. If bottom depth is shallow, on bottom depth.
real function znight_funLF(this,bottom, the_mode, in_diapause) result(zval)
   class(Late_finmarchicus)     :: this
   real               :: bottom, coeff_lo, intercept_lo
   integer            :: the_mode
   integer            :: in_diapause

   if (in_diapause == 0) then
      ! 0 = in diapause
      if (the_mode == 1) then
         ! upper mode
         ! optimal depth is XXm except if the bottom depth is shallower
         intercept_lo = 37
         zval = min(bottom, intercept_lo)
      else
         ! lower mode
         intercept_lo = 34.4
         coeff_lo     = 0.56
         if (bottom <= 150) then
            zval = (intercept_lo + 150 * coeff_lo) * bottom / 150
         else
            zval = intercept_lo + coeff_lo * bottom
         endif
      endif
   else
      ! 1 = active
      if (the_mode == 1) then
         ! upper mode
         ! optimal depth is XXm except if the bottom depth is shallower
         intercept_lo = 41
         zval = min(bottom, intercept_lo)
      else
         ! lower mode
         intercept_lo = 1.6
         coeff_lo     = 0.69
         if (bottom <= 150) then
            zval = (intercept_lo + 150 * coeff_lo) * bottom / 150
         else
            zval = intercept_lo + coeff_lo * bottom
         endif
      endif
   endif

end function znight_funLF



real function sigma_funLF(this,bottom, the_mode) result(sigma)
   class(Late_finmarchicus)     :: this
   real               :: bottom
   integer            :: the_mode
   integer            :: in_diapause

   ! If in diapause
   if (in_diapause == 0) then
      ! If it is the upper mode
      if (the_mode == 1) then
        sigma = 0.90
      ! if it is the lower mode
      else
        sigma = 0.0021 * bottom !0.0057
      endif
   else
      if (the_mode == 1) then
        sigma = 0.80
      else
        sigma = 0.0022 * bottom !0.0061
      endif
   endif

end function sigma_funLF



! proportion of importance of upper mode cmpared to lower mode
real function mode_prop_funLF(this,part_prop, bottom, in_diapause) result(propval)
   class(Late_finmarchicus)     :: this
   real               :: bottom
   real               :: part_prop
   integer            :: in_diapause

   if (in_diapause == 0) then
      propval = 0.30
   else
      propval = 0.44
   endif

end function mode_prop_funLF





!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!72
!! Young_finmarchicus class procedures
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!72


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!72
!! Getters

! Environmental variables
real function get_tpYF(this)
  class(Young_finmarchicus) :: this
  get_tpYF = this%tp
end function get_tpYF

! Size
real function get_massYF(this)
  class(Young_finmarchicus) :: this
  get_massYF = this%mass
end function get_massYF


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!72
!! Setters

! Environmental variables
subroutine set_tpYF(this, tp)
  class(Young_finmarchicus) :: this
  real         :: tp
  this%tp = tp
end subroutine set_tpYF

! Size
subroutine set_massYF(this, mass)
  class(Young_finmarchicus) :: this
  real         :: mass
  this%mass = mass
end subroutine set_massYF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!72
!! Young_finmarchicus class procedures

! Daytime preferred depth
! few meters above bottom depth. If bottom depth is shallow, on bottom depth.
real function zday_funYF(this,bottom, the_mode) result(zval)
   class(Young_finmarchicus)     :: this
   real               :: bottom, coeff_lo, intercept_lo
   integer             :: the_mode

   if (the_mode == 1) then
      ! optimal depth is XXm except if the bottom depth is shallower
      intercept_lo = 28.2
      zval = min(bottom, intercept_lo)
   else
      intercept_lo = -15.9
      coeff_lo     = 0.74
      if (bottom <= 150) then
         zval = (intercept_lo + 150 * coeff_lo) * bottom / 150
      else
         zval = intercept_lo + coeff_lo * bottom
      endif
   endif

end function zday_funYF


! Nightime preferred depth
! few meters above bottom depth. If bottom depth is shallow, on bottom depth.
real function znight_funYF(this,bottom, the_mode) result(zval)
   class(Young_finmarchicus)     :: this
   real               :: bottom, coeff_lo, intercept_lo
   integer             :: the_mode

   if (the_mode == 1) then
      ! optimal depth is XXm except if the bottom depth is shallower
      intercept_lo = 28.2
      zval = min(bottom, intercept_lo)
   else
      intercept_lo = -15.9
      coeff_lo     = 0.74
      if (bottom <= 150) then
         zval = (intercept_lo + 150 * coeff_lo) * bottom / 150
      else
         zval = intercept_lo + coeff_lo * bottom
      endif
   endif

end function znight_funYF



real function sigma_funYF(this,bottom, the_mode) result(sigma)
   class(Young_finmarchicus)     :: this
   real               :: bottom
   integer            :: the_mode

   if (the_mode == 1) then
     sigma = 0.537
   else
     sigma = 0.00096 * bottom !0.0026
   endif

end function sigma_funYF



! proportion of importance of upper mode cmpared to lower mode
real function mode_prop_funYF(this,part_prop, bottom) result(propval)
   class(Young_finmarchicus)     :: this
   real               :: bottom
   real               :: part_prop

   propval = 0.86

end function mode_prop_funYF


!--------------------------------------------------------------------72

end module partSubClass_mod
