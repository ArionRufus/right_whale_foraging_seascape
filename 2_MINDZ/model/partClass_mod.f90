module partClass_mod

  use rng_par_zig ! random generator module
  use abstList_mod
  
  implicit none

  ! Constants used within classes
  real, parameter                      :: EPS = 10e-6

  !--------------------------------------------------------------------72
  ! Grid dimensions
  ! WARNING: here only interested in the first 20 layers (~550 m) !!!
  ! ilo = 46 in NEMO


#ifdef CIOPS
  integer,   parameter                 :: m = 801, n = 701, ilo = 69 !m = 1410, n = 945, ilo = 69
#else
  integer,   parameter                 :: m = 197, n = 234, ilo = 20
#endif

  real,      parameter                 :: pi = 3.141592653589

  ! Topography/domain
  
  integer                              :: bathy_id              ! bathymetry values id

  real,      dimension(ilo)            :: dd, dz, up_depth, mid_depth

  integer*2, dimension(m,n)            :: nlayer

  integer,   dimension(m,n,ilo)        :: bathy

  real,      dimension(m,n)            :: dlx, dly, lon, lat

  ! Time steps and dates

  integer                              :: initial, final, kstep, dt, dtnc, kd_ini

  integer                              :: istart, iend, ystart, mstart, dstart, yend, mend, dend

  integer,   dimension(3)              :: today, now, tnc_origin

  integer,   dimension(6)              :: time, tzero

  real                                 :: time_sunrise, time_sunset, p_night !time_go_down, time_go_up

  logical                              :: light_day, nc_repeat
  
  real                                 :: init_xpo, init_ypo

  character(len=100)                   :: tnc_origin_out, tnc_units

  ! Input file

  character(len=150)                   :: infile, infile2, infile_loc

  !--------------------------------------------------------------------72
  ! Particle tracking

  integer,   parameter                 :: max_part=10000000

  integer                              :: istart_part, iend_part

  integer                              :: part_spec, part_fac, part_freq, part_duration, init_part, in_diapause

  integer                              :: npart, npart0, tpart, isteps_part, output_freq, output_freq2
  
  integer                              :: y_Astart, m_Astart, d_Astart, y_Aend, m_Aend, d_Aend
  
  integer                              :: next_available_part

  integer,   dimension(m,n)            :: nxy

  integer,   dimension(:), allocatable :: ipart, jpart, part_mater, part_serial, part_birthdate

  real,      dimension(m,n)            :: part_grid_init, salt2, phyto, parsurf, temp2

  real,      dimension(m,n,2)          :: salt

  real,      dimension(m,n,ilo)        :: u, v, u1, v1, u2, v2, du, dv, temp
#ifdef VERT

  real,      dimension(m,n,ilo)        :: w, w1, w2, dw
#endif
#ifdef YUMY

  real,      dimension(m,n,ilo)        :: diatom, flag, microzoo, food, phyto3d
#endif

  ! Position inner variables & parameters
  real                                 :: xpo, ypo, zpo, bdep

  integer                                 :: ipo, jpo

  integer                                 :: mode  !upper (1) or lower (2) mode in the vertical distribution

  real                                   :: mode_prop  !importance of yhe upper mde compared to the lower mode

  real                                 :: xk1, xk2, xk3, xk4, &
                                          yk1, yk2, yk3, yk4, &
                                          zk1, zk2, zk3, zk4

  real                                 :: theta, sigma

  logical                              :: dvm, cell_out

  ! I/O
  logical                              :: part_traj_out, part_ftle_out, part_nc_out

  character(len=400)                   :: partfile, parsfile, outfile

  !--------------------------------------------------------------------72
  ! I/O

  namelist /run_nml/  ystart, mstart, dstart,             &
                      yend,   mend,   dend,               &
                      infile, infile2, infile_loc,  nc_repeat,         &
                      tnc_origin, tnc_origin_out,         &
                      tnc_units        

  namelist /part_nml/ cell_out,                           &
                      part_spec, in_diapause, part_serial, part_mater, &
                      part_birthdate, dvm, theta, sigma,  &
                      part_fac, part_freq, part_duration, &
                      output_freq, init_part,             &
                      y_Astart, m_Astart, d_Astart,       &
                      y_Aend,   m_Aend,   d_Aend,         &
                      partfile, parsfile, outfile,        &
                      init_xpo, init_ypo, time         

  !--------------------------------------------------------------------72
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! BELOW IS THE GENERIC LIST OF PARTICLES SECTION !!
  !! Users should NOT modify this module            !!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  public              :: partList
  type, extends(list) :: partList 

   contains

     ! FUNCTIONS

     ! Get number of particles in the list
     procedure         :: n_part
     
     ! Get particle pointed by iterator
     procedure, public :: current => currentPart

     ! SUBROUTINES

     ! Allow acces to a specific particle
     procedure         :: thisPart
    
     ! Add particle in list
     procedure, public :: addPart

  end type partList

  !--------------------------------------------------------------------72
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! BELOW IS THE GENERIC Zoo CLASS SECTION         !!
  !! Define individual advected particle properties !!
  !!                                                !!
  !! Users should NOT modify this module...         !!
  !! ...                                            !!
  !! ... unless you REALLY know what you're doing   !! 
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Allows acces to the Zoo class
  public :: Zoo

  type Zoo

     ! "Species", or rather particles category
     integer :: part_spec

     ! In diapause or not
     integer :: in_diapause

     ! Time
     integer :: istart
     integer :: isteps

     ! Swimming behaviour
     logical :: dvm, light_day

     real    :: zday, znight, sigma_local
 
     ! Particle status
     integer :: part_status

     ! Transport process
     logical :: cell_out

     real    :: xpo
     real    :: ypo
     real    :: zpo
     real    :: bdep
     integer :: ipo, jpo
     integer :: mode
     real    :: mode_prop
     real    :: xk1, xk2, xk3, xk4
     real    :: yk1, yk2, yk3, yk4
     real    :: zk1, zk2, zk3, zk4

   contains

     ! Getters
     procedure, public :: get_species
     procedure, public :: get_diapause

     procedure, public :: get_istart
     procedure, public :: get_isteps

     procedure, public :: get_dvm
     procedure, public :: get_cell_out
     procedure, public :: get_part_status
     procedure, public :: get_zday
     procedure, public :: get_sigma_local
     procedure, public :: get_znight
     procedure, public :: get_light_day

     procedure, public :: get_xpo
     procedure, public :: get_xk1
     procedure, public :: get_xk2
     procedure, public :: get_xk3
     procedure, public :: get_xk4

     procedure, public :: get_ypo
     procedure, public :: get_yk1
     procedure, public :: get_yk2
     procedure, public :: get_yk3
     procedure, public :: get_yk4

     procedure, public :: get_zpo
     procedure, public :: get_bdep
     procedure, public :: get_mode
     procedure, public :: get_mode_prop
     procedure, public :: get_zk1
     procedure, public :: get_zk2
     procedure, public :: get_zk3
     procedure, public :: get_zk4

     procedure, public :: get_ipo
     procedure, public :: get_jpo

     ! Setters
     procedure, public :: set_species
     procedure, public :: set_diapause

     procedure, public :: set_istart
     procedure, public :: set_isteps

     procedure, public :: set_dvm
     procedure, public :: set_cell_out
     procedure, public :: set_part_status
     procedure, public :: set_zday
     procedure, public :: set_sigma_local
     procedure, public :: set_znight
     procedure, public :: set_light_day

     procedure, public :: set_xpo
     procedure, public :: set_xk1
     procedure, public :: set_xk2
     procedure, public :: set_xk3
     procedure, public :: set_xk4

     procedure, public :: set_ypo
     procedure, public :: set_yk1
     procedure, public :: set_yk2
     procedure, public :: set_yk3
     procedure, public :: set_yk4

     procedure, public :: set_zpo
     procedure, public :: set_bdep
     procedure, public :: set_mode
     procedure, public :: set_mode_prop
     procedure, public :: set_zk1
     procedure, public :: set_zk2
     procedure, public :: set_zk3
     procedure, public :: set_zk4

     procedure, public :: set_ipo
     procedure, public :: set_jpo

     ! Methods
     procedure, public :: zday_funZ
     procedure, public :: znight_funZ
     !procedure, public :: sigma_funZ

  end type Zoo

  !--------------------------------------------------------------------72

contains

  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!72
  !! LIST procedures
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!72

  integer function n_part(this)
    class(partList) :: this 
    integer         :: n 

    call this%reset()

    n = 0
    do while( this%moreValues() )
       n = n + 1
       call this%next()
    enddo

    n_part = n

  end function n_part

  function currentPart(this)
    class(partList)   :: this
    class(*), pointer :: v
    type(Zoo)         :: currentPart

    v => this%currentValue()

    select type(v)
    type is (Zoo)
       currentPart = v
    end select

  end function currentPart

  subroutine thisPart(this, partPosition)
    class(partList) :: this
    integer         :: partPosition, i

    call this%reset()

    do i = 1, ( partPosition - 1 )
       call this%next()
    end do

  end subroutine thisPart

  subroutine addPart(this, value)
    class(partList)       :: this
    class(Zoo)            :: value ! user defined class of particle
    class(*), allocatable :: v

    allocate( v, source = value )

    call this%addValue(v)

  end subroutine addPart

  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!72
  !! Getters
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!72

  ! Particle's category
  
  integer function get_species(this)
    class(Zoo) :: this
    get_species = this%part_spec
  end function get_species

  ! In diapause or not

  integer function get_diapause(this)
    class(Zoo) :: this
    get_diapause = this%in_diapause
  end function get_diapause


  ! Time variables
  
  integer function get_istart(this)
    class(Zoo) :: this
    get_istart = this%istart
  end function get_istart

  integer function get_isteps(this)
    class(Zoo) :: this
    get_isteps = this%isteps
  end function get_isteps

  ! Swimming behaviour

  logical function get_dvm(this)
    class(Zoo) :: this
    get_dvm = this%dvm
  end function get_dvm

  logical function get_light_day(this)
    class(Zoo) :: this
    get_light_day = this%light_day
  end function get_light_day
  
  real function get_zday(this)
    class(Zoo) :: this
    get_zday = this%zday
  end function get_zday
  
  real function get_sigma_local(this)
    class(Zoo) :: this
    get_sigma_local = this%sigma_local
  end function get_sigma_local
  
  integer function get_part_status(this)
    class(Zoo) :: this
    get_part_status = this%part_status
  end function get_part_status

  real function get_znight(this)
    class(Zoo) :: this
    get_znight = this%znight
  end function get_znight

  ! Advection scheme
  
  logical function get_cell_out(this)
    class(Zoo) :: this
    get_cell_out = this%cell_out
  end function get_cell_out

  real function get_xpo(this)
    class(Zoo) :: this
    get_xpo = this%xpo
  end function get_xpo

  real function get_xk1(this)
    class(Zoo) :: this
    get_xk1 = this%xk1
  end function get_xk1

  real function get_xk2(this)
    class(Zoo) :: this
    get_xk2 = this%xk2
  end function get_xk2

  real function get_xk3(this)
    class(Zoo) :: this
    get_xk3 = this%xk3
  end function get_xk3

  real function get_xk4(this)
    class(Zoo) :: this
    get_xk4 = this%xk4
  end function get_xk4

  real function get_ypo(this)
    class(Zoo) :: this
    get_ypo = this%ypo
  end function get_ypo

  real function get_yk1(this)
    class(Zoo) :: this
    get_yk1 = this%yk1
  end function get_yk1

  real function get_yk2(this)
    class(Zoo) :: this
    get_yk2 = this%yk2
  end function get_yk2

  real function get_yk3(this)
    class(Zoo) :: this
    get_yk3 = this%yk3
  end function get_yk3

  real function get_yk4(this)
    class(Zoo) :: this
    get_yk4 = this%yk4
  end function get_yk4

  real function get_zpo(this)
    class(Zoo) :: this
    get_zpo = this%zpo
  end function get_zpo

  real function get_bdep(this)
    class(Zoo) :: this
    get_bdep = this%bdep
  end function get_bdep

  real function get_ipo(this)
    class(Zoo) :: this
    get_ipo = this%ipo
  end function get_ipo

  real function get_jpo(this)
    class(Zoo) :: this
    get_jpo = this%jpo
  end function get_jpo

  real function get_mode(this)
    class(Zoo) :: this
    get_mode = this%mode
  end function get_mode

  real function get_mode_prop(this)
    class(Zoo) :: this
    get_mode_prop = this%mode_prop
  end function get_mode_prop


  real function get_zk1(this)
    class(Zoo) :: this
    get_zk1 = this%zk1
  end function get_zk1

  real function get_zk2(this)
    class(Zoo) :: this
    get_zk2 = this%zk2
  end function get_zk2

  real function get_zk3(this)
    class(Zoo) :: this
    get_zk3 = this%zk3
  end function get_zk3

  real function get_zk4(this)
    class(Zoo) :: this
    get_zk4 = this%zk4
  end function get_zk4

  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!72
  !! Setters
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!72

  ! Particle's category
  
  integer function set_species(this)
    class(Zoo) :: this
    this%part_spec = set_species
  end function set_species

  ! Particle's category

  integer function set_diapause(this)
    class(Zoo) :: this
    this%in_diapause = set_diapause
  end function set_diapause

  ! Time variables
  
  subroutine set_istart(this, istart)
    class(Zoo) :: this
    integer    :: istart
    this%istart = istart
  end subroutine set_istart

  subroutine set_isteps(this, isteps)
    class(Zoo) :: this
    integer    :: isteps
    this%isteps = isteps
  end subroutine set_isteps

  ! Swimming behaviour

  subroutine set_dvm(this, dvm)
    class(Zoo) :: this
    logical    :: dvm
    this%dvm = dvm
  end subroutine set_dvm

  subroutine set_light_day(this, light_day)
    class(Zoo) :: this
    logical    :: light_day
    this%light_day = light_day
  end subroutine set_light_day
  
  subroutine set_zday(this, zday)
    class(Zoo) :: this
    real       :: zday
    this%zday = zday
  end subroutine set_zday
  
  subroutine set_sigma_local(this, sigma_local)
    class(Zoo) :: this
    real       :: sigma_local
    this%sigma_local = sigma_local
  end subroutine set_sigma_local

  subroutine set_part_status(this, part_status)
    class(Zoo) :: this
    integer       :: part_status
    this%part_status = part_status
  end subroutine set_part_status

  subroutine set_znight(this, znight)
    class(Zoo) :: this
    real       :: znight
    this%znight = znight
  end subroutine set_znight
 
  ! Advection scheme
  
  subroutine set_cell_out(this, cell_out)
    class(Zoo) :: this
    logical    :: cell_out
    this%cell_out = cell_out
  end subroutine set_cell_out

    subroutine set_xpo(this, x)
    class(Zoo) :: this
    real       :: x
    this%xpo = x
  end subroutine set_xpo

  subroutine set_xk1(this, xk1)
    class(Zoo) :: this
    real       :: xk1
    this%xk1 = xk1
  end subroutine set_xk1

  subroutine set_xk2(this, xk2)
    class(Zoo) :: this
    real       :: xk2
    this%xk2 = xk2
  end subroutine set_xk2

  subroutine set_xk3(this, xk3)
    class(Zoo) :: this
    real       :: xk3
    this%xk3 = xk3
  end subroutine set_xk3

  subroutine set_xk4(this, xk4)
    class(Zoo) :: this
    real       :: xk4
    this%xk4 = xk4
  end subroutine set_xk4

  subroutine set_ypo(this, y)
    class(Zoo) :: this
    real       :: y
    this%ypo = y
  end subroutine set_ypo

  subroutine set_yk1(this, yk1)
    class(Zoo) :: this
    real       :: yk1
    this%yk1 = yk1
  end subroutine set_yk1

  subroutine set_yk2(this, yk2)
    class(Zoo) :: this
    real       :: yk2
    this%yk2 = yk2
  end subroutine set_yk2

  subroutine set_yk3(this, yk3)
    class(Zoo) :: this
    real       :: yk3
    this%yk3 = yk3
  end subroutine set_yk3

  subroutine set_yk4(this, yk4)
    class(Zoo) :: this
    real       :: yk4
    this%yk4 = yk4
  end subroutine set_yk4

  subroutine set_zpo(this, z)
    class(Zoo) :: this
    real       :: z
    this%zpo = z
  end subroutine set_zpo

  subroutine set_bdep(this, bdep)
    class(Zoo) :: this
    real       :: bdep
    this%bdep = bdep
  end subroutine set_bdep

  subroutine set_ipo(this, i)
    class(Zoo) :: this
    integer    :: i
    this%ipo = i
  end subroutine set_ipo

  subroutine set_jpo(this, j)
    class(Zoo) :: this
    integer    :: j
    this%jpo = j
  end subroutine set_jpo

  subroutine set_mode(this, mode)
    class(Zoo) :: this
    integer    :: mode
    this%mode = mode
  end subroutine set_mode

  subroutine set_mode_prop(this, mode_prop)
    class(Zoo) :: this
    real       :: mode_prop
    this%mode_prop = mode_prop
  end subroutine set_mode_prop

  subroutine set_zk1(this, zk1)
    class(Zoo) :: this
    real       :: zk1
    this%zk1 = zk1
  end subroutine set_zk1

  subroutine set_zk2(this, zk2)
    class(Zoo) :: this
    real       :: zk2
    this%zk2 = zk2
  end subroutine set_zk2

  subroutine set_zk3(this, zk3)
    class(Zoo) :: this
    real       :: zk3
    this%zk3 = zk3
  end subroutine set_zk3

  subroutine set_zk4(this, zk4)
    class(Zoo) :: this
    real       :: zk4
    this%zk4 = zk4
  end subroutine set_zk4

  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!72
  !! ZOO class procedures
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!72

  ! Daytime preferred depth
  real function zday_funZ(this, surf) result(zval)
    class(Zoo)     :: this
    real, optional :: surf

    if ( present(surf) ) then
       zval = 7.53 * max(24.,surf) - 64.2 ! M. norvegica
    else
       zval = 10.
    endif
     
  end function zday_funZ

  ! Nighttime preferred depth
  real function znight_funZ(this) result(zval)
    class(Zoo) :: this
    integer    :: k
    real       :: r1

    call rng_uni(r1)
    k = max( 1, ceiling( r1 * 10. ) )
    zval = mid_depth(k)
     
  end function znight_funZ



  !--------------------------------------------------------------------72

end module partClass_mod
