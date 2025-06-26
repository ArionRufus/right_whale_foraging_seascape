
subroutine zeros

  use partSubClass_mod
  implicit none

  !---------------------------------------------------------------------
  ! 1-D arrays

  dd        = 0.

  dz        = 0.

  up_depth  = 0.

  mid_depth = 0.

  !---------------------------------------------------------------------
  ! 2-D arrays

  dlx       = 0.

  dly       = 0.

  nlayer    = 0

  !---------------------------------------------------------------------
  ! 3-D arrays

  salt      = 0.0

  temp      = 0.0

  u         = 0.0

  v         = 0.0

  u1        = 0.0

  v1        = 0.0

  u2        = 0.0

  v2        = 0.0

#ifdef VERT
  w         = 0.0
  
  w1        = 0.0
  
  w2        = 0.0
#endif

#ifdef YUMY
  phyto     = 0.0
  
  phyto3d   = 0.0

  diatom    = 0.0
  
  flag      = 0.0
  
  microzoo  = 0.0
  
  food      = 0.0
  
  parsurf   = 0.0
#endif

  !---------------------------------------------------------------------
  nc_repeat = .false.
  

  !---------------------------------------------------------------------
  ! Particules data

  xpo       = 0.

  ypo       = 0.

  zpo       = 0.

  xk1       = 0.

  yk1       = 0.

  zk1       = 0.

  xk2       = 0.

  yk2       = 0.

  zk2       = 0.

  xk3       = 0.

  yk3       = 0.

  zk3       = 0.

  xk4       = 0.

  yk4       = 0.

  zk4       = 0.

  mode      = 1 !at initialisation, all the modes are surface modes

  mode_prop = 1

  dvm       = .true.

  !---------------------------------------------------------------------
  ! Default particules parameter: see part_nml namelist

  part_spec     = 0

  part_fac      = 1

  output_freq   = 1

  part_nc_out   = .false.

  part_traj_out = .false.

  mature        = .false.

  partfile      = './data/part_domain_unif.dat'

  outfile       = '../NETCDF/output'

  !---------------------------------------------------------------------
  ! Default phyisiological parameters: see physio_nml namelist

  !---------------------------------------------------------------------

  return

end subroutine zeros

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine init_random_seed
  
  use iso_fortran_env, only: int32, int64, real32, real64
  use rng_par_zig
  
#ifdef PAR
  use omp_lib
#endif

  implicit none

  integer :: i, nt

  integer(int64) :: seed(2)

  ! Just some initial entropy for a reproducible random sequence
  integer(int64), parameter :: base(2) = [int(Z'1B2F65A43216D531', int64), &
                                          int(Z'264665753F152E1C', int64)]

  !---------------------------------------------------------------------

  nt   = 1
  seed = base

#ifdef PAR
  !$OMP PARALLEL
  !nt = omp_get_num_threads()
  !$OMP END PARALLEL

  ! Jump nthread-times for each image so that each image
  ! has a disjunct set of random sequences for its threads.
  !call rng_jump(seed, nt*(myim-1))
#endif

  call rng_init(nt, seed)

end subroutine init_random_seed

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
