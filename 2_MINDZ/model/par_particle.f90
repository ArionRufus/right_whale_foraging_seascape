subroutine par_particle

  use partSubClass_mod

  implicit none

  !--------------------------------------------------------------------72
  ! Parameters values for each species 

  sex       = sex_ini(part_spec) 
  activity  = activity_ini(part_spec) 
  mean_mass = mean_mass_ini(part_spec)
  min_mass  = min_mass_ini(part_spec)
  max_mass  = max_mass_ini(part_spec)
  a_molt    = a_molt_ini(part_spec)
  b_molt    = b_molt_ini(part_spec)
  k0        = k0_ini(part_spec)
  h0        = h0_ini(part_spec)
  A         = A_ini(part_spec)
  r0        = r0_ini(part_spec)
  rb        = rb_ini(part_spec)
  er        = er_ini(part_spec)
  t_lim     = t_lim_ini(part_spec)

#ifdef DEBUG
  write(*,'(/a,/a,i12,13(/a,d12.6))') ' List of particles parameters :', &
                                      ' sex         = ', sex,            &
                                      ' mean_mass   = ', mean_mass,      &
                                      ' min_mass    = ', min_mass,       &
                                      ' max_mass    = ', max_mass,       &
                                      ' a_molt      = ', a_molt,         &
                                      ' b_molt      = ', b_molt,         &
                                      ' k0          = ', k0,             &
                                      ' h0          = ', h0,             &
                                      ' A           = ', A,              &
                                      ' r0          = ', r0,             &
                                      ' rb          = ', rb,             &
                                      ' er          = ', er,             &
                                      ' t_lim       = ', t_lim
#endif

end subroutine par_particle
