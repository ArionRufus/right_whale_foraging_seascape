subroutine ll_to_xy_NEMO(rlon,rlat,x,y)

  USE OPA_GET_XY_LL

  real rlon,rlat,x,y

  lon00 = rlon
  lat00 = rlat

  call opa_proj

  x     = ri00
  y     = rj00

  return

end subroutine ll_to_xy_NEMO

!-----------------------------------------------------------------------

subroutine xy_to_ll_NEMO(x,y,rlon,rlat)

  USE OPA_GET_XY_LL

  real rlon,rlat,x,y

  ri00 = x
  rj00 = y

  call opa_proj_inv

  rlon = lon00
  rlat = lat00

  return

end subroutine xy_to_ll_NEMO
