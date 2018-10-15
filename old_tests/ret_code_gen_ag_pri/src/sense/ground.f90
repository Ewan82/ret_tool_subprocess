!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!> \file ground.f90
!
!> PROJECT : Sentinel Synergy Study
!
!> DESCRIPTION: Fortran90 implemenation of Ground class ported from Sentinel Simulator (sense/model.py)
!
!> \authors The Inversion Lab (Michael Vossbeck) 
!
!> \date  February 2018
!                               
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!***********************************************************
!     calc_rho
!
!> @brief calculate coherent p-polarized reflectivity  (Ref: Eq. 11.11 (Ulaby, 2014))
!
!---------
! Note that the specular reflectivity is corrected by a roughness term
!         if ks>0.2
!
!         however, a sensitivity analysis showed that even for ks==0.2
!         deviations can be up to 15% for typical incidence angles
!         Only in case that ks << 0.1, the correction can be neglected.
!         We therefore always use the roughness correction factor!
!         TODO: unclear so far how this relates to surface (soil) scattering models
!---------
!
!> @param[in]   theta  incidence angle [rad]
!> @param[in]   ks     (surface) roughness parameter
!> @param[in]   v      reflectivity ver. pol.
!> @param[in]   h      reflectivity hor. pol.
!> @param[out]  rho_v  coherent p-polarized reflectivity ver. pol.
!> @param[out   rho_h  coherent p-polarized reflectivity hor. pol.
!
! Ref: sense/model.py, l177 ff.
!
subroutine calc_rho(theta, ks, v, h, rho_v, rho_h)
  implicit none
  ! arguments
  real(kind=8), intent(in) :: theta
  real(kind=8), intent(in) :: ks
  real(kind=8), intent(in) :: v,h
  real(kind=8), intent(out) :: rho_v, rho_h

  !--
  rho_v = v * exp(-4._8 * cos(theta)**2 * ks**2)
  rho_h = h * exp(-4._8 * cos(theta)**2 * ks**2)
end subroutine calc_rho


!***********************************************************
!     sigma_gcg
!
!> @brief compute ground-canopy-ground interaction
!
!> @param[in]  sigma_vol_back
!> @param[in]  theta
!> @param[in]  rho_v
!> @param[in]  rho_h
!> @param[in]  t_h
!> @param[in]  t_v
!> @param[in]  ke_v
!> @param[in]  ke_h
!> @param[out] s0gcg
!
! Ref: sense/model.py, l205 ff.
!
subroutine sigma_gcg(sigma_vol_back, theta, rho_v, rho_h, t_v, t_h, ke_v, ke_h, s0gcg)
  use mo_sensimul_s1, only:sense_fv
  implicit none
  ! arguments
  real(kind=8), intent(in) :: sigma_vol_back(3) !hh,vv,hv
  real(kind=8), intent(in) :: theta, rho_v, rho_h, t_v, t_h, ke_v, ke_h
  real(kind=8), intent(out) :: s0gcg(3)         !hh,vv,hv

  !-- HH
  s0gcg(1) = sigma_vol_back(2)*cos(theta)*rho_h*rho_h*(t_h**2-t_h**4)/(ke_h + ke_h)

  !-- VV
  s0gcg(2) = sigma_vol_back(1)*cos(theta)*rho_v*rho_v*(t_v**2-t_v**4)/(ke_v + ke_v)

  !-- HV
  if( sigma_vol_back(3).eq.sense_fv ) then
     s0gcg(3) = sense_fv
  else
     s0gcg(3) = sigma_vol_back(3)*cos(theta)*rho_h*rho_v*(t_h*t_v-t_h**2*t_v**2)/(ke_h + ke_v)
  endif
end subroutine sigma_gcg


!***********************************************************
!     sigma_cg
!
!> @brief calculate canopy ground scattering coefficient This is based on Eq. 11.17 (last term) in Ulaby (2014) and 11.14 in Ulaby (2014)
!  (for co-pol, coherent addition can be made as an option)
!
!> @param[in]  sigma_vol_bistatic
!> @param[in]  d
!> @param[in]  coherent
!> @param[in]  rho_v
!> @param[in]  rho_h
!> @param[in]  t_h
!> @param[in]  t_v
!> @param[out] s0cg
!
! Ref: sense/model.py, l215 ff.
!
subroutine sigma_cg(sigma_vol_bistatic, d, coherent, rho_v, rho_h, t_v, t_h, s0cg)
  use mo_sensimul_s1, only:sense_fv
  implicit none
  ! arguments
  real(kind=8), intent(in) :: sigma_vol_bistatic(3) !hh,vv,hv
  real(kind=8), intent(in) :: d, rho_v, rho_h, t_v, t_h
  logical, intent(in) :: coherent
  real(kind=8), intent(out) :: s0cg(3)              !hh,vv,hv
  ! local decls
  real(kind=8) :: n

  if( coherent ) then
     n = 2._8
  else
     n = 1._8
  endif

  !-- HH
  s0cg(1) = n*sigma_vol_bistatic(1)*d*(rho_h + rho_h)*t_h*t_h

  !-- VV
  s0cg(2) = n*sigma_vol_bistatic(2)*d*(rho_v + rho_v)*t_v*t_v

  !-- HV
  if( sigma_vol_bistatic(3).eq.sense_fv ) then
     s0cg(3) = sense_fv
  else
     s0cg(3) = n*sigma_vol_bistatic(3)*d*(rho_v + rho_h)*t_h*t_v
  endif

end subroutine sigma_cg


!***********************************************************
!     ground_sigma
!
!> @brief calculate backscattering coefficient (Eq. 11.4, p.463 Ulaby (2014))
!
!> @param[in]  soil_backscatter
!> @param[in]  t_v
!> @param[in]  t_h
!> @param[out] gnd_sigma
!
! Ref: sense/model.py, l243 ff.
!
subroutine ground_sigma(soil_backscatter, t_v, t_h, gnd_sigma)
  use mo_sensimul_s1, only:sense_fv
  implicit none
  ! arguments
  real(kind=8), intent(in) :: soil_backscatter(3) !hh,vv,hv
  real(kind=8), intent(in) :: t_v, t_h
  real(kind=8), intent(out) :: gnd_sigma(3)

  !-- HH
  gnd_sigma(1) = soil_backscatter(1)*t_h*t_h

  !-- VV
  gnd_sigma(2) = soil_backscatter(2)*t_v*t_v

  !-- HV
  if( soil_backscatter(3).eq.sense_fv ) then
     gnd_sigma(3) = sense_fv
  else
     gnd_sigma(3) = soil_backscatter(3)*t_v*t_h
  endif

end subroutine ground_sigma
