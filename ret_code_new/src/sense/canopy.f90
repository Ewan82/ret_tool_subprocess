!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!> \file canopy.f90
!
!> PROJECT : Sentinel Synergy Study
!
!> DESCRIPTION: Fortran90 implemenation of canopy class ported from Sentinel Simulator (sense/canopy.py)
!
!> \authors The Inversion Lab (Michael Vossbeck) 
!
!> \date  February 2018
!                               
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!***********************************************************
!     canopy_extinction_coeffs
!
!> @brief compute extinction coefficients from LAI (microwave) and single-scattering albedo
!
!> @param[in]  omega   single scattering albedo
!> @param[in]  lai_mw  Leaf-Area Index (MW-compliant)
!> @param[out] ke_h
!> @param[out] ke_v
!> @param[out] ks_h
!> @param[out] ks_v
!
! Ref: sense/canopy.py, l21 ff.
!      simulator.py, l100 ff.
!
subroutine canopy_extinction_coeffs( omega, lai_mw, ke_v, ke_h, ks_v, ks_h )
  implicit none
  ! arguments
  real(kind=8), intent(in) :: omega, lai_mw
  real(kind=8), intent(out) :: ke_h, ke_v, ks_h, ks_v

  ke_h = lai_mw
  ke_v = lai_mw
  ks_h = omega*lai_mw
  ks_v = omega*lai_mw
end subroutine canopy_extinction_coeffs


!***********************************************************
!     canopy_transmissivity
!
!> @brief compute canopy transmissivity
!
!> @param[in] ke_h  volume extinction coefficient (hor. pol) [Np/m]
!> @param[in] ke_v  volume extinction coefficient (ver. pol) [Np/m]
!> @param[in] d     height of canopy layer
!> @param[in] theta incidence angle
!> @param[out] t_h  transmissivity, h-polarised
!> @param[out] t_v  transmissivity, v-polarised
!
! Ref: sense/canopy.py, l299 ff.
!      sense/canopy.py, l356 ff
!
subroutine canopy_transmissivity( ke_v, ke_h, d, theta, t_v, t_h )
  implicit none
  ! arguments
  real(kind=8), intent(in) :: ke_h, ke_v, d, theta
  real(kind=8), intent(out) :: t_h, t_v
  ! local decls
  real(kind=8) :: tau_h, tau_v

  tau_h = tau_fct(ke_h)
  tau_v = tau_fct(ke_v)

  t_h = exp(-tau_h)
  t_v = exp(-tau_v)

contains
  real(kind=8) function tau_fct(k) !assumption: extinction is isotropic
    implicit none
    real(kind=8), intent(in) :: k
    tau_fct = k*d/cos(theta)
  end function tau_fct

end subroutine canopy_transmissivity


!***********************************************************
!     canopy_transmissivity
!
!> @brief calculate canopy volume contribution only (Eq. 11.10 + 11.16 as seen in 11.17, Ulaby (2014))
!
! Ref: sense/model.py, l372 ff.
!
subroutine canopy_sigma_c(sigma_vol_back, theta, ke_v, ke_h, t_v, t_h, s0c)
  use mo_sensimul_s1, only:sense_fv
  implicit none
  ! arguments
  real(kind=8), intent(in) :: sigma_vol_back(3) !hh,vv,hv
  real(kind=8), intent(in) :: theta, ke_v, ke_h, t_v, t_h
  real(kind=8), intent(out) :: s0c(3)

  !-- HH
  s0c(1) = (1._8 - t_h*t_h)*(sigma_vol_back(1)*cos(theta))/(ke_h + ke_h)

  !-- VV
  s0c(2) = (1._8 - t_v*t_v)*(sigma_vol_back(2)*cos(theta))/(ke_v + ke_v)
  
  if( sigma_vol_back(3).eq.sense_fv ) then
     s0c(3) = sense_fv
  else
     s0c(3) = (1._8 - t_h*t_v)*(sigma_vol_back(3)*cos(theta))/(ke_h + ke_v)
  endif

end subroutine canopy_sigma_c
