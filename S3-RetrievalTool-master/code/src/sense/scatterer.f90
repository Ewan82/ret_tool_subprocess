!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!> \file scatterer.f90
!
!> PROJECT : Sentinel Synergy Study
!
!> DESCRIPTION: Fortran90 implemenation of scatterer class ported from Sentinel Simulator (sense/scatterer.py)
!
!> \authors The Inversion Lab (Michael Vossbeck) 
!
!> \date  February 2018
!                               
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!***********************************************************
!     ScatRayleigh_sigma_v_back
!
!> @brief compute backscattering coefficient from volume extinction coefficient(s) for Rayleigh scatterer
!
!> @param[in]  sigma_s_hh   volume extinction coefficient
!> @param[in]  sigma_s_vv   volume extinction coefficient
!> @param[out] backscatter  volume backscattering coefficient
!
! Ref: sense/scatterer.py, l56 ff.
!
subroutine ScatRayleigh_sigma_v_back(sigma_s_vv, sigma_s_hh, backscatter)
  use mo_sensimul_s1, only:sense_fv
  implicit none
  real(kind=8), intent(in) :: sigma_s_vv, sigma_s_hh
  real(kind=8), intent(out) :: backscatter(3) !'hh','vv','hv'

  backscatter(1) = 1.5_8*sigma_s_hh
  backscatter(2) = 1.5_8*sigma_s_vv
  backscatter(3) = sense_fv
end subroutine ScatRayleigh_sigma_v_back


!***********************************************************
!     ScatRayleigh_sigma_v_bistatic
!
!>  @brief compute backscattering coefficient from volume extinction coefficient(s) for Rayleigh scatterer (Ulceby (Eq. 11.22))
!
!> @param[in]  sigma_s_hh   volume extinction coefficient
!> @param[in]  sigma_s_vv   volume extinction coefficient
!> @param[out] backscatter  volume backscattering coefficient
!
! Ref: sense/scatterer.py, l60 ff.
!
subroutine ScatRayleigh_sigma_v_bistatic(sigma_s_vv, sigma_s_hh, backscatter)
  implicit none
  real(kind=8), intent(in) :: sigma_s_vv, sigma_s_hh
  real(kind=8), intent(out) :: backscatter(3) !'hh','vv','hv'

  call ScatRayleigh_sigma_v_back(sigma_s_vv, sigma_s_hh, backscatter)
end subroutine ScatRayleigh_sigma_v_bistatic
