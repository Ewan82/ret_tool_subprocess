!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!> \file core.f90
!
!> PROJECT : Sentinel Synergy Study
!
!> DESCRIPTION: Fortran90 implemenation of core routines ported from Sentinel Simulator (sense/core.py)
!
!> \authors The Inversion Lab (Michael Vossbeck) 
!
!> \date  February 2018
!                               
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!***********************************************************
!     fresnel0
!
!> @brief calculate the Nadir Fresnel reflectivity (e.g. Ulaby (2014), eq. 10.36)
!
!> @param[in]  e  complex relative dielectric permitivity
!
!> \return reflectivity
!
! Ref: sense/core.py, l4 ff.
!
real(kind=8) function fresnel0(e)
  implicit none
  complex(kind=8), intent(in) :: e
  complex(kind=8) :: hlp
  ! interfaces
  real(kind=8), external :: csq2

  ! fresnel0 = abs( (1._8 - sqrt(e))/(1._8+sqrt(e)) )**2
  hlp = (1._8 - sqrt(e))/(1._8+sqrt(e))
  fresnel0 = csq2(hlp)
end function fresnel0


!***********************************************************
!     calc_reflection_coefficients
!
!> @brief calculate reflection coefficients (Woodhouse, 2006; Eq. 5.54, 5.55)
!
!> @param[in]   theta  incidence angle [rad]
!> @param[in]   eps    relative dielectric permitivity (complex)
!> @param[out]  rho_v  reflection coefficent, vertical pol. (complex)
!> @param[out   rho_h  reflection coefficient, horizontal pol. (complex)
!
! Ref: sense/core.py, l47 ff.
!
subroutine calc_reflection_coefficients(theta, eps, rho_v, rho_h)
  implicit none
  ! arguments
  real(kind=8),intent(in) :: theta
  complex(kind=8), intent(in) :: eps
  complex(kind=8), intent(out) :: rho_v, rho_h
  ! local variables
  real(kind=8) :: co, si2
  complex(kind=8) :: hlp

  co = cos(theta)
  si2 = sin(theta)**2
  hlp = sqrt(eps-si2)

  rho_v = (eps*co-hlp)/(eps*co+hlp)
  rho_h = (co-hlp)/(co+hlp)

end subroutine calc_reflection_coefficients


!***********************************************************
!     reflectivity
!
!> @brief calculate reflectivity for H and V polarisation (table 2.5 Ulaby (2014), assumes specular surface)
!
!> @param[in]   theta  incidence angle [rad]
!> @param[in]   eps    relative dielectric permitivity (complex)
!> @param[out]  v      reflectivity ver. pol.
!> @param[out   h      reflectivity hor. pol.
!
! Ref: sense/core.py, l26 ff.
!
subroutine reflectivity(theta, eps, v, h)
  implicit none
  ! arguments
  real(kind=8),intent(in) :: theta
  complex(kind=8), intent(in) :: eps
  real(kind=8), intent(out) :: v, h
  ! interfaces
  real(kind=8), external :: csq2
  ! local vars
  complex(kind=8) :: rho_v, rho_h
  call calc_reflection_coefficients(theta, eps, rho_v, rho_h)

  ! v = abs(rho_v)**2
  ! h = abs(rho_h)**2
  ! v = real(rho_v)**2 + aimag(rho_v)**2
  ! h = real(rho_h)**2 + aimag(rho_h)**2
  v = csq2(rho_v)
  h = csq2(rho_h)
end subroutine reflectivity
