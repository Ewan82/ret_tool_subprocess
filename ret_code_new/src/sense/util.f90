!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!> \file util.f90
!
!> PROJECT : Sentinel Synergy Study
!
!> DESCRIPTION: Fortran90 implemenation of utilities from Sentinel Simulator (sense/dielectric/util.py)
!
!> \authors The Inversion Lab (Michael Vossbeck) 
!
!> \date  February 2018
!                               
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!***********************************************************
!     f2lam
!
!> @brief for given frequency in GHz compute wavelength [m]
!
!> @param[in]  f   frequencey [GHz]
!> \return wavelengh [m]
!
real(kind=8) function f2lam(f)
  implicit none
  ! arguments
  real(kind=8), intent(in) :: f
  real(kind=8), parameter :: c0 = 299792458._8 !-- speed of light [m/s]

  f2lam = c0/(f*1.e9_8)
end function f2lam


!***********************************************************
!     csq2
!
!> @brief computes square of absolute value of complex number 
!
!> @param[in]  z  complex number
!> \return     |z|^2
!
real(kind=8) function csq2(z)
  implicit none
  complex(kind=8), intent(in) :: z

  !  csq2 = abs(z)**2
  csq2 = real(z)**2 + aimag(z)**2
end function csq2
