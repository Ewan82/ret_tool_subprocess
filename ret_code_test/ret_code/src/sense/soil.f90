!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!> \file soil.f90
!
!> PROJECT : Sentinel Synergy Study
!
!> DESCRIPTION: Fortran90 implemenation ported from Sentinel Simulator (sense/soil.py)
!
!> \authors The Inversion Lab (Michael Vossbeck) 
!
!> \date  February 2018
!                               
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!***********************************************************
!     convert_mv2eps
!
!> @brief calculate dielectric permittivity from (volumetric) soil moisture using Dobson85
!
!> @param[in]  mv     volumetric soil moisture [m**3/m**3]
!> @param[in]  bulk   bulk density [g/cm**3]
!> @param[in]  alpha  as computed by dobson85_init
!> @param[in]  beta1  as computed by dobson85_init
!> @param[in]  beta2  as computed by dobson85_init
!> @param[in]  ew     dielectric permittivity of free water
!> @param[out] eps    relative dielectric permittivity
!
! Ref: sense/soil.py, l61 ff. (but only Dobson85 is supported)
!
subroutine convert_mv2eps(mv, bulk, alpha, beta1, beta2, ew, eps)
  implicit none
  ! arguments
  real(kind=8), intent(in) :: mv, bulk, alpha, beta1, beta2
  complex(kind=8), intent(in) :: ew
  complex(kind=8), intent(out) :: eps
  ! local variables
  external dobson85_eps

  call dobson85_eps(mv, bulk, alpha, beta1, beta2, ew, eps)

end subroutine convert_mv2eps


!***********************************************************
!     soil_ks
!
!> @brief calculate surface roughness from surface height and frequency
!
!> @param[in]  s   surface rms height [m]
!> @param[in]  f   frequency [GHz]
!> @param[out] ks  roughness parameter
!
! Ref: sense/soil.py, l51 ff.
!
subroutine soil_ks(s, f, ks)
  implicit none
  ! arguments
  real(kind=8), intent(in) :: s, f
  real(kind=8), intent(out) :: ks
  ! local variables
  real(kind=8), external :: f2lam
  real(kind=8), parameter :: pi = 3.1415926535879
  real(kind=8) :: k

  k = (2._8*pi)/f2lam(f)
  ks = s*k
end subroutine soil_ks
