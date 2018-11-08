!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!> \file dobson85.f90
!
!> PROJECT : Sentinel Synergy Study
!
!> DESCRIPTION: Fortran90 implemenation of Dobson85 dielectric model ported from Sentinel Simulator (sense/dielectric/dobson85.py)
!
!> \authors The Inversion Lab (Michael Vossbeck) 
!
!> \date  February 2018
!                               
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!***********************************************************
!     dobson85_init
!
!> @brief initialisation of Dobson85 model parameters (eq. 4.68, Ulaby (2014))
!
!> @param[in]  sand  sand content as fractional volume
!> @param[in]  clay  clay content as fractional volume
!> @param[in]  bulk  bulk density [g/cm**3]
!> @param[out] alpha
!> @param[out] beta1
!> @param[out] beta2
!> @param[out] sigma
!
subroutine dobson85_init(sand, clay, bulk, alpha, beta1, beta2, sigma)
  implicit none
  ! arguments
  real(kind=8), intent(in) :: sand, clay, bulk
  real(kind=8), intent(out) :: alpha, beta1, beta2, sigma

  alpha = 0.65_8
  beta1 = 1.27_8 - 0.519_8*sand - 0.152_8*clay
  beta2 = 2.06_8 - 0.928_8*sand - 0.255_8*clay
  sigma = -1.645_8 + 1.939_8*bulk - 2.256_8*sand + 1.594_8*clay
end subroutine dobson85_init


!***********************************************************
!     dobson85_ew
!
!> @brief compute dielectric permittivity of free water using a simplistic approach (eq. 4.69, Ulaby (2014))
!
!> @param[in]  f      frequency
!> @param[in]  sigma  as computed by dobson85_init
!> @param[out] ew     dielectric permittivity of free water
!
subroutine dobson85_ew(f, sigma, ew)
  implicit none
  ! arguments
  real(kind=8), intent(in) :: f
  real(kind=8), intent(in) :: sigma
  complex(kind=8), intent(out) :: ew
  ! local vars
  real(kind=8), parameter :: f0 = 18.64_8 ! relaxation frequency [GHz]
  real(kind=8) :: hlp, e1, e2

  hlp = f/f0
  e1 = 4.9_8 + 74.1_8/(1._8 + hlp**2)
  e2 = (74.1_8*hlp)/(1._8 + hlp**2) + 6.46_8*sigma/f
  ew = cmplx(e1,e2,kind=8)

end subroutine dobson85_ew


!***********************************************************
!     dobson85_eps
!
!> @brief calculate dielectric permittivity (Eq. 4.66 (Ulaby et al., 2014))
!
!> @param[in]  mv     volumetric soil moisture [m**3/m**3]
!> @param[in]  bulk   bulk density [g/cm**3]
!> @param[in]  alpha  as computed by dobson85_init
!> @param[in]  beta1  as computed by dobson85_init
!> @param[in]  beta2  as computed by dobson85_init
!> @param[in]  ew     dielectric permittivity of free water
!> @param[out] eps    relative dielectric permittivity
!
subroutine dobson85_eps(mv, bulk, alpha, beta1, beta2, ew, eps)
  implicit none
  ! arguments
  real(kind=8), intent(in) :: mv, bulk, alpha, beta1, beta2
  complex(kind=8), intent(in) :: ew
  complex(kind=8), intent(out) :: eps
  ! local variables
  real(kind=8) :: e1,e2

  e1 = (1._8 + 0.66_8*bulk + mv**beta1*real(ew)**alpha - mv)**(1._8/alpha)
  e2 = aimag(ew)*mv**beta2

  eps = cmplx(e1,e2,kind=8)
end subroutine dobson85_eps
