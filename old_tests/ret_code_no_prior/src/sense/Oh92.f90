!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!> \file Oh92.f90
!
!> PROJECT : Sentinel Synergy Study
!
!> DESCRIPTION: Fortran90 implemenation of ...
!
!> \authors The Inversion Lab (Michael Vossbeck) 
!
!> \date  February 2018
!                               
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine Oh92_backscatter(eps, ks, theta, backscatter )
  use mo_sensimul_s1, only:pi
  implicit none
  ! arguments
  complex(kind=8), intent(in) :: eps
  real(kind=8), intent(in) :: ks, theta
  real(kind=8), intent(out) :: backscatter(3) !hh,vv,hv
  ! local declarations
  real(kind=8) :: v,h
  real(kind=8) :: frsn0
  real(kind=8) :: p, q, vv0
  real(kind=8), external :: fresnel0
  external reflectivity

  !-- Nadir Fresnel reflectivity
  frsn0 = fresnel0(eps)

  !--
  call reflectivity(theta, eps, v, h)

  !
  p = calc_p()
  q = calc_q()
  vv0 = calc_vv()

  backscatter(1) = p*vv0
  backscatter(2) = vv0
  backscatter(3) = q*vv0
  
contains
  real(kind=8) function calc_p()
    implicit none
    real(kind=8) :: a

    a = 1._8/(3._8*frsn0)
    calc_p = (1._8 - (2._8*theta/pi)**a * exp(-ks))**2
  end function calc_p

  real(kind=8) function calc_q()
    implicit none
    calc_q = 0.23_8*frsn0**0.5_8 * (1._8 -exp(-ks))
  end function calc_q

  real(kind=8) function calc_vv()
    implicit none
    real(kind=8) :: a,b

    a = 0.7_8*(1._8 - exp(-0.65_8*ks**1.8))
    b = cos(theta)**3._8 * (v + h) / sqrt(p)
    calc_vv = a*b
  end function calc_vv
end subroutine Oh92_backscatter
