!***********************************************************
!> \file  cost.f90
!> \brief evaluates total cost function
!> \authors MV/TK, The Inversion Lab
!> \date  January 2018
!***********************************************************
!***********************************************************
! SUBROUTINE COST()
!> \brief evaluates total cost function that combines 
!>        terms for observations, model, and prior
!> \authors The Inversion Lab
!> \date  January 2018
subroutine cost(n, x, m, f)
  use mo_retrieval, only:retr_use_prior_term, retr_use_state_term
  implicit none
  ! arguments
  integer, intent(in) :: n         !< length of control
  integer, intent(in) :: m         !< length of output vector
  real(kind=8), intent(in) :: x(n) !< control vector
  real(kind=8), intent(out) :: f   !< function value
  ! local
  real(kind=8) :: obsdiff(m), costobs
  real(kind=8) :: priordiff(n), costprior
  real(kind=8) :: modeldiff(n), costmodel
  real(kind=8) :: xphys(n) 
  logical :: ldebug
  ! externals
  external x2p, state_model, residual_prior

  !-- set flags
  ldebug = .true.

  !-- control vector in physical units
  call x2p(n, x, xphys)

  !-- compute misfit
  call misfit (n, xphys, m, obsdiff)
  costobs = 0.5_8 * (sum(obsdiff**2))
  if( ldebug ) then
     write(*,'(a,e25.16)') ' DIAG::cost:cost_obs=  ',costobs
  endif

  !-- compute state model differences
  if( retr_use_state_term ) then
     call H_m(n, xphys, modeldiff)
     costmodel = 0.5_8 * (sum(modeldiff**2))
  else
     costmodel = 0._8
  endif
  if( ldebug ) then
     write(*,'(a,e25.16)') ' DIAG::cost:cost_model=',costmodel
  endif

  !-- compute prior diff
  if( retr_use_prior_term ) then
     call residual_prior (n,x,priordiff)
     costprior = 0.5_8 * (sum(priordiff**2))
  else
     costprior = 0._8
  endif
  if( ldebug ) then
     write(*,'(a,e25.16)') ' DIAG::cost:cost_prior=',costprior
  endif

  !-- sum-up total cost value
  f = costobs + costmodel + costprior
end subroutine cost
