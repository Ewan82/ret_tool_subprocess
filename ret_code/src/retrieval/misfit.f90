!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!> \file misfit.f90
!
!> PROJECT : Sentinel Synergy Study
!
!> DESCRIPTION: computation of misfit between observations and simulated equivalents
!
!> \authors The Inversion Lab (Michael Vossbeck, Thomas Kaminski) 
!
!> \date  April 2018
!                               
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!***********************************************************
!     misfit
!
!> @brief determines the misfit vector between observations and simulated
!         observation equivalents
!
!> @details TBD
!
!> @param[in]  n       overall (1D) length of control vector
!> @param[in]  x       complete (normalised) control vector
!> @param[in]  m       overall (1D) dimension of misfit vector
!                      (must equal dimension of simulation vector)
!> @param[out] obsdiff misfit vector
!
subroutine misfit(n, x, m, obsdiff)
  use mo_sensimul, only:get_m_s1,get_m_s2,sim_fill_value
  implicit none
  ! arguments
  integer, intent(in) :: n, m
  real(kind=8), intent(in)  :: x(n)
  real(kind=8), intent(out) :: obsdiff(m)
  ! local
  real(kind=8) :: y(m), yobs(m), syobs(m)
  logical :: succeed
  integer :: nobs_s1,nobs_s2
  integer :: m_s1,m_s2
  integer :: i
  ! externals
  external simulate_s1s2,getobs


  obsdiff(1:m) = 0._8

  !-- read obs
  call getobs(m, yobs, syobs, nobs_s1, nobs_s2, succeed)
  if( .not.succeed ) then
     write(*, '(a)') ' FATAL::misfit:getobs failed, cannot continue.'
     stop
  endif

  !-- run S1+S2 simulation
  call simulate_s1s2(n, x, m, y)

  m_s1 = get_m_s1()
  m_s2 = get_m_s2()
  if( nobs_s1.gt.0 ) then
     do i=1,m_s1
        if( yobs(i).ne.sim_fill_value.and.y(i).ne.sim_fill_value ) then
           obsdiff(i) = (y(i)-yobs(i))/syobs(i)
        endif
     enddo
  endif

  if( nobs_s2.gt.0 )  then
     do i=m_s1+1,m_s1+m_s2
        if( yobs(i).ne.sim_fill_value.and.y(i).ne.sim_fill_value ) then
           obsdiff(i) = (y(i)-yobs(i))/syobs(i)
        endif
     enddo
  endif

end subroutine misfit

