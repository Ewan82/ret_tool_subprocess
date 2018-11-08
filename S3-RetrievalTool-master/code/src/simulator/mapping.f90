!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> \file mapping.f90
!> \brief mapping of 1D control vector to physical variables
!> \authors The Inversion Lab (Michael Vossbeck, Thomas Kaminski) 
!> \date  April 2018
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine x2p(n,x,xp)
  implicit none
  ! arguments
  integer, intent(in) :: n
  real(kind=8), intent(in) :: x(n)
  real(kind=8), intent(out) :: xp(n)
  ! local decls
  real(kind=8) :: sx(n),pr(n)

  call getprior(n, pr, sx)

  xp = x*sx

end subroutine x2p

subroutine control_vector_split(n, x, n_s1, x_s1, n_s2, x_s2)
  use mo_sensimul
  implicit none
  ! arguments
  integer, intent(in) :: n, n_s1, n_s2
  real(kind=8), intent(in) :: x(n)
  real(kind=8), intent(out) :: x_s1(n_s1), x_s2(n_s2)
  ! local decls
  integer :: nexpect
  integer:: i,j,i_s1,i_s2

  !-- dimensional consistency
  nexpect = nparam_s1+npts*nsc
  if( n.ne.nexpect ) then
     write(*, '(a,i4,1x,a,i4)') ' FATAL::control_vector_split:n=',&
          n,'does not equal expected size=',nexpect
     stop
  else if( n_s1.ne.nparam_s1+npts_s1*nsc ) then
     write(*, '(a,i4,1x,a,i4)') ' FATAL::control_vector_split:n_s1=',&
          n_s1,'does not equal expected size=',nparam_s1+npts_s1*nsc
     stop
  else if( n_s2.ne.npts_s2*nsc ) then
     write(*, '(a,i4,1x,a,i4)') ' FATAL::control_vector_split:n_s2=',&
          n_s2,'does not equal expected size=',npts_s2*nsc
     stop
  endif

  !-- map S1 parameter
  x_s1(1:nparam_s1) = x(1:nparam_s1)

  j    = nparam_s1+1 ! position in full control vector
  i_s1 = nparam_s1+1 ! position in S1 restricted control vector
  i_s2 = 1           ! position in S2 restrictred control vector
  do i=1,npts
     !-- point with S1 simulation
     if( idx_is_s1(i) ) then
        x_s1(i_s1:i_s1+nsc-1) = x(j:j+nsc-1)
        i_s1 = i_s1 + nsc
     endif
     !-- point with S2 simulation
     if( idx_is_s2(i) ) then
        x_s2(i_s2:i_s2+nsc-1) = x(j:j+nsc-1)
        i_s2 = i_s2 + nsc
     endif
     !-- increment position in control vector
     j = j + nsc
  enddo

end subroutine control_vector_split
