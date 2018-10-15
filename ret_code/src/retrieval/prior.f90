subroutine residual_prior (n, x, priordiff)
  implicit none
  ! arguments
  integer, intent(in) :: n
  real(kind=8), intent(in) :: x(n)
  real(kind=8), intent(out) :: priordiff(n)
  ! local
  real(kind=8) :: sx(n), x0(n) 
  ! externals
  external getprior

  !-- get prior
  call getprior(n, x0, sx)

  priordiff = (x-x0/sx)

end subroutine residual_prior
