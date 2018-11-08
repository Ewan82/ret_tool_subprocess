subroutine ua (m, n, x, c, debug, diag, ierr)
  implicit none
  ! externals
  external cost_jacobian
  ! arguments
  integer, intent(in)    :: n, m           ! jacobian dimensions
  integer, intent(inout) :: ierr           ! error flag
  real(kind=8), intent(in)       :: x(n)   ! optimum
  real(kind=8), intent(out)      :: c(n,n) ! covar of uncertainty
  logical, intent(in) :: debug, diag
  ! local  
  integer :: iflag, i, j, l, k
  real(kind=8)    :: fvec(m), fjac(m,n), w(n), u(m,n), v(n,n)
  logical :: matu = .true., matv = .true.

  iflag = 2
  call cost_jacobian ( m, n, x, fvec, fjac, m, iflag )

  ierr = 0
  call svd ( m, n, fjac, w, matu, u, matv, v, ierr )
  if (debug.or.diag) write ( *, '(a)' ) ' '
  if (debug.or.diag) write ( *, '(a,i10)' ) ' INFO::ua:SVD done with ierr = ', ierr
  if (debug.or.diag) write ( *, '(a)' ) ' '

  c = 0.
  do k = 1, n
     do l = 1, n
        do i = 1, n
           c(k,l) = c(k,l) + v(k,i) * v(l,i)* 1/(w(i)**2) ! prior is in data
        enddo
     enddo
  enddo

end subroutine ua
