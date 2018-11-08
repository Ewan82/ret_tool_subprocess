!***********************************************************
!     fapar_target
!
!> @brief Implementation of FAPAR target processor
!
!> @param[in]   n     length of control vector
!> @param[in]   x     control vector normalised by prior uncertainty
!                     (expected ordering is: S1 related parameter(s) followed by
!                     by state variables LAI,HC,SM per 'simulation point')
!> @param[in]   xcov  covariance matrix of control vector ('Cx.b')
!                     (in normalised control  vector units)
!
!> \authors MV/TK, The Inversion Lab
!> \date    February 2018
!
subroutine fapar_target(n, x, xcov)
  use mo_sensimul
  use mo_util, only:print_matrix
  implicit none
  ! arguments
  integer, intent(in)      :: n          !< length of control vector
  real(kind=8), intent(in) :: x(n)       !< control vector (normalised, no physical units)
  real(kind=8), intent(in) :: xcov(n,n)  !< posterior covariance matrix
  ! local decls
  real(kind=8) :: x0(n), sxpr(n)
  real(kind=8), allocatable :: fapar_pr(:), faparunc_pr(:), fapar_po(:), faparunc_po(:)
  real(kind=8), allocatable :: tjac(:,:) !-- FAPAR Jacobian
  real(kind=8), allocatable :: ct_pr(:,:), ct_po(:,:)
  logical :: succeed
  integer :: i,j,k,l,m
  character(len=256) :: target_fname

  !-- get the prior
  call initx(n, x0, sxpr)

  !-- number of time-points (FAPAR can only be simulated at those where geometry is available)
  m = get_npts()

  !--
  allocate( fapar_pr(m), fapar_po(m), faparunc_pr(m), faparunc_po(m) )
  allocate( tjac(m,n), ct_pr(m,m), ct_po(m,m) )

  !-- FAPAR: ***prior,posterior***
  write(*, '(a)') ' INFO::fapar_target:'//&
       'start computation of FAPAR at prior...'
  call simulate_fapar(n, x0, m, fapar_pr)
  write(*, '(a)') ' INFO::fapar_target:'//&
       '...DONE'
  write(*, '(a)') ' INFO::fapar_target:'//&
       'start computation of FAPAR at x...'
  call simulate_fapar(n, x, m, fapar_po)
  write(*, '(a)') ' INFO::fapar_target:'//&
       '...DONE'

  !-- FAPAR Jacobian
  write(*, '(a)') ' INFO::fapar_target:'//&
       'start computation of FAPAR Jacobian at x...'
  call fapar_jacobian(m, n, x, tjac)
  write(*, '(a)') ' INFO::fapar_target:'//&
       '...DONE'

  !-- uncertainty propagation for prior/post FAPAR
  !-- compute: T'*Cx*T'^t (D5v2, Eq. 1.11)
  ! NOTE::since we are operating 'x' in normalised (i.e. scaled by prior sigma)
  !       coordinates, Ct_pr is just the identity
  do j=1,m
     do i=1,m
        ct_pr(i,j) = 0._8
        ct_po(i,j)  = 0._8
        do k=1,n
           !-- tjac^t(k,j) = tjac(j,k)
           ct_pr(i,j) = ct_pr(i,j) + tjac(i,k)*tjac(j,k)
           do l=1,n
              !-- tjac^t(l,j) = tjac(j,l)
              ct_po(i,j) = ct_po(i,j) + tjac(i,k)*xcov(k,l)*tjac(j,l)
           enddo
        enddo
     enddo
  enddo

  !-- extract uncertainties only (diagonal)
  !   but only when FAPAR was simulated
  faparunc_pr = sim_fill_value
  faparunc_po = sim_fill_value
  do i=1,m
     if( fapar_pr(i).ne.sim_fill_value ) then
        faparunc_pr(i)  = sqrt(ct_pr(i,i))
        faparunc_po(i)  = sqrt(ct_po(i,i))
     endif
  enddo
     
  !-- write prior FAPAR
  target_fname = 'fapar_prior.nc'
  call ncwrt_target_vector(target_fname, 'fapar', m, fapar_pr, faparunc_pr, succeed)
  if( succeed ) then
     write(*, '(a)') ' INFO::fapar_target:'//&
          'generated target data file ***'//trim(target_fname)//'***'
  else
     write(*, '(a)') ' ERROR::fapar_target:'//&
          'error occurred when writing target data file ***'//trim(target_fname)//'***'
  endif

  !-- write post FAPAR
  target_fname = 'fapar_post.nc'
  call ncwrt_target_vector(target_fname, 'fapar', m, fapar_po, faparunc_po, succeed)
  if( succeed ) then
     write(*, '(a)') ' INFO::fapar_target:'//&
          'generated target data file ***'//trim(target_fname)//'***'
  else
     write(*, '(a)') ' ERROR::fapar_target:'//&
          'error occurred when writing target data file ***'//trim(target_fname)//'***'
  endif

  !-- dispose memory
  deallocate( fapar_pr, fapar_po, faparunc_pr, faparunc_po )
  deallocate( tjac, ct_pr, ct_po )

end subroutine fapar_target
