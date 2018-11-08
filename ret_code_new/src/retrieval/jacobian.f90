subroutine cost_jacobian ( mm, n, x, y, y_jac, ldfjac, iflag )
  use diffsizes
  use mo_retrieval, only:ngrad,costp,costm,costo
  use mo_retrieval, only:retr_use_prior_term,retr_use_model_term
  implicit none
  ! arguments
  integer, intent(in) :: mm, n, ldfjac, iflag
  real(kind=8), intent(in)  :: x(n)
  real(kind=8), intent(out) :: y(mm),  y_jac(ldfjac,n)
  ! externals
  external misfit, H_m, residual_prior
  ! local
  real(kind=8) xd(nbdirsmax,n)
  real(kind=8) obsdiff(mm-2*n),obsdiffd(nbdirsmax,mm-2*n)
  real(kind=8) modeldiff(n), modeldiffd(nbdirsmax,n)
  real(kind=8) priordiff(n), priordiffd(nbdirsmax,n)
  integer  :: i,m,nbdirs
  logical :: ldebug = .false.
  real(kind=8) :: costobsd(n), costmodeld(n), costpriord(n) !-- for gradient diagnostic

  !-- size of simulation vector
  m = mm-2*n


  if (iflag == 1) then ! evaluate function
     obsdiff = 0.
     call misfit(n, x, m, obsdiff)
     modeldiff = 0.
     if (retr_use_model_term) call H_m   (n,x,modeldiff)
     priordiff = 0.
     if (retr_use_prior_term) call residual_prior(n,x,priordiff)
     y = 0.
     y(1:m)=obsdiff
     y(m+1:m+n)=modeldiff
     y(m+n+1:m+2*n)=priordiff
     if (ldebug) then
        print '(9xe20.6,a20,5e20.6)', sum(obsdiff**2)+sum(modeldiff**2)+sum(priordiff**2), '-', &
             sum(obsdiff**2), sum(modeldiff**2), sum(priordiff**2), x(1:min(n,2))
     endif
  else if (iflag == 2) then ! evaluate function + jacobian
     y = 0.
     y_jac = 0.
     
     !-- compute all n directional derivatives
     nbdirs = n

     !-- set directional seed matrix
     xd = 0.
     do i = 1, nbdirs
        xd(i,i) = 1
     enddo

     !-- observational part of cost-jacobian
     obsdiff = 0.
     obsdiffd = 0.
     call MISFIT_FWV(n, x, xd, m, obsdiff, obsdiffd, nbdirs)
     y(1:m)=obsdiff
     y_jac(1:m,:) = transpose(obsdiffd(1:n,1:m))

     !-- model part of cost-jacobian
     modeldiff = 0.
     modeldiffd = 0.
     if (retr_use_model_term) call H_m_fwv( n, x, xd, modeldiff, modeldiffd, nbdirs)
     y(m+1:m+n)=modeldiff
     y_jac(m+1:m+n,:) = transpose(modeldiffd(1:n,1:n))

     !-- prior part of cost-jacobian
     priordiff = 0.
     priordiffd = 0.
     if (retr_use_prior_term) call residual_prior_fwv(n, x, xd, priordiff, priordiffd, nbdirs)
     y(m+n+1:m+2*n)=priordiff
     y_jac(m+n+1:m+2*n,:) = transpose(priordiffd(1:n,1:n))


     !gradient for diagnostics:
     do i=1,n
        ! cost obs
        costobsd(i) = sum(obsdiff*obsdiffd(i, :))
        ! cost model
        costmodeld(i) = sum(modeldiff*modeldiffd(i, :))
        ! cost prior
        costpriord(i) = sum(priordiff*priordiffd(i, :))
        ! cost prior
     end do
     ngrad = sqrt(sum((costobsd+costmodeld+costpriord)**2))
     costp = 0.5_8*sum(priordiff**2)
     costo = 0.5_8*sum(obsdiff**2)
     costm = 0.5_8*sum(modeldiff**2)

     if( ldebug ) then
        print '(7a20)' ,'cost_all','ngrad','costo','costm','costp','x(1)','x(2)'
        print '(7e20.6)', costp+costo+costm, &
             ngrad, costo, costm, costp, x(1:min(n,2))
     endif
  else
     print*, 'wrong value of iflag = ',iflag
     stop
  endif
end subroutine cost_jacobian


subroutine fapar_jacobian( m, n, x, tjac )
  use diffsizes
  implicit none
  ! arguments
  integer, intent(in) :: m, n
  real(kind=8), intent(in) :: x(n)
  real(kind=8), intent(out) :: tjac(m,n)
  ! local declarations
  external simulate_fapar_fwv
  real(kind=8) :: y(m), yd(nbdirsmax,m), xd(nbdirsmax,n)
  integer :: i

  !-- set directions
  xd = 0._8
  do i=1,n
     xd(i,i) = 1._8
  enddo

  !-- directional derivatives from AD code
  call simulate_fapar_fwv(n, x, xd, m, y, yd, n)

  !-- need to transpose
  tjac = transpose( yd(1:n,1:m) )

end subroutine fapar_jacobian
