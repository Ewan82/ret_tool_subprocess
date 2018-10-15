subroutine jacobian ( mm, n, x, y, y_jac, ldfjac, iflag )
  use DIFFSIZES
  use mo_retrieval
  implicit none
  ! arguments
  integer, intent(in) :: mm, n, ldfjac, iflag
  real(kind=8), intent(in)  :: x(n)
  real(kind=8), intent(out) :: y(mm),  y_jac(ldfjac,n)
  ! local
  real(kind=8) xd(nbdirsmax,n)
  real(kind=8) obsdiff(mm-2*n),obsdiffd(nbdirsmax,mm-2*n)
  real(kind=8) modeldiff(n), modeldiffd(nbdirsmax,n)
  real(kind=8) priordiff(n), priordiffd(nbdirsmax,n)
  integer  :: i,j
  ! just for printing gradient
  real(kind=8) :: argod(n,mm-2*n), argmd(n, n), argpd(n,n)
  real(kind=8) :: costobsd(n), costmodeld(n), costpriord(n)

  real(kind=8) :: xphys(n)

  integer :: nbdirs ! Variable to be declared without kind statement 
                    ! to conform with generated derivative code
  integer :: m
  ! externals
  external misfit, H_m, residual_prior
  m = mm-2*n

  call x2p(n, x, xphys)


  ! print*, 'MVMV::mm=',mm,'m=',m,'n=',n
  ! print*, 'MVMV::shape(y_jac)=',shape(y_jac)

  if (iflag == 1) then ! evaluate function
! could insert call to fcn here
     obsdiff = 0.
     call misfit(n, xphys, m, obsdiff)
     modeldiff = 0.
     if (retr_use_state_term) call H_m   (n,xphys,modeldiff)
     priordiff = 0.
     if (retr_use_prior_term) call residual_prior(n,x,priordiff)
     y = 0.
     y(1:m)=obsdiff
     y(m+1:m+n)=modeldiff
     y(m+n+1:m+2*n)=priordiff
     if (.true.) print '(9xe20.6,a20,5e20.6)', sum(obsdiff**2)+sum(modeldiff**2)+sum(priordiff**2), '-', &
          sum(obsdiff**2), sum(modeldiff**2), sum(priordiff**2), x(1:min(n,2))
  else if (iflag == 2) then ! evaluate function + jacobian
     y = 0.
     y_jac = 0.
     nbdirs = n

     xd = 0.
     do i = 1, n
        xd(i,i) = 1
     enddo
     obsdiff = 0.
     obsdiffd = 0.
     call MISFIT_FWV(n, xphys, xd, m, obsdiff, obsdiffd, nbdirs)
     ! print*, 'MVMV::jacobian:shape(ranspose(obsdiffd)=',shape(transpose(obsdiffd))
     y_jac(1:m,:) = transpose(obsdiffd(1:n,1:m))

     xd = 0.
     do i = 1, n
        xd(i,i) = 1
     enddo
     modeldiff = 0.
     modeldiffd = 0.
     if (retr_use_state_term) call H_m_fwv( n, xphys, xd, modeldiff, modeldiffd, nbdirs)
     ! print*, 'MVMV::jacobian:shape(ranspose(modeldiffd)=',shape(transpose(modeldiffd))
     y_jac(m+1:m+n,:) = transpose(modeldiffd(1:n,1:n))

     xd = 0.
     do i = 1, n
        xd(i,i) = 1
     enddo
     priordiff = 0.
     priordiffd = 0.
     if (retr_use_prior_term) call residual_prior_fwv(n, x, xd, priordiff, priordiffd, nbdirs)
     ! print*, 'MVMV::jacobian:shape(ranspose(priordiffd)=',shape(transpose(priordiffd))
     y_jac(m+n+1:m+2*n,:) = transpose(priordiffd(1:n,1:n))

     !MVDEBUG
     ! write(*,'(a)') ' MVDEBUG::y_jac:'
     ! do i=m+n+1,m+2*n
     !    do j=1,n
     !       write(*,'(e12.6,1x)',advance='no') y_jac(i,j)
     !    enddo
     !    print*
     ! enddo

     y(1:m)=obsdiff
     y(m+1:m+n)=modeldiff
     y(m+n+1:m+2*n)=priordiff

     !gradient for diagnostics:
     DO i=1,n
        ! cost obs
        argod(i, :) = 2*obsdiff*obsdiffd(i, :)
        costobsd(i) = 0.5*SUM(argod(i, :))
        ! cost model
        argmd(i, :) = 2*modeldiff*modeldiffd(i, :)
        costmodeld(i) = 0.5*SUM(argmd(i, :))
        ! cost prior
        argpd(i, :) = 2*priordiff*priordiffd(i, :)
        costpriord(i) = 0.5*SUM(argpd(i, :))
     END DO
     ngrad = sqrt(sum((costobsd+costmodeld+costpriord)**2))
     costp = sum(priordiff**2)
     costo = sum(obsdiff**2)
     costm = sum(modeldiff**2)
     if (.true.) then
        print '(7a20)' ,'cost_all','ngrad','costo','costm','costp','x(1)','x(2)'
        print '(7e20.6)', costp+costo+costm, &
          ngrad, costo, costm, costp, x(1:min(n,2))
     endif

     ! icall = icall + 1
  else
     print*, 'wrong value of iflag = ',iflag
     stop
  endif
end subroutine jacobian
