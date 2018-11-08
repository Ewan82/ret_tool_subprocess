!***********************************************************
!> \file  runretrieval.f90
!> \brief includes main program for minimisation
!> \authors MV/TK, The Inversion Lab
!> \date  September 2018
!***********************************************************
! We make doxygen ignore the following routine
!> \cond
                             
module sizes
  integer :: nn !< length of control vector
  integer :: mm !< length of output vector
end module sizes
!> \endcond

!***********************************************************
!  PROGRAM: runretrieval
!
!> \brief Driver for minimisation routine LBFGSB with reverse mode derivative
!> \authors MV/TK, The Inversion Lab
!> \date  March 2018
!***********************************************************
program runretrieval
  use sizes
  use mo_sensimul, only:get_npts,finishdf
  implicit none
  !-- constants
  character(len=*), parameter ::  xfname = 'x.b'
  character(len=*), parameter :: sxfname = 'sx.b'
  character(len=*), parameter :: cxfname = 'Cx.b'
  !--
  integer :: n, mf, mjac
  real(kind=8), allocatable :: x(:), x0(:), sx(:), g(:), ysim(:)
  real(kind=8), allocatable :: cx(:,:) !-- posterior control vector covariance
  real(kind=8) :: pert
  !-- interface to lbfgsb optimiser
  character(len=*), parameter :: opti_ok1 = 'CONVERGENCE: NORM_OF_PROJECTED_GRADIENT_<=_PGTOL'
  character(len=*), parameter :: opti_ok2 = 'CONVERGENCE: REL_REDUCTION_OF_F_<=_FACTR*EPSMCH'
  character(len=60) :: opti_desc
  logical :: opti_succeed    !-- true iff after optimisation opti_desc is in [opti_ok1,opti_ok2]
  real(kind=8) :: pgtol
  !--
  logical :: run_fapar_target !-- whether to apply FAPAR target processor
                             !   can be disabled on command line, since
                             !   this is the most time-consuming part
  integer :: ldebug
  logical :: tstflg
  integer :: ierr, iostat
  logical :: diag
  !-- commandline
  integer :: iarg,narg, nconsumed
  character(len=32) :: argval, argname

  !-- initialise
  run_fapar_target = .true.
  opti_succeed = .false.
  ldebug = 1
  diag = .true.

  !-- command line handling
  narg=command_argument_count()
  if(narg>0)then
     !loop across options
     nconsumed=0
     arg: do iarg=1,narg
        if( nconsumed.gt.0 ) then
           nconsumed = nconsumed - 1
           cycle arg
        endif
        call get_command_argument(iarg,argname)
        select case(adjustl(argname))
        case("--help","-h")
           write(*,'(/)')
           write(*,'(a)') "INFO: This is program 'runretrieval', look into code to know options"
           stop 0
           exit arg
        case("--dbglev")
           call get_command_argument(iarg+1, argval)
           read( argval, * ) ldebug
           nconsumed = 1
        case('--no_fapar')
           run_fapar_target = .false.
        case default
           write(*,'(a)') "ERROR::Option ",adjustl(argname),"unknown"
           stop 0
        end select
     end do arg
  end if

  !-------------------
  !-- initialise retrieval procedure
  !
  call retrieval_read_ctl()
  call retrieval_dump_ctl()
  !-- get gradient tolerance condition
  call retr_get_gradient_tol(pgtol)
  !-- get perturbation of prior for initial guess
  call retr_get_prior_pert(pert)

  !-------------------
  !-- init model
  write(*, '(a)') ' INFO::runretrieval:calling initf...'
  call initf(n, mf)
  write(*, '(a)') ' INFO::runretrieval:...done.'
  if( ldebug.gt.0 ) then
     write(*, '(a,2(a,i3,1x))') ' DEBUG::runretrieval:initf yields:',&
          'n=',n,'mf=',mf
  endif

  !-- set sizes in module
  nn = n
  mm = mf

  ! init unknowns
  allocate( x(n),x0(n),sx(n),g(n),ysim(mf) )
  allocate( cx(n,n) )
  write(*, '(a)') ' INFO::runretrieval:calling initx...'
  call initx(n, x0, sx)
  write(*, '(a)') ' INFO::runretrieval:...done.'

  !-- init state-vector regularisation model
  call optim_is_state_term_enabled(tstflg)
  if( tstflg ) then
     call state_model_set()
  endif

  !-- perturb initial parameter values
  x = x0
  x = x*(1+pert)

  !-- run optimiser
  write(*, '(a)') ' INFO::runretrieval:start optimise...'
  call optimise(n, x, g, pgtol, opti_desc, ldebug)
  write(*, '(a)') ' INFO::runretrieval:...optimise DONE'
  opti_succeed = (opti_desc.eq.opti_ok1).or.(opti_desc.eq.opti_ok2)
  write(*, '(a,1x,l1)') ' INFO::runretrieval::opti_succeed:',opti_succeed
  
  !-- save optimal state-vector (normalised units), logging
  call printoptimum(n,x0,sx,x,xfname,sxfname)

  !-- uncertainty propagation on state
  write(*, '(a)') ' INFO::runretrieval:start uncertainty analysis on state...'
  mjac = mf+2*n !-- simulation + model + prior
  call ua(mjac, n, x, cx, ldebug.ge.2, diag, ierr)
  write(*, '(a)') ' INFO::runretrieval:...uncertainty analysis on state DONE'

  !-- save posterior uncertainty matrix 'Cx.b' (normalisd units)
  open(unit=1, file=cxfname, iostat=iostat, form='unformatted')
  write(1) cx
  close(1)

  !-- logging
  call printua(n, cx, ldebug.ge.2)

  !-- saving prior/posterior control vector
  !   (with uncertainties) in physical units
  call finishdf(n, x, cx, opti_succeed)

  !-- target operator: FAPAR
  if( run_fapar_target ) then
     call fapar_target(n, x, cx)
  endif

  !-- dispose memory
  deallocate( x, x0, sx, cx, g, ysim )
end program runretrieval


!***********************************************************
!  subroutine fg
!
!> @brief evaluation of function and gradient,
!>        calls subroutine cost() and its derivative routine
!> \authors MV/TK, The Inversion Lab
!> \date  January 2018
!
subroutine optimise(n, x, g, pgtol, task, dbglev)
  implicit none
  ! arguments
  integer, intent(in) :: n
  real(kind=8), intent(inout) :: x(n)
  real(kind=8), intent(out) :: g(n)
  character(len=60), intent(out) :: task !-- to interface with lbfgsb optimisation routine
  real(kind=8), intent(in) :: pgtol
  integer, intent(in) :: dbglev
  ! local decls
  integer,  parameter     :: m = 5, iprint = 1
  real(kind=8), parameter :: factr  = 1.0d+1   !-- extremely high accuracy in lbfgsb
  character(len=60)      :: csave
  logical                :: lsave(4)
  integer                :: isave(44)
  real(kind=8)               :: f
  real(kind=8)               :: dsave(29)
  integer,  allocatable  :: nbd(:), iwa(:)
  real(kind=8), allocatable  :: l(:), u(:), wa(:)
  integer :: iloop
  character(len=4) :: iter_str

  !-- allocate work space
  allocate ( nbd(n), l(n), u(n))
  allocate ( iwa(3*n) )
  allocate ( wa(2*m*n + 5*n + 11*m*m + 8*m) )

  !-------------------
  ! control vector boundaries
  !
  write(*, '(a)') ' INFO::optimise:load control vector boundary settings by calling initb...'
  ! ! declaring no bounds
  ! nbd = 0 ! no bounds
  ! alternative: setting upper and lower bounds
  call initb(n, nbd, l, u)
  write(*, '(a)') ' INFO::optimise:...done.'
  write(*, '(a)') ' INFO::optimise:start optimiser...'
  !     We start the iteration by initializing task.
  task = 'START'
  !     The beginning of the loop
  iloop = 0
  do while(task(1:2).eq.'FG'.or.task.eq.'NEW_X'.or. &
       task.eq.'START') 
     iloop = iloop + 1
     write(iter_str,'(i4.4)') iloop
     !     This is the call to the L-BFGS-B code.
     if( dbglev.ge.2 ) then
        print*, 'DEBUG::iloop=',iloop,'task(before)=',trim(task)
     endif
     call setulb ( n, m, x, l, u, nbd, f, g, factr, pgtol, &
          wa, iwa, task, iprint,&
          csave, lsave, isave, dsave )
     if( dbglev.ge.2 ) then
        write(*,'(a)') ' INFO::optimise:'//&
             'iloop='//iter_str//' task='//trim(task)
     endif
     if (task(1:2) .eq. 'FG') then
        call fg (n,x,f,g,dbglev) ! function and gradient evaluation
     end if
  end do
  !     end of loop do while
  write(*, '(a)') ' INFO::optimise:...optimiser finished.'

  !-- deallocate work space
  deallocate( nbd, l, u)
  deallocate( iwa )
  deallocate( wa )

end subroutine optimise


!***********************************************************
!  subroutine fg
!
!> @brief evaluation of function and gradient,
!>        calls subroutine cost() and its derivative routine
!> \authors MV/TK, The Inversion Lab
!> \date  January 2018
!
subroutine fg (n,x,f,g,ldebug)
  use sizes
  implicit none
  ! arguments
  integer, intent(in)              :: n       !< length of control vector
  real (kind=8), intent(in)    :: x(n)        !< control vector
  real (kind=8), intent(out)   :: f           !< function value
  real (kind=8), intent(out)   :: g(n)        !< gradient value
  integer, intent(in) :: ldebug
  ! local
  real (kind=8) :: fb
  integer :: i

  f = -999._8
  g = 0.
  fb = 1.
  if(ldebug.ge.3 ) then
     print*, 'fg:x------------------------------'
     do i=1,n
        print*, i,x(i)
     enddo
     print*, 'fg:end-x------------------------------'
  endif
  call cost_bw(nn, x, g, mm, f, fb) ! cost_b does not return f
  if(ldebug.ge.3) then
     print*, 'fg:gradient------------------------------'
     do i=1,n
        print*, i,g(i)
     enddo
     print*, 'fg:end-gradient------------------------------'
  endif
  call cost(nn,x,mm, f)
  if(ldebug.ge.3) then
     print*, 'fg:cost------------------------------'
     print*, f
     print*, 'fg:end-cost------------------------------'
  endif
end subroutine fg
   

! We make doxygen ignore the following routine
!> \cond
subroutine printoptimum(n,x0,sx,x,xfname,sxfname)
  implicit none
  ! arguments
  integer, intent(in)          :: n
  real(kind=8),    intent(in)  :: x(n), x0(n), sx(n)
  character(len=*), intent(in) :: xfname, sxfname
  ! local
  integer ( kind = 4 ) :: i

  write ( *, '(a)' ) ' '
  write ( *, '(a9x,3a20)' ) 'i', 'Prior', 'Posterior', 'Change [%]'
  write ( *, '(a)' ) ' '
  do i = 1, n
     print '(i9x,2e20.6,f20.5)', i, sx(i)*x0(i), sx(i)*x(i), 100*(x(i)-x0(i))/x0(i)
  enddo

  open (unit=1, file=xfname, form='unformatted')
  write (1) x
  close (1)
  write(*, '(a)') ' INFO::writing final control vector to '//trim(xfname)

  open (unit=1, file=sxfname, form='unformatted')
  write (1) sx
  close (1)
  write(*, '(a)') ' INFO::writing prior uncertainty (physical units) to '//trim(sxfname)
end subroutine printoptimum
!> \endcond
