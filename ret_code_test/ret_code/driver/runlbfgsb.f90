!***********************************************************
!> \file  runlbfgsb.f90
!> \brief includes main program for minimisation
!> \authors MV/TK, The Inversion Lab
!> \date  Febuary 2018
!***********************************************************
! We make doxygen ignore the following routine
!> \cond
                             
module sizes
  integer :: nn !< length of control vector
  integer :: mm !< length of output vector
end module sizes
!> \endcond

!***********************************************************
!  PROGRAM: runopti
!
!> \brief Driver for minimisation routine LBFGSB with reverse mode derivative
!> \authors MV/TK, The Inversion Lab
!> \date  March 2018
!***********************************************************
program runopti
  use sizes
  implicit none
  integer :: iter, n, mf
  real (kind=8), allocatable :: x(:), x0(:), sx(:), g(:), ysim(:)

  integer,  parameter    :: m = 5, iprint = 1
  real(kind=8), parameter    :: factr  = 1.0d+1
  real(kind=8) :: pgtol, pert
  !
  character(len=60)      :: task, csave
  logical                :: lsave(4)
  integer                :: isave(44)
  real(kind=8)               :: f
  real(kind=8)               :: dsave(29)
  integer,  allocatable  :: nbd(:), iwa(:)
  real(kind=8), allocatable  :: l(:), u(:), wa(:)
  integer :: ldebug
  logical :: tstflg
  integer :: i
  integer :: iloop
  ! commandline
  integer :: narg, nconsumed
  character(len=32) :: argval, argname

  !-- initialise
  ldebug = 1


  !Check if any arguments are found
  narg=command_argument_count()
  !Loop over the arguments
  if(narg>0)then
     !loop across options
     nconsumed=0
     arg: do i=1,narg
        if( nconsumed.gt.0 ) then
           nconsumed = nconsumed - 1
           cycle arg
        endif
        call get_command_argument(i,argname)
        select case(adjustl(argname))
        case("--help","-h")
           write(*,'(/)')
           write(*,'(a)') "INFO: This is program 'runopti', look into code to know options"
           stop 0
           exit arg
        case("--dbglev")
           call get_command_argument(i+1, argval)
           read( argval, * ) ldebug
           nconsumed = 1
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
  write(*, '(a)') ' INFO::runopti:calling initf...'
  call initf(n, mf)
  write(*, '(a)') ' INFO::runopti:...done.'
  if( ldebug.gt.0 ) then
     write(*, '(a,2(a,i3,1x))') ' DEBUG::runopti:initf yields:',&
          'n=',n,'mf=',mf
  endif

  !-- set sizes in module
  nn = n
  mm = mf


  !-------------------
  ! init unknowns
  allocate (x(n),x0(n),sx(n),g(n),ysim(mf))
  write(*, '(a)') ' INFO::runopti:calling initx...'
  call initx(n, x0, sx)
  write(*, '(a)') ' INFO::runopti:...done.'

  !-- init state-vector regularisation model
  call optim_is_state_term_enabled(tstflg)
  if( tstflg ) then
     call state_model_set()
  endif


  !-- perturb initial parameter values
  x = x0
  x = x*(1+pert)

  ! optimiser settings
  allocate ( nbd(n), l(n), u(n))
  allocate ( iwa(3*n) )
  allocate ( wa(2*m*n + 5*n + 11*m*m + 8*m) )

  !-------------------
  ! control vector boundaries
  !
  write(*, '(a)') ' INFO::runopti:load control vector boundary settings by calling initb...'
  ! ! declaring no bounds
  ! nbd = 0 ! no bounds
  ! alternative: setting upper and lower bounds
  call initb(n, nbd, l, u)
  write(*, '(a)') ' INFO::runopti:...done.'
  write(*, '(a)') ' INFO::runopti:start optimiser...'
  !     We start the iteration by initializing task.
  task = 'START'
  !     The beginning of the loop
  iloop = 0
  do while(task(1:2).eq.'FG'.or.task.eq.'NEW_X'.or. &
       task.eq.'START') 
     iloop = iloop + 1
     !     This is the call to the L-BFGS-B code.
     if( ldebug.gt.2 ) then
        print*, 'DEBUG::iloop=',iloop,'task(before)=',trim(task)
     endif
     call setulb ( n, m, x, l, u, nbd, f, g, factr, pgtol, &
          wa, iwa, task, iprint,&
          csave, lsave, isave, dsave )
     if( ldebug.gt.2 ) then
        print*, 'DEBUG::task(after)=',trim(task)
     endif
     if (task(1:2) .eq. 'FG') then
        call fg (n,x,f,g,ldebug) ! function and gradient evaluation
     end if
  end do
  !     end of loop do while
  write(*, '(a)') ' INFO::runopti:...optimiser finished.'

  call printoptimum(n,x0,sx,x)

  ! !-- run combined S1+S2 simulation
  ! write(*,'(a)') ' INFO::runsim:calling combinded S1+S2 simulation with simulate_s1s2...'
  ! call simulate_s1s2(n, x, mf, ysim)
  ! write(*,'(a)') ' INFO::runsim:...done.'
  ! !-- create simulation output
  ! call finishf(n, x, mf, ysim)

end program runopti

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
  integer, intent(in)              :: n           !< length of control vector
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
  if(ldebug.ge.2 ) then
     print*, 'fg:x------------------------------'
     do i=1,n
        print*, i,x(i)
     enddo
     print*, 'fg:end-x------------------------------'
  endif
  call cost_bw(nn, x, g, mm, f, fb) ! cost_b does not return f
  if(ldebug.ge.2) then
     print*, 'fg:gradient------------------------------'
     do i=1,n
        print*, i,g(i)
     enddo
     print*, 'fg:end-gradient------------------------------'
  endif
  call cost(nn,x,mm, f)
  if(ldebug.ge.2) then
     print*, 'fg:cost------------------------------'
     print*, f
     print*, 'fg:end-cost------------------------------'
  endif
end subroutine fg
   
! We make doxygen ignore the following routine
!> \cond
subroutine printoptimum(n,x0,sx,x)
  implicit none
  ! arguments
  integer ( kind = 4 ), intent(in) ::  n
  real (kind=8),    intent(in) :: x(n), x0(n), sx(n)
  ! local
  integer ( kind = 4 ) :: i

  write ( *, '(a)' ) ' '
  write ( *, '(a9x,3a20)' ) 'i', 'Prior', 'Posterior', 'Change [%]'
  write ( *, '(a)' ) ' '
  do i = 1, n
     print '(i9x,2e20.6,f20.5)', i, sx(i)*x0(i), sx(i)*x(i), 100*(x(i)-x0(i))/x0(i)
  enddo

  open (unit=1, file='x.b', form='unformatted')
  write (1) x
  close (1)
  print*, 'writing final control vector to x.b'

  open (unit=1, file='sx.b', form='unformatted')
  write (1) sx
  close (1)
  print*, 'writing prior uncertainty to sx.b'
end subroutine printoptimum
!> \endcond
