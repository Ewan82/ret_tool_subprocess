program runua
  use mo_sensimul, only:finishdf
  implicit none
  integer :: n, mf, mjac
  integer :: ierr
  real(kind=8), allocatable :: x(:), sx(:), c(:,:)
  integer :: iostat
  logical :: exist, tstflg
  logical :: ldebug = .false., diag = .true.
  logical :: run_fapar_target !-- whether to apply FAPAR target processor
                             !   can be disabled on command line, since
                             !   this is the most time-consuming part
  ! commandline
  integer :: iarg,narg, nconsumed
  character(len=32) :: argval, argname

  !-- initialise
  run_fapar_target = .true.
  ldebug = .false.
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
           write(*,'(a)') "INFO: This is program 'runua', look into code to know options"
           stop 0
           exit arg
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

  !-------------------
  !-- init model
  write(*, '(a)') ' INFO::runua:calling initf...'
  call initf(n, mf)
  mjac = mf + 2*n
  write(*, '(a)') ' INFO::runua:...done.'
  if( ldebug ) then
     write(*, '(a,6(a,i3,1x))') ' DEBUG::runua:initf yields:',&
          'n=',n,'mf=',mf
  endif

  !-------------------
  ! init unknowns
  allocate (x(n),sx(n),c(n,n))
  write(*, '(a)') ' INFO::runua:calling initx...'
  call initx(n, x, sx)
  write(*, '(a)') ' INFO::runua:...done.'

  !-- control vector might be read from external file
  !   resulting from a retrieval procedure
  inquire(FILE='x.b',EXIST=exist)
  if(exist) then
     write(*, '(a)') ' INFO::runua:reading control vector from file ***'//'x.b'//'***'
     open (unit=1, file='x.b', form='unformatted', iostat=iostat)
     if( iostat.ne.0 ) then
        write(*, '(a)') ' FATAL::runua:'//&
             'failed opening existing file ***x.b***'
        stop
     endif
     read(1, iostat=iostat) x
     if( iostat.ne.0 ) then
        write(*, '(a)') ' FATAL::runua:'//&
             'failed reading control vector from file ***x.b***'
        stop
     endif
     close (1)
     write(*,'(a)') ' INFO::runua:...done'
  endif

  !-- init dynamic state-vector model
  call optim_is_state_term_enabled(tstflg)
  if( tstflg ) then
     call state_model_set()
  endif

  !-- uncertainty analysis
  call ua(mjac, n, x, c, ldebug, diag, ierr)

  !-- logging/output
  call wsigma(6,n,sx,c)
  call wsigma(1,n,sx,c)

  !-- post processing ("targets")
  call finishdf(n, x, c)

  !-- informational output
  call printua(n, c, ldebug)

  !-- target operator: FAPAR
  if( run_fapar_target ) then
     call fapar_target(n, x, c)
  endif

  !-- dispose memory
  deallocate( x, sx, c )

end program runua
