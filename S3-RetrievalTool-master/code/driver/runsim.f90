program runsim
  use mo_sensimul
  implicit none

  ! commandline controlled
  logical :: silent = .false.!true should yield less console output
  logical :: s2_loss_visnir = .false.   !-- set S2 simulations in the visible bands to missing value
  logical :: s1_loss_vh      = .false.
  logical :: s1a_loss_vh     = .false.
  logical :: s1b_loss_vh     = .false.
  logical :: s1_loss_vv      = .false.
  logical :: s1a_loss_vv     = .false.
  logical :: s1b_loss_vv     = .false.
  logical :: s2_loss_swnir   = .false.     !-- set S2 simulations in the NIR bands to missing value
  logical :: s2a_loss_visnir = .false.     !-- set S2A simulations in the visible bands to missing value
  logical :: s2a_loss_swnir  = .false.     !-- set S2A simulations in the NIR bands to missing value
  logical :: s2b_loss_visnir = .false.     !-- set S2B simulations in the visible bands to missing value
  logical :: s2b_loss_swnir  = .false.     !-- set S2B simulations in the NIR bands to missing value

  ! other
  integer             :: j, n, m
  real(kind=8), allocatable :: x(:), sx(:), ysim(:), xphys(:)
  logical :: exist
  logical :: ldebug = .true.

  !-------------------
  ! command line
  integer :: narg !#of arg & counter of arg
  integer :: nskip
  character(len=64) :: argname
  character(len=64) :: argval


  !-------------------
  !     command line
  !
  narg=command_argument_count()
  !Loop over the arguments
  if(narg>0)then
     !loop across options
     nskip=0
     arg: do j=1,narg
        if( nskip.gt.0 ) then
           nskip = nskip - 1
           cycle arg
        endif
        call get_command_argument(j,argname)
        select case(adjustl(argname))
        case("--help","-h")
           write(*,'(/)')
           write(*,'(a)') "INFO: This is program 'runsim', look into code to see options"
           stop 0
           exit arg
        case('--s1_loss_vh')
           s1_loss_vh = .true.
        case('--s1a_loss_vh')
           s1a_loss_vh = .true.
        case('--s1b_loss_vh')
           s1b_loss_vh = .true.
        case('--s1_loss_vv')
           s1_loss_vv = .true.
        case('--s1a_loss_vv')
           s1a_loss_vv = .true.
        case('--s1b_loss_vv')
           s1b_loss_vv = .true.
        case('--s2_loss_visnir')
           s2_loss_visnir = .true.
        case('--s2_loss_swnir')
           s2_loss_swnir = .true.
        case('--s2a_loss_visnir')
           s2a_loss_visnir = .true.
        case('--s2a_loss_swnir')
           s2a_loss_swnir = .true.
        case('--s2b_loss_visnir')
           s2b_loss_visnir = .true.
        case('--s2b_loss_swnir')
           s2b_loss_swnir = .true.
        case("-s")
           silent = .true.
        case default
           write(*,'(a)') "ERROR::Option ",adjustl(argname),"unknown"
           stop 0
        end select
     end do arg
  end if


  !-------------------
  !-- init model
  write(*, '(a)') ' INFO::runsim:calling initf...'
  call initf(n, m)
  write(*, '(a)') ' INFO::runsim:...done.'
  if( ldebug ) then
     write(*, '(a,2(a,i3,1x))') ' DEBUG::runsim:initf yields:',&
          'n=',n,'m=',m
  endif


  !-------------------
  ! allocate arrays
  allocate(x(n),xphys(n),sx(n),ysim(m))

  !-------------------
  ! init unknowns !!!physical values!!!
  write(*, '(a)') ' INFO::runsim:calling initx...'
  call initx(n, x, sx)
  write(*, '(a)') ' INFO::runsim:...DONE.'

  !-------------------
  ! possibly overwrite
  !
  inquire(FILE='x-external.b',EXIST=exist)
  if (exist) then
     open (unit=1, file='x-external.b', form='unformatted')
     read (1) x
     close (1)
     write(*, '(a)') ' INFO::runsim:have read control vector from x-external.b'//&
          ' !!!expecting scaled control vector!!!'
  endif


  !--
  call x2p(n, x, xphys)

  if( ldebug ) then
     write(*,'(a)') ' DEBUG::runsim::calling simulate_s1s2 at x ...'
     write(*, '(a3,3(a15))' ) 'j', 'x-physical', 'x-scaled', 'x-sigma'
     do j=1,n
        write(*, '(i3,3(f15.8))' ) j, xphys(j), x(j), sx(j)
     enddo
  endif

  !-- run combined S1+S2 simulation
  write(*,'(a)') ' INFO::runsim:calling combinded S1+S2 simulation with simulate_s1s2...'
  call simulate_s1s2(n, x, m, ysim)
  write(*,'(a)') ' INFO::runsim:...done.'

  !-- apply filtering (if selected by user)
  !-- S1 failure scenarios
  if( s1_loss_vh ) then
     s1_failure_satellite = 'S1'
     s1_failure_pol       = 'vh'
     write(*, '(a)')' INFO::runsim:apply failure scenario "S1 loss VH"'
  else if( s1a_loss_vh ) then
     s1_failure_satellite = 'S1A'
     s1_failure_pol       = 'vh'
     write(*, '(a)')' INFO::runsim:apply failure scenario "S1 loss VH"'
  else if( s1B_loss_vh ) then
     s1_failure_satellite = 'S1B'
     s1_failure_pol       = 'vh'
     write(*, '(a)')' INFO::runsim:apply failure scenario "S1 loss VH"'
  else if( s1_loss_vv ) then
     s1_failure_satellite = 'S1'
     s1_failure_pol       = 'vv'
     write(*, '(a)')' INFO::runsim:apply failure scenario "S1 loss VV"'
  else if( s1a_loss_vv ) then
     s1_failure_satellite = 'S1A'
     s1_failure_pol       = 'vv'
     write(*, '(a)')' INFO::runsim:apply failure scenario "S1 loss VV"'
  else if( s1B_loss_vv ) then
     s1_failure_satellite = 'S1B'
     s1_failure_pol       = 'vv'
     write(*, '(a)')' INFO::runsim:apply failure scenario "S1 loss VV"'
  endif
  !-- S2 failure scenarios
  if( s2_loss_visnir ) then
     s2_failure_satellite = 'S2'
     write(*, '(a)')' INFO::runsim:apply failure scenario "S2 loss vis/nir bands" (1-9)'
     s2_failure_bands(1:10) = .true.
  else if( s2_loss_swnir ) then
     s2_failure_satellite = 'S2'
     write(*, '(a)')' INFO::runsim:apply failure scenario "S2 loss swnir bands" (10-12)'
     s2_failure_bands(11:13) = .true.
  else if( s2a_loss_visnir ) then
     s2_failure_satellite = 'S2A'
     write(*, '(a)')' INFO::runsim:apply failure scenario "S2A loss vis/nir bands" (1-9)'
     s2_failure_bands(1:10) = .true.
  else if( s2a_loss_swnir ) then
     s2_failure_satellite = 'S2A'
     write(*, '(a)')' INFO::runsim:apply failure scenario "S2A loss swnir bands" (10-12)'
     s2_failure_bands(11:13) = .true.
  else if( s2b_loss_visnir ) then
     s2_failure_satellite = 'S2B'
     write(*, '(a)')' INFO::runsim:apply failure scenario "S2B loss vis/nir bands" (1-9)'
     s2_failure_bands(1:10) = .true.
  else if( s2b_loss_swnir ) then
     s2_failure_satellite = 'S2B'
     write(*, '(a)')' INFO::runsim:apply failure scenario "S2B loss swnir bands" (10-12)'
     s2_failure_bands(11:13) = .true.
  endif

  !-- 
  call simulation_failure_filter(m, ysim)

  !-- create simulation output
  call finishf(n, x, m, ysim)



  ! dispose memory
  call dispose()

  contains

    subroutine dispose()
      implicit none
      if( allocated(x) )    deallocate(x)
      if( allocated(xphys) )    deallocate(xphys)
      if( allocated(sx) )   deallocate(sx)
      if( allocated(ysim) ) deallocate(ysim)
    end subroutine dispose

end program runsim
