program runua
  implicit none
  integer :: n, mf, mjac
  integer :: ierr, i, l, k
  real(kind=8), allocatable :: x(:), sx(:), c(:,:)
  integer :: iostat
  logical :: exist
  logical :: ldebug = .false., diag = .true.
  logical :: tstflg

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
  call wsigma(6,n,sx,c)
  call wsigma(1,n,sx,c)

  !-- post processing ("targets")
  call finishdf(n, x, c)


  do l = 1, n
     do k = 1, n
        c(k,l) = c(k,l)/sqrt(c(k,k))/sqrt(c(l,l))
     enddo
  enddo
  if( ldebug ) then
     write ( *, '(a)' ) ' '
     write ( *, '(a40)' ) 'Correlations of Posterior Uncertainty'
     write ( *, '(9x,18a5)' ) '1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18'
     do i = 1, n
        print '(i9x,18(f4.2x))', i, c(i,:)
     enddo
  endif

  write ( *, '(a)' ) ' '
  write ( *, '(a40)' ) 'Extreme off diagonal correlations'
  write ( *, '(a20x,a20,2a5)' )   'Min/Max ', 'value', 'i', 'j'
  do i = 1, n
     c(i,i) = 0.
  enddo
  write ( *, '(a20x,f20.5,2i5)' )   'Min ', minval(c), minloc(c)
  write ( *, '(a20x,f20.5,2i5)' )   'Max ', maxval(c), maxloc(c)
  write ( *, '(a)' ) ' '



  !-- dispose memory
  deallocate( x, sx, c )

end program runua


subroutine wsigma(unit,n,sx,c)
  implicit none
  ! argruments
  integer ( kind = 4 ), intent(in) :: n, unit
  real ( kind = 8 ), intent(out)   :: sx(n), c(n,n)
  ! local
  integer ( kind = 4 )             :: i
  if (unit.ne.6) open (unit=unit,file='sigma.dat', form='formatted', action='write')
  write (unit, '(a1)' ) '#'
  write (unit, '(a1,a8x,3a20)' ) '#','i', 'Prior Sigma', 'Post Sigma', 'Unc. Red. [%]'
  write (unit, '(a1)' ) '#'
  do i = 1, n
     write (unit,'(i9x,2e20.6,f20.5)') i, sx(i), sqrt(c(i,i))*sx(i), 100*(1-sqrt(c(i,i)))
  enddo
  if (unit.ne.6) close (unit=unit)
end subroutine wsigma
