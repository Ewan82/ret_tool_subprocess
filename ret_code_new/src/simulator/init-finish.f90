!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> \file init-finish.f90
!> \brief provides interfaces to control and run generic S1/S2 observation operators
!>               suitable for use within optimisation framework.
!> \authors The Inversion Lab (Michael Vossbeck, Thomas Kaminski) 
!> \date  April 2018
!                               
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!***********************************************************
!  MODULE mo_prior
!> \brief stores prior control vector
!> \authors The Inversion Lab (Michael Vossbeck, Thomas Kaminski) 
!> \date April 2018
module mo_prior
  implicit none

  !-- the 1D prior control vector will consist of
  !   param_s1 + param_s2 + state
  real(kind=8), allocatable :: xpr(:) !< prior control vector
  real(kind=8), allocatable :: sxpr(:)!< prior control vector uncertainty

end module mo_prior


!***********************************************************
!  MODULE mo_obs
!> \brief stores observations and their uncertainties
!> \authors The Inversion Lab (Michael Vossbeck, Thomas Kaminski) 
!> \date April 2018
module mo_obs
  implicit none
  real(kind=8), allocatable :: obs_vector_s1(:)    !< s1 observations
  real(kind=8), allocatable :: obsunc_vector_s1(:) !< s1 observations uncertainty
  real(kind=8), allocatable :: obs_vector_s2(:)    !< s2 observations
  real(kind=8), allocatable :: obsunc_vector_s2(:) !< s2 observations uncertainty
  logical :: obs_are_set = .false.

contains
  integer function get_nobs_s1()
    implicit none
    get_nobs_s1 = 0
    if( allocated(obs_vector_s1) ) get_nobs_s1 = size(obs_vector_s1)
  end function get_nobs_s1
  integer function get_nobs_s2()
    implicit none
    get_nobs_s2 = 0
    if( allocated(obs_vector_s2) ) get_nobs_s2 = size(obs_vector_s2)
  end function get_nobs_s2

end module mo_obs

!***********************************************************
!     initf
!
!> @brief determines sizes of control vector (n) and overall size of
!>        simulation vector (m)
!>         In addition necessary initialisations to actually run the observational operator(s)
!>         are performed.
!
!> @details This is a preparation routine for any downstream client application
!>          in order to prepare the run time environment for the
!>          Sentinel Simulator retrieval system.
!>          In addition to returning the sizes of the control and the simulation
!>          vector, this method organises all preparatory allocations and initialisations.
!>          It therefore relies on reading the external files 'retrconfig.nc' and
!>          'retrprior.nc' (as described in D5v1, chapter 2.2.1)
!
!> @param[out]  n  overall length of control vector
!> @param[out]  m  overall length of simulation vector
!
subroutine initf(n, m)
  use mo_obs
  use mo_sensimul
  implicit none
  ! arguments
  integer, intent(out) :: n,m
  ! local
  logical :: ldebug = .true.
  logical :: succeed

  !-- basic initialisation
  call simulator_init()

  !-- read retrieval configuration
  call ncrd_retrieval_config(config_file, succeed)
  if( .not.succeed ) then
     write(*, '(a)') ' FATAL::initf:reading retrieval configuration from '//&
          '***'//trim(config_file)//'*** failed, cannot continue.'
     stop
  else
     write(*, '(A)') ' INFO::initf:'//&
          'retrieval configuration was read from ***'//trim(config_file)//'***'
     if( ldebug ) then
        write(*, '(a,3(a,i4,1x))') ' DEBUG::initf:ncrd_retrieval_config has finished, yields:',&
             'npts=',npts,'npts_s1=',npts_s1,'npts_s2=',npts_s2
     endif
  endif

  !-- load site-specific settings
  call load_sitenml(site_nml)

  !-- n: number of 'time-points' times number of states per point
  n    = get_n()

  !-- m: dimension of simulation vector
  m    = get_m()

  !-- set prior control vector
  call setprior()

contains
  subroutine setprior()
    use mo_prior
    use mo_sensimul, only:prior_file
    implicit none
    ! local decls
    logical :: io_succeed
    logical :: ldebug = .true.

    call ncrd_retrieval_prior( prior_file, io_succeed )

    if( .not. io_succeed ) then
       write(*, '(a)') ' FATAL::setprior:'//&
            'reading prior control vector from ***'//trim(prior_file)//'*** failed, cannot continue.'
       stop
    else
       write(*, '(a)') ' INFO::setprior:'//&
            'prior control vector was read from ***'//trim(prior_file)//'***'
       if( ldebug ) then
          write(*, '(a)') ' DEBUG::setprior:ncrd_retrieval_prior terminated successfully'
       endif
    endif

  end subroutine setprior

  !***********************************************************
  !     load_sitenml
  !
  !> @brief reads site description namelist file and assigns values to
  !         global variables used within the retrieval system.
  !
  !> @details The expected format of the namelist file is shown below:
  !>------------------------------
  !> &site_params
  !>     !location
  !>     lat = 12.88  !latitude of site
  !>     lon = 48.684  !longitude of site
  
  !>     !soil characteristics
  !>     clay = 0.23  !fraction soil texture clay
  !>     sand = 0.27  !fraction soil texture sand
  !>     bulk = 1.65  !soil bulk density (g cm-3)
  
  !>     !microwave sense params
  !>     freq = 5.405  !frequency (GHz)
  !>     s = 0.015  !surface rms height (m)
  !>     lai_coeff = 0.1  !coefficient of lai for which to calculate extinction and volume scattering coefficients
  !>     omega = 0.1  !coefficient for calculation of volume scattering coefficients
  
  !>     !optical semi discrete params
  !>     mode = 'fast'  !mode to run semiDiscrete in 'fast' or 'slow'
  !>     rsl1 = 0.2  !weight of the first soil vector, default: 0.2
  !>     sm_coeff = 0.5  !weighting of soil moisture impact, bound between (0,1)
  !>     cab = 75.0  !leaf chlorophyl concentration, default: 75.0
  !>     cw = 0.01  !equivelant leaf water thickness, default: 0.01
  !> /
  !>------------------------------
  !
  !> @param[in]   fname  fully qualified name of namelist file
  !
  !
  subroutine load_sitenml(fname)
    use mo_sensimul_s1, only:freq, sand, clay, bulk, sfc_rms_height, omega
    use mo_sensimul_s1, only:rt_surface, rt_canopy
    use mo_sensimul_s2, only:cab,cw,rsl1_default,sm_coeff
    implicit none
    ! arguments
    character(len=*), intent(in) :: fname
    ! local decls
    integer, parameter :: iounit=18
    integer :: iostat
    logical :: exist
    real(kind=8) :: lon, lat
    real(kind=8) :: s, lai_coeff
    real :: rsl1
    character(len=4) :: mode
    namelist /site_params/ &
         lat,&
         lon,&
         clay,&
         sand,&
         bulk,&
         freq,&
         s,&
         lai_coeff,&
         omega,&
         mode,&
         rsl1,&
         sm_coeff,&
         cab,&
         cw
         
    !--
    inquire(file=fname, exist=exist)
    if( .not. exist ) then
       write(*, '(a)') ' FATAL::load_sitenml:'//&
            'non-existent expected file ***'//trim(fname)//'***'
       stop
    else
       open(unit=iounit, file=fname)
       read(unit=iounit, nml=site_params, iostat=iostat)
       if( iostat.gt.0 ) then
          write(*, '(a)') ' FATAL::load_sitenml:error reading namelist file '//&
               '***'//trim(fname)//'***'
          stop
       else
          write(*, '(a)' ) ' INFO::load_sitenml:'//&
               'successfully read site namelist file ***'//trim(fname)//'***'
       endif
       close(iounit)

       !-- mapping to retrieval system variables
       !   (i.e. map variable names as used in namelist to
       !    those defined in the retrieval tool)
       sfc_rms_height = s
       rsl1_default   = rsl1

    endif

  end subroutine load_sitenml


end subroutine initf


!***********************************************************
!     initx
!
!> @brief method to access the prior control vector in normalised units
!>        and the prior uncertainty in physical units
!>
!> @param[in]   n   overall length of control vector
!> @param[out]  x   control vector (in normalised units)
!> @param[out]  sx  uncertainty of control vector elements (physical units)
!
subroutine initx(n, x, sx)
  implicit none
  ! arguments
  integer, intent(in) :: n
  real(kind=8), intent(out) :: x(n), sx(n)

  !-- access the prior
  call getprior(n, x, sx)

  !-- scale to normalised control vector
  x = x/sx

end subroutine initx


!***********************************************************
!     getprior
!
!> @brief method to access the prior control vector in phyiscal units
!>        and the prior uncertainty in physical units
!>
!> @param[in]   n   overall length of control vector
!> @param[out]  x   control vector (in physical units)
!> @param[out]  sx  uncertainty of control vector elements (physical units)
!
subroutine getprior(n, x, sx)
  use mo_prior
  implicit none
  ! arguments
  integer, intent(in) :: n
  real(kind=8), intent(out) :: x(n), sx(n)

  if( size(xpr).ne.n ) then
     write(*, '(a,i4,1x,a,i4)') ' FATAL::getprior:n=',n,'does not equal size of prior=',size(xpr)
     stop
  endif

  !-- prior
  x(1:n) = xpr(1:n)

  !-- prior uncertainty
  sx(1:n) = sxpr(1:n)

end subroutine getprior


!***********************************************************
!     getprior_s1
!
!> @brief method to access the S1 relevant slice of the prior control vector
!>        in phyiscal units and the prior uncertainty in physical units
!>
!> @param[in]   n   length of S1 related elements in control vector
!> @param[out]  x   control vector (in physical units)
!> @param[out]  sx  uncertainty of control vector elements (physical units)
!
subroutine getprior_s1(n, x, sx)
  use mo_prior
  use mo_sensimul
  ! arguments
  integer, intent(in) :: n
  real(kind=8), intent(out) :: x(n), sx(n)
  ! local decls
  integer :: i,i_s1,j

  !-- dimensional consistency
  if( n.ne.get_nc_s1() ) then
     write(*, '(a,2(a,i4,1x))') ' FATAL::getprior_s1:',&
          'expected=',get_nc_s1(),'got=',n
     stop
  endif

  x(1:nparam_s1) = xpr(1:nparam_s1)
  sx(1:nparam_s1) = sxpr(1:nparam_s1)

  j = nparam_s1+1    ! pos in full control vector
  i_s1 = nparam_s1+1 ! pos in S1 restricted control vector
  do i=1,npts
     if( idx_is_s1(i) ) then
        x (i_s1:i_s1+nsc-1) = xpr (j:j+nsc-1)
        sx(i_s1:i_s1+nsc-1) = sxpr(j:j+nsc-1)
        i_s1 = i_s1 + nsc ! increment pos in restricted control vector
     endif
     j = j + nsc ! increment pos in full prior control vector
  enddo

end subroutine getprior_s1


!***********************************************************
!     getprior_s2
!
!> @brief method to access the S2 relevant slice of the prior control vector
!>        in phyiscal units and the prior uncertainty in physical units
!>
!> @param[in]   n   length of S1 related elements in control vector
!> @param[out]  x   control vector (in physical units)
!> @param[out]  sx  uncertainty of control vector elements (physical units)
!
subroutine getprior_s2(n, x, sx)
  use mo_prior
  use mo_sensimul
  ! arguments
  integer, intent(in) :: n
  real(kind=8), intent(out) :: x(n), sx(n)
  ! local decls
  integer :: i,i_s2,j

  !-- dimensional consistency
  if( n.ne.get_nc_s2() ) then
     write(*, '(a,2(a,i4,1x))') ' FATAL::getprior_s2:',&
          'expected=',get_nc_s2(),'got=',n
     stop
  endif

  j = nparam_s1+1    ! pos in full control vector
  i_s2 = 1           ! pos in S2 restricted control vector
  do i=1,npts
     if( idx_is_s2(i) ) then
        x (i_s2:i_s2+nsc-1) = xpr (j:j+nsc-1)
        sx(i_s2:i_s2+nsc-1) = sxpr(j:j+nsc-1)
        i_s2 = i_s2 + nsc ! increment pos in control vector restricted to S2
     endif
     j = j + nsc !-- increment pos in full prior control vector
  enddo

end subroutine getprior_s2


!***********************************************************
!     getobs
!
!> @brief Method to load observations and their respective uncertainties
!>        into the retrieval system.
!>
!> @details observations will be read from files 'obs_s1.nc' and 'obs_s2.nc'
!>          (as denoted in D5v2, Section 2.2.1).
!>          The retrieval system is capable of handling the absence of either
!>          of these files. Internally, the I/O will be done only once and
!>          the observations are cached to speedup repeated access.
!
!> @param[in]   m        length of simulation vector (S1 and S2)
!> @param[out]  yobs     observation vector
!> @param[out]  syobs    uncertainties of observation vector elements
!> @param[out]  nobs_s1  actual number of S1 observations in observation vector
!> @param[out]  nobs_s2  acutal number of S2 observations in observation vector
!> @param[out]  succeed  flag indicating 
!
subroutine getobs(m, yobs, syobs, nobs_s1, nobs_s2, succeed)
  use mo_sensimul
  use mo_obs
  implicit none
  ! arguments
  integer, intent(in) :: m                       ! # of observations
  real(kind=8), intent(out) :: yobs(m), syobs(m) ! observations and sigma
  integer, intent(out) :: nobs_s1, nobs_s2       ! #S1/S2 observations
  logical, intent(out) :: succeed
  ! local variables
  logical :: exist, io_succeed
  integer :: m_s1,m_s2
  character(len=256) :: fname

  !-- initialise return values
  succeed = .false.
  nobs_s1 = 0
  nobs_s2 = 0
  yobs(1:m)  = sim_fill_value
  syobs(1:m) = sim_fill_value

  !-- get sizes of simulation-vectors
  m_s1 = get_m_s1()
  m_s2 = get_m_s2()

  if( .not. obs_are_set ) then

     write(*, '(a)') ' INFO::getobs:initial loading of observations...'

     !-- loading procedure is now entered...
     obs_are_set = .true.

     !-- S1 related observations
     fname = obs_file_s1
     inquire(FILE=fname, exist=exist)
     if( .not.exist ) then
        write(*,'(a)') ' INFO::getobs:no S1 observations ***'//trim(fname)//'*** available.'
     else
        write(*,'(a)') ' INFO::getobs:reading S1 observations from ***'//trim(fname)//'***'
        allocate( obs_vector_s1(m_s1), obsunc_vector_s1(m_s1) )
        call ncrd_obs_vector(fname, 's1', m_s1, obs_vector_s1, obsunc_vector_s1, io_succeed)
        if( .not. io_succeed ) then
           write(*, '(a)') ' FATAL::getobs:reading S1 observatsions from file '//&
                '***'//trim(fname)//'*** failed.'
           return
        endif
        write(*, '(a)') ' INFO::getobs:...loading S1 observations done.'
     endif

     !-- S2 related observations
     fname = obs_file_s2
     inquire(FILE=fname, exist=exist)
     if( .not.exist ) then
        write(*,'(a)') ' INFO::getobs:no S2 observations ***'//trim(fname)//'*** available.'
     else
        write(*,'(a)') ' INFO::getobs:reading S2 observations from ***'//trim(fname)//'***'
        allocate( obs_vector_s2(m_s2), obsunc_vector_s2(m_s2) )
        call ncrd_obs_vector(fname, 's2', m_s2, obs_vector_s2, obsunc_vector_s2, io_succeed)
        if( .not. io_succeed ) then
           write(*, '(a)') ' FATAL::getobs:reading S2 observatsions from file '//&
                '***'//trim(fname)//'*** failed.'
           return
        endif
        write(*, '(a)') ' INFO::getobs:...loading observations done.'
     end if
     write(*, '(a)') ' INFO::getobs:...initial loading done'

  endif

  !-- map observations into 1D simulation vector
  nobs_s1 = get_nobs_s1()
  if( nobs_s1.gt.0 ) then
     yobs( 1     :m_s1     ) = obs_vector_s1
     syobs(1     :m_s1     ) = obsunc_vector_s1
  endif
  nobs_s2 = get_nobs_s2()
  if( nobs_s2.gt.0 ) then
     yobs( m_s1+1:m_s1+m_s2) = obs_vector_s2
     syobs(m_s1+1:m_s1+m_s2) = obsunc_vector_s2
  endif

  !-- successful return
  succeed = .true.
end subroutine getobs


!***********************************************************
!     simulation_failure_filter
!
!> @brief Postprocessing method to filter out specific polarisations
!>        or particular bands from the simulation.
!
!> \authors MV/TK, The Inversion Lab
!> \date    August 2018
subroutine simulation_failure_filter(m, y)
  use mo_sensimul
  use mo_sensimul_s2,only:nb_s2
  implicit none
  ! arguments
  integer, intent(in)      :: m        !< length of output vector
  real(kind=8), intent(inout) :: y(m)  !< output vector
  ! local
  integer :: m_s1,m_s2
  integer :: i,ib,i_s1,i_s2,simtyp
  logical :: s1_match,s2_match

  !-- get simulation sizes
  m_s1 = get_m_s1()
  m_s2 = get_m_s2()

  !-- S1
  if( npts_s1.gt.0 ) then
     do i_s1=1,npts_s1
        i = timept_idxs_s1(i_s1)  !-- index in full timepts array
        simtyp = timept_simtyp(i) !-- simulation type for this index
        !-- S1A and S1B failure
        if( s1_failure_satellite.eq.'S1' ) then !-- any S1 observation
           s1_match = .true.
        !-- S1A or S1B failure
        else
           !-- S1A failure
           s1_match = s1_failure_satellite.eq.'S1A'.and.iand(simtyp,2**2).eq.2**2
           !-- S2B failure
           s1_match = s1_match.or.(s1_failure_satellite.eq.'S1B'.and.iand(simtyp,2**3).eq.2**3)
        endif
        if( s1_match ) then
           if( s1_failure_pol.eq.'vh' ) then
              y((i_s1-1)*2+1) = sim_fill_value !-- 'VH' is first polarisation at each time-point
           else if( s1_failure_pol.eq.'vv' ) then
              y((i_s1-1)*2+2) = sim_fill_value !-- 'VV' is first polarisation at each time-point
           endif
        endif       
     enddo
  end if
  
  !-- S2
  if( npts_s2.gt.0 ) then
     do i_s2=1,npts_s2
        i = timept_idxs_s2(i_s2)  !-- index in full timepts array
        simtyp = timept_simtyp(i) !-- simulation type for this index
        !-- S2A *and* S2B failure
        if( s2_failure_satellite.eq.'S2' ) then !-- any S2 observation
           s2_match = .true.
        !-- S2A *or* S2B failure
        else
           !-- S2A
           s2_match = s2_failure_satellite.eq.'S2A'.and.iand(simtyp,2**4).eq.2**4
           !-- S2B
           s2_match = s2_match.or.(s2_failure_satellite.eq.'S2B'.and.iand(simtyp,2**5).eq.2**5)
        endif
        if( s2_match ) then
           do ib=1,nb_s2
              if( s2_failure_bands(ib) ) then
                 y(m_s1+(i_s2-1)*nb_s2+ib) = sim_fill_value
              endif
           enddo
        endif
     enddo
  endif

end subroutine simulation_failure_filter


subroutine simulation_synthetic_unc(m, y, yunc)
  use mo_sensimul
  implicit none
  ! arguments
  integer, intent(in)      :: m        !< length of simulation vector
  real(kind=8), intent(in) :: y(m)     !< simulation vector
  real(kind=8), intent(out) :: yunc(m) !< synthetic simulation uncertainty vector
  ! local
  integer :: m_s1
  real(kind=8) :: s1unc_floor !-- in linear units
  real(kind=8) :: s1runc_lu !-- absolute uncertainty in [dB] yields relative uncertainty
                            !   in linear units

  !-- get simulation sizes
  m_s1 = get_m_s1()

  !-- Microwave: equations used to derive backscatter uncertainty:
  !   for XX in [VH,VV] we have
  !   XX_db(XX)  = 10*log10(XX)    !-- backscatter [dB]
  !   XX(XX_db)  = 10**(XX_db/10)  !-- backscatter, linear unit
  !
  !   sXX_db: a constant uncertainty in [dB] units with value 's1unc_db' is assumed
  !
  !   sXX is determined conservatively by:
  !   2*sXX = [ XX(XX_db+sXX_db) - XX(XX_db-sXX_db) ]
  !         = XX(XX_db) * [ 10**(sXX_db/10.) - 10**(-sXX_db/10.) ]
  !         = XX * [ 10**(sXX_db/10.) - 10**(-sXX_db/10.) ]
  !         = XX * [ 10**(sXX_db/10.) - 1./10**(sXX_db/10.) ]
  !
  !   NOTE: constant uncertainty in [dB] turns into relative uncertainty
  !         in linear units!
  !
  !   In Fortan code below
  !   XX(XX_db) are just the backscatter values we have computed (linear units),
  !   and we are using s1runc_lu := [ 10**(sXX_db/10.) - 1./10**(sXX_db/10.) ]

  !-- S1
  if( npts_s1.gt.0 ) then
     write(*,'(a,1(a,f5.2,1x))') ' INFO::simulation_synthetic_unc:'//&
          'S1 synthetic observational unceratainties will be generated with:',&
          'uncertainty[dB]=',s1unc_db

     !-- backscatter uncertainty [dB] converted as backscatter value
     s1runc_lu = 0.5_8 * ( 10**(s1unc_db/10._8) -1._8/(10**(s1unc_db/10._8)) )

     !-- constant uncertainty in [dB] unit turns into relative uncertainty in linear units
     yunc(1:m_s1) = y(1:m_s1) * s1runc_lu
     !-- apply floor value
     s1unc_floor  = 10**(s1floor_db/10._8) * s1runc_lu
     yunc(1:m_s1) = max(yunc(1:m_s1), s1unc_floor)

     !-- discard where 'y' is alread fill value
     where(y(1:m_s1).eq.sim_fill_value)
        yunc(1:m_s1) = sim_fill_value
     end where
  endif

  !-- S2
  if( npts_s2.gt.0 ) then
     write(*,'(a,2(a,f10.5,1x))') ' INFO::simulation_synthetic_unc:'//&
          'S2 synthetic observational unceratainties will be generated with:',&
          'relunc=',syn_s2unc_rel,'unc_floor=',syn_s2unc_floor
     yunc(m_s1+1:m) = max(y(m_s1+1:m)*syn_s2unc_rel, syn_s2unc_floor)
     where(y(m_s1+1:m).eq.sim_fill_value)
        yunc(m_s1+1:m) = sim_fill_value
     end where
  endif

end subroutine simulation_synthetic_unc


!***********************************************************
!     finishf
!
!> @brief Postprocessing method to write the simulation vector
!>        to binary and/or NetCDF file(s).
!
!> \authors MV/TK, The Inversion Lab
!> \date    February 2018
subroutine finishf(n, x, m, y)
  use mo_sensimul
  implicit none
  ! arguments
  integer, intent(in)      :: n     !< length of control vector
  integer, intent(in)      :: m     !< length of output vector
  real(kind=8), intent(in) :: x(n)  !< control vector
  real(kind=8), intent(in) :: y(m)  !< output vector
  ! local
  character(len=256) :: ncfname
  logical :: succeed
  integer :: m_s1,m_s2,iounit
  real(kind=8) :: xpr(n),sx(n),yunc(m)
  logical :: ldebug = .false.

  !-- get simulation sizes
  m_s1 = get_m_s1()
  m_s2 = get_m_s2()


  !-- binary output
  if( ldebug ) then
     iounit = 1
     call wrt_binary('sensimul.b', iounit, m, y)
     if( npts_s1.gt.0 ) then
        call wrt_binary('sensimul_s1.b', iounit, m_s1, y(1:m_s1))
     endif
     if( npts_s2.gt.0 ) then
        call wrt_binary('sensimul_s2.b', iounit, m_s2, y(m_s1+1:m_s1+m_s2))
     endif
  endif

  !-- for the synthetic data files we need to add uncertainties
  call simulation_synthetic_unc(m, y, yunc)

  !-- NetCDF output for S1 simulation
  if( npts_s1.gt.0 ) then
     ncfname = 'sensimul_s1.nc'
     call ncwrt_simulation_vector(ncfname, 's1', m_s1, y(1:m_s1), yunc(1:m_s1), succeed)
     if( succeed ) then
        write(*,'(a)') ' INFO::finishf:simulation for S1 domain has been written to file '//&
             '***'//trim(ncfname)//'***'
     else
        write(*,'(a)') ' ERROR::finishf:writing simulation for S1 domain failed'
     endif
  endif

  !-- NetCDF output for S2 simulation
  if( npts_s2.gt.0 ) then
     ncfname = 'sensimul_s2.nc'
     call ncwrt_simulation_vector(ncfname, 's2', m_s2, y(m_s1+1:m), yunc(m_s1+1:m), succeed)
     if( succeed ) then
        write(*,'(a)') ' INFO::finishf:simulation for S2 domain has been written to file '//&
             '***'//trim(ncfname)//'***'
     else
        write(*,'(a)') ' ERROR::finishf:writing simulation for S2 domain failed'
     endif
  endif

  !-- combined output
  if( ldebug ) then
     call getprior(n, xpr, sx)
     ncfname = 'simulation.nc'
     call ncwrt_state_and_sim(ncfname, n, x, sx, m, y, succeed)
     if( succeed ) then
        write(*,'(a)') ' INFO::finishf:'//&
             'state and simulation vector have been written to file '//&
             '***'//trim(ncfname)//'***'
     else
        write(*,'(a)') ' ERROR::finishf:'//&
             'errors occurred while writing file ***'//trim(ncfname)//'***'
     endif
  endif

contains
  subroutine wrt_binary(fname, unit, dsiz, data)
    implicit none
    ! arguments
    character(len=*), intent(in) :: fname
    integer, intent(in) :: unit, dsiz
    real(kind=8), intent(in) :: data(dsiz)
    ! local decls
    integer :: iostat

    open (unit=unit, file=fname, form='unformatted', iostat=iostat)
    if( iostat.ne.0 ) then
       write(*,'(a)') ' ERROR::finishf:wrt_binary:failed to open '//&
            '***'//trim(fname)//'*** for writing.'
    else
       write(unit, iostat=iostat) data
       if( iostat.ne.0 ) then
          write(*,'(a)') ' ERROR::finishf:wrt_binary:failed to open '//&
               '***'//trim(fname)//'*** for writing.'
       else
          write(*,'(a)') ' INFO::finishf:wrt_binary:'//&
               'file ***'//trim(fname)//'*** has been written.'
       endif
       close(unit)
    endif

  end subroutine wrt_binary

end subroutine finishf


!***********************************************************
!     finishdf
!
!> @brief Postprocessing method to be called after retrieval or mba
!>        in order to write all prototype tool output files.
!
!> @param[in]   n  length of control vector for multiple simulations in MW and optical domain
!> @param[in]   x  control vector normalised by prior uncertainty
!                  (expected ordering is: S1 related parameter(s) followed by
!                   by state variables LAI,HC,SM per 'simulation point')
!> @param[in]   c  posterior covariance matrix 'Cx.b' with respect to normalised control
!                  vector units.
!
!> \authors MV/TK, The Inversion Lab
!> \date    February 2018
!
subroutine finishdf(n, x, cov, optiflag)
  ! arguments
  integer, intent(in)           :: n          !< length of control vector
  real(kind=8), intent(in)      :: x(n)       !< control vector (normalised, no physical units)
  real(kind=8), intent(in)      :: cov(n,n)   !< posterior covariance matrix
  logical, intent(in), optional :: optiflag   !< boolean flag characterising success of optimisaton
  ! local
  character(len=256) :: ncfname
  real(kind=8) :: xphys(n), sx(n), xpr(n), sxpr(n)
  logical :: succeed
  integer :: i
  ! interface
  interface
     logical function ncwrt_controlvector(fname, n, x, sx, optiflag) result(succeed)
       ! arguments
       character(len=*), intent(in)  :: fname
       integer, intent(in)           :: n
       real(kind=8), intent(in)      :: x(n)
       real(kind=8), intent(in)      :: sx(n)
       logical, intent(in), optional :: optiflag
     end function ncwrt_controlvector
  end interface

  !-- get the prior
  call getprior(n, xpr, sxpr)

  !-- compute posterior uncertainties
  do i=1,n
     sx(i)    = sxpr(i)*sqrt(cov(i,i))
     xphys(i) = sxpr(i)*x(i)
  enddo

  !-- write prior control to NetCDF file
  ncfname = 'controlvector_prior.nc'
  succeed = ncwrt_controlvector(ncfname, n, xpr, sxpr)
  if( succeed ) then
     write(*, '(a)') ' INFO::finishdf:'//&
          'target file ***'//trim(ncfname)//'*** successfully written.'
  else
     write(*, '(a)') ' ERROR::finishdf:'//&
          'problems occurred while writing target file '//&
          '***'//trim(ncfname)//'***.'
  endif

  !-- write posterior control to NetCDF file
  !-- NOTE: using *physical* units in output
  ncfname = 'controlvector_post.nc'
  succeed = ncwrt_controlvector(ncfname, n, x*sxpr, sx, optiflag)
  if( succeed ) then
     write(*, '(a)') ' INFO::finishdf:'//&
          'target file ***'//trim(ncfname)//'*** successfully written.'
  else
     write(*, '(a)') ' ERROR::finishdf:'//&
          'problems occurred while writing target file '//&
          '***'//trim(ncfname)//'***.'
  endif

end subroutine finishdf


! We make doxygen ignore the following routine
!> \cond
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
!> \endcond


! We make doxygen ignore the following routine
!> \cond
subroutine printua(n, c, debug)
  implicit none
  ! arguments
  integer, intent(in) :: n
  real(kind=8), intent(in) :: c(n,n)
  logical, intent(in) :: debug
  ! local decls
  real(kind=8) :: ccopy(n,n)
  integer :: i,k,l

  ccopy = c

  do l = 1, n
     do k = 1, n
        ccopy(k,l) = ccopy(k,l)/sqrt(ccopy(k,k))/sqrt(ccopy(l,l))
     enddo
  enddo

  if( debug ) then
     write ( *, '(a)' ) ' '
     write ( *, '(a40)' ) 'Correlations of Posterior Uncertainty'
     write ( *, '(9x,18a5)' ) '1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17', '18'
     do i = 1, n
        print '(i9x,18(f4.2x))', i, ccopy(i,:)
     enddo
  endif

  write ( *, '(a)' ) ' '
  write ( *, '(a40)' ) 'Extreme off diagonal correlations'
  write ( *, '(a20x,a20,2a5)' )   'Min/Max ', 'value', 'i', 'j'
  do i = 1, n
     ccopy(i,i) = 0.
  enddo
  write ( *, '(a20x,f20.5,2i5)' )   'Min ', minval(ccopy), minloc(ccopy)
  write ( *, '(a20x,f20.5,2i5)' )   'Max ', maxval(ccopy), maxloc(ccopy)
  write ( *, '(a)' ) ' '



end subroutine printua
!> \endcond
