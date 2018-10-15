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

  !-- n: number of 'points' to be multiplied with number of states per point
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
  integer :: i

  !-- get simulation sizes
  m_s1 = get_m_s1()
  m_s2 = get_m_s2()

  !-- S1
  if( npts_s1.gt.0 ) then
     if( s1_failure_pol.eq.'vh' ) then
        y(1:m_s1:2) = sim_fill_value
     else if( s1_failure_pol.eq.'vv' ) then
        y(2:m_s1:2) = sim_fill_value
     endif
  end if
  
  if( npts_s2.gt.0 ) then
     filterloop:do i=1,nb_s2
        if( s2_failure_bands(i) ) then
           y(m_s1+i:m:nb_s2) = sim_fill_value
        else
           cycle filterloop
        endif
     enddo filterloop
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
  real(kind=8) :: s1unc_factor

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
  !   and we are using s1unc_factor := 10**(sXX_db/10.)

  !-- S1
  if( npts_s1.gt.0 ) then
     write(*,'(a,1(a,f5.2,1x))') ' INFO::simulation_synthetic_unc:'//&
          'S1 synthetic observational unceratainties will be generated with:',&
          'uncertainty[dB]=',s1unc_db

     !-- backscatter uncertainty [dB] converted as backscatter value
     s1unc_factor = 10** ( s1unc_db/10._8)

     !-- constant uncertainty in [dB] unit turns into relative uncertainty in linear units
     yunc(1:m_s1) = y(1:m_s1) * 0.5_8 * (s1unc_factor - 1._8/s1unc_factor)
 
  endif

  !-- S2
  if( npts_s2.gt.0 ) then
     write(*,'(a,2(a,f10.5,1x))') ' INFO::simulation_synthetic_unc:'//&
          'S2 synthetic observational unceratainties will be generated with:',&
          'relunc=',syn_s2unc_rel,'unc_floor=',syn_s2unc_floor
     yunc(m_s1+1:m) = max(y(m_s1+1:m)*syn_s2unc_rel, syn_s2unc_floor)
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

  !-- get simulation sizes
  m_s1 = get_m_s1()
  m_s2 = get_m_s2()


  !-- binary output
  iounit = 1
  call wrt_binary('sensimul.b', iounit, m, y)
  if( npts_s1.gt.0 ) then
     call wrt_binary('sensimul_s1.b', iounit, m_s1, y(1:m_s1))
  endif
  if( npts_s2.gt.0 ) then
     call wrt_binary('sensimul_s2.b', iounit, m_s2, y(m_s1+1:m_s1+m_s2))
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
!> @brief Postprocessing method to be called after mba
!>        to binary and/or NetCDF file(s).
!
!> \authors MV/TK, The Inversion Lab
!> \date    February 2018
!
subroutine finishdf(n, x, cov)
  use mo_sensimul
  implicit none
  ! arguments
  integer, intent(in)      :: n         !< length of control vector
  real(kind=8), intent(in) :: x(n)      !< control vector (normalised, no physical units)
  real(kind=8), intent(in) :: cov(n,n)  !< posterior covariance matrix
  ! local
  character(len=256) :: ncfname
  real(kind=8) :: xphys(n), sx(n), xpr(n), sxpr(n)
  logical :: succeed
  integer :: i

  !-- get the prior
  call getprior(n, xpr, sxpr)

  !-- compute posterior uncertainties
  do i=1,n
     sx(i)    = sxpr(i)*sqrt(cov(i,i))
     xphys(i) = sxpr(i)*x(i)
  enddo

  !-- write prior control to NetCDF file
  ncfname = 'controlvector_prior.nc'
  call ncwrt_controlvector(ncfname, n, xpr, sxpr, succeed)
  if( succeed ) then
     write(*, '(a)') ' INFO::finishdf:'//&
          'target file ***'//trim(ncfname)//'*** successfully written.'
  else
     write(*, '(a)') ' ERROR::finishdf:'//&
          'problems occurred while writing target file '//&
          '***'//trim(ncfname)//'***.'
  endif

  !-- write posterior control to NetCDF file
  ncfname = 'controlvector_post.nc'
  call ncwrt_controlvector(ncfname, n, x*sxpr, sx, succeed)
  if( succeed ) then
     write(*, '(a)') ' INFO::finishdf:'//&
          'target file ***'//trim(ncfname)//'*** successfully written.'
  else
     write(*, '(a)') ' ERROR::finishdf:'//&
          'problems occurred while writing target file '//&
          '***'//trim(ncfname)//'***.'
  endif

end subroutine finishdf
