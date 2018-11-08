!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!> \file mo_sensimul.f90
!> \brief provides environment to run the Fortran implementation 
!>        of the Sentinel Simulator within the prototype retrieval system
!> \authors The Inversion Lab (Michael Vossbeck, Thomas Kaminski) 
!> \date  February 2018
!                               
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  MODULE mo_sensimul
module mo_sensimul
!> \brief provides environment to run the Fortran implementation 
!>        of the Sentinel Simulator in the assimilation system.
!> \authors The Inversion Lab (Michael Vossbeck, Thomas Kaminski) 
!> \date  February-August 2018
  implicit none

  private :: cnt_lines


  !-- constants
  integer, parameter :: nsc   = 3 !number of state-vector components (per 'time-point')
                                  !(LAI,HC,SM) for now

  integer, parameter :: nparam_s1 = 1      ! #parameters involved in S1 simulation
                                           !  (lai-coefficient relating optical to MW LAI)
  integer, parameter :: nparam_s2 = 0      ! #parameters involved in S2 simulation


  !-- simulation/retrieval configuration settings
  !   (to be determined by reading external file 'retrconfig.nc')
  integer :: npts    = -1
  integer :: npts_s1 = -1
  integer :: npts_s2 = -1
  real(kind=8), allocatable :: iv_geom(:,:)       !     4 x npts     (sza,saa,vza,vaa)
  real(kind=8), allocatable :: timepts(:)         !         npts
  integer, allocatable      :: timept_simtyp(:)   !         npts     (0:other, bit0:S1, bit1:S2)
  integer, allocatable      :: timept_idxs_s1(:)  !         npts_s1  !position(s) of S1 in schedule
  integer, allocatable      :: timept_idxs_s2(:)  !         npts_s2  !position(s) of S2 in schedule
  character(256)            :: time_unit = ''

  !-- external files to set-up and initialise simulation/retrieval system
  !-- files below should have been made available by the rs_setup.py utility
  character(len=256), parameter :: s2_srf_file      = 'input/s2a.srf'
  character(len=256), parameter :: solar_spect_file = 'input/ASTMG173.csv'
  !-- NOTE:files below *must* reside in the current working directory,
  !--      and are typically created by invoking rs_pre.py (preprocessing tool)
  character(len=256), parameter :: config_file = 'retrconfig.nc'
  character(len=256), parameter :: prior_file  = 'retrprior.nc'
  character(len=256), parameter :: model_file  = 'retrmodel.nc'
  character(len=256), parameter :: site_nml    = 'site.nml'
  character(len=256), parameter :: obs_file_s1 = 'obs_s1.nc'
  character(len=256), parameter :: obs_file_s2 = 'obs_s2.nc'

  !--
  real(kind=8), parameter :: sim_fill_value = -99999._8

  !-- p r e p a r a t i o n   o f   s y n t h e t i c   d a t a
  !-- S1: we apply an uncertainty of 0.4dB (as suggested by BR, ref. PM7 notes)
  real(kind=8), parameter :: s1unc_db = 0.4_8
  !-- S1: we apply an uncertainty floor value (in linear units!)
  !       which corresponds to 0.4db of a value of -22db converted to linear units
  real(kind=8), parameter :: s1floor_db = -22._8
  !-- for S2 Philip did look-up a value for the relative uncertainty
  !   (which however is likely for the sensor itself and does not include
  !    errors from the atmospheric correction)
  real(kind=8), parameter :: syn_s2unc_rel   = 0.05_8
  !-- using 5% relative uncertainty of BRF=0.01 as floor value:
  real(kind=8), parameter :: syn_s2unc_floor = 0.01_8*syn_s2unc_rel
  !-- preparation of failure scenarios
  character(len=2)     :: s1_failure_pol = '' !-- of ['','vh','vv']
  character(len=3)     :: s1_failure_satellite = '' !-- 'S1','S1A','S1B'
  logical, allocatable :: s2_failure_bands(:) !-- nb_s2
  character(len=3)     :: s2_failure_satellite = '' !-- 'S2','S2A','S2B'


  !-----------------------------
  !    i n t e r f a c e s
  !
  interface
     subroutine finishdf(n, x, cov, optiflag)
       implicit none
       ! arguments
       integer, intent(in)           :: n
       real(kind=8), intent(in)      :: x(n)
       real(kind=8), intent(in)      :: cov(n,n)
       logical, intent(in), optional :: optiflag
     end subroutine finishdf
  end interface
contains

  subroutine simulator_init()
    use mo_sensimul_s1, only:sense_init
    use mo_sensimul_s2, only:semid_init,nb_s2
    implicit none

    !-- set filter
    s1_failure_pol = ''
    if( allocated(s2_failure_bands) ) then
       write(*, '(a)') ' FATAL::simulator_init:s2_failure_bands already allocated within initialisation!'
       stop
    else
       allocate(s2_failure_bands(nb_s2))
       s2_failure_bands = .false.
    endif

    !-- initialise MW operator
    call sense_init()

    !-- initialise operator for optical domain
    call semid_init(s2_srf_file, solar_spect_file)

  end subroutine simulator_init


  !****************************************
  !    get_npol
  !
  !> @brief returns number of Sentinel1 polarisations
  !
  !> \return  number of Sentinel1 polarisations
  !
  pure integer function get_npol()
    use mo_sensimul_s1, only:npol
    get_npol = npol
  end function get_npol


  !****************************************
  !    get_nbs2
  !
  !> @brief returns number of Sentinel2 bands
  !
  !> \return  number of Sentinel2 bands
  !
  pure integer function get_nbs2()
    use mo_sensimul_s2, only:nb_s2
    get_nbs2 = nb_s2
  end function get_nbs2


  !****************************************
  !    get_nc_s1
  !
  !> @brief returns (1D) overall length S1 part in control vector (including parameter(s))
  !
  !> \return  length of S1 part in control vector
  !
  pure integer function get_nc_s1()
    implicit none
    get_nc_s1 = nparam_s1+npts_s1*nsc
  end function get_nc_s1

  !****************************************
  !    get_nc_s2
  !
  !> @brief returns (1D) overall length S2 part in control vector
  !
  !> \return  length of S2 part in control vector
  !
  pure integer function get_nc_s2()
    implicit none
    get_nc_s2 = nparam_s2+npts_s2*nsc
  end function get_nc_s2

  !****************************************
  !    get_n
  !
  !> @brief returns (1D) overall length of control vector (including parameter(s))
  !
  !> \return  length of control vector
  !
  pure integer function get_n()
    implicit none
    get_n = nparam_s1+npts*nsc
  end function get_n

  !****************************************
  !    get_npts
  !
  !> @brief returns (1D) overall number of time-points
  !
  !> \return  number of time-points
  !
  pure integer function get_npts()
    implicit none
    get_npts = npts
  end function get_npts

  !****************************************
  !    get_ns
  !
  !> @brief returns (1D) overall length of state variables in control vector
  !
  !> \return  length of state vector
  !
  pure integer function get_ns()
    implicit none
    get_ns = npts*nsc
  end function get_ns

  !****************************************
  !    get_np
  !
  !> @brief returns (1D) overall length of parameter part in control vector
  !
  !> \return  length of parameter vector
  !
  pure integer function get_np()
    implicit none
    get_np = nparam_s1 + nparam_s2
  end function get_np

  !****************************************
  !    get_m_s1
  !
  !> @brief returns (1D) length of simulation vector in MW domain
  !
  !> \return  length of simulation vector
  !
  pure integer function get_m_s1()
    use mo_sensimul_s1, only:npol
    implicit none
    get_m_s1 = npol*npts_s1
  end function get_m_s1

  !****************************************
  !    get_m_s2
  !
  !> @brief returns (1D) length of simulation vector in optical domain
  !
  !> \return  length of simulation vector
  !
  pure integer function get_m_s2()
    use mo_sensimul_s2, only:nb_s2
    implicit none
    get_m_s2 = nb_s2*npts_s2
  end function get_m_s2

  !****************************************
  !    get_m
  !
  !> @brief returns (1D) overall length of simulation vector
  !
  !> \return  length of simulation vector
  !
  pure integer function get_m()
    use mo_sensimul_s1, only:npol
    implicit none
    get_m = get_m_s1() + get_m_s2()
  end function get_m

  !****************************************
  !    idx_is_s1
  !
  !> @brief returns whether timeindex 'i' references a S1 acquisition time
  !
  !> \return  flag
  !
  pure logical function idx_is_s1(i)
    implicit none
    integer, intent(in) :: i !< index in time-series
    idx_is_s1 = iand(timept_simtyp(i),1).eq.1
  end function idx_is_s1
  
  !****************************************
  !    idx_is_s2
  !
  !> @brief returns whether timeindex 'i' references a S2 acquisition time
  !
  !> \return  flag
  !
  pure logical function idx_is_s2(i)
    implicit none
    integer, intent(in) :: i !< index in time-series
    idx_is_s2 = iand(timept_simtyp(i),2).eq.2
  end function idx_is_s2

  !****************************************
  !    idx_is_other
  !
  !> @brief returns whether timeindex 'i' references an extra target time
  !
  !> \return  flag
  !
  pure logical function idx_is_other(i)
    implicit none
    integer, intent(in) :: i !< index in time-series
!    idx_is_other = .not.(idx_is_s1(i).or.idx_is_s2(i))
    idx_is_other = (timept_simtyp(i).eq.0)
  end function idx_is_other


  !****************************************
  !    sensimul_finalise
  !
  !> @brief cleanup stuff to be called after retrieval system has finished.
  !
  subroutine sensimul_finalise()
    implicit none
    if( allocated(iv_geom) )          deallocate(iv_geom)
    if( allocated(timepts) )          deallocate(timepts)
    if( allocated(timept_simtyp) )    deallocate(timept_simtyp)
    if( allocated(timept_idxs_s1) )   deallocate(timept_idxs_s1)
    if( allocated(timept_idxs_s2) )   deallocate(timept_idxs_s2)
    if( allocated(s2_failure_bands) ) deallocate(s2_failure_bands)
  end subroutine sensimul_finalise

  integer function cnt_lines(fname)
    implicit none
    character(len=*), intent(in) :: fname
    integer, parameter :: lunit = 18

    cnt_lines = 0
    open(lunit, file=fname)
    do
       read (lunit,*, end=10)
       cnt_lines = cnt_lines + 1
    end do
10  close (lunit)
  end function cnt_lines


end module mo_sensimul
