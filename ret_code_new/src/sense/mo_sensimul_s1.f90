!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!> \file mo_sensimul_s1.f90
!
!> PROJECT : Sentinel Synergy Study
!
!> DESCRIPTION: provision of variables used to run the microwave and/or optical model of the Sentinel Simulator (which are not regarded as "active" in the context of the retrieval tool.
!
!> \authors The Inversion Lab (Michael Vossbeck) 
!
!> \date  February 2018
!                               
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module mo_sensimul_s1
  implicit none

  ! number of polarisations for back-scatter values
  integer, parameter :: npol = 2   ! vh,vv

  !-- miscellaneous
  real(kind=8), parameter :: pi = 3.1415926535879
  !-- fill-value
  real(kind=8), parameter :: sense_fv = -99999. ! python code of SentinelSimulator uses
                                                ! NaN to indicate not computed sigma

  !-- dielectric soil model
  character(len=*), parameter :: soilm_dielectric = 'dobson85'
  !-- surface model
  character(len=*), parameter :: rt_surface = 'Oh92'
  !-- canopy model
  character(len=*), parameter :: rt_canopy  = 'turbid_rayleigh'
  !-- whether backscatter is computed in coherent/non-coherent manner
  !   Ref. sense/model.py, l76
  logical, parameter      :: coherent = .true.


  !-- S1 frequency
  real(kind=8) :: freq = 5.405_8
  !-- S1 incidence angle (SENSE package operates internally operates in [rad])
  real(kind=8), parameter :: theta_deg_wallerfing = 37._8
  real(kind=8) :: theta = (theta_deg_wallerfing*pi)/180._8


  !-- site-specific settings (by default the Wallering settings are applied.)
  !-- sand content (fractional volume)
  real(kind=8) :: sand = 0.27_8
  !-- clay content (fractional volume)
  real(kind=8) :: clay = 0.23_8
  !-- bulk density [g cm-3]
  real(kind=8) :: bulk = 1.65_8
  !-- surface rms height [m]
  real(kind=8) :: sfc_rms_height = 0.015_8
  !-- single scattering albedo (calculation of volume scattering coefficients)
  real(kind=8) :: omega = 0.1_8


  !-------------------
  !     run-time values
  !     (ATTENTION: should not become active in terms of AD)
  !
  !-- Dobson85 specific?
  real(kind=8) :: alpha, beta1, beta2, sigma
  !-- surface roughness
  real(kind=8) :: ks
  !-- permittivity of free water
  complex(kind=8) :: ew

contains
  !***********************************************************
  !     sense_init
  !
  !> @brief initialisation of non-active parameter(s) to run the SENSE RT model
  !
  !> @details comprises computation of
  !           - surface roughness (for fixed frequency and surface rms height)
  !           - dielectric permittivity of free water
  !
  subroutine sense_init()
    implicit none
    external soil_ks, dobson85_init, dobson85_ew

    !-- surface roughness
    call soil_ks(sfc_rms_height, freq, ks)

    !-- dielectric permittivity of free water
    if( soilm_dielectric.eq.'dobson85' ) then
       call dobson85_init(sand, clay, bulk, alpha, beta1, beta2, sigma)
       call dobson85_ew(freq, sigma, ew)
    else
       write(*,'(a)') ' FATAL::sense_init: dielectric soil model -->'//trim(soilm_dielectric)//&
            '<-- not (yet) supported.'
       stop 1
    endif

  end subroutine sense_init

end module mo_sensimul_s1
