!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> \file target_operator.f90
!> \brief provides implementation of the target operators in the
!>        optical and microwave domain.
!> \authors The Inversion Lab (Michael Vossbeck, Thomas Kaminski) 
!> \date  August 2018
!                               
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!***********************************************************
!     allbackscatter_operator_1geom
!
!> @brief Implementation of the Sentinel Synergy Study target operator
!>        in the microwave domain (denoted H_1 in document D5) for a single state and
!>        incidence angle.
!
!> @details For the given input state (LAI_coeff, LAI, HC, and SM) the observation operator
!>          derives backscatter values in VH, VV and HH polarisation.
!>          The implementation is based on the Fortran port of the SENSE model.
!
!> @param[in]   statev   state components required to run S1 simulation
!>                       consisting of LAI conversion coefficient, (optical) LAI, HC, and SM
!> @param[in]   iv_geom  illumination-view geomtry in order sza,saa,vza,vaa (in degrees)
!> @param[out]  bscat    back-scatter values in VH/VV/HH polarisation
!
subroutine allbackscatter_operator_1geom( statev, iv_geom, bscat )
  use mo_sensimul, only:nsc,nparam_s1
  use mo_sensimul_s1, only:theta,pi !-- viewing zenith-angle for microwave simulation
  implicit none
  ! arguments
  real(kind=8), intent(in)  :: statev(nparam_s1+nsc)
  real(kind=8), intent(in)  :: iv_geom(4) !sza,saa,vza,vaa
  real(kind=8), intent(out) :: bscat(3)   !VH,VV,HH
  ! externals
  external sigma0_vhvvhh !core SENSE routine
  ! local declarations
  real(kind=8) :: lai_coeff, lai, hc, sm
  real(kind=8) :: s0vh, s0vv, s0hh
  real(kind=8) :: vza_deg  !-- S1 viewing zenith angle

  !-- map state variables
  lai_coeff = statev(1)
  lai       = statev(2)
  hc        = statev(3)
  sm        = statev(4)
  !--
  vza_deg   = iv_geom(3)
  !-- convert to RAD
  theta = vza_deg*(pi/180._8)

  !-- call SENSE core routine
  call sigma0_vhvvhh(lai, lai_coeff, hc, sm, s0vh, s0vv, s0hh)

  !-- map backscatter values
  bscat(1) = s0vh
  bscat(2) = s0vv
  bscat(3) = s0hh

end subroutine allbackscatter_operator_1geom


!***********************************************************
!     fapar_operator_1geom
!
!> @brief Implementation of the FAPAR target operator.
!
!> @details For the given input state (LAI_coeff, LAI, HC, and SM) the observation operator
!>          derives FAPAR, i.e. the Fraction of Absorbed Photosynthetically Active radiation in the canopy in the visible domain (400nm - 700nm)
!
!> @param[in]   statev   state components required to run S1 simulation
!>                       consisting of LAI conversion coefficient, (optical) LAI, HC, and SM
!> @param[in]   iv_geom  illumination-view geomtry in order sza,saa,vza,vaa (in degrees)
!> @param[out]  fapar    FAPAR
!
subroutine fapar_operator_1geom( statev, iv_geom, fapar )
  use mo_sensimul, only:nsc
  use mo_sensimul_s2
  ! arguments
  real(kind=8), intent(in)  :: statev(nsc)
  real(kind=8), intent(in)  :: iv_geom(4) !sza,saa,vza,vaa
  real(kind=8), intent(out) :: fapar
  ! local decls
  real :: lai, hc, sm
  real :: rsl1
  real :: theta_i, phi_i, theta_v, phi_v
  real :: fpar_nadim(nw1nm) !'nadime_prospect_price_fast_1geom' does not expect kind=8 !
  ! externals
  real, external :: sm_to_rsl1
  external nadime_prospect_price_full_1geom

  !-- map state variables
  !   (might potentially cast from r8 to r4)
  lai = statev(1)
  hc  = statev(2)
  sm  = statev(3)
  !-- 1st spectral vector of soil reflectance (of Price's EOF)
  !   is modulated by the soil moisture state
  rsl1 = sm_to_rsl1(sm, rsl1_default, sm_coeff)

  !-- geometry
  theta_i = iv_geom(1)
  phi_i   = iv_geom(2)
  theta_v = iv_geom(3)
  phi_v   = iv_geom(4)

  !-- call coupled optical model for full wavelength-spectrum
  call nadime_prospect_price_full_1geom(nw1nm, &
       theta_i, phi_i, theta_v, phi_v, &
       lad, rpl, lai, hc, &
       vai, cab, cw, cp, cc, &
       rsl1, rsl2, rsl3, rsl4, &
       fpar_nadim)

  !-- ! ! ! IMPL-MISSING::integration of 1nm fpar to FAPAR ! ! !
  write(*, '(a)') ' ERROR::fapar_operator_1geom:FAPAR computation incomplete, just setting -9999'
  fapar = -9999._8
end subroutine fapar_operator_1geom
