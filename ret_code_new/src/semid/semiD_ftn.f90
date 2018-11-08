!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!> \file semiD_ftn.f90
!
!> PROJECT : Sentinel Synergy Study
!
!> DESCRIPTION: main routine of the coupled Semidiscrite-Prospect-Price model provided by
!               Tristan Quaife and implemented as mixed Fortran/C code being ported
!               to pure Fortran implementation.
!               This aims for smooth application of an AD tool and may also yield small performance
!               gain.
!
!> \authors The Inversion Lab (Michael Vossbeck) 
!
!> \date  February/April/August 2018
!                               
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!***********************************************************
!     nadim_prospect_price_fast_1geom
!
!> @brief runs coupled NADIM-PROSPECT-PRICE model and computes  BRFs for 'nb' (broad) bands and a single illumination/view geometry
!
!> @details TBD
!
!> @param[in]  nb           number of wave-bands
!> @param[in]  nw1nm        number of 1nm wavelengths in range [400nm,2500nm]
!> @param[in]  resp_mat     responses for every wavelength in [400nm,2500nm] and all 'nb' bands
!> @param[in]  theta_i      illumination zenith angle [deg]
!> @param[in]  phi_i        illumination azimuth angle [deg]
!> @param[in]  theta_v      viewing zenith angle [deg]
!> @param[in]  phi_v        viewing azimuth angle [deg]
!> @param[in]  lad          Leaf angle distribution (1:Planophile, 2:Erectophile, 3:Plagiophile, 4:Extremophile, 5:Uniforme
!> @param[in]  rpl          radius of single leaf
!> @param[in]  lai          leaf area index [m^2/m^2]
!> @param[in]  hc           canopy height [m]
!> @param[in]  vai          leaf structure (PROSPECT)
!> @param[in]  cab          leaf chlorophyll (PROSPECT)
!> @param[in]  cw           leaf water equivelent thickness (PROSPECT)
!> @param[in]  cp           protein concentration (PROSPECT)
!> @param[in]  cc           cellulose and lignin (PROSPECT)
!> @param[in]  rsl1         weight of first spectral vector of the soil reflectance (PRICE)
!> @param[in]  rsl2         weight of second spectral vector of the soil reflectance (PRICE)
!> @param[in]  rsl3         weight of third spectral vector of the soil reflectance (PRICE)
!> @param[in]  rsl4         weight of fourth spectral vector of the soil reflectance (PRICE)
!> @param[out] brf          Bidirectional reflectance factor
!
subroutine nadim_prospect_price_fast_1geom( nb, nw1nm, resp_mat, &
     theta_i, phi_i, theta_v, phi_v, &
     lad, rpl, lai, hc, &
     vai, cab, cw, cp, cc, &
     rsl1, rsl2, rsl3, rsl4, &
     brf )
  implicit none
  ! constants
  integer, parameter :: niv = 1 !number of view--illumination geometries
  ! arguments
  integer, intent(in) :: nb, nw1nm
  real, intent(in) :: theta_i, phi_i, theta_v, phi_v
  integer, intent(in) :: lad
  real, intent(in) :: rpl, lai, hc
  real(kind=8), intent(in) :: vai, cab, cw, cp, cc
  real, intent(in) :: rsl1, rsl2, rsl3, rsl4
  real(kind=8), intent(in) :: resp_mat(nw1nm, nb)
  real, intent(out) :: brf(nb)
  ! externals
  external price_soil_fullSpectrum_interp1nm_ftn
  external prospect_fullSpectrum_interp1nm_ftn
  real(kind=8), external :: convolve
  ! local variables
  integer :: ib
  real(kind=8) :: p_rl(nw1nm), p_tl(nw1nm) ! leaf reflectance/transmittance
  real(kind=8) :: p_rs(nw1nm)              ! soil reflectance
  real(kind=8) :: rl_8(nb), tl_8(nb), rs_8(nb)         ! values for individual S2 bands
  real :: rl(nb), tl(nb), rs(nb)                       ! ...as normal 'real'
  real :: theta_v_asarr(niv), phi_v_asarr(niv)
  real :: brf_nadim(1)

  ! ! D E B U G
  ! integer, save :: cnt = 0
  ! character(len=1) :: cbuf
  ! character(len=2) :: ccbuf
  ! cnt = cnt + 1
  ! write(cbuf, '(i1)') cnt

  !-- scalar to 1D array
  theta_v_asarr(1) = theta_v
  phi_v_asarr(1)   = phi_v

  !-- get leaf optical properties: [400nm,2500nm] at 1nm resolution
  call prospect_fullSpectrum_interp1nm_ftn(vai, cab, cw, cp, cc, p_rl, p_tl)

  ! get soil reflectance: [400nm,2500nm] at 1nm resolution
  call price_soil_fullSpectrum_interp1nm_ftn(rsl1, rsl2, rsl3, rsl4, p_rs)



  !-- loop over S2 bands
!!! !$AD II-LOOP
  bndloop:do ib=1,nb
     rs_8(ib) = convolve(nw1nm, p_rs, resp_mat(1:nw1nm,ib))
     rl_8(ib) = convolve(nw1nm, p_rl, resp_mat(1:nw1nm,ib))
     tl_8(ib) = convolve(nw1nm, p_tl, resp_mat(1:nw1nm,ib))
     rs(ib) = rs_8(ib)
     rl(ib) = rl_8(ib)
     tl(ib) = tl_8(ib)

     !-- BRF computation by NADIM
     call nadimbrf( theta_i, phi_i, niv, theta_v_asarr, phi_v_asarr, &
          lad, rs(ib), hc, lai, rpl, rl(ib), tl(ib), &
          brf_nadim )

     !-- save NADIM BRF
     brf(ib) = brf_nadim(1)
  enddo bndloop

end subroutine nadim_prospect_price_fast_1geom


!***********************************************************
!     nadime_prospect_price_full_1geom
!
!> @brief runs the extended coupled NADIM-PROSPECT-PRICE model for every wavelength at 1nm resolution and derives the bihemisphic reflectance factor (BRF), the fraction absorbed in the canopy (FPAR), the albedo, and the transmission.
!
!> @details TBD
!
!> @param[in]  nw1nm        number of 1nm wavelengths in range [400nm,2500nm]
!> @param[in]  nw           number of wave-bands actuall used (starting from 400nm)
!> @param[in]  resp_mat     responses for every wavelength in [400nm,2500nm] and all 'nb' bands
!> @param[in]  theta_i      illumination zenith angle [deg]
!> @param[in]  phi_i        illumination azimuth angle [deg]
!> @param[in]  theta_v      viewing zenith angle [deg]
!> @param[in]  phi_v        viewing azimuth angle [deg]
!> @param[in]  lad          Leaf angle distribution (1:Planophile, 2:Erectophile, 3:Plagiophile, 4:Extremophile, 5:Uniforme
!> @param[in]  rpl          radius of single leaf
!> @param[in]  lai          leaf area index [m^2/m^2]
!> @param[in]  hc           canopy height [m]
!> @param[in]  vai          leaf structure (PROSPECT)
!> @param[in]  cab          leaf chlorophyll (PROSPECT)
!> @param[in]  cw           leaf water equivelent thickness (PROSPECT)
!> @param[in]  cp           protein concentration (PROSPECT)
!> @param[in]  cc           cellulose and lignin (PROSPECT)
!> @param[in]  rsl1         weight of first spectral vector of the soil reflectance (PRICE)
!> @param[in]  rsl2         weight of second spectral vector of the soil reflectance (PRICE)
!> @param[in]  rsl3         weight of third spectral vector of the soil reflectance (PRICE)
!> @param[in]  rsl4         weight of fourth spectral vector of the soil reflectance (PRICE)
!> @param[out] brf          Bidirectional reflectance factor
!
subroutine nadime_prospect_price_full_1geom( nw1nm, nw, &
     theta_i, phi_i, theta_v, phi_v, &
     lad, rpl, lai, hc, &
     vai, cab, cw, cp, cc, &
     rsl1, rsl2, rsl3, rsl4, &
     brf, fpar, albedo, trans )
  implicit none
  ! constants
  integer, parameter :: niv = 1 !number of view--illumination geometries
  ! arguments
  integer, intent(in) :: nw    !-- actual  number of 1nm wave-lengths [400nm,400nm+nw-1]
  integer, intent(in) :: nw1nm !-- maximal number of 1nm wave-lenghts [400nm,2500nm]
  real, intent(in) :: theta_i, phi_i, theta_v, phi_v
  integer, intent(in) :: lad
  real, intent(in) :: rpl, lai, hc
  real(kind=8), intent(in) :: vai, cab, cw, cp, cc
  real, intent(in) :: rsl1, rsl2, rsl3, rsl4
  real, intent(out) :: brf(nw), fpar(nw), albedo(nw), trans(nw)
  ! externals
  external price_soil_fullSpectrum_interp1nm_ftn
  external prospect_fullSpectrum_interp1nm_ftn
  ! local variables
  integer :: iw
  real(kind=8) :: p_rl(nw1nm), p_tl(nw1nm) ! leaf reflectance/transmittance
  real(kind=8) :: p_rs(nw1nm)              ! soil reflectance
  real :: rl, tl, rs                       ! ...as normal 'real'
  !-- Note: single NADIM call is for single illumination/view geometry
  real :: theta_v_asarr(niv), phi_v_asarr(niv)
  real :: brf_nadim(niv), fpar_nadim(niv), albedo_nadim(niv), trans_nadim(niv)

  !-- scalar to 1D array
  theta_v_asarr(1) = theta_v
  phi_v_asarr(1)   = phi_v

  !-- get leaf optical properties: [400nm,2500nm] at 1nm resolution
  call prospect_fullSpectrum_interp1nm_ftn(vai, cab, cw, cp, cc, p_rl, p_tl)

  ! get soil reflectance: [400nm,2500nm] at 1nm resolution
  call price_soil_fullSpectrum_interp1nm_ftn(rsl1, rsl2, rsl3, rsl4, p_rs)


  !-- loop over all wavelengthts
  bndloop:do iw=1,nw

     !-- set soil-reflectance, leaf reflectance/transmittance
     rs = p_rs(iw)
     rl = p_rl(iw)
     tl = p_tl(iw)

     !-- BRF computation by NADIM
     call nadimbrfe( theta_i, phi_i, niv, theta_v_asarr, phi_v_asarr, &
          lad, rs, hc, lai, rpl, rl, tl, &
          brf_nadim, fpar_nadim, albedo_nadim, trans_nadim )

     !-- save NADIM output
     brf(iw)    = brf_nadim(1)
     fpar(iw)   = fpar_nadim(1)
     albedo(iw) = albedo_nadim(1)
     trans(iw)  = trans_nadim(1)
  enddo bndloop

end subroutine nadime_prospect_price_full_1geom


!***********************************************************
!     convolve
!
!> @brief convolves to spectra
!
!> @details TBD
!
!> @param[in]  nw      length/size of spectra
!> @param[in]  spect1  first spectrum
!> @param[in]  spect2  second spectrum
!
real(kind=8) function convolve(nw, spect1, spect2)
  implicit none
  ! arguments
  integer, intent(in) :: nw
  real(kind=8), intent(in) :: spect1(nw), spect2(nw)
  ! local decls
  real(kind=8) :: spect1_dot_spect2, spect2_sum

  spect2_sum = sum(spect2)
  spect1_dot_spect2 = sum(spect1*spect2)
  convolve = spect1_dot_spect2/spect2_sum

end function convolve


!***********************************************************
!     sm_to_rsl1
!
!> @brief modulates the weight of first empirical orthogonal function of soil reflectance model (Price 1990) by the soil moisture
!
!> @details TBD
!
!> @param[in]  sm        volumetric soilmoisture (in range [0,1])
!> @param[in]  rsl1      1st spectral vector of soil reflectance (of Price's EOF)
!> @param[in]  sm_coeff  weighting of soil moisture impact, bound between (0,1)
!> \return  weight of first empirical orthogonal function
!
real function sm_to_rsl1(sm, rsl1, sm_coeff)
  implicit none
  real, intent(in) :: sm
  real, intent(in) :: rsl1
  real, intent(in) :: sm_coeff
  sm_to_rsl1 = rsl1*(1.-sm_coeff*sm) ! SM measured volumetric,but not in percent (sm/100))
end function sm_to_rsl1
