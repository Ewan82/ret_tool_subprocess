!***********************************************************
!     sigma0_hhvv
!
!> @brief compute total single-scattering backscattering coefficents for HH and VV polarisations
!
!> @details 
!
!> @param[in]  lai        leaf-area index [m2/m2] as used in optical domain
!> @param[in]  lai_coeff  conversion coefficient for lai optical to microwave
!> @param[in]  canht      canopy height [m]
!> @param[in]  mv         volumetric soilmoisture [m3/m3]
!> @param[out] s0hh       backscattering coefficient, HH polarisation
!> @param[out] s0vv       backscattering coefficient, VV polarisation
!
subroutine sigma0_hhvv(lai, lai_coeff, canht, mv, s0hh, s0vv)
  use mo_sensimul_s1
  implicit none
  ! arguments
  real(kind=8), intent(in) :: lai, lai_coeff, canht, mv
  real(kind=8), intent(out) :: s0hh, s0vv
  ! externals
  external dobson85_eps
  external reflectivity, calc_rho
  external canopy_extinction_coeffs, canopy_transmissivity
  external ScatRayleigh_sigma_v_back, ScatRayleigh_sigma_v_bistatic
  external Oh92_backscatter
  external canopy_sigma_c, ground_sigma, sigma_cg, sigma_gcg
  ! local declarations
  complex(kind=8) :: eps
  real(kind=8) :: v,h,rho_v,rho_h, t_v, t_h, ke_v, ke_h, ks_h, ks_v
  real(kind=8) :: lai_mw
  real(kind=8) :: soil_backscatter(3)                      !hh,vv,hv
  real(kind=8) :: sigma_vol_bistatic(3), sigma_vol_back(3) !hh,vv,hv
  real(kind=8) :: s0c(3), s0g(3), s0cgt(3), s0gcg(3)       !hh,vv,hv

  !-- microwave LAI
  lai_mw = lai*lai_coeff

  !-- soil-moisture ==> eps (dielectric permittivity)
  call dobson85_eps(mv, bulk, alpha, beta1, beta2, ew, eps)

  !-- reflectivity
  call reflectivity(theta, eps, v, h)

  !-- coherent p-polarized reflectivity
  call calc_rho(theta, ks, v, h, rho_v, rho_h)

  !-- extinction coefficients
  call canopy_extinction_coeffs( omega, lai_mw, ke_v, ke_h, ks_v, ks_h )

  !-- canopy transmissivity
  call canopy_transmissivity( ke_v, ke_h, canht, theta, t_v, t_h )

  !-- canopy contribution
  !-- (Ulaby, eq. 11.22)
  call ScatRayleigh_sigma_v_back(ks_v, ks_h, sigma_vol_back)
  call canopy_sigma_c(sigma_vol_back, theta, ke_v, ke_h, t_v, t_h, s0c)

  !-- ground contribution
  call Oh92_backscatter(eps, ks, theta, soil_backscatter)
  call ground_sigma(soil_backscatter, t_v, t_h, s0g)

  !-- total canopy-ground
  call ScatRayleigh_sigma_v_bistatic(ks_v, ks_h, sigma_vol_bistatic)
  call sigma_cg(sigma_vol_bistatic, canht, coherent, rho_v, rho_h, t_v, t_h, s0cgt)  

  !-- ground-canopy-ground
  call sigma_gcg(sigma_vol_back, theta, rho_v, rho_h, t_v, t_h, ke_v, ke_h, s0gcg)

  !-- sum-up contributions
  s0hh = s0c(1) + s0g(1) + s0cgt(1) + s0gcg(1)
  s0vv = s0c(2) + s0g(2) + s0cgt(2) + s0gcg(2)

end subroutine sigma0_hhvv


!***********************************************************
!     sigma0_vhvv
!
!> @brief compute total single-scattering backscattering coefficents for VH and VV polarisations
!
!> @details 
!
!> @param[in]  lai        leaf-area index [m2/m2] as used in optical domain
!> @param[in]  lai_coeff  conversion coefficient for lai optical to microwave
!> @param[in]  canht      canopy height [m]
!> @param[in]  mv         volumetric soilmoisture [m3/m3]
!> @param[out] s0vh       backscattering coefficient, VH polarisation
!> @param[out] s0vv       backscattering coefficient, VV polarisation
!
subroutine sigma0_vhvv(lai, lai_coeff, canht, mv, s0vh, s0vv)
  use mo_sensimul_s1
  implicit none
  ! arguments
  real(kind=8), intent(in) :: lai, lai_coeff, canht, mv
  real(kind=8), intent(out) :: s0vh, s0vv
  ! externals
  external dobson85_eps
  external reflectivity, calc_rho
  external canopy_extinction_coeffs, canopy_transmissivity
  external ScatRayleigh_sigma_v_back, ScatRayleigh_sigma_v_bistatic
  external Oh92_backscatter
  external canopy_sigma_c, ground_sigma, sigma_cg, sigma_gcg
  ! local declarations
  complex(kind=8) :: eps
  real(kind=8) :: v,h,rho_v,rho_h, t_v, t_h, ke_v, ke_h, ks_h, ks_v
  real(kind=8) :: lai_mw
  real(kind=8) :: soil_backscatter(3)                      !hh,vv,hv
  real(kind=8) :: sigma_vol_bistatic(3), sigma_vol_back(3) !hh,vv,hv
  real(kind=8) :: s0c(3), s0g(3), s0cgt(3), s0gcg(3)       !hh,vv,hv

  !-- microwave LAI
  lai_mw = lai*lai_coeff

  !-- soil-moisture ==> eps (dielectric permittivity)
  call dobson85_eps(mv, bulk, alpha, beta1, beta2, ew, eps)

  !-- reflectivity
  call reflectivity(theta, eps, v, h)

  !-- coherent p-polarized reflectivity
  call calc_rho(theta, ks, v, h, rho_v, rho_h)

  !-- extinction coefficients
  call canopy_extinction_coeffs( omega, lai_mw, ke_v, ke_h, ks_v, ks_h )

  !-- canopy transmissivity
  call canopy_transmissivity( ke_v, ke_h, canht, theta, t_v, t_h )

  !-- canopy contribution
  !-- (Ulaby, eq. 11.22)
  call ScatRayleigh_sigma_v_back(ks_v, ks_h, sigma_vol_back)
  call canopy_sigma_c(sigma_vol_back, theta, ke_v, ke_h, t_v, t_h, s0c)

  !-- ground contribution
  call Oh92_backscatter(eps, ks, theta, soil_backscatter)
  call ground_sigma(soil_backscatter, t_v, t_h, s0g)

  !-- total canopy-ground
  call ScatRayleigh_sigma_v_bistatic(ks_v, ks_h, sigma_vol_bistatic)
  call sigma_cg(sigma_vol_bistatic, canht, coherent, rho_v, rho_h, t_v, t_h, s0cgt)  

  !-- ground-canopy-ground
  call sigma_gcg(sigma_vol_back, theta, rho_v, rho_h, t_v, t_h, ke_v, ke_h, s0gcg)

  !-- sum-up contributions
  !-- for S1 we have sender=receiver, thus VH equals HV (reciprocal theorem)
  s0vv = s0c(2) + s0g(2) + s0cgt(2) + s0gcg(2)
  s0vh = 0._8
  if( s0c(3).ne.sense_fv )   s0vh = s0vh + s0c(3)
  if( s0g(3).ne.sense_fv )   s0vh = s0vh + s0g(3)
  if( s0cgt(3).ne.sense_fv ) s0vh = s0vh + s0cgt(3)
  if( s0gcg(3).ne.sense_fv ) s0vh = s0vh + s0gcg(3)


end subroutine sigma0_vhvv


!***********************************************************
!     sigma0_vhvvhh
!
!> @brief compute total single-scattering backscattering coefficents for HH and VV polarisations
!
!> @details 
!
!> @param[in]  lai        leaf-area index [m2/m2] as used in optical domain
!> @param[in]  lai_coeff  conversion coefficient for lai optical to microwave
!> @param[in]  canht      canopy height [m]
!> @param[in]  mv         volumetric soilmoisture [m3/m3]
!> @param[out] s0vh       backscattering coefficient, VH polarisation
!> @param[out] s0vv       backscattering coefficient, VV polarisation
!> @param[out] s0hh       backscattering coefficient, HH polarisation!
!
subroutine sigma0_vhvvhh(lai, lai_coeff, canht, mv, s0vh, s0vv, s0hh)
  use mo_sensimul_s1
  implicit none
  ! arguments
  real(kind=8), intent(in) :: lai, lai_coeff, canht, mv
  real(kind=8), intent(out) :: s0vh, s0vv, s0hh
  ! externals
  external dobson85_eps
  external reflectivity, calc_rho
  external canopy_extinction_coeffs, canopy_transmissivity
  external ScatRayleigh_sigma_v_back, ScatRayleigh_sigma_v_bistatic
  external Oh92_backscatter
  external canopy_sigma_c, ground_sigma, sigma_cg, sigma_gcg
  ! local declarations
  complex(kind=8) :: eps
  real(kind=8) :: v,h,rho_v,rho_h, t_v, t_h, ke_v, ke_h, ks_h, ks_v
  real(kind=8) :: lai_mw
  real(kind=8) :: soil_backscatter(3)                      !hh,vv,vh
  real(kind=8) :: sigma_vol_bistatic(3), sigma_vol_back(3) !hh,vv,vh
  real(kind=8) :: s0c(3), s0g(3), s0cgt(3), s0gcg(3)       !hh,vv,vh

  !-- microwave LAI
  lai_mw = lai*lai_coeff

  !-- soil-moisture ==> eps (dielectric permittivity)
  call dobson85_eps(mv, bulk, alpha, beta1, beta2, ew, eps)

  !-- reflectivity
  call reflectivity(theta, eps, v, h)

  !-- coherent p-polarized reflectivity
  call calc_rho(theta, ks, v, h, rho_v, rho_h)

  !-- extinction coefficients
  call canopy_extinction_coeffs( omega, lai_mw, ke_v, ke_h, ks_v, ks_h )

  !-- canopy transmissivity
  call canopy_transmissivity( ke_v, ke_h, canht, theta, t_v, t_h )

  !-- canopy contribution
  !-- (Ulaby, eq. 11.22)
  call ScatRayleigh_sigma_v_back(ks_v, ks_h, sigma_vol_back)
  call canopy_sigma_c(sigma_vol_back, theta, ke_v, ke_h, t_v, t_h, s0c)

  !-- ground contribution
  call Oh92_backscatter(eps, ks, theta, soil_backscatter)
  call ground_sigma(soil_backscatter, t_v, t_h, s0g)

  !-- total canopy-ground
  call ScatRayleigh_sigma_v_bistatic(ks_v, ks_h, sigma_vol_bistatic)
  call sigma_cg(sigma_vol_bistatic, canht, coherent, rho_v, rho_h, t_v, t_h, s0cgt)  

  !-- ground-canopy-ground
  call sigma_gcg(sigma_vol_back, theta, rho_v, rho_h, t_v, t_h, ke_v, ke_h, s0gcg)

  !-- sum-up contributions
  s0hh = s0c(1) + s0g(1) + s0cgt(1) + s0gcg(1)
  s0vv = s0c(2) + s0g(2) + s0cgt(2) + s0gcg(2)
  s0vh = 0._8
  if( s0c(3).ne.sense_fv )   s0vh = s0vh + s0c(3)
  if( s0g(3).ne.sense_fv )   s0vh = s0vh + s0g(3)
  if( s0cgt(3).ne.sense_fv ) s0vh = s0vh + s0cgt(3)
  if( s0gcg(3).ne.sense_fv ) s0vh = s0vh + s0gcg(3)
end subroutine sigma0_vhvvhh
