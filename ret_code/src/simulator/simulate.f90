!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> \file simulate.f90
!> \brief provides implementation of the observation operators in the
!>              optical and microwave domain.mapping of 1D control vector to physical variables
!> \authors The Inversion Lab (Michael Vossbeck, Thomas Kaminski) 
!> \date  April 2018
!                               
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!***********************************************************
!     simulate_s1_1geom
!
!> @brief Implementation of the Sentinel Synergy Study observation operator
!>        in the microwave domain (denoted H_1 in document D5) for a single state.
!
!> @details For the given input state (LAI_coeff, LAI,HC, and SM) the observation operator
!>          derives backscatter values in VH and VV polarisation.
!>          The implementation is based on the Fortran port of the SENSE model.
!
!> @param[in]   statev   state components (physical units) required to run S1 simulation
!>                       consisting of LAI conversion coefficient, (optical) LAI, HC, and SM
!> @param[out]  bscat    back-scatter values in VH/VV polarisation
!
subroutine simulate_s1_1geom( statev, iv_geom, bscat )
  use mo_sensimul, only:nsc,nparam_s1
  use mo_sensimul_s1, only:theta,pi !-- viewing zenith-angle for microwave simulation
  implicit none
  ! arguments
  real(kind=8), intent(in)  :: statev(nparam_s1+nsc)
  real(kind=8), intent(in)  :: iv_geom(4) !sza,saa,vza,vaa
  real(kind=8), intent(out) :: bscat(2) !VH,VV
  ! externals
  external sigma0_vhvv !core SENSE routine
  ! local declarations
  real(kind=8) :: lai_coeff, lai, hc, sm
  real(kind=8) :: s0vh, s0vv
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
  call sigma0_vhvv(lai, lai_coeff, hc, sm, s0vh, s0vv)

  !-- map backscatter values
  bscat(1) = s0vh
  bscat(2) = s0vv

end subroutine simulate_s1_1geom


!***********************************************************
!     simulate_s2_1geom
!
!> @brief Implementation of the Sentinel Synergy Study observation operator
!>        in the optical domain (denoted H_2 in document D5) for a single state.
!         
!
!> @details For the given input state (LAI,HC, and SM) the observation operator
!>          derives TOC BRFs for all 13 S2 wave-bands. The implementation is based
!>          on the coupled SemiDiscrete-PROSPECT-PRICE model.
!>          Running the SemiDiscrete model requires to specify the illumination- and
!>          view geometry.
!
!> @param[in]   statev   state-vector consisting of LAI, HC, and SM
!> @param[in]   iv_geom  illumination-view geomtry in order sza,saa,vza,vaa (in degrees)
!> @param[out]  brf      TOC BRFs for 13 S2 wave-bands
!
subroutine simulate_s2_1geom( statev, iv_geom, brf )
  use mo_sensimul, only:nsc
  use mo_sensimul_s2
  implicit none
  ! arguments
  real(kind=8), intent(in)  :: statev(nsc)
  real(kind=8), intent(in)  :: iv_geom(4) !sza,saa,vza,vaa
  real(kind=8), intent(out) :: brf(nb_s2)
  ! local decls
  real :: lai, hc, sm
  real :: rsl1
  real :: theta_i, phi_i, theta_v, phi_v
  real :: brf_(nb_s2) !'nadim_prospect_price_fast_1geom' does not expect kind=8 !
  ! externals
  real, external :: sm_to_rsl1
  external nadim_prospect_price_fast_1geom


  !-- map state variables
  !   (potential cast from r8 to r4)
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

  !-- call coupled optical model
  call nadim_prospect_price_fast_1geom(nb_s2, nw1nm, s2_response_mat,&
       theta_i, phi_i, theta_v, phi_v, &
       lad, rpl, lai, hc, &
       vai, cab, cw, cp, cc, &
       rsl1, rsl2, rsl3, rsl4, &
       brf_)

  !-- re-cast buffer
  brf(1:nb_s2) = real(brf_(1:nb_s2), kind=8)
end subroutine simulate_s2_1geom


!***********************************************************
!     simulate_s1
!
!> @brief Implementation of S1 simulation.
!>        Backscatter values in VH and VV polarisation will be computed for the
!>        given control vector.
!
!> @details This routine provides the implementation of the observation operator H_1
!>          as denoted in Figure 1.2 in document D5v1.
!
!> @param[in]   n  length of control-vector for multiple simulations in the MW domain
!> @param[in]   x  S1 relevant part of full control vector in physical units
!                  (expected ordering is: S1 related parameter(s) followed by
!                   by state-variables LAI,HC,SM per state)
!> @param[in]   m  length of output vector
!> @param[out]  y  simulated backscatter values (VH,VV) per state
!
subroutine simulate_s1(n, x, m, y)
  use mo_sensimul, only:nparam_s1, nsc, npts_s1, nparam_s1, get_nc_s1, get_m_s1
  use mo_sensimul, only:iv_geom, timept_idxs_s1
  use mo_sensimul_s1, only:npol !VH,VV
  implicit none
  ! arguments
  integer, intent(in) :: n, m
  real(kind=8), intent(in) :: x(n)
  real(kind=8), intent(out) :: y(m)
  ! externals
  external simulate_s1_1geom
  ! local decls
  integer :: ipt, ipt_s1
  integer :: i0,i1,j0,j1
  real(kind=8) :: statev(nparam_s1+nsc)
  real(kind=8) :: single_ivgeom(4)
  !-- dimension consistency
  if( n.ne.get_nc_s1() ) then
     write(*, '(a,2(a,i3,1x))') ' FATAL::simulate_s1:inconsistent length of state vector.',&
          'expected=',get_nc_s1(),'got=',n
     stop
  endif
  if( m.ne.get_m_s1() ) then
     write(*, '(a,2(a,i3,1x))') ' FATAL::simulate_s1:inconsistent length of result vector.',&
          'expected=',get_m_s1(),'got=',m
  endif

  !-- map S1 *parameter part* into model input vector
  !   NOTE: S1 parameter is part of control vector
  statev(1:nparam_s1) = x(1:nparam_s1)

  !-- loop over time-points
  i0 = nparam_s1 + 1 !-- initial position in 'x' vector
  j0 = 1             !-- initial position in 'y' vector
  simloop:do ipt=1,npts_s1
     i1 = i0 + nsc - 1
     j1 = j0 + npol - 1
     statev(nparam_s1+1:nparam_s1+nsc) = x(i0:i1)
     ipt_s1 = timept_idxs_s1(ipt)
     single_ivgeom(1:4) = iv_geom(1:4,ipt_s1)
     call simulate_s1_1geom( statev, single_ivgeom, y(j0:j1) )
     !-- increment start position
     i0 = i1 + 1
     j0 = j1 + 1
  enddo simloop

end subroutine simulate_s1


!***********************************************************
!     simulate_s2
!
!> @brief Implementation of the S2 observation operator.
!>        Computes top-of-canopy BRFs in all 13 Sentinel2 bands for every
!>        state in the control vector.
!
!
!> @param[in]   n  length of control vector for multiple simulations in the optical domain
!> @param[in]   x  S2 relevant part of full control vector in physical units.
!                  (expected ordering is: state variables LAI,HC,SM per state)
!> @param[in]   m  length of output vector
!> @param[out]  y  simulated TOC BRFs (for all 13 S2 wave-bands) and multiple observations
!                  (ordering will be BRF1-BRF12 per state)
!
subroutine simulate_s2(n, x, m, y)
  use mo_sensimul, only:nsc, npts_s2, iv_geom, get_nc_s2, get_m_s2, timept_idxs_s2
  use mo_sensimul_s2, only:nb_s2
  implicit none
  ! arguments
  integer, intent(in) :: n, m
  real(kind=8), intent(in) :: x(n)
  real(kind=8), intent(out) :: y(m)
  ! external
  external simulate_s2_1geom
  ! local decls
  integer :: ipt, ipt_s2
  integer :: ii0,ii1,j0,j1
  real(kind=8) :: statev(nsc)
  real(kind=8) :: single_ivgeom(4)


  !-- dimension consistency
  if( n.ne.get_nc_s2() ) then
     write(*, '(a,2(a,i3,1x))') ' FATAL::simulate_s2:inconsistent length of state vector!',&
          'expected=',get_nc_s2(),'got=',n
     stop
  endif
  if( m.ne.get_m_s2() ) then
     write(*, '(a,2(a,i3,1x))') ' FATAL::simulate_s2:inconsistent length of result vector.',&
          'expected=',get_m_s2(),'got=',m
  endif

  !-- loop over all time-points
  ii0 = 1
  j0 = 1
!!! !$AD II-LOOP
  simloop:do ipt=1,npts_s2
     ii1 = ii0 + nsc - 1
     j1  =  j0 + nb_s2 - 1
     statev(1:nsc) = x(ii0:ii1)
     ipt_s2 = timept_idxs_s2(ipt)
     single_ivgeom(1:4) = iv_geom(1:4,ipt_s2)
     call simulate_s2_1geom( statev, single_ivgeom, y(j0:j1) )
     !-- increment start position
     ii0 = ii1 + 1
     j0  =  j1 + 1
  enddo simloop

end subroutine simulate_s2


!***********************************************************
!     simulate_s1s2
!
!> @brief Implementation of the combined S1 and S2 observation operator.
!>        Computes backscatter values (VH,VV) and top-of-canopy BRFs in 13 Sentinel2 bands
!>        for the respective states in the control vector.
!
!
!> @param[in]   n  length of control vector for multiple simulations in MW and optical domain
!> @param[in]   x  control vector normalised by prior uncertainty
!                  (expected ordering is: S1 related parameter(s) followed by
!                   by state variables LAI,HC,SM per 'simulation point')
!> @param[in]   m  length of output vector
!> @param[out]  y  simulated backscatter values (VH,VV) and simulated TOC BRFs
!                  (for all 13 S2 wave-bands)
!                  (Ordering is: simulated backscatter values per 'simulation point'
!                                followed by BRF values per 'simulation point')
!
subroutine simulate_s1s2(n, x, m, y)
  use mo_sensimul, only:get_nc_s1,get_nc_s2,get_m_s1,get_m_s2,get_n,get_m
  use mo_sensimul, only:sim_fill_value
  implicit none
  ! arguments
  integer, intent(in) :: n, m
  real(kind=8), intent(in) :: x(n)
  real(kind=8), intent(out) :: y(m)
  ! local decls
  real(kind=8) :: xp_s1(get_nc_s1()),xp_s2(get_nc_s2())
  integer :: n_s1,n_s2,m_s1,m_s2
  real(kind=8) :: xphys(n)

  !-- dimension consistency
  if( n.ne.get_n() ) then
     write(*, '(a,2(a,i3,1x))') ' FATAL::simulate_s1s2:inconsistent length of state vector!',&
          'expected=',get_n(),'got=',n
     stop
  endif
  if( m.ne.get_m() ) then
     write(*, '(a,2(a,i3,1x))') ' FATAL::simulate_s1s2:inconsistent length of simulation vector.',&
          'expected=',get_m(),'got=',m
     stop
  endif

  !-- initialise output
  y = sim_fill_value

  !-- convert normalised to physical control vector
  call x2p(n, x, xphys)

  !-- extract specific control vectors
  n_s1 = get_nc_s1()
  n_s2 = get_nc_s2()
  m_s1 = get_m_s1()
  m_s2 = get_m_s2()
  call control_vector_split(n, xphys, n_s1, xp_s1, n_s2, xp_s2)

  !-- S1 simulation
  if( n_s1.gt.0 ) then
     call simulate_s1(n_s1, xp_s1, m_s1, y(1:m_s1))
  endif

  !-- S2 simulation
  if( n_s2.gt.0 ) then
     call simulate_s2(n_s2, xp_s2, m_s2, y(m_s1+1:m_s1+m_s2))
  endif

end subroutine simulate_s1s2
