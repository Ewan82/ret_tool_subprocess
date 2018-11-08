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
  !-- NOTE:'nadime_prospect_price_fast_1geom' uses plain 'real' as type declarator
  real :: brf_(nwvis1nm), fpar_(nwvis1nm), albedo_(nwvis1nm), trans_(nwvis1nm)
  real(kind=8) :: fpar_r8(nwvis1nm)
  ! externals
  real, external :: sm_to_rsl1
  external nadime_prospect_price_full_1geom
  real(kind=8), external :: convolve

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

  !-- call coupled optical model for nwvis1nm wavelenghts starting at 400nm
  call nadime_prospect_price_full_1geom(nw1nm, nwvis1nm, &
       theta_i, phi_i, theta_v, phi_v, &
       lad, rpl, lai, hc, &
       vai, cab, cw, cp, cc, &
       rsl1, rsl2, rsl3, rsl4, &
       brf_, fpar_, albedo_, trans_)

  !-- ensure correct precision
  fpar_r8 = real(fpar_, kind=8)

  !-- integration over wave-lenghts (convolve)
  fapar = convolve(nwvis1nm, fpar_r8, solspect_response)
end subroutine fapar_operator_1geom


!***********************************************************
!     simulate_fapar
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
subroutine simulate_fapar(n, x, m, y)
  use mo_sensimul, only:nsc, iv_geom
  use mo_sensimul, only:get_np,get_n,get_npts
  use mo_sensimul, only:idx_is_s2, sim_fill_value
  implicit none
  ! arguments
  integer, intent(in) :: n, m
  real(kind=8), intent(in) :: x(n)
  real(kind=8), intent(out) :: y(m)
  ! externals
  external fapar_operator_1geom
  ! local decls
  real(kind=8) :: xphys(n)
  integer :: ipt
  integer :: ii0,ii1
  real(kind=8) :: statev(nsc)
  real(kind=8) :: single_ivgeom(4)
  !-- debugging
  real(kind=8) :: init,finish
  logical :: ldebug

  !-- dimension consistency
  if( n.ne.get_n() ) then
     write(*, '(a,2(a,i3,1x))') ' FATAL::simulate_fapar:inconsistent length of state vector!',&
          'expected=',get_n(),'got=',n
     stop
  endif
  if( m.ne.get_npts() ) then
     write(*, '(a,2(a,i3,1x))') ' FATAL::simulate_fapar:inconsistent length of simulation vector.',&
          'expected=',get_npts(),'got=',m
     stop
  endif

  !-- initialise
  ldebug = .false.
  y(1:m) = sim_fill_value

  !-- get physical control vector
  call x2p(n, x, xphys)

  !-- loop over all time-points
  ii0 = get_np() + 1 !-- state-vector components start *after* parameter
  simloop:do ipt=1,m
     if( .not.idx_is_s2(ipt) ) then !-- no S2 simulation at 'ipt', continue to next point
        ii0 = ii0 + nsc
     else
        if( ldebug ) then
           write(*, '(a,i2)') ' DEBUG::simulate_fapar:starting ipt=',ipt
           call cpu_time(init)
        endif
        !-- extract 'nsc' state vector components
        !   for current time-point
        ii1 = ii0 + nsc - 1
        statev(1:nsc) = xphys(ii0:ii1)
        !-- get geometry for this time-point
        single_ivgeom(1:4) = iv_geom(1:4,ipt)
        call fapar_operator_1geom( statev, single_ivgeom, y(ipt) )
        !-- increment start position
        ii0 = ii1 + 1
        if( ldebug ) then
           call cpu_time(finish)
           write(*, '(a,i2,1x,a,f6.2)') ' DEBUG::simulate_fapar:finished ipt=',ipt,&
                'elapsed=',(finish-init)
        endif
     endif
  enddo simloop

end subroutine simulate_fapar
