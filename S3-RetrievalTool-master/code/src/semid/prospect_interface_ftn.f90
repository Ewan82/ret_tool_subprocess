subroutine prospect_fullSpectrum_interp1nm_ftn( vai, cab, cw, cp, cc, refl, tran )
  implicit none
  ! constants
  integer, parameter :: wvl_min = 400
  integer, parameter :: wvl_max = 2500
  integer, parameter :: nw = wvl_max - wvl_min +1
  ! arguments
  real(kind=8), intent(in) :: vai, cab, cw, cp, cc
  real(kind=8), intent(out) :: refl(nw), tran(nw)
  ! local variables
  external leaftwo
  !-- MVO-ATTENTION::dimension=2500 compliant with hard-coded setting in routine
  !                  leaftwo (prospect.redux.f)
  ! double precision :: main_r(2500), main_t(2500)
  double precision :: main_r(0:2499), main_t(0:2499)
  integer :: upper, lower
  real :: fraction
  integer :: w

  !-- 
  call leaftwo(main_r, main_t, vai, cab, cw, cp, cc)

  !-- 
  do w=1,nw
     lower = int((w-1)/5.0) !NOTE::w-1 since main_t,main_r use 0-indexing
     upper = lower + 1
     fraction = real(w-1)/5.0 - real(lower)
     !--MVO:: attention type (real vs. double precision)
     !--MVO:: ATTENTION:for the last wavelength, the upper index is
     !                  beyond the initialised range in main_t,main_r,
     !                  and the value equals the one at the lower index
     if( w.eq.nw ) then
        tran(w) = main_t(lower)
        refl(w) = main_r(lower)
     else
        tran(w) = main_t(lower)*(1.-fraction) + main_t(upper)*fraction
        refl(w) = main_r(lower)*(1.-fraction) + main_r(upper)*fraction
     endif
  enddo

end subroutine prospect_fullSpectrum_interp1nm_ftn


subroutine prospect_monochromatic_ftn( wavelength, vai, cab, cw, cp, cc, refl, tran )
  implicit none
  ! arguments
  real, intent(in) :: wavelength
  real(kind=8), intent(in) :: vai, cab, cw, cp, cc
  real, intent(out) :: refl, tran
  ! local variables
  external leaftwo
  !-- MVO-ATTENTION::dimension=2500 compliant with hard-coded setting in routine
  !                  leaftwo (prospect.redux.f)
  double precision :: main_r(0:2499), main_t(0:2499)
  integer :: upper, lower
  real :: fracw, fraction

  !-- check wavelength
  if( wavelength.lt.400. .or. wavelength.gt.2500. ) then
     write(*, '(a)') 'FATAL::price_soil_ftn::wavelength beyond range [400,2500]. Cannot continue!'
     stop
  endif

  !--
  fracw = (wavelength-400.)/5.0
  lower = int(fracw)
  upper = lower + 1
  fraction = fracw - lower

  call leaftwo(main_r, main_t, vai, cab, cw, cp, cc)

  !--MVO:: attention type (real vs. double precision)
  tran = main_t(lower)*(1.-fraction) + main_t(upper)*fraction
  refl = main_r(lower)*(1.-fraction) + main_r(upper)*fraction

end subroutine prospect_monochromatic_ftn
