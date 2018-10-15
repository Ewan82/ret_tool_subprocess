!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!> \file mo_retrieval.f90
!
!> PROJECT : Sentinel Synergy Study
!
!> DESCRIPTION: provides environment to control Optimisation settings
!
!> \authors The Inversion Lab (Michael Vossbeck, Thomas Kaminski) 
!
!> \date  April-August 2018
!                               
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module mo_retrieval
  implicit none

  !-- bound specification of parameter and state-components
  !   NOTE:below DEFAULT values suitable for an agricultural site are set,
  !        but these *may* be overwritten by reading from external bounds_file
  !        if this is present in the working directory
  character(len=256) :: bounds_file  = 'retrctlvecbounds.nml'
  real(kind=8) :: lai_coeff_lobound = 1.e-8_8 !-- CRITICAL: LAI_mw must not become zero
  real(kind=8) :: lai_coeff_hibound = 1._8    !
  real(kind=8) :: lai_lobound   = 1.e-8_8  !-- CRITICAL: LAI_mw must not become zero
                                         !   (otherwise we are getting zero extinction coefficients)
  real(kind=8) :: lai_hibound   = 10._8    !-- maybe 8.
  real(kind=8) :: canht_lobound = 0._8     !
  real(kind=8) :: canht_hibound = 4.5_8    !-- maybe 3.5
  real(kind=8) :: sm_lobound    = 0._8     !-- maybe 0.1
  real(kind=8) :: sm_hibound    = 0.5_8    !-- maybe 0.4
  namelist /ctlvector_bounds/ &
       lai_coeff_lobound, lai_coeff_hibound, &
       lai_lobound, lai_hibound, &
       canht_lobound, canht_hibound, &
       sm_lobound, sm_hibound

     
  logical :: retr_use_prior_term = .true.
  logical :: retr_use_model_term = .true.
  logical :: retr_prior_diag     = .false.
  logical :: retr_state_diag     = .false.

  !--
  real(kind=8) :: prior_pert   = 0._8
  real(kind=8) :: gradient_tol = 1.e-5_8

  !-- settings here may be read from external bounds_file
  character(len=256) :: retrctl_file = 'retrctl.nml'
  namelist /retrctl/ retr_use_prior_term,retr_use_model_term,gradient_tol,prior_pert

  !-- diagnostics
  real(kind=8) :: ngrad, costp, costm, costo

end module mo_retrieval
