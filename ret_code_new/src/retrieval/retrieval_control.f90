!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!> \file optim_control.f90
!
!> PROJECT : Sentinel Synergy Study
!
!> DESCRIPTION: provides interfaces to control Optimisation settings,
!               interfaces reside in global namespace such that downstream
!               application do not need to include the underlying module
!
!> \authors The Inversion Lab (Michael Vossbeck) 
!
!> \date  April 2018
!                               
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine optim_enable_prior_term()
  use mo_retrieval, only:retr_use_prior_term
  retr_use_prior_term = .true.
end subroutine optim_enable_prior_term


subroutine optim_disable_prior_term()
  use mo_retrieval, only:retr_use_prior_term
  retr_use_prior_term = .false.
end subroutine optim_disable_prior_term


subroutine optim_is_prior_term_enabled(flag)
  use mo_retrieval, only:retr_use_prior_term
  implicit none
  logical, intent(out) :: flag
  flag = retr_use_prior_term
end subroutine optim_is_prior_term_enabled


subroutine optim_enable_state_term()
  use mo_retrieval, only:retr_use_model_term
  retr_use_model_term = .true.
end subroutine optim_enable_state_term


subroutine optim_disable_state_term()
  use mo_retrieval, only:retr_use_model_term
  retr_use_model_term = .false.
end subroutine optim_disable_state_term


subroutine optim_is_state_term_enabled(flag)
  use mo_retrieval, only:retr_use_model_term
  implicit none
  logical, intent(out) :: flag
  flag = retr_use_model_term
end subroutine optim_is_state_term_enabled


subroutine retr_get_gradient_tol(tol)
  use mo_retrieval, only:gradient_tol
  implicit none
  real(kind=8), intent(out) :: tol
  tol = gradient_tol
end subroutine retr_get_gradient_tol


subroutine retr_get_prior_pert(pert)
  use mo_retrieval, only:prior_pert
  implicit none
  real(kind=8), intent(out) :: pert
  pert = prior_pert
end subroutine retr_get_prior_pert


subroutine retrieval_read_ctl()
  use mo_retrieval
  implicit none
  ! local decls
  integer :: iounit = 18
  logical :: exist

  inquire(FILE=retrctl_file,exist=exist)
  if( exist ) then
     write(*,'(a)') ' INFO::retrieval_read_ctl:'//&
          'read settings from file ***'//trim(retrctl_file)//'***...'
     open (UNIT=iounit, FILE=retrctl_file )
     read (UNIT=iounit, nml=retrctl)
     close(iounit)
     write(*,'(a)') ' INFO::retrieval_read_ctl:...reading done.'
  else
     write (*, '(a)') ' WARN::retrieval_read_ctl:'//&
          'control file ***'//trim(retrctl_file)//'*** not present,'//&
          ' will rely on internal defaults!'
  end if

end subroutine retrieval_read_ctl


!***********************************************************
!     getbounds
!
!> @brief utility method to be called before running lbfgsb minimiser
!>        which allows to retieve the bounds of the control vector
!>        components.
!>
!> @param[in]   n    length control vector
!> @param[out]  nbd  number of bounds per control vector element (upper,lower)
!> @param[out]  lx   lower bound of control vector elements (normalised units)
!> @param[out]  ux   upper bound of control vector elements (normalised units)
!
subroutine getbounds(n, nbd, lx, ux)
  use mo_retrieval
  use mo_sensimul, only:nparam_s1,nsc
  implicit none
  ! arguments
  integer, intent(in) :: n
  integer, intent(out) :: nbd(n)
  real(kind=8), intent(out) :: lx(n), ux(n)
  ! local decls
  logical :: exist

  !--   NOTE: we will have lower/upper bound for each component
  nbd = 2     !lower/upper

  !-- check for an external bounds file
  inquire(FILE=bounds_file, EXIST=exist)
  if( exist ) then
     call load_boundsnml(bounds_file)
     write(*, '(a)') ' INFO::initb:bounds were read from file ***'//trim(bounds_file)//'***'
  else
     write(*, '(a)') ' INFO::initb:using default bounds'
  endif

  !-- report bounds
  write(*, '(a,2(a,1x,f8.3,1x))') ' INFO::initb:applying bounds: ',&
       'lai_coeff_lobound=',lai_coeff_lobound,'lai_coeff_hibound=',lai_coeff_hibound
  write(*, '(a,2(a,1x,f8.3,1x))') ' INFO::initb:applying bounds: ',&
       'lai_lobound=',lai_lobound,'lai_hibound=',lai_hibound
  write(*, '(a,2(a,1x,f8.3,1x))') ' INFO::initb:applying bounds: ',&
       'canht_lobound=',canht_lobound,'canht_hibound=',canht_hibound
  write(*, '(a,2(a,1x,f8.3,1x))') ' INFO::initb:applying bounds: ',&
       'sm_lobound=',sm_lobound,'sm_hibound=',sm_hibound

  !-- bounds on lai_coeff
  lx(1) = lai_coeff_lobound
  ux(1) = lai_coeff_hibound

  !-- LAI bounds, 1st component
  lx(nparam_s1+1:n:nsc) = lai_lobound
  ux(nparam_s1+1:n:nsc) = lai_hibound

  !-- canopy-height bounds, 2nd component
  lx(nparam_s1+2:n:nsc) = canht_lobound
  ux(nparam_s1+2:n:nsc) = canht_hibound

  !-- SM bounds
  lx(nparam_s1+3:n:nsc) = sm_lobound
  ux(nparam_s1+3:n:nsc) = sm_hibound

contains
  subroutine load_boundsnml(fname)
    implicit none
    ! arguments
    character(len=*), intent(in) :: fname
    ! local decls
    integer, parameter :: iounit=18
    integer :: iostat
    logical :: exist

    !--
    inquire(file=fname, exist=exist)
    if( .not. exist ) then
       write(*, '(a)') ' FATAL::load_boundsnml:'//&
            'non-existent expected file ***'//trim(fname)//'***'
       stop
    else
       open(unit=iounit, file=fname)
       read(unit=iounit, nml=ctlvector_bounds, iostat=iostat)
       if( iostat.gt.0 ) then
          write(*, '(a)') ' FATAL::load_boundsnml:error reading namelist file '//&
               '***'//trim(fname)//'***'
          stop
       endif
       close(iounit)
    endif
  end subroutine load_boundsnml

end subroutine getbounds


!***********************************************************
!     initb
!
!> @brief utility method to be called before running lbfgsb minimiser
!>        which allows to retieve the bounds of the control vector
!>        components.
!>
!> @param[in]   n    length control vector
!> @param[out]  nbd  number of bounds per control vector element (upper,lower)
!> @param[out]  lx   lower bound of control vector elements (normalised units)
!> @param[out]  ux   upper bound of control vector elements (normalised units)
!
subroutine initb(n, nbd, lx, ux)
  implicit none
  ! arguments
  integer, intent(in) :: n
  integer, intent(out) :: nbd(n)
  real(kind=8), intent(out) :: lx(n), ux(n)
  ! local decls
  real(kind=8) :: x(n), sx(n)

  !-- get boundary values
  call getbounds(n, nbd, lx, ux)

  !-- get prior for scaling
  call getprior(n, x, sx)

  !-- apply scaling
  lx = lx/sx
  ux = ux/sx

end subroutine initb


subroutine retrieval_dump_ctl()
  use mo_retrieval
  implicit none
  write(*,'(a,1x,l1)') ' INFO::use_prior_term:',retr_use_prior_term
  write(*,'(a,1x,l1)') ' INFO::use_state_term:',retr_use_model_term
  write(*,'(a,1x,e10.5)') ' INFO::gradient_tolerance:',gradient_tol
  write(*,'(a,1x,e10.5)') ' INFO::priopr_pert:',prior_pert
end subroutine retrieval_dump_ctl
