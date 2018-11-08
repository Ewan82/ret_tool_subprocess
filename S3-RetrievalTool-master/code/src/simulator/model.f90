!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!> \file state.f90
!> \brief provides interfaces to control and run generic S1/S2 observation operators
!>               suitable for use within optimisation framework.
!> \authors The Inversion Lab (Michael Vossbeck, Thomas Kaminski) 
!> \date  April 2018
!                               
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!***********************************************************
!  MODULE mo_model
!> \brief stores the dynamical model parameters
!> \authors The Inversion Lab (Michael Vossbeck, Thomas Kaminski) 
!> \date April 2018
module mo_model
  integer:: nstates = -1
  real(kind=8), allocatable :: mat_model(:,:)   !< matrix applied to state-vector
  real(kind=8), allocatable :: offset_model(:)  !< this is the 'b'
  real(kind=8), allocatable :: unc_model(:)     !< this is the 'sigma_model'
  real(kind=8), allocatable :: a_model(:)       !< the 'a' for the simple model
contains
  logical function is_set()
    is_set = (nstates.ge.0)
  end function is_set
  subroutine state_init(n)
    implicit none
    integer, intent(in) :: n
    nstates = n
    allocate( mat_model(nstates,nstates), offset_model(nstates), unc_model(nstates) )
    allocate( a_model(nstates) )
  end subroutine state_init
  subroutine state_dispose()
    implicit none
    if( allocated(mat_model) )    deallocate(mat_model)
    if( allocated(offset_model) ) deallocate(offset_model)
    if( allocated(unc_model) )    deallocate(unc_model)
    if( allocated(a_model) )      deallocate(a_model)
  end subroutine state_dispose
end module mo_model


subroutine state_model_dispose()
  use mo_model, only: state_dispose
  implicit none
  call state_dispose()
end subroutine state_model_dispose


subroutine state_model_set()
  use mo_sensimul, only:npts,nsc,get_ns,model_file
  use mo_model
  implicit none
  ! local decls
  real(kind=8) :: a_dft=1._8, b_dft=0._8, munc_dft=1._8
  integer :: ns
  integer :: ipt,isppt
  integer :: irow,icol
  logical :: exist,succeed
  ! externals
  external ncrd_retrieval_model

  if( is_set() ) then
     write(*,'(a)') ' FATAL::state_model_set:should actually only be called for not '//&
          'initialised state model.'
     stop
  else
     !-- length of state-vector
     ns = get_ns()
     !-- allocate structures
     call state_init(ns)

     write(*,'(a)') ' INFO::state_model_set:'//&
          'will be setting-up simple state-model of type "mx(i+1)=a(i+1)*x(i) + b(i+1)"'
     !-- eventually we would start with a simple first model:
     !   - simulates the components of the state vector independent of each other
     !   - each state vector component is equal to that component at the "previous" point,
     !     multiplied by a factor (a) plus a constant (b)
     !   - a and b depend on state vector component and point
     !-- 
     !   m_lai(1)   =        lai(1) + b(1)
     !   m_lai(i+1) = a(i+1)*lai(i) + b(i+1) i=1,...,npts-1
     !   ...and similar for canopy-height, soil moisture
     inquire(file=model_file,exist=exist)
     if( .not.exist ) then
        write(*,'(a,3(a,f5.2,1x))') ' WARN::state_model_set:'//&
             'model file ***'//trim(model_file)//'*** not present, will set simple defaults ',&
             'a=',a_dft,'b=',b_dft,'munc=',munc_dft
        a_model      = a_dft
        offset_model = b_dft
        unc_model    = munc_dft
     else
        write(*,'(a)') ' INFO::state_model_set:'//&
             'reading model configuration from file ***'//trim(model_file)//'***...'
        call ncrd_retrieval_model(model_file, succeed)
        if( .not.succeed ) then
           write(*,'(a)') ' FATAL::state_model_set:'//&
                'reading model configuration failed.'
        else
           write(*,'(a)') ' INFO::state_model_set:...reading done.'
        endif
     endif

     !-- set full matrix
     mat_model = 0._8
     do ipt=1,npts
        do isppt=1,nsc
           irow = isppt + (ipt-1)*nsc
           if( ipt.eq.1 ) then
              icol = isppt
           else
              icol = isppt + (ipt-2)*nsc
           endif
           mat_model(irow,icol) = a_model(irow)
        enddo
     enddo

  endif

end subroutine state_model_set


!***********************************************************
!     H_m
!
!> @brief implements Eq. (1.2) in the prototype tool description, the model 'M'determines overall size of control vector (n) and overall size of
!         simulation vector (m)
!         In addition necessary initialisations to actually run the observational operator(s)
!         are performed.
!
!> @param[in]  n  overall (1D-)length of control vector
!> @param[in]  x  overall (1D) control vector (normlalised)
!> @param[out] statediff differences between control vector and the mapped control vector
!
subroutine H_m(n, x, statediff)
  use mo_sensimul, only:get_n,get_ns,get_np
  use mo_model
  implicit none
  ! arguments
  integer, intent(in) :: n
  real(kind=8), intent(in) :: x(n)
  real(kind=8), intent(out) :: statediff(n)
  ! local decls
  integer :: nc,ns,np
  real(kind=8) :: xm(n), xphys(n)
  integer :: i,j

  nc = get_n()  !-- #control vector
  ns = get_ns() !-- number of states
  np = get_np() !-- number of parameter, first element of state-vector starts at np+1

  !-- dimensional consistency
  if( n.ne.nc ) then
     write(*,'(a(2a,i3,1x))') ' FATAL::H_m:'//&
          'inconsistent dimension of state vector:','expected=',nc,'got=',n
     stop
  endif

  !-- parameter do not contribute here
  statediff(1:np) = 0._8
  xm(1:np)        = 0._8 !these are *not* used below

  !-- control vector in physical units
  call x2p(n, x, xphys)

  !-- xm = A*xs + b
  do i=1,ns
     xm(np+i) = offset_model(i)
     do j=1,ns
        xm(np+i) = xm(np+i) + mat_model(i,j)*xphys(np+j)
     enddo
  end do

  !--
  do i=1,ns
     statediff(np+i) = (xm(np+i)-xphys(np+i))/unc_model(i)
  enddo
  !iLab::AD-Problem
  ! statediff(np+1:np+ns) = (xm(np+1:np+ns)-x(np+1:np+ns))/unc_model(1:ns)
end subroutine H_m

