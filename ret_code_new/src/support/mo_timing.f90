!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!> \file mo_timing.f90
!
!> DESCRIPTION: support module for timing purpose(s)
!
!> \authors The Inversion Lab (Michael Vossbeck) 
!
!> \date  March 2017
!                               
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module mo_timing
  implicit none

  !-----------------------------
  !     timer environment
  !
  integer, private, parameter :: dp = kind(1.0d0)
  integer, private :: crate, csecmax
  integer, private :: dt
  integer, private :: iclock0(8), iclock(8)
  integer, private :: isec0, isec
  logical, private :: tm_ok

  public :: btiming, etiming

contains

  subroutine btiming()
    implicit none
    if( .not. tm_ok ) then
       call system_clock( count_rate=crate, count_max=csecmax )
       if( crate.eq.0 ) then
          print '(a,a)', &
               '  ERROR::btiming::system colock not available, ',&
               &'no timing possible.'
          print*
          stop
       else
          tm_ok = .true.
       endif
    end if
    call date_and_time( values=iclock0 )
    call system_clock( count=isec0 )
  end subroutine btiming


  subroutine etiming(deltat)
    implicit none
    real(kind=dp), intent(out) :: deltat
    integer :: idt
    integer :: iday,iday0
    integer :: imon,imon0

    if( .not. tm_ok ) then
       print '(2x,a,a)', &
            'ERROR::etiming:: should not have been called, ',&
            '                 because no system clock is available!'
       print*
       stop 1
    end if
    call system_clock( count=isec )
    call date_and_time( values=iclock )

    imon  = iclock(2)
    imon0 = iclock0(2)
    iday  = iclock(3)  !
    iday0 = iclock0(3) !
    if( imon .ne. imon0 ) then
       write(*,'(2x,a)') 'WARNING::etiming::timing that exceeds monthly bounds is not yet supported!'
       write(*,'(2x,a)') '                  erroneously returning 0.0 to the caller.'
       deltat = 0.0D0
    elseif( iday.eq.iday0 ) then
       idt = (isec-isec0) - dt
       deltat = real(idt,dp)/crate
    elseif( iday.gt.iday0 ) then
       isec=isec + (iday-iday0)*csecmax
       deltat=float((isec-isec0) - dt)/crate   
    else
       print*
       write(*,'(2x,a)') 'SEVERE-ERROR::etiming::should have never reached here. Will abort!'
       stop
    endif
  end subroutine etiming

  !-----------------------------
  !     time_now (ISO8601 compliant)
  !
  function time_now()
    character(len=23) :: time_now !YYYYMMDDThh:mm:ss+xx:yy
    character(len=8) :: date
    character(len=10) :: time
    character(len=5) :: zone
    character(len=2) :: hh,mm,ss
    integer, dimension(8) :: values
    call date_and_time(date, time, zone, values)
    write(hh,'(I0.2)') values(5)
    write(mm,'(I0.2)') values(6)
    write(ss,'(I0.2)') values(7)
    time_now = date//'T'//hh//':'//mm//':'//ss//zone(1:3)//':'//zone(4:5)
    return
  end function time_now

end module mo_timing
