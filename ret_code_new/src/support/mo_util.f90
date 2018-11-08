module mo_util

  contains
    !==============================
    !
    !          time_now
    !
    ! Purpose:
    ! - return ISO6801 compliant string representation
    !   (!YYYYMMDDThh:mm:ss+xx:yy) of current system time
    !
    function time_now()
      implicit none
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

    !==============================
    !
    !          split_string
    !
    ! Purpose:
    ! - split a string into 2 either side of a delimiter token
    !
    subroutine split_string(instring, mxnsplit, delim, ntok, splitted_string)
      implicit none
      ! arguments
      character(len=*), intent(in) :: instring
      character(len=1), intent(in) :: delim
      integer, intent(in) :: mxnsplit
      integer, intent(out) :: ntok
      character(len=*), intent(out) :: splitted_string(mxnsplit)
      ! local variables
      character(len=256) :: act_string
      integer :: index
      integer :: iscan

      ntok = -1
      splitted_string(1:mxnsplit) = ''
      iscan = 0
      act_string = trim(instring)
      index = scan(act_string,delim)
      delimloop:do while( index.gt.0 )
         iscan = iscan + 1
         splitted_string(iscan) = trim(act_string(1:index-1))
         act_string = act_string(index+1:)
         index = scan(trim(act_string),delim)
      end do delimloop

      ntok = iscan+1
      if( iscan.eq.0 ) then
         splitted_string = trim(instring)
      else
         splitted_string(iscan+1) = trim(act_string)
      endif

    end subroutine split_string


    logical function file_exists(fname)
      implicit none
      character(len=*), intent(in) :: fname
      inquire(file=fname,exist=file_exists)
    end function file_exists


    integer function get_file_unit(lu_max)
      implicit none
      integer, intent(in) :: lu_max
      !   get_file_unit returns a unit number that is not in use
      integer :: lu, m, iostat
      logical :: opened

      get_file_unit = -1
      if( lu_max.lt.0 )  then
         write(*,'(a)') ' FATAL::get_file_unit::lu_max is negative'
         stop
      endif

      m = lu_max
      do lu = m,1,-1
         inquire (unit=lu, opened=opened, iostat=iostat)
         if (iostat.ne.0) then
            cycle
         else if (.not.opened) then
            get_file_unit = lu
            exit
         endif
      end do

      return
    end function get_file_unit

    subroutine print_vector(n, vec, float_fmt)
      implicit none
      integer, intent(in) :: n
      real(kind=8), intent(in) :: vec(n)
      character(len=*)  :: float_fmt
      integer :: i

      do i=1,n
         write(*, float_fmt, advance='no') vec(i)
      enddo
      write(*, '(/)')

    end subroutine print_vector


    subroutine print_matrix(m, n, mat, float_fmt)
      implicit none
      integer, intent(in) :: m, n
      real(kind=8), intent(in) :: mat(m,n)
      character(len=*)  :: float_fmt
      character(len=32) :: fmt
      character(len=3)  :: istr
      integer :: i

      if( n<10 ) then
         write(istr, '(i1)') n
      else if( n<100 ) then
         write(istr, '(i2)') n
      else if( n<1000 ) then
         write(istr, '(i3)') n
      else
         write(*,'(a)') ' FATAL::print_matrix:n to large...'
         stop
      endif

      fmt = '(' // trim(istr) // trim(float_fmt) // ')'

      do i=1,m
         print fmt, mat(i,1:n)
      enddo

    end subroutine print_matrix

end module mo_util
