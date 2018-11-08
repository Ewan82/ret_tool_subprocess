subroutine real_fw(x, xw, y, yw)
  implicit none
  ! arguments
  complex(kind=8), intent(in) :: x, xw
  real(kind=8), intent(out) :: y, yw
  y = real(x, kind=8)
  yw = real(xw, kind=8)
end subroutine real_fw

subroutine real_fwv(x, xw, y, yw, nbdirs)
  use DIFFSIZES
  implicit none
  ! arguments
  integer, intent(in) :: nbdirs
  complex(kind=8), intent(in) :: x, xw(nbdirsmax)
  real(kind=8), intent(out) :: y, yw(nbdirsmax)
  y = real(x, kind=8)
  yw = real(xw, kind=8)
end subroutine real_fwv

subroutine real_bw(x, xb, yb)
  implicit none
  ! arguments
  complex(kind=8), intent(in) :: x
  complex(kind=8), intent(inout) :: xb
  real(kind=8), intent(inout) :: yb
  xb = xb + cmplx(yb, 0._8)
end subroutine real_bw

subroutine aimag_fw(x, xw, y, yw)
  implicit none
  ! arguments
  complex(kind=8), intent(in) :: x, xw
  real(kind=8), intent(out) :: y, yw
  y = aimag(x)
  yw = aimag(xw)
end subroutine aimag_fw

subroutine aimag_fwv(x, xw, y, yw, nbdirs)
  use DIFFSIZES
  implicit none
  ! arguments
  integer, intent(in) :: nbdirs
  complex(kind=8), intent(in) :: x, xw(nbdirsmax)
  real(kind=8), intent(out) :: y, yw(nbdirsmax)
  y = aimag(x)
  yw = aimag(xw)
end subroutine aimag_fwv

subroutine aimag_bw(x, xb, yb)
  implicit none
  ! arguments
  complex(kind=8), intent(in) :: x
  complex(kind=8), intent(inout) :: xb
  real(kind=8), intent(inout) :: yb
  xb = xb + cmplx(0._8, yb, kind=8)
end subroutine aimag_bw


real(kind=8) function csq2_fw1(z, z_fw, csq2)
  implicit none
  complex(kind=8), intent(in) :: z
  complex(kind=8), intent(in) :: z_fw
  real(kind=8), intent(out) :: csq2
  intrinsic real
  intrinsic aimag

  csq2_fw1 = 2*real(z)*real(z_fw) + 2*aimag(z)*aimag(z_fw)
  csq2 = real(z)**2 + aimag(z)**2
end function csq2_fw1


subroutine csq2_bw1(z, z_bw, csq2_bw)
  implicit none
  complex(kind=8), intent(in) :: z
  complex(kind=8) :: z_bw
  real(kind=8) :: csq2_bw
  intrinsic real
  intrinsic aimag

  z_bw = conjg(cmplx(2*real(z),2*aimag(z),kind=8))*csq2_bw
end subroutine csq2_bw1
