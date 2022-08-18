module by03
  ! Bodek-Yang neutrino structure functions from hep-ex/0203009
  ! Nucl. Phys. Proc. Suppl. 112 (2002) 70-76
  use iso_fortran_env, only: dp=>real64
  use grv94
  implicit none
  public :: F2_emu, F2_neutrino, F2_antineutrino, corr_fac_2xF1_F2, &
    xF3_neutrino, xF3_antineutrino
  private
  real(dp), parameter :: A = 1.753_dp, B = 0.624_dp, C = 0.188_dp, &
    M = 0.938272046_dp, mc2 = 1.3_dp**2
  integer, parameter :: nf = 4 ! u, d, s, c
  real(dp), dimension(nf), parameter :: ei2 = ([2, 1, 1, 2]/3.0_dp)**2
contains
  subroutine calc_xw(x, Q2, xw)
    real(dp), intent(in) :: x, Q2
    real(dp), intent(out) :: xw

    xw = x*(Q2 + B)/(Q2 + A*x)
  end subroutine calc_xw

  function R_world(x, Q2)
    real(dp) :: R_world
    real(dp), intent(in) :: x, Q2

    real(dp) :: theta

    theta = 1 + 12*Q2/(Q2 + 1)*0.125_dp**2/(0.125_dp**2 + x**2)
    R_world = 0.0635_dp/log(Q2/0.04_dp)*theta &
      + 0.5747_dp/Q2 - 0.3534_dp/(Q2**2 + 0.09_dp)
  end function R_world

  function R(x, Q2)
    real(dp) :: R
    real(dp), intent(in) :: x, Q2

    if (Q2 < 0.35_dp) then
      R = 3.207_dp*Q2/(Q2**2 + 1)*R_world(x, 0.35_dp)
    else
      R = R_world(x, Q2)
    end if
  end function R

  subroutine calc_pdf(x, Q2, xq, xqbar)
    real(dp), intent(in) :: x, Q2
    real(dp), dimension(nf), intent(out) :: xq, xqbar

    real(dp) :: xw, Q2w, xi_c, uv, dv, us, ds, delta_du

    call calc_xw(x, Q2, xw)
    Q2w = max(Q2, 0.24_dp)

    uv = xuv(xw, Q2w)
    dv = xdv(xw, Q2w)
    us = (xud(xw, Q2w) - xDelta(xw, Q2w))/2
    ds = (xud(xw, Q2w) + xDelta(xw, Q2w))/2
    delta_du = -0.0161_dp + 0.0549_dp*x + 0.355_dp*x**2 - 0.193_dp*x**3

    xqbar(1) = us/(1 + delta_du*us/(us + ds)) 
    xq(1) = uv/(1 + delta_du*uv/(uv + dv)) + xqbar(1)

    xqbar(2) = (ds + us*delta_du)/(1 + delta_du*us/(us + ds))
    xq(2) = (dv + uv*delta_du)/(1 + delta_du*uv/(uv + dv)) + xqbar(2)

    xq(3) = xs(xw, Q2w)
    xqbar(3) = xq(3)

    xi_c = x*(1 + mc2/Q2)
    if (xi_c <= 1.0_dp) then
      xq(4) = xs(xi_c, Q2) ! FIXME?
    else
      xq(4) = 0.0_dp
    end if
    xqbar(4) = xq(4)

    xq = Q2/(Q2 + C)*xq
    xqbar = Q2/(Q2 + C)*xqbar
  end subroutine calc_pdf

  function F2_emu(x, Q2) result(F2)
    real(dp) :: F2
    real(dp), intent(in) :: x, Q2

    real(dp), dimension(nf) :: xq, xqbar

    call calc_pdf(x, Q2, xq, xqbar)
    F2 = sum(ei2*(xq + xqbar))
  end function F2_emu

  function F2_antineutrino(x, Q2) result(F2)
    real(dp) :: F2
    real(dp), intent(in) :: x, Q2

    real(dp), dimension(nf) :: xq, xqbar

    call calc_pdf(x, Q2, xq, xqbar)
    F2 = 2*(xq(1) + xqbar(2) + xqbar(3) + xq(4))
  end function F2_antineutrino

  function F2_neutrino(x, Q2) result(F2)
    real(dp) :: F2
    real(dp), intent(in) :: x, Q2

    real(dp), dimension(nf) :: xq, xqbar

    call calc_pdf(x, Q2, xq, xqbar)
    F2 = 2*(xq(2) + xqbar(1) + xqbar(4) + xq(3))
  end function F2_neutrino

  function corr_fac_2xF1_F2(x, Q2) result(corr_fac)
    real(dp) :: corr_fac != 2 x F1/F2
    real(dp), intent(in) :: x, Q2

    corr_fac = (1 + 4*(M*x)**2/Q2)/(1 + R(x, Q2))
  end function corr_fac_2xF1_F2

  function xF3_antineutrino(x, Q2) result(xF3)
    real(dp) :: xF3
    real(dp), intent(in) :: x, Q2

    real(dp), dimension(nf) :: xq, xqbar

    call calc_pdf(x, Q2, xq, xqbar)
    xF3 = 2*(xq(1) - xqbar(2) - xqbar(3) + xq(4))
  end function xF3_antineutrino

  function xF3_neutrino(x, Q2) result(xF3)
    real(dp) :: xF3
    real(dp), intent(in) :: x, Q2

    real(dp), dimension(nf) :: xq, xqbar

    call calc_pdf(x, Q2, xq, xqbar)
    xF3 = 2*(xq(2) - xqbar(1) - xqbar(4) + xq(3))
  end function xF3_neutrino
end module by03
