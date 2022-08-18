module bodek2021
  ! Update Bodek-Yang model from arXiv:2108.09240, first iteration
  use iso_fortran_env, only: dp=>real64
  use grv98
  implicit none
  private
  public :: F2_emu, F2_neutrino, F2_antineutrino, corr_fac_2xF1_F2, &
    xF3_neutrino, xF3_antineutrino
  real(dp), parameter :: M = 0.938272046_dp, mc2 = 1.32_dp**2, &
    sin2_thc = 0.2249_dp, cos2_thc = 1 - sin2_thc
  real(dp), parameter :: A = 0.419_dp, B = 0.223_dp, C_v1 = 0.544_dp, &
    C_v2 = 0.431_dp, C_sea = 0.380_dp, N = 1.011_dp
  integer, parameter :: nf = 3 ! u, d, s
  real(dp), parameter :: ei2(nf) = ([2, 1, 1]/3.0_dp)**2
contains
  subroutine calc_xiw(x, Q2, xiw)
    real(dp), intent(in) :: x, Q2
    real(dp), intent(out) :: xiw

    xiw = 2*x*(Q2 + B)/(Q2*(1 + sqrt(1 + 4*(M*x)**2/Q2)) + 2*A*x)
  end subroutine calc_xiw

  function GD(Q2)
    real(dp) :: GD, Q2

    GD = 1/(1 + Q2/0.71_dp)**2
  end function GD

  subroutine calc_pdf(xi_w, Q2, xq, xqbar)
    real(dp), intent(in) :: xi_w, Q2
    real(dp), dimension(nf), intent(out) :: xq, xqbar

    real(dp) :: Q2_08, K_vector_sea, K_vector_valence
    real(dp) :: u_valence, d_valence, u_sea, ubar_sea, d_sea, dbar_sea, &
      s_sea, sbar_sea
    real(dp) :: d_sea_grv98, dbar_sea_grv98, u_sea_grv98, ubar_sea_grv98, &
      d_valence_grv98, u_valence_grv98, u_d, Delta

    Q2_08 = max(Q2, 0.8_dp)

    u_d = xud(xi_w, Q2_08)
    Delta = xDelta(xi_w, Q2_08)
    dbar_sea_grv98 = (u_d - Delta)/2
    d_sea_grv98 = dbar_sea_grv98
    ubar_sea_grv98 = (u_d + Delta)/2
    u_sea_grv98 = ubar_sea_grv98
    d_valence_grv98 = xdv(xi_w, Q2_08)
    u_valence_grv98 = xuv(xi_w, Q2_08)

    d_sea = 1.05_dp*d_sea_grv98
    dbar_sea = 1.05_dp*dbar_sea_grv98
    u_sea = 1.05_dp*u_sea_grv98
    ubar_sea = 1.05_dp*ubar_sea_grv98
    d_valence = d_valence_grv98 - 0.05_dp*(d_sea_grv98 + dbar_sea_grv98)
    u_valence = u_valence_grv98 - 0.05_dp*(u_sea_grv98 + ubar_sea_grv98)

    s_sea = xs(xi_w, Q2_08)
    sbar_sea = s_sea

    K_vector_sea = Q2/(Q2 + C_sea)
    K_vector_valence = (1 - GD(Q2)**2)*((Q2 + C_v2)/(Q2 + C_v1))

    xq(1) = N*(K_vector_valence*u_valence + K_vector_sea*u_sea)
    xq(2) = N*(K_vector_valence*d_valence + K_vector_sea*d_sea)
    xq(3) = N*K_vector_sea*s_sea

    xqbar(1) = N*K_vector_sea*ubar_sea
    xqbar(2) = N*K_vector_sea*dbar_sea
    xqbar(3) = N*K_vector_sea*sbar_sea
  end subroutine calc_pdf

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

  function F2_emu(x, Q2) result(F2)
    real(dp) :: F2
    real(dp), intent(in) :: x, Q2

    real(dp), dimension(nf) :: xq, xqbar
    real(dp) :: xi_w

    call calc_xiw(x, Q2, xi_w)
    call calc_pdf(xi_w, Q2, xq, xqbar)
    F2 = sum(ei2*(xq + xqbar))
  end function F2_emu

  function F2_antineutrino(x, Q2) result(F2)
    real(dp) :: F2
    real(dp), intent(in) :: x, Q2

    real(dp), dimension(nf) :: xq, xqbar
    real(dp) :: F2_ncp, F2_cp, xi_w

    call calc_pdf(x, Q2, xq, xqbar)
    
    call calc_xiw(x, Q2, xi_w)
    call calc_pdf(xi_w, Q2, xq, xqbar)
    F2_ncp = 2*(xq(1) + xqbar(2) + xqbar(3))

    xi_w = 2*x*(Q2 + mc2 + B)/(Q2*(1 + sqrt(1 + 4*(M*x)**2/Q2)) + 2*A*x)
    if (xi_w < 1.0_dp) then
      call calc_pdf(x, Q2, xq, xqbar)
    else
      xq = 0; xqbar = 0
    end if
    F2_cp = 2*(xq(1) + xqbar(2) + xqbar(3))
    F2 = F2_ncp + F2_cp
  end function F2_antineutrino

  function F2_neutrino(x, Q2) result(F2)
    real(dp) :: F2
    real(dp), intent(in) :: x, Q2

    real(dp), dimension(nf) :: xq, xqbar
    real(dp) :: F2_ncp, F2_cp, xi_w

    call calc_xiw(x, Q2, xi_w)
    call calc_pdf(xi_w, Q2, xq, xqbar)
    F2_ncp = 2*((xq(1) + xq(2))/2*cos2_thc + xq(3)*sin2_thc &
      + (xqbar(1) + xqbar(2))/2)

    xi_w = 2*x*(Q2 + mc2 + B)/(Q2*(1 + sqrt(1 + 4*(M*x)**2/Q2)) + 2*A*x)
    if (xi_w < 1.0_dp) then
      call calc_pdf(x, Q2, xq, xqbar)
    else
      xq = 0; xqbar = 0
    end if
    F2_cp = 2*((xq(1) + xq(2))/2*sin2_thc + xq(3)*cos2_thc)
    F2 = F2_ncp + F2_cp
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
    real(dp) :: xF3_ncp, xF3_cp, xi_w, K_charm, H

    call calc_xiw(x, Q2, xi_w)
    call calc_pdf(x, Q2, xq, xqbar)
    xF3_ncp = 2*(xq(1) - xqbar(2) - xqbar(3))

    xi_w = 2*x*(Q2 + mc2 + B)/(Q2*(1 + sqrt(1 + 4*(M*x)**2/Q2)) + 2*A*x)
    K_charm = Q2/(Q2 + mc2)
    if (xi_w < 1.0_dp) then
      call calc_pdf(x, Q2, xq, xqbar)
    else
      xq = 0; xqbar = 0
    end if
    H = 0.914_dp + 0.296_dp*x - 0.374_dp*x**2 + 0.165_dp*x**3
    xF3_cp = H*K_charm*2*(xq(1) - xqbar(2) - xqbar(3))
    xF3 = xF3_ncp + xF3_cp
  end function xF3_antineutrino

  function xF3_neutrino(x, Q2) result(xF3)
    real(dp) :: xF3
    real(dp), intent(in) :: x, Q2

    real(dp), dimension(nf) :: xq, xqbar
    real(dp) :: xF3_ncp, xF3_cp, xi_w, K_charm, H

    call calc_xiw(x, Q2, xi_w)
    call calc_pdf(x, Q2, xq, xqbar)
    xF3_ncp = 2*((xq(1) + xq(2))/2*cos2_thc + xq(3)*sin2_thc &
      - (xqbar(1) + xqbar(2))/2)

    xi_w = 2*x*(Q2 + mc2 + B)/(Q2*(1 + sqrt(1 + 4*(M*x)**2/Q2)) + 2*A*x)
    K_charm = Q2/(Q2 + mc2)
    if (xi_w < 1.0_dp) then
      call calc_pdf(x, Q2, xq, xqbar)
    else
      xq = 0; xqbar = 0
    end if
    H = 0.914_dp + 0.296_dp*x - 0.374_dp*x**2 + 0.165_dp*x**3
    xF3_cp = H*K_charm*2*((xq(1) + xq(2))/2*sin2_thc + xq(3)*cos2_thc)
    xF3 = xF3_ncp + xF3_cp
  end function xF3_neutrino
end module bodek2021
