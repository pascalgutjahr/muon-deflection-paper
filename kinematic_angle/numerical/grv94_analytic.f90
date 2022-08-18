module grv94
  ! Glück, Reya, Vogt, Dynamical parton distributions nof the proton and
  ! small-x physics. Z. Phys. C 67, 433-447 (1995)
  use iso_fortran_env, only: dp=>real64
  implicit none
  private
  public :: xuv, xdv, xDelta, xg, xud, xs
  real(dp), parameter :: mu2_LO = 0.23_dp
contains
  subroutine calc_s(Q2, s)
    real(dp), intent(in) :: Q2
    real(dp), intent(out) :: s

    real(dp), parameter :: aux = 0.232_dp**2

    s = log(log(Q2/aux)/log(mu2_LO/aux))
  end subroutine calc_s

  function xuv(x, Q2)
    real(dp) :: xuv
    real(dp), intent(in) :: x, Q2

    real(dp) :: N, A, B, C, D, aa, bb, s

    call calc_s(Q2, s)
    aa = 0.590_dp - 0.024_dp*s
    bb = 0.131_dp + 0.063_dp*s
    N = 2.284_dp + 0.802_dp*s + 0.055_dp*s**2
    A = -0.449_dp - 0.138_dp*s - 0.076_dp*s**2
    B = 0.213_dp + 2.669_dp*s - 0.728_dp*s**2
    C = 8.854_dp - 9.135_dp*s + 1.979_dp*s**2
    D = 2.997_dp + 0.753_dp*s - 0.076_dp*s**2
    xuv = N*x**aa*(1 + A*x**bb + B*x + C*x**1.5_dp)*(1 - x)*D
  end function xuv

  function xdv(x, Q2)
    real(dp) :: xdv
    real(dp), intent(in) :: x, Q2

    real(dp) :: N, A, B, C, D, aa, bb, s

    call calc_s(Q2, s)
    aa = 0.376_dp
    bb = 0.486_dp + 0.062_dp*s
    N = 0.371_dp + 0.083_dp*s + 0.039_dp*s**2
    A = -0.509_dp + 3.310_dp*s - 1.248_dp*s**2
    B = 12.41_dp - 10.52_dp*s + 2.267_dp*s**2
    C = 6.373_dp - 6.208_dp*s + 1.418_dp*s**2
    D = 3.691_dp + 0.799_dp*s - 0.071_dp*s**2
    xdv = N*x**aa*(1 + A*x**bb + B*x + C*x**1.5_dp)*(1 - x)*D
  end function xdv

  function xDelta(x, Q2)
    real(dp) :: xDelta
    real(dp), intent(in) :: x, Q2

    real(dp) :: N, A, B, C, D, aa, bb, s

    call calc_s(Q2, s)
    aa = 0.409_dp - 0.005_dp*s
    bb = 0.799_dp + 0.071_dp*s
    N = 0.082_dp + 0.014_dp*s + 0.008_dp*s**2
    A = -38.07_dp + 36.13_dp*s - 0.656_dp*s**2
    B = 90.31_dp - 74.15_dp*s + 7.645_dp*s**2
    C = 0.0_dp
    D = 7.486_dp + 1.217_dp*s - 0.159_dp*s**2
    xDelta = N*x**aa*(1 + A*x**bb + B*x + C*x**1.5_dp)*(1 - x)*D
  end function xDelta

  function xg(x, Q2) result(xw)
    real(dp) :: xw
    real(dp), intent(in) :: x, Q2

    real(dp) :: aa, A, B, C, bb, alpha, E, Ep, beta, D, s

    call calc_s(Q2, s)
    alpha = 0.524_dp; beta = 1.088_dp
    aa = 1.742_dp - 0.930_dp*s; bb = -0.399_dp*s**2
    A = 7.486_dp - 2.185_dp*s
    B = 16.69_dp - 22.74_dp*s + 5.779_dp*s**2
    C = -25.59_dp + 29.71_dp*s - 7.296_dp*s**2
    D = 2.792_dp + 2.215_dp*s + 0.422_dp*s**2 - 0.104_dp*s**3
    E = 0.807_dp + 2.005_dp*s; Ep = 3.841_dp + 0.316_dp*s
    xw = (x*aa*(A + B*x + C*x**2)*log(1/x)**bb + s**alpha &
      *exp(-E + sqrt(Ep*s**beta*log(1/x))))*(1 - x)**D
  end function xg

  function xud(x, Q2) result(xw)
    real(dp) :: xw
    real(dp), intent(in) :: x, Q2

    real(dp) :: aa, A, B, C, bb, alpha, E, Ep, beta, D, s

    call calc_s(Q2, s)
    alpha = 1.451_dp; beta = 0.271_dp
    aa = 0.410_dp - 0.232_dp*s; bb = 0.534_dp - 0.457*s
    A = 0.890_dp - 0.140_dp*s; B = -0.981_dp
    C = 0.320_dp + 0.683_dp*s
    D = 4.752_dp + 1.164_dp*s + 0.286_dp*s**2
    E = 4.119_dp + 1.713_dp*s; Ep = 0.682_dp + 2.978_dp*s
    xw = (x*aa*(A + B*x + C*x**2)*log(1/x)**bb + s**alpha &
      *exp(-E + sqrt(Ep*s**beta*log(1/x))))*(1 - x)**D
  end function xud

  function xs(x, Q2) result(xwp)
    real(dp) :: xwp
    real(dp), intent(in) :: x, Q2

    real(dp) :: alpha, aa, A, B, D, E, Ep, beta, s

    call calc_s(Q2, s)
    alpha = 0.914_dp; beta = 0.577_dp
    aa = 1.798_dp - 0.596_dp*s
    A = -5.548_dp + 3.669_dp*sqrt(s) - 0.616_dp*s
    B = 18.92_dp - 16.73_dp*sqrt(s) + 5.168_dp*s
    D = 6.379_dp - 0.350_dp*s + 0.142_dp*s**2
    E = 3.981_dp + 1.638_dp*s; Ep = 6.402_dp
    xwp = s**alpha/log(1/x)**aa*(1 + A*sqrt(x) + B*x)*(1 - x)**D &
      *exp(-E + sqrt(Ep*s**beta*log(1/x)))
  end function xs
end module grv94
