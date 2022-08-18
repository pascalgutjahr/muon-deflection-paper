module block
  use iso_fortran_env, only: dp=>real64
  implicit none
  private
  public :: F2

  real(dp), parameter :: M2 = 0.753_dp, mu2 = 2.82_dp, &
    a0 = 8.205e-4_dp, a1 = -5.148e-2_dp, a2 = -4.725e-3_dp, &
    b0 = 2.217e-3_dp, b1 = 1.244e-2_dp, b2 = 5.958e-4_dp, &
    c0 = 0.255_dp, c1 = 1.475e-1_dp, &
    n = 11.49_dp, lambda = 2.430_dp
contains
  function F2(x, Q2)
    real(dp) :: F2
    real(dp), intent(in) :: x, Q2

    real(dp) :: A, B, C, D, L_Q2, L_x

    L_Q2 = log(1 + Q2/mu2)
    A = a0 + a1*L_Q2 + a2*L_Q2**2
    B = b0 + b1*L_Q2 + b2*L_Q2**2
    C = c0 + c1*L_Q2

    D = Q2*(Q2 + lambda*M2)/(Q2 + M2)**2

    L_x = log(Q2/x/(Q2 + mu2))
    F2 = D*(1 - x)**n*(C + A*L_x + B*L_x**2)
  end function F2
end module block
