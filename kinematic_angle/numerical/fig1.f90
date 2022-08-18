program fig1
  use by03
  use ogpf
  use iso_fortran_env, only: dp=>real64
  implicit none
  real(dp) :: x
  real(dp), dimension(12), parameter :: xval = [0.070_dp, 0.100_dp, &
    0.140_dp, 0.180_dp, 0.225_dp, 0.275_dp, 0.350_dp, 0.450_dp, 0.550_dp, &
    0.650_dp, 0.750_dp, 0.850_dp]
  real(dp), parameter :: Q2_min = 0.3_dp, Q2_max = 300.0_dp
  integer :: j, k
  integer, parameter :: n_pts = 20
  real(dp), dimension(n_pts) :: Q2, F2
  type(gpf) :: gp
  character(len=5) :: x_str

  do j = 1, 1!size(xval)
    x = xval(j)
    write (x_str, '(F5.3)') x
    call gp%options('set title "x = ' // x_str // '"')
    call gp%options('set xlabel "Q^2"; set ylabel "F_2"')
    do k = 1, n_pts
      Q2(k) = Q2_min*(Q2_max/Q2_min)**((k - 1)/real(n_pts - 1, dp))
      F2(k) = F2_emu(x, Q2(k))
    end do
    call gp%semilogx(Q2, F2)
  end do
end program fig1
