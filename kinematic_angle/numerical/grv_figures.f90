program grv_figures
  use iso_fortran_env, only: dp=>real64
  use grv94
  use ogpf
  implicit none
  real(dp), parameter :: Q2_val(3) = [5, 15, 35], &
    x_min = 1e-4_dp, x_max = 1e-2_dp
  integer, parameter :: n_pts = 20, nf = 3
  real(dp) :: Q2
  real(dp), dimension(n_pts) :: x, F2
  real(dp), dimension(nf), parameter :: ei2 = ([2, 1, 1]/3.0_dp)**2
  integer :: j, k
  type(gpf) :: gp
  character(len=2) :: q2_str

  do j = 1, size(Q2_val)
    Q2 = Q2_val(j)
    write (q2_str, '(I2)') int(Q2)
    call gp%options('set title "Q2 = ' // q2_str // '"')
    call gp%options('set xlabel "x"; set ylabel "F_2"')
    do k = 1, n_pts
      x(k) = x_min*(x_max/x_min)**((k - 1)/real(n_pts - 1, dp))
      F2(k) = F2_emu(x(k), Q2)
    end do
    call gp%semilogx(x, F2)
  end do
contains
  subroutine calc_pdf(x, Q2, xq, xqbar)
    real(dp), intent(in) :: x, Q2
    real(dp), dimension(nf), intent(out) :: xq, xqbar

    real(dp) :: uv, dv, us, ds

    uv = xuv(x, Q2)
    dv = xdv(x, Q2)
    us = (xud(x, Q2) - xDelta(x, Q2))/2
    ds = (xud(x, Q2) + xDelta(x, Q2))/2

    xqbar(1) = us
    xq(1) = uv + xqbar(1)

    xqbar(2) = ds
    xq(2) = dv + xqbar(2)

    xq(3) = xs(x, Q2)
    xqbar(3) = xq(3)
  end subroutine calc_pdf

  function F2_emu(x, Q2) result(F2)
    real(dp) :: F2
    real(dp), intent(in) :: x, Q2

    real(dp), dimension(nf) :: xq, xqbar

    call calc_pdf(x, Q2, xq, xqbar)
    F2 = sum(ei2*(xq + xqbar))
  end function F2_emu
end program grv_figures
