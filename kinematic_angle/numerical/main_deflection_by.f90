program main_deflection
  use iso_fortran_env, only: dp=>real64
  use multidim_integrate, only: dcuhre
  use kinematic
  use by03
  implicit none
  real(dp), parameter :: E_min = 10.0_dp, E_max = 10e6_dp, &
    ml = 0.105683715_dp, GF = 1.1663787e-5_dp, MW = 80.385_dp, &
    pi = 3.141592653589793238_dp
  integer, parameter :: n_pts = 7
  real(dp) :: E
  integer :: j

  integer, parameter :: ndim = 2, numfun = 3, minpts = 10, maxpts = int(1e6), &
    key = 0, restar = 0
  real(dp), parameter :: epsrel = 1e-3_dp, epsabs = 1e-30_dp
  real(dp), dimension(ndim), parameter :: a = 0.0_dp, b = 1.0_dp
  real(dp), dimension(numfun) :: integral, abserr
  integer :: neval, ifail, nsub

  do j = 1, n_pts
    E = E_min*(E_max/E_min)**((j - 1)/real(n_pts - 1, dp))
    call dcuhre(ndim, numfun, a, b, minpts, maxpts, integrand_neutrino, &
      epsabs, epsrel, key, restar, integral, abserr, neval, ifail, nsub)
    !write (*,*) E, integral
    write (*,*) E, &
      integral(2)/integral(1), &
      sqrt(integral(3)/integral(1) - (integral(2)/integral(1))**2), &
      ifail
  end do
contains
  subroutine integrand_neutrino(numdim, x, nfun, f)
    integer, intent(in) :: numdim, nfun
    real(dp), intent(in) :: x(numdim)
    real(dp), intent(out) :: f(nfun)

    real(dp) :: y, Q2, costh, min_y, max_y, min_Q2, max_Q2, jacobian

    min_y = y_min(E, ml)
    max_y = y_max(E, ml)
    y = min_y*(max_y/min_y)**x(1)
    jacobian = y*log(max_y/min_y)

    min_Q2 = Q2_min(E, y, ml)
    max_Q2 = Q2_max(E, y)
    Q2 = min_Q2*(max_Q2/min_Q2)**x(2)
    jacobian = jacobian * Q2*log(max_Q2/min_Q2)

    costh = acos(cos_th(E, y, Q2, ml))*180/pi

    f(1) = dsigma_neutrino(E, y, Q2)*jacobian
    f(2) = costh*f(1)
    f(3) = costh*f(2)
  end subroutine integrand_neutrino

  function dsigma_neutrino(E, y, Q2) result(dsigma)
    real(dp) :: dsigma
    real(dp), intent(in) :: E, y, Q2

    real(dp) :: x, xF1, F2, xF3

    x = xBj(E, y, Q2)
    F2 = F2_neutrino(x, Q2)
    xF1 = F2*corr_fac_2xF1_F2(x, Q2)/2
    xF3 = xF3_neutrino(x, Q2)
    dsigma = GF**2/(4*pi)*(MW**2/(Q2 + MW**2))**2 &
      /y*((1 - y - (x*y*mp)**2/Q2)*F2 + y**2*xF1 + (y - y**2/2)*xF3)
  end function dsigma_neutrino

  function dsigma_antineutrino(E, y, Q2) result(dsigma)
    real(dp) :: dsigma
    real(dp), intent(in) :: E, y, Q2

    real(dp) :: x, xF1, F2, xF3

    x = xBj(E, y, Q2)
    F2 = F2_antineutrino(x, Q2)
    xF1 = F2*corr_fac_2xF1_F2(x, Q2)/2
    xF3 = xF3_antineutrino(x, Q2)
    dsigma = GF**2/(4*pi)*(MW**2/(Q2 + MW**2))**2 &
      /y*((1 - y - (x*y*mp)**2/Q2)*F2 + y**2*xF1 - (y - y**2/2)*xF3)
  end function dsigma_antineutrino
end program main_deflection
