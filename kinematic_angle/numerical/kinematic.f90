module kinematic
  ! Kinematic relations
  use iso_fortran_env, only: dp=>real64
  implicit none
  real(dp), parameter :: mp = 0.938272046_dp
contains
  function xBj(E, y, Q2)
    real(dp) :: xBj
    real(dp), intent(in) :: E, y, Q2

    xBj = Q2/(2*mp*E*y)
  end function xBj

  function cos_th(E, y, Q2, ml)
    real(dp) :: cos_th
    real(dp), intent(in) :: E, y, Q2, ml

    !cos_th = 1 - Q2/(2*E**2*(1 - y))
    cos_th = (2*E**2*(1 - y) - Q2 - ml**2)/(2*E*sqrt(E**2 - ml**2))
  end function cos_th

  function Q2_max(E, y)
    real(dp) :: Q2_max
    real(dp), intent(in) :: E, y

    Q2_max = 2*mp*E*y
  end function Q2_max

  function Q2_min(E, y, ml)
    real(dp) :: Q2_min
    real(dp), intent(in) :: E, y, ml

    Q2_min = ml**2*y/(1 - y)
  end function Q2_min

  function y_min(E, ml)
    real(dp) :: y_min
    real(dp), intent(in) :: E, ml

    real(dp), parameter :: m_pion = 134.9766e-3_dp

    y_min = (m_pion + m_pion**2/mp)/E
  end function y_min

  function y_max(E, ml)
    real(dp) :: y_max
    real(dp), intent(in) :: E, ml

    y_max = 1 - ml/E
  end function y_max
end module kinematic
