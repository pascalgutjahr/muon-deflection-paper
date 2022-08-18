module grv98
  ! Glück, Reya, Vogt, 1998 LO PDFs via PDFLIB
  ! EPJC 5 (1998) 461
  use iso_fortran_env, only: dp=>real64
  implicit none
  private
  public :: xuv, xdv, xDelta, xg, xud, xs
  logical, save :: initialized = .false.
contains
  subroutine init()
    character(len=20), dimension(20) :: parm
    real(dp), dimension(20) :: val

    parm(1) = 'Init0'; val(1) = 0.0_dp
    call pdfset(parm, val)

    parm(1) = 'Nptype'; val(1) = 1
    parm(2) = 'Ngroup'; val(2) = 5
    parm(3) = 'Nset';   val(3) = 12
    call pdfset(parm, val)

    initialized = .true.
  end subroutine

  function xuv(x, Q2)
    real(dp) :: xuv
    real(dp), intent(in) :: x, Q2

    real(dp) :: xdv, xubar, xdbar, xs_, xc, xb, xt, xg_

    if (.not. initialized) call init
    call structm(x, Q2, xuv, xdv, xubar, xdbar, xs_, xc, xb, xt, xg_)
  end function xuv

  function xdv(x, Q2)
    real(dp) :: xdv
    real(dp), intent(in) :: x, Q2

    real(dp) :: xuv_, xubar, xdbar, xs_, xc, xb, xt, xg_

    if (.not. initialized) call init
    call structm(x, Q2, xuv_, xdv, xubar, xdbar, xs_, xc, xb, xt, xg_)
  end function xdv

  function xDelta(x, Q2)
    real(dp) :: xDelta
    real(dp), intent(in) :: x, Q2

    real(dp) :: xuv_, xdv, xubar, xdbar, xs_, xc, xb, xt, xg_

    if (.not. initialized) call init
    call structm(x, Q2, xuv_, xdv, xubar, xdbar, xs_, xc, xb, xt, xg_)
    xDelta = xubar - xdbar
  end function xDelta

  function xg(x, Q2) result(xw)
    real(dp) :: xw
    real(dp), intent(in) :: x, Q2

    real(dp) :: xuv_, xdv, xubar, xdbar, xs_, xc, xb, xt

    if (.not. initialized) call init
    call structm(x, Q2, xuv_, xdv, xubar, xdbar, xs_, xc, xb, xt, xw)
  end function xg

  function xud(x, Q2) result(xw)
    real(dp) :: xw
    real(dp), intent(in) :: x, Q2

    real(dp) :: xuv_, xdv, xubar, xdbar, xs_, xc, xb, xt, xg_

    if (.not. initialized) call init
    call structm(x, Q2, xuv_, xdv, xubar, xdbar, xs_, xc, xb, xt, xg_)
    xw = xubar + xdbar
  end function xud

  function xs(x, Q2) result(xwp)
    real(dp) :: xwp
    real(dp), intent(in) :: x, Q2

    real(dp) :: xuv_, xdv, xubar, xdbar, xc, xb, xt, xg_

    call structm(x, Q2, xuv_, xdv, xubar, xdbar, xwp, xc, xb, xt, xg_)
  end function xs
end module grv98
