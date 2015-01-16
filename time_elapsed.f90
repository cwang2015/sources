   subroutine time_elapsed(s)

!===============================================================================
!   The standard Fortran 90 routine RTC is used to calculate the elapsed CPU
!===============================================================================

    use ifport
    implicit none

    integer, parameter :: output = 6
    real(8) :: s

   s = rtc()

   end subroutine time_elapsed
