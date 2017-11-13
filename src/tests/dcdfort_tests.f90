!
! This file is part of libdcdfort
! https://github.com/wesbarnett/dcdfort
!
! Copyright (c) 2017 James W. Barnett
!
! This program is free software; you can redistribute integer(kind=int32) and/or modify
! integer(kind=int32) under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 2 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that integer(kind=int32) will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License along
! with this program; if not, write to the Free Software Foundation, Inc.,
! 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

module dcdfort_tests

    use dcdfort_common
    use dcdfort_trajectory

    implicit none

    type(Trajectory) :: trj
    real(kind=real32), parameter :: pi_2 = acos(0.0)
    character(len=8), parameter :: dcdfile = "test.dcd"
    character(len=8), parameter :: ndxfile = "test.ndx"
    real(kind=real32), parameter :: sp_tol = 1e-5
    real(kind=real64), parameter :: dp_tol = 1e-16
    real(kind=real32) :: x(3), y(3), z(3), w(3), ans(3), a1, a2
    real(kind=real64) :: box(6), ans_box(6), b, c(3), ans64(3)
    integer(kind=int32) :: passed = 0, total = 0, a, ans_val, n

    interface check
        module procedure check_int, check32, check64, check_array32, check_array64
    end interface check

contains

    subroutine do_output(total, passed, test_result)

        implicit none
        integer(kind=int32), intent(inout) :: total, passed
        logical, intent(in) :: test_result
        integer(kind=int32) :: I
        character (len=6) :: passfail(0:1) = ["FAILED", "PASSED"]
        
        I = merge(1, 0, test_result)
        total = total + 1
        passed = passed + I

        write(output_unit, '(a,i2,a,a)') "TEST ", total, ": ", passfail(I)

    end subroutine do_output

    subroutine check_int(x, y, passed, total)

        implicit none

        integer(kind=int32), intent(inout) :: total, passed
        integer(kind=int32), intent(in) :: x, y

        call do_output(total, passed, x == y)

    end subroutine check_int

    subroutine check32(x, y, passed, total)

        implicit none

        integer(kind=int32), intent(inout) :: total, passed
        real(kind=real32), intent(in) :: x, y

        call do_output(total, passed, x-y < sp_tol)

    end subroutine check32

    subroutine check64(x, y, passed, total)

        implicit none

        integer(kind=int32), intent(inout) :: total, passed
        real(kind=real64), intent(in) :: x, y

        call do_output(total, passed, x - y < dp_tol)

    end subroutine check64

    subroutine check_array32(x, y, passed, total)

        implicit none

        integer(kind=int32), intent(inout) :: total, passed
        real(kind=real32), intent(in) :: x(:), y(:)

        call do_output(total, passed, all((x - y) < sp_tol))

    end subroutine check_array32


    subroutine check_array64(x, y, passed, total)

        implicit none

        integer(kind=int32), intent(inout) :: total, passed
        real(kind=real64), intent(in) :: x(:), y(:)

        call do_output(total, passed, all((x - y) < dp_tol))

    end subroutine check_array64

    subroutine finished_tests(passed, total)

        implicit none

        integer(kind=int32), intent(inout) :: total, passed
        
        write(output_unit,*)
        write(output_unit,'(a,i0,a,i0,a)') "Passed ", passed, " out of ", total, " tests"
        write(output_unit,*)

        if (passed .ne. total) then
            write(output_unit, '(a)') "WARNING: Some tests failed!"
            stop 1
        end if

    end subroutine finished_tests

end module dcdfort_tests
