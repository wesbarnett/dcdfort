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
    real(kind=real64), parameter :: PI = 2.0d0*dacos(0.0d0)
    character(len=8), parameter :: dcdfile = "test.dcd"
    character(len=8), parameter :: ndxfile = "test.ndx"
    real(kind=real32) :: x(3), y(3), z(3), w(3), ans(3), box(6), ans_box(6), b, c
    integer(kind=int32) :: passed = 0, total = 0, a, ans_val, n

    interface check
        module procedure check_int, check_real, check_array, check_array_2d
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

        call do_output(total, passed, x .eq. y)

    end subroutine check_int

    subroutine check_real(x, y, passed, total)

        implicit none

        integer(kind=int32), intent(inout) :: total, passed
        real(kind=real32), intent(in) :: x, y
        real(kind=real32) :: tol = 1e-4

        call do_output(total, passed, abs(x-y) .le. tol)

    end subroutine check_real

    subroutine check_array(x, y, passed, total)

        implicit none

        integer(kind=int32), intent(inout) :: total, passed
        real(kind=real32), intent(in) :: x(:), y(:)
        real(kind=real32) :: tol = 1e-4

        call do_output(total, passed, all(abs(x - y) .le. tol))

    end subroutine check_array

    subroutine check_array_2d(x, y, passed, total)

        implicit none

        integer(kind=int32), intent(inout) :: total, passed
        real(kind=real32), intent(in) :: x(:,:), y(:,:)
        real(kind=real32) :: tol = 1e-6

        call do_output(total, passed, all(abs(x - y) .le. tol))

    end subroutine check_array_2d

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
