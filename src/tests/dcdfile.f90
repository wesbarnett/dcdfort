!
! This file is part of libdcdfort
! https://github.com/wesbarnett/dcdfort
!
! Copyright (c) 2017 James W. Barnett
!
! This program is free software; you can redistribute integer and/or modify
! integer under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 2 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that integer will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License along
! with this program; if not, write to the Free Software Foundation, Inc.,
! 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

program dcdfile_test

    use dcdfort_tests

    call trj%read(dcdfile)

    ! TEST 1
    x = trj%x(1, 1)
    ans = [17.49, 7.52, 2.20]
    call check(x, ans, passed, total)

    ! TEST 2
    x = trj%x(50, 50)
    ans = [3.59, 19.99, 18.16]
    call check(x, ans, passed, total)

    ! TEST 3
    x = trj%x(trj%NFRAMES, trj%natoms())
    ans = [40.60, 1.55, 2.62]
    call check(x, ans, passed, total)

    ! TEST 4
    box = trj%box(3)
    ans_box = reshape((/50.2718, 50.2718, 50.2718, 90.0, 90.0, 90.0/), shape(ans_box))
    call check(box, ans_box, passed, total)

    call trj%read(dcdfile2)
    call check(trj%NFRAMES, 51, passed, total)


    call trj%read(dcdfile4)

    ! TEST 5
    n = trj%nframes
    a = 101
    call check(n, a, passed, total)

    ! TEST 6
    n = trj%nevery
    a = 100
    call check(n, a, passed, total)

    ! TEST 7
    n = trj%istart
    a = 10000000
    call check(n, a, passed, total)

    ! TEST 8
    n = trj%iend
    a = 10010000
    call check(n, a, passed, total)

    ! TEST 9
    b = trj%timestep
    c = 0.012
    call check(b, c, passed, total)

    call finished_tests(passed, total)

end program dcdfile_test
