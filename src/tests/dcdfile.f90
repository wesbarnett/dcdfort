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
    ans = [10.547903, -6.632888, 8.799228]
    call check(x, ans, passed, total)

    ! TEST 2
    x = trj%x(5, 100)
    ans = [-2.448210, 7.568568, -0.057638]
    call check(x, ans, passed, total)

    ! TEST 3
    x = trj%x(trj%NFRAMES, trj%natoms())
    ans = [-4.320391, -5.902228, -8.670344]
    call check(x, ans, passed, total)

    ! TEST 4
    box = trj%box(1)
    ans_box = reshape((/22.7344379, 22.7344379, 22.7344379, 90.0000000, 90.0000000, 90.0000000/), shape(ans_box))
    call check(box, ans_box, passed, total)

    ! TEST 5
    box = trj%box(5)
    ans_box = reshape((/22.7259502, 22.4160500, 23.0503407, 90.0732346, 89.1723709, 90.7712402/), shape(ans_box))
    call check(box, ans_box, passed, total)

    ! TEST 5
    box = trj%box(trj%NFRAMES)
    ans_box = reshape((/22.8251362, 22.3129768, 23.0636845, 90.1219635, 90.0261078, 90.6143799/), shape(ans_box))
    call check(box, ans_box, passed, total)

    ! TEST 6
    n = trj%nframes
    a = 11
    call check(n, a, passed, total)

    ! TEST 7
    n = trj%nevery
    a = 1000
    call check(n, a, passed, total)

    ! TEST 8
    n = trj%istart
    a = 0
    call check(n, a, passed, total)

    ! TEST 9
    n = trj%iend
    a = 10000
    call check(n, a, passed, total)

    ! TEST 10
    b = trj%timestep
    c = 0.012
    call check(b, c, passed, total)

    call finished_tests(passed, total)

end program dcdfile_test
