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

    implicit none

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
    ans_box = reshape((/22.734438235370995d0, 22.734438235370995d0, 22.734438235370995d0, 0.0d0, 0.0d0, 0.0d0/), shape(ans_box))
    call check(box, ans_box, passed, total)

    ! TEST 5
    box = trj%box(5)
    ans_box = reshape((/22.725949956792412d0, 22.416050319868443d0, 23.050339930644768d0, &
        -1.2781211362092404d-003, 1.4444318634328334d-002, -1.3460336751738845d-002/), &
        shape(ans_box))
    call check(box, ans_box, passed, total)

    ! TEST 5
    box = trj%box(trj%NFRAMES)
    ans_box = reshape((/22.825135714633518d0, 22.312977215060247d0, 23.063685023095289d0, &
        -2.1286804723239386d-003, -4.5563434429383672d-004,  -1.0722788675108733d-002 /), &
        shape(ans_box))
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
    a1 = trj%timestep
    a2 = 0.012
    call check(a1, a2, passed, total)

    ! TEST 11
    call trj%read(dcdfile, every=2)
    x = trj%x(5, 5)
    ans = [-11.111647, -9.451539, 8.782396]
    call check(x, ans, passed, total)

    ! TEST 12
    call trj%read(dcdfile, skip=4)
    x = trj%x(6, 5)
    ans = [-11.111647, -9.451539, 8.782396]
    call check(x, ans, passed, total)

    ! TEST 13
    call trj%read(dcdfile, skip=4, every=2)
    x = trj%x(3, 5)
    ans = [-11.111647, -9.451539, 8.782396]
    call check(x, ans, passed, total)

    ! TEST 13
    call trj%read(dcdfile, skip=5, every=5)
    x = trj%x(1, 5)
    ans = [-11.111647, -9.451539, 8.782396]
    call check(x, ans, passed, total)

    call finished_tests(passed, total)

end program dcdfile_test
