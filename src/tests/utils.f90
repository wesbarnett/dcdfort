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

program utils_test 

    use dcdfort_utils, only: distance, magnitude, bond_angle, dihedral_angle, pbc
    use dcdfort_tests

    implicit none

    ! TEST 1
    x = [5.5, 5.5, 3.5]
    y = [3.6, 4.7, 5.0]
    box = [ 3.5d0, 4.5d0, 4.0d0, 0.0d0, 0.0d0, 0.0d0 ]
    b = distance(dble(x), dble(y), box)
    call check(real(b), 2.33452, passed, total)

    ! TEST 2
    b = magnitude(dble(x))
    call check(real(b), 8.52936, passed, total)

    ! TEST 3
    x = [0.0, 0.0, 0.0]
    y = [0.0, 1.0, 0.0]
    z = [1.0, 1.0, 0.0]
    b = bond_angle(dble(x), dble(y), dble(z), box)
    call check(real(b), pi_2, passed, total)

    ! TEST 4
    w = [1.0, 1.0, 1.0]
    b = dihedral_angle(dble(x), dble(y), dble(z), dble(w), box)
    call check(real(b), pi_2, passed, total)

    ! TEST 5
    w = [1.0, 1.0, -1.0]
    b = dihedral_angle(dble(x), dble(y), dble(z), dble(w), box)
    call check(real(b), -pi_2, passed, total)

    ! TEST 6
    x = [0.0, 0.0, 0.0]
    y = [1.1, 1.1, 1.1]
    box = [ 1.0d0, 1.0d0, 1.0d0, 0.0d0, 0.0d0, 0.0d0 ]
    b = distance(dble(x), dble(y), box)
    call check(real(b), sqrt(3.0*0.1**2), passed, total)

    ! TEST 7
    x = [0.0, 0.0, 0.0]
    y = [1.1, 1.1, 1.1]
    box = [ 1.0d0, 1.0d0, 1.0d0, 0.0d0, 0.0d0, 0.0d0 ]
    b = distance(dble(x), dble(y), box)
    call check(real(b), sqrt(3.0*0.1**2), passed, total)

    ! TEST 8
    x = [0.0, 0.0, 0.0]
    y = [1.1, 1.0, 1.0]
    box = [ 1.0d0, 1.0d0, 1.0d0, 0.0d0, 0.0d0, 0.0d0 ]
    b = distance(dble(x), dble(y), box)
    call check(real(b), 0.1, passed, total)

    ! TEST 9
    x = [0.0, 0.0, 0.0]
    y = [1.0, 1.1, 1.0]
    box = [ 1.0d0, 1.0d0, 1.0d0, 0.0d0, 0.0d0, 0.0d0 ]
    b = distance(dble(x), dble(y), box)
    call check(real(b), 0.1, passed, total)

    ! TEST 10
    x = [0.0, 0.0, 0.0]
    y = [0.0, 0.7, 0.0]
    box = [ 1.0d0, 1.0d0, 1.0d0, 0.5d0, 0.0d0, 0.0d0 ]
    b = distance(dble(x), dble(y), box)
    call check(real(b), 0.3, passed, total)

    ! TEST 11
    x = [0.0, 0.0, 0.0]
    y = [0.6, 0.966025, 0.0]
    box = [ 1.0d0, 1.0d0, 1.0d0, 0.0d0, 0.0d0, 0.5d0 ]
    c = pbc(dble(y)-dble(x), box)
    ans64 = [0.1d0, 0.1d0, 0.0d0]
    call check(real(c), real(ans64), passed, total) 

    ! TEST 12
    x = [0.0, 0.0, 0.0]
    y = [0.4, 0.533, 0.0]
    box = [ 1.0d0, 0.5d0, 1.0d0, 0.0d0, 0.0d0, 0.5d0 ]
    c = pbc(dble(y)-dble(x), box)
    ans64 = [0.15d0, 0.1d0, 0.0d0]
    call check(real(c), real(ans64), passed, total) 

    ! TEST 13
    x = [0.0, 0.0, 0.0]
    y = [0.6, 0.6, 0.0]
    box = [ 1.0d0, 1.0d0, 1.0d0, 0.0d0, 0.5d0, 0.0d0 ]
    c = pbc(dble(y)-dble(x), box)
    ans64 = [-0.4d0, -0.4d0, 0.0d0]
    call check(real(c), real(ans64), passed, total) 

    ! TEST 14
    x = [0.0, 0.0, 0.0]
    y = [0.6, 0.6, 0.566]
    box = [ 1.0d0, 1.0d0, 1.0d0, 0.0d0, 0.5d0, 0.0d0 ]
    c = pbc(dble(y)-dble(x), box)
    ans64 = [0.1d0, -0.4d0, -0.3d0]
    call check(real(c), real(ans64), passed, total) 

    ! TEST 15
    x = [0.0, 0.0, 0.0]
    y = [0.0, 0.0, 0.566]
    box = [ 1.0d0, 1.0d0, 1.0d0, 0.5d0, 0.0d0, 0.0d0 ]
    c = pbc(dble(y)-dble(x), box)
    ans64 = [0.0d0, 0.5d0, -0.3d0]
    call check(real(c), real(ans64), passed, total) 

    ! TEST 16
    x = [0.0, 0.0, 0.0]
    y = [0.1, 0.5, -0.6]
    box = [ 0.95d0, 0.9d0, 0.8d0, 0.3d0, -0.6d0, -0.3d0 ]
    c = pbc(dble(y)-dble(x), box)
    ans64 = [-0.11d0, -0.2579d0, 3.2038d-02]
    call check(real(c), real(ans64), passed, total) 

    call finished_tests(passed, total)

end program utils_test
