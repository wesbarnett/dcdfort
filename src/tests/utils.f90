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

    use dcdfort_utils
    use dcdfort_tests

    ! TEST 1
    x = [5.5d0, 5.5d0, 3.5d0]
    y = [3.6d0, 4.7d0, 5.0d0]
    box = [ 3.5d0, 4.5d0, 4.0d0, 0.0d0, 0.0d0, 0.0d0 ]
    b = distance(x, y, box)
    call check(b, 2.33452d0, passed, total)

    ! TEST 2
    b = magnitude(x)
    call check(b, 8.52936d0, passed, total)

    ! TEST 3
    x = [0.0d0, 0.0d0, 0.0d0]
    y = [0.0d0, 1.0d0, 0.0d0]
    z = [1.0d0, 1.0d0, 0.0d0]
    b = bond_angle(x, y, z, box)
    call check(b, PI/2.0d0, passed, total)

    ! TEST 4
    w = [1.0d0, 1.0d0, 1.0d0]
    b = dihedral_angle(x, y, z, w, box)
    call check(b, PI/2.0d0, passed, total)

    ! TEST 5
    w = [1.0d0, 1.0d0, -1.0d0]
    b = dihedral_angle(x, y, dble(z), dble(w), box)
    call check(b, -PI/2.0d0, passed, total)

    ! TEST 6
    x = [0.0d0, 0.0d0, 0.0d0]
    y = [1.1d0, 1.1d0, 1.1d0]
    box = [ 1.0d0, 1.0d0, 1.0d0, 0.0d0, 0.0d0, 0.0d0 ]
    b = distance(x, y, box)
    call check(b, dsqrt(3.0d0*0.1d0**2), passed, total)

    ! TEST 7
    x = [0.0d0, 0.0d0, 0.0d0]
    y = [1.1d0, 1.1d0, 1.1d0]
    box = [ 1.0d0, 1.0d0, 1.0d0, 0.0d0, 0.0d0, 0.0d0 ]
    b = distance(x, y, box)
    call check(b, sqrt(3.0d0*0.1d0**2), passed, total)

    ! TEST 8
    x = [0.0d0, 0.0d0, 0.0d0]
    y = [1.1d0, 1.0d0, 1.0d0]
    box = [ 1.0d0, 1.0d0, 1.0d0, 0.0d0, 0.0d0, 0.0d0 ]
    b = distance(x, y, box)
    call check(b, 0.1d0, passed, total)

    ! TEST 9
    x = [0.0d0, 0.0d0, 0.0d0]
    y = [1.0d0, 1.1d0, 1.0d0]
    box = [ 1.0d0, 1.0d0, 1.0d0, 0.0d0, 0.0d0, 0.0d0 ]
    b = distance(x, y, box)
    call check(b, 0.1d0, passed, total)

    ! TEST 10
    x = [0.0d0, 0.0d0, 0.0d0]
    y = [0.0d0, 0.7d0, 0.0d0]
    box = [ 1.0d0, 1.0d0, 1.0d0, 0.5d0, 0.0d0, 0.0d0 ]
    b = distance(x, y, box)
    call check(b, 0.3d0, passed, total)

    ! TEST 11
    x = [0.0d0, 0.0d0, 0.0d0]
    y = [0.6d0, 0.966025d0, 0.0d0]
    box = [ 1.0d0, 1.0d0, 1.0d0, 0.0d0, 0.0d0, 0.5d0 ]
    z = pbc(y-x, box)
    ans = [0.1d0, 0.1d0, 0.0d0]
    call check(z, ans, passed, total) 

    ! TEST 12
    x = [0.0d0, 0.0d0, 0.0d0]
    y = [0.4d0, 0.533d0, 0.0d0]
    box = [ 1.0d0, 0.5d0, 1.0d0, 0.0d0, 0.0d0, 0.5d0 ]
    z = pbc(y-x, box)
    ans = [0.15d0, 0.1d0, 0.0d0]
    call check(z, ans, passed, total) 

    ! TEST 13
    x = [0.0d0, 0.0d0, 0.0d0]
    y = [0.6d0, 0.6d0, 0.0d0]
    box = [ 1.0d0, 1.0d0, 1.0d0, 0.0d0, 0.5d0, 0.0d0 ]
    z = pbc(y-x, box)
    ans = [-0.4d0, -0.4d0, 0.0d0]
    call check(z, ans, passed, total) 

    ! TEST 14
    x = [0.0d0, 0.0d0, 0.0d0]
    y = [0.6d0, 0.6d0, 0.566d0]
    box = [ 1.0d0, 1.0d0, 1.0d0, 0.0d0, 0.5d0, 0.0d0 ]
    z = pbc(y-x, box)
    ans = [0.1d0, -0.4d0, -0.3d0]
    call check(z, ans, passed, total) 

    ! TEST 15
    x = [0.0d0, 0.0d0, 0.0d0]
    y = [0.0d0, 0.0d0, 0.566d0]
    box = [ 1.0d0, 1.0d0, 1.0d0, 0.5d0, 0.0d0, 0.0d0 ]
    z = pbc(y-x, box)
    ans = [0.0d0, 0.5d0, -0.3d0]
    call check(z, ans, passed, total) 

    ! TEST 16
    x = [0.0d0, 0.0d0, 0.0d0]
    y = [0.1d0, 0.5d0, -0.6d0]
    box = [ 0.95d0, 0.9d0, 0.8d0, 0.3d0, -0.6d0, -0.3d0 ]
    z = pbc(y-x, box)
    ans = [-0.11d0, -0.2579d0, 3.2038d-02]
    call check(z, ans, passed, total) 

    call finished_tests(passed, total)

end program utils_test
