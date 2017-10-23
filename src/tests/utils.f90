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

    ! TODO with alpha, beta, gamma
    ! TEST 1
!   x = [5.5, 5.5, 3.5]
!   box = reshape((/5.0, 0.0, 0.0, &
!                   0.0, 5.0, 0.0, &
!                   2.5, 2.5, 3.5/), shape(box))
!   x = pbc(dble(x), dble(box))
!   ans = [-2.000, -2.000, 0.000]
!   call check(x, ans, passed, total) 

    ! TODO with alpha, beta, gamma
    ! TEST 2
!   x = [5.5, 5.5, 3.5]
!   b = distance(dble(x), dble(ans), dble(box))
!   call check(b, 0.0, passed, total)

    ! TEST 1
    x = [5.5, 5.5, 3.5]
    y = [3.6, 4.7, 5.0]
    box = [ 3.5, 4.5, 4.0, 90.0, 90.0, 90.0 ]
    b = distance(dble(x), dble(y), dble(box))
    call check(b, 2.33452, passed, total)

    ! TEST 2
    b = magnitude(dble(x))
    call check(b, 8.52936, passed, total)

    ! TEST 3
    x = [0.0, 0.0, 0.0]
    y = [0.0, 1.0, 0.0]
    z = [1.0, 1.0, 0.0]
    b = bond_angle(dble(x), dble(y), dble(z), dble(box))
    call check(b, real(PI/2.0d0), passed, total)

    ! TEST 4
    w = [1.0, 1.0, 1.0]
    b = dihedral_angle(dble(x), dble(y), dble(z), dble(w), dble(box))
    call check(b, real(PI/2.0d0), passed, total)

    ! TEST 5
    w = [1.0, 1.0, -1.0]
    b = dihedral_angle(dble(x), dble(y), dble(z), dble(w), dble(box))
    call check(b, real(-PI/2.0d0), passed, total)

    ! TEST 6
    x = [0.0, 0.0, 0.0]
    y = [1.1, 1.1, 1.1]
    box = [ 1.0, 1.0, 1.0, 90.0, 90.0, 90.0 ]
    b = distance(dble(x), dble(y), dble(box))
    call check(b, sqrt(3.0*0.1**2), passed, total)

    call finished_tests(passed, total)

end program utils_test
