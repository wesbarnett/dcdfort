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

program read_next_test 

    use dcdfort_tests

    call trj%open(dcdfile)

    ! TEST 1
    a = trj%skip_next(4)
    a = trj%read_next(6)
    x = trj%x(6, 5)
    ans = [-11.111647, -9.451539, 8.782396]
    call check(x, ans, passed, total)

    ! TEST 2
    box = trj%box(1)
    ans_box = reshape((/22.7259502, 22.4160500, 23.0503407, 90.0732346, 89.1723709, 90.7712402/), shape(ans_box))
    call check(box, ans_box, passed, total)

    ! TEST 2
    ans_val = 6
    call check(a, ans_val, passed, total) 

    ! TEST 3
    a = trj%read_next()
    ans_val = 1
    call check(a, ans_val, passed, total) 

    a = trj%skip_next(200)
    ans_val = 0
    call check(a, ans_val, passed, total) 

    ! TEST 4
    a = trj%read_next()
    ans_val = 0
    call check(a, ans_val, passed, total) 
    call trj%close()

    call finished_tests(passed, total)

end program read_next_test
