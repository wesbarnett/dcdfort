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

program ndxfile_grp_test 

    use dcdfort_tests

    implicit none

    ! TEST 1
    call trj%read(dcdfile, ndxfile, "gasGroup")
    x = trj%x(10, 10)
    ans = [-1.284096, -5.972089, 2.997065]
    call check(x, ans, passed, total)

    ! TEST 2
    a = trj%natoms()
    ans_val = 100
    call check(a, ans_val, passed, total) 

    call finished_tests(passed, total)

end program ndxfile_grp_test
