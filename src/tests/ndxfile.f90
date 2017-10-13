!
! This file is part of libdcdfort
! https://github.com/wesbarnett/dcdfort
!
! Copyright (c) 2016,2017 James W. Barnett
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

program ndxfile_test 

    use dcdfort_tests
    use dcdfort_index

    type(IndexFile) :: ndx

    call trj%read(dcdfile, ndxfile)

    ! TEST 1
    x = trj%x(1, 10, "gasGroup")
    ans = [0.128691, -8.314034, 6.497755]
    call check(x, ans, passed, total)

    x = trj%x(5, trj%natoms("monGroup")-5, "monGroup")
    ans = [-10.029089, -1.459442, 5.690279]
    call check(x, ans, passed, total)

    x = trj%x(trj%NFRAMES, trj%natoms("gasGroup"), "gasGroup")
    ans = [-4.320391, -5.902228, -8.670344]
    call check(x, ans, passed, total)

    ! TEST 2
    a = trj%natoms()
    ans_val = 10100
    call check(a, ans_val, passed, total) 

    ! TEST 3
    a = trj%natoms("System")
    ans_val = 10100
    call check(a, ans_val, passed, total) 

    ! TEST 4
    a = trj%natoms("gasGroup")
    ans_val = 100
    call check(a, ans_val, passed, total) 

    ! TEST 5
    a = trj%natoms("FOOBAR")
    ans_val = 0
    call check(a, ans_val, passed, total) 

    ! TEST 6
    a = trj%natoms("TEST")
    ans_val = 1
    call check(a, ans_val, passed, total) 

    ! TEST 7
    a = trj%nframes
    ans_val = 11
    call check(a, ans_val, passed, total) 

    call ndx%read(ndxfile, trj%natoms())

    ! TEST 8
    a = ndx%get("gasGroup",99)
    ans_val = 10099
    call check(a, ans_val, passed, total) 

    ! TEST 9
    a = ndx%get("gasGroup",100)
    ans_val = 10100
    call check(a, ans_val, passed, total) 

    ! TEST 10
    a = ndx%get("gasGroup",trj%natoms("gasGroup"))
    ans_val = 10100
    call check(a, ans_val, passed, total) 

    ! TEST 11
    a = ndx%get("gasGroup",trj%natoms("gasGroup")-1)
    ans_val = 10099
    call check(a, ans_val, passed, total) 

    call finished_tests(passed, total)

end program ndxfile_test
