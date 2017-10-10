!
! This file is part of libdcdfort
! https://github.com/wesbarnett/dcdfort
!
! Copyright (c) 2017 by James W. Barnett
!
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 2 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License along
! with this program; if not, write to the Free Software Foundation, Inc.,
! 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

module dcdfort_writer

    use dcdfort_common
    use iso_c_binding, only: C_NULL_CHAR

    implicit none

    real(8), parameter :: pi = 2.0d0*dacos(0.0d0)

    type, public :: dcdwriter
        integer :: u
        integer :: nframes_pos, iend_pos, curr_pos, nframes, iend, nevery
    contains
        procedure :: open => dcdwriter_open
        procedure :: write_header => dcdwriter_write_header
        procedure :: close => dcdwriter_close
        procedure :: write_next => dcdwriter_write_next
    end type dcdwriter

contains

    subroutine dcdwriter_open(this, filename)

        implicit none
        character (len=*) :: filename
        class(dcdwriter), intent(inout) :: this
        integer :: filesize
        logical :: ex

        open(newunit=this%u, file=trim(filename), form="unformatted", access="stream", status="new")

    end subroutine dcdwriter_open

    subroutine dcdwriter_write_header(this, istart, nevery, timestep, natoms)

        implicit none

        integer :: dummy, istart, nevery, natoms, i, n
        real :: timestep
        class(dcdwriter), intent(inout) :: this
        character (len=79) :: remarks1, remarks2
        character (len=8) :: date
        character (len=10) :: time

        call date_and_time(date=date,time=time)
        remarks1 = "Created by libdcdfort"
        remarks2 = "REMARK Created on "//date//" "//time

        this%nframes = 0
        this%iend = istart
        this%nevery = nevery

        write(this%u) 84
        write(this%u) "CORD"

        inquire(unit=this%u, pos=this%nframes_pos)
        ! Number of snapshots in file
        write(this%u) this%nframes

        ! Timestep of first snapshot
        write(this%u) istart

        ! Save snapshots every this many steps
        write(this%u) nevery

        inquire(unit=this%u, pos=this%iend_pos)
        ! Timestep of last snapshot
        write(this%u) this%iend

        do i = 1, 4
            write(this%u) 0
        end do

        ! Has unit cell
        write(this%u) 1

        ! Simulation timestep
        write(this%u) timestep

        do i = 1, 9
            write(this%u) 0
        end do

        ! Pretend to be CHARMM version 24
        write(this%u) 24
        write(this%u) 84
        write(this%u) 164

        write(this%u) 2
        write(this%u) remarks1//C_NULL_CHAR
        write(this%u) remarks2//C_NULL_CHAR

        write(this%u) 164
        write(this%u) 4

        ! Number of atoms in each snapshot
        write(this%u) natoms

        write(this%u) 4

    end subroutine dcdwriter_write_header

    subroutine dcdwriter_close(this)

        implicit none
        class(dcdwriter), intent(inout) :: this

        close(this%u)

    end subroutine dcdwriter_close

    subroutine dcdwriter_write_next(this, xyz, box_in)

        implicit none
        real, intent(in) :: xyz(:,:)
        real(8), intent(in) :: box_in(6)
        real(8) :: box(6)
        integer :: dummy
        class(dcdwriter), intent(inout) :: this
        integer :: coord_size

        coord_size = size(xyz,2)*4
        box = box_in
    
        ! Should be 48 (6 double precision floats)
        write(this%u) 48

        box(4) = (90.0 - box(4)) * dsin(pi/180.0)
        box(5) = (90.0 - box(5)) * dsin(pi/180.0)
        box(6) = (90.0 - box(6)) * dsin(pi/180.0)
        write(this%u) box(1) ! A
        write(this%u) box(6) ! gamma
        write(this%u) box(2) ! B
        write(this%u) box(5) ! beta
        write(this%u) box(4) ! alpha
        write(this%u) box(3) ! C

        write(this%u) 48
        write(this%u) coord_size

        write(this%u) xyz(1,:)

        write(this%u) 48
        write(this%u) coord_size

        write(this%u) xyz(2,:)

        write(this%u) 48
        write(this%u) coord_size

        write(this%u) xyz(3,:)

        write(this%u) 48

        inquire(unit=this%u, pos=this%curr_pos)

        this%nframes = this%nframes+1
        this%iend = this%iend + this%nevery

        write(this%u, pos=this%nframes_pos) this%nframes
        write(this%u, pos=this%iend_pos) this%iend
        write(this%u, pos=this%curr_pos)

    end subroutine dcdwriter_write_next

end module dcdfort_writer
