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

!> @file
!> @author James W. Barnett, Columbia University
!
!> @brief Module that contains dcdwriter class

module dcdfort_writer

    use dcdfort_common
    use iso_c_binding, only: C_NULL_CHAR

    implicit none

    real(8), parameter :: pi = 2.0d0*dacos(0.0d0)

    !> @brief dcdwriter class
    type, public :: dcdwriter
        integer, private :: u
        integer, private :: nframes_pos, iend_pos, curr_pos
        integer, private :: nframes
        integer, private :: iend
        integer, private :: nevery
    contains
        !> Opens new file to write to
        procedure :: open => dcdwriter_open
        !> Writes header to new DCD file
        procedure :: write_header => dcdwriter_write_header
        !> Closes DCD file
        procedure :: close => dcdwriter_close
        !> Writes header to new DCD file
        procedure :: write_next => dcdwriter_write_next
    end type dcdwriter

contains

    !> @brief Opens new file to write to
    !> @param[inout] this dcdwriter object
    !> @param[in] filename name of new DCD file to write to
    subroutine dcdwriter_open(this, filename)

        implicit none
        character (len=*), intent(in) :: filename
        class(dcdwriter), intent(inout) :: this
        integer :: filesize
        logical :: ex

        open(newunit=this%u, file=trim(filename), form="unformatted", access="stream", status="new")

    end subroutine dcdwriter_open

    !> @brief Writes header to new DCD file
    !> @param[inout] this dcdwriter object
    !> @param[in] istart first timestep in file 
    !> @param[in] nevery how often snapshots written (in timesteps)
    !> @param[in] timestep simulation timestep
    !> @param[in] natoms number of atoms in each snapshot
    subroutine dcdwriter_write_header(this, istart, nevery, timestep, natoms)

        implicit none

        integer :: dummy, i, n
        integer, intent(in) :: istart, nevery, natoms
        real, intent(in) :: timestep
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

    !> @brief Closes DCD file
    !> @param[inout] this dcdwriter object
    subroutine dcdwriter_close(this)

        implicit none
        class(dcdwriter), intent(inout) :: this

        close(this%u)

    end subroutine dcdwriter_close

    !> @brief Writes snapshot to an open DCD file
    !> @details Writes a new snapshot to a DCD file. Header should have already been written.
    !> @param[inout] this dcdwriter object
    !> @param[in] xyz coordinates of all atoms in this snapshot
    !> @param[in] box_in box dimensions for this snapshot
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

        ! Go back and update header
        write(this%u, pos=this%nframes_pos) this%nframes
        write(this%u, pos=this%iend_pos) this%iend
        write(this%u, pos=this%curr_pos)

    end subroutine dcdwriter_write_next

end module dcdfort_writer
