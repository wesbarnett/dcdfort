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

module dcdfort_reader

    implicit none

    real(8), parameter :: pi = 2.0d0*dacos(0.0d0)

    type, public :: dcdfile
        integer :: u
    contains
        procedure :: open => dcdfile_open
        procedure :: read_header => dcdfile_read_header
        procedure :: close => dcdfile_close
        procedure :: read_next => dcdfile_read_next
        procedure :: skip_next => dcdfile_skip_next
    end type dcdfile

contains

    subroutine dcdfile_open(this, filename)

        implicit none
        character (len=*) :: filename
        class(dcdfile), intent(inout) :: this

        open(newunit=this%u, file=trim(filename), form="unformatted", access="stream")

    end subroutine dcdfile_open

    subroutine dcdfile_read_header(this, nframes, istart, nevery, iend, timestep, natoms)

        implicit none

        integer :: dummy, nframes, istart, nevery, iend, natoms, i, ntitle
        character (len=4) :: cord_string
        character (len=80) :: title_string
        real :: timestep
        class(dcdfile), intent(inout) :: this

        ! Should be 84
        read(this%u) dummy

        ! Should be 'CORD'
        read(this%u) cord_string

        ! Number of snapshots in file
        read(this%u) nframes

        ! Timestep of first snapshot
        read(this%u) istart

        ! Save snapshots every this many steps
        read(this%u) nevery

        ! Timestep of last snapshot
        read(this%u) iend

        do i = 1, 5
            read(this%u) dummy
        end do

        ! Simulation timestep
        read(this%u) timestep

        do i = 1, 12
            read(this%u) dummy
        end do

        read(this%u) ntitle
        do i = 1, ntitle
            read(this%u) title_string
            write(*,'(a)') trim(title_string)
        end do

        read(this%u) dummy
        read(this%u) dummy

        ! Number of atoms in each snapshot
        read(this%u) natoms

        read(this%u) dummy

    end subroutine dcdfile_read_header

    subroutine dcdfile_close(this)

        implicit none
        class(dcdfile), intent(inout) :: this

        close(this%u)

    end subroutine dcdfile_close

    subroutine dcdfile_read_next(this, xyz, box)

        implicit none
        real, allocatable, intent(inout) :: xyz(:,:)
        real(8), intent(inout) :: box(6)
        integer :: dummy
        class(dcdfile), intent(inout) :: this
    
        ! Should be 48
        read(this%u) dummy

        read(this%u) box(1) ! A
        read(this%u) box(6) ! gamma
        read(this%u) box(2) ! B
        read(this%u) box(5) ! beta
        read(this%u) box(4) ! alpha
        read(this%u) box(3) ! C
        if (box(4) >= -1.0 .and. box(4) <= 1.0 .and. box(5) >= -1.0 .and. box(5) <= 1.0 .and. &
            box(6) >= -1.0 .and. box(6) <= 1.0) then
            box(4) = 90.0 - asin(box(4)) * 90.0 / pi
            box(5) = 90.0 - asin(box(5)) * 90.0 / pi
            box(6) = 90.0 - asin(box(6)) * 90.0 / pi
        end if

        ! 48 again
        read(this%u) dummy
        read(this%u) dummy

        read(this%u) xyz(1,:)

        ! 48 again
        read(this%u) dummy
        read(this%u) dummy

        read(this%u) xyz(2,:)

        ! 48 again
        read(this%u) dummy
        read(this%u) dummy

        read(this%u) xyz(3,:)

        read(this%u) dummy

    end subroutine dcdfile_read_next

    subroutine dcdfile_skip_next(this)

        implicit none
        real :: real_dummy
        integer :: dummy, i
        real(8) :: box_dummy(6)
        class(dcdfile), intent(inout) :: this
    
        ! Should be 48
        read(this%u) dummy

        read(this%u) box_dummy

        ! 48 again
        read(this%u) dummy
        read(this%u) dummy

        do i = 1, dummy/4
            read(this%u) real_dummy
        end do

        ! 48 again
        read(this%u) dummy
        read(this%u) dummy

        do i = 1, dummy/4
            read(this%u) real_dummy
        end do

        ! 48 again
        read(this%u) dummy
        read(this%u) dummy

        do i = 1, dummy/4
            read(this%u) real_dummy
        end do

        read(this%u) dummy

    end subroutine dcdfile_skip_next

end module dcdfort_reader
