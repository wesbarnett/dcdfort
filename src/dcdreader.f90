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
!> @brief Module that contains dcdreader class

module dcdfort_reader

    use dcdfort_common
    use iso_c_binding, only: C_NULL_CHAR

    implicit none

    !> @brief dcdwriter class
    type, public :: dcdfile
        integer, private :: u
        integer(8) :: filesize, framesize
    contains
        !> Opens file to read from
        procedure :: open => dcdfile_open
        !> Reads header of open DCD file
        procedure :: read_header => dcdfile_read_header
        !> Closes DCD file
        procedure :: close => dcdfile_close
        !> Reads next frame into memory
        procedure :: read_next => dcdfile_read_next
        !> Skips reading this frame into memory
        procedure :: skip_next => dcdfile_skip_next
    end type dcdfile

contains

    !> @brief Opens file to read from
    !> @param[inout] this dcdreader object
    !> @param[in] filename name of of DCD file to read from
    subroutine dcdfile_open(this, filename)

        implicit none
        character (len=*), parameter :: magic_string = "CORD"
        integer, parameter :: magic_number = 84
        character (len=*), intent(in) :: filename
        class(dcdfile), intent(inout) :: this
        integer :: line1, charmm_version, has_extra_block, four_dimensions
        character (len=4) :: line2
        logical :: ex

        ! Does file exist?
        inquire(file=trim(filename), exist=ex, size=this%filesize)
        if (ex .eqv. .false.) then
            call error_stop_program(trim(filename)//" does not exist.")
        end if

        ! Open file in native endinness
        open(newunit=this%u, file=trim(filename), form="unformatted", access="stream")

        ! Read in magic number of magic string
        read(this%u,pos=1) line1
        read(this%u) line2

        ! Ensure the magic number and string are correct, if not we'll swap the endinness
        if (line1 .ne. magic_number .or. line2 .ne. magic_string) then

            ! Try converting to the reverse endianness
            close(this%u)
            open(newunit=this%u, file=trim(filename), form="unformatted", access="stream", convert="swap")

            read(this%u,pos=1) line2
            read(this%u) line2

            ! We tried both native and reverse endiness and didn't have magic number or string
            if (line1 .ne. magic_number .or. line2 .ne. magic_string) then
                call error_stop_program("This DCD file format is not supported, or the file header is corrupt.")
            end if

        end if

        ! Check if the file identifies as CHARMM (LAMMPS pretends to be CHARMM v. 24)
        read(this%u, pos=85) charmm_version
        if (charmm_version .eq. 0) then
            call error_stop_program("DCD file indicates it is not CHARMM. Only CHARMM-style DCD files are supported.")
        end if

        ! We only support files with the extra unitcell block
        read(this%u, pos=49) has_extra_block
        if (has_extra_block .ne. 1) then
            call error_stop_program("DCD file indicates it does not have unit cell information. Only DCD files with&
                & unit cell information are supported.")
        end if

        ! We don't support files with four dimensions
        read(this%u) four_dimensions
        if (four_dimensions .eq. 1) then
            call error_stop_program("DCD file indicates it has four dimensions. Only DCD files with three dimensions&
                & are supported.")
        end if

    end subroutine dcdfile_open

    !> @brief Reads header of open DCD file
    !> @details Should be called after open()
    !> @param[inout] this dcdreader object
    !> @param[in] nframes rnumber of frames (snapshots) in file
    !> @param[out] istart first timestep of trajectory file
    !> @param[out] nevery how often (in timesteps) file was written to
    !> @param[out] iend last timestep of trajectory file
    !> @param[out] timestep timestep of simulation
    !> @param[out] natoms number of atoms in each snapshot
    subroutine dcdfile_read_header(this, nframes, istart, nevery, iend, timestep, natoms)

        implicit none

        integer, intent(out) :: nframes, istart, nevery, iend, natoms
        integer :: i, ntitle, n, framesize, nframes2, dummy
        character (len=80) :: title_string
        real, intent(out) :: timestep
        class(dcdfile), intent(inout) :: this

        read(this%u, pos=9) nframes, istart, nevery, iend

        read(this%u, pos=45) timestep

        read(this%u, pos=97) ntitle
        if (ntitle > 0) then
            write(error_unit,'(a)') prompt//"The following titles were found:"
        end if
        do i = 1, ntitle
            read(this%u) title_string

            n = 1
            do while (n .le. 80 .and. title_string(n:n) .ne. C_NULL_CHAR)
                n = n + 1
            end do

            write(error_unit,'(a)') prompt//"  "//trim(title_string(1:n))
        end do

        read(this%u) dummy, dummy
        if (dummy .ne. 4) then
            call error_stop_program("This DCD file format is not supported, or the file header is corrupt.")
        end if

        ! Number of atoms in each snapshot
        read(this%u) natoms, dummy
        if (dummy .ne. 4) then
            call error_stop_program("This DCD file format is not supported, or the file header is corrupt.")
        end if

        ! Each frame has natoms*3 (4 bytes each) = natoms*12
        ! plus 6 box dimensions (8 bytes each) = 48
        ! Additionally there are 32 bytes of file information in each frame
        this%framesize = natoms*12 + 80
        ! Header is 276 bytes
        nframes2 = (this%filesize-276)/this%framesize
        if ( nframes2 .ne. nframes) then
            write(error_unit,'(a,i0,a,i0,a)') prompt//"WARNING: Header indicates ", nframes, &
                &" frames, but file size indicates ", nframes2, "." 
            nframes = nframes2
        end if

    end subroutine dcdfile_read_header

    !> @brief Closes DCD file
    !> @param[inout] this dcdreader object
    subroutine dcdfile_close(this)

        implicit none
        class(dcdfile), intent(inout) :: this

        close(this%u)

    end subroutine dcdfile_close

    !> @brief Reads next frame into memory
    !> @param[inout] this dcdreader object
    !> @param[inout] xyz coordinates of all atoms in snapshot
    !> @param[inout] box dimensions of snapshot
    subroutine dcdfile_read_next(this, xyz, box)

        implicit none
        real, allocatable, intent(inout) :: xyz(:,:)
        real(8), intent(inout) :: box(6)
        integer :: dummy(2)
        class(dcdfile), intent(inout) :: this
    
        ! Should be 48
        read(this%u) dummy(1)

        !            A       gamma   B       beta    alpha   C
        read(this%u) box(1), box(6), box(2), box(5), box(4), box(3)

        ! 48, then no. of bytes for x coordinates, x coordinates (repeat for y and z coordinates)
        read(this%u) dummy, xyz(:,1), dummy, xyz(:,2), dummy, xyz(:,3)

        read(this%u) dummy(1)

    end subroutine dcdfile_read_next

    !> @brief Skips reading this frame into memory
    !> @param[inout] this dcdreader object
    !> @param[in] n number of frames to skip
    subroutine dcdfile_skip_next(this, n)

        implicit none
        real :: real_dummy
        integer :: dummy, i
        integer(8) :: pos, newpos
        integer, intent(in), optional :: n
        real(8) :: box_dummy(6)
        class(dcdfile), intent(inout) :: this
   
        inquire(unit=this%u, pos=pos)

        ! We subtract 4 bytes so that the next read of the 4-byte integer will line things up properly for the next read
        if (.not. present(n)) then
            newpos = pos + this%framesize - 4
        else
            newpos = pos + this%framesize*n - 4
        end if
        
        read(this%u, pos=newpos) dummy

    end subroutine dcdfile_skip_next

end module dcdfort_reader
