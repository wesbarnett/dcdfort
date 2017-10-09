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

! TODO: groups/types, lammps file
module dcdfort_trajectory

    use, intrinsic :: iso_fortran_env
    use dcdfort_common
    use dcdfort_index
    use dcdfort_reader
    implicit none
    private

    type :: Frame
        real, allocatable :: xyz(:,:)
        real(8) :: box(6)
    end type

    type, public :: Trajectory
        type(dcdfile) :: dcd
        type(Frame), allocatable :: frameArray(:)
        type(IndexFile) :: ndx
        integer :: NFRAMES
        integer :: NUMATOMS, N, frames_read
        integer :: FRAMES_REMAINING
        integer :: ISTART ! timestep of first frame in trajectory file
        integer :: IEND  ! timestep of last frame in trajectory file
        integer :: NEVERY ! frequency trajectory was saved (in time steps)
        real :: timestep ! simulation time step
        logical :: read_only_index_group
    contains
        procedure :: open => trajectory_open
        procedure :: natoms => trajectory_get_natoms
        procedure :: read_next => trajectory_read_next
        procedure :: x => trajectory_get_xyz
        procedure :: box => trajectory_get_box
        procedure :: read => trajectory_read
        procedure :: close => trajectory_close
        procedure :: skip_next => trajectory_skip_next
    end type

contains

    subroutine trajectory_open(this, filename_in, ndxfile)

        implicit none
        class(Trajectory), intent(inout) :: this
        character (len=*), intent(in) :: filename_in
        character (len=*), intent(in), optional :: ndxfile
        character (len=206) :: filetype
        logical :: ex
        integer :: i, j

        call this%dcd%open(trim(filename_in))
        call this%dcd%read_header(this%nframes, this%istart, this%nevery, this%iend, this%timestep, this%NUMATOMS)

        this%FRAMES_REMAINING = this%NFRAMES

        this%N = this%NUMATOMS ! Save for use when user selects just one group
        if (present(ndxfile)) call this%ndx%read(ndxfile, this%NUMATOMS)

        write(error_unit,'(a)') "Opened "//trim(filename_in)//" for reading."
        write(error_unit,'(i0,a)') this%NUMATOMS, " atoms present in system."
        write(error_unit,'(i0,a)') this%NFRAMES, " frames present in trajectory file."
        write(error_unit,'(a,i0)') "First timestep saved: ", this%ISTART
        write(error_unit,'(a,i0)') "Last timestep saved: ", this%IEND
        write(error_unit,'(a,f12.6)') "Simulation timestep: ", this%timestep
        write(error_unit,'(a,i0,a)') "Trajectory written every ", this%nevery, " timesteps."

    end subroutine trajectory_open

    function trajectory_get_natoms(this, group)

        implicit none
        integer :: trajectory_get_natoms
        class(Trajectory), intent(inout) :: this
        character (len=*), intent(in), optional :: group

        if (this%read_only_index_group .and. present(group)) then
            call error_stop_program("Do not specify an index group in natoms() when already specifying an & 
                &index group with read() or read_next().")
        end if

        trajectory_get_natoms = merge(this%ndx%get_natoms(group), this%NUMATOMS, present(group))

    end function trajectory_get_natoms

    function trajectory_skip_next(this, F)

        class(Trajectory), intent(inout) :: this
        integer, intent(in), optional :: F
        integer :: trajectory_skip_next, i, stat, N

        ! If the user specified how many frames to read and it is greater than one, use it
        N = merge(F, 1, present(F))

        ! Are we near the end of the file?
        N = min(this%FRAMES_REMAINING, N)
        this%FRAMES_REMAINING = this%FRAMES_REMAINING - N

        do i = 1, N
            call this%dcd%skip_next()
            if (stat .ne. 0) then
                write(error_unit,'(a,i0,a)') "Skipped ", i-1, " frames."
                trajectory_skip_next = i-1
                exit
            end if
        end do
        write(error_unit,'(a,i0,a)') "Skipped ", N, " frames."
        trajectory_skip_next = N

    end function trajectory_skip_next

    function trajectory_read_next(this, F, ndxgrp)

        implicit none
        integer :: trajectory_read_next
        class(Trajectory), intent(inout) :: this
        integer, intent(in), optional :: F
        character (len=*), optional :: ndxgrp
        integer :: I, J, N, STAT
        real, allocatable :: xyz(:,:)

        ! If the user specified how many frames to read and it is greater than one, use it
        N = merge(F, 1, present(F))

        ! Are we near the end of the file?
        N = min(this%FRAMES_REMAINING, N)
        this%FRAMES_REMAINING = this%FRAMES_REMAINING - N

        if (allocated(this%frameArray)) deallocate(this%frameArray)
        allocate(this%frameArray(N))

        write(error_unit,*)

        this%read_only_index_group = .false.

        if (present(ndxgrp)) then

            allocate(xyz(3,this%N))
            this%NUMATOMS = this%natoms(trim(ndxgrp))
            do I = 1, N

                if (modulo(I, 1000) .eq. 0) call print_frames_saved(I)

                allocate(this%frameArray(I)%xyz(3,this%NUMATOMS))
                call this%dcd%read_next(xyz, this%frameArray(I)%box)

                do J = 1, size(this%ndx%group)
                    if (trim(this%ndx%group(J)%title) .eq. trim(ndxgrp)) then
                        this%frameArray(I)%xyz = xyz(:,this%ndx%group(J)%LOC)
                        exit
                    end if
                end do

            end do
            deallocate(xyz)

            this%read_only_index_group = .true.

        else

            do I = 1, N

                if (modulo(I, 1000) .eq. 0) call print_frames_saved(I)

                allocate(this%frameArray(I)%xyz(3,this%NUMATOMS))
                call this%dcd%read_next(this%frameArray(I)%xyz, this%frameArray(I)%box)

            end do

        end if

        this%frames_read = N
        trajectory_read_next = N
        call print_frames_saved(N)

    end function trajectory_read_next

    subroutine print_frames_saved(I)

        implicit none
        integer, intent(in) :: I
        write(error_unit,'(a,i0)') achar(27)//"[1A"//achar(27)//"[K"//"Frames saved: ", I

    end subroutine print_frames_saved

    function trajectory_get_xyz(this, frame, atom, group)

        implicit none
        real :: trajectory_get_xyz(3)
        integer, intent(in) :: frame, atom
        integer :: atom_tmp, natoms
        class(Trajectory), intent(inout) :: this
        character (len=1024) :: message
        character (len=*), intent(in), optional :: group

        call trajectory_check_frame(this, frame)

        if (this%read_only_index_group .and. present(group)) then
            call error_stop_program("Do not specify an index group in x() when already specifying an & 
                &index group with read() or read_next().")
        end if

        atom_tmp = merge(this%ndx%get(group, atom), atom, present(group))
        natoms = merge(this%natoms(group), this%natoms(), present(group))

        if (atom > natoms .or. atom < 1) then
            write(message, "(a,i0,a,i0,a)") "Tried to access atom number ", atom_tmp, " when there are ", &
                natoms, ". Note that Fortran uses one-based indexing."
            call error_stop_program(trim(message))
        end if

        trajectory_get_xyz = this%frameArray(frame)%xyz(:,atom_tmp)

    end function trajectory_get_xyz

    subroutine trajectory_read(this, dcdfile, ndxfile, ndxgrp)

        implicit none
        class(Trajectory), intent(inout) :: this
        character (len=*), optional :: ndxfile, ndxgrp
        character (len=*) :: dcdfile
        integer :: N

        call this%open(dcdfile, ndxfile)

        N = this%read_next(this%NFRAMES, ndxgrp)

        call this%close()

    end subroutine trajectory_read

    function trajectory_get_box(this, frame)

        implicit none
        real :: trajectory_get_box(6)
        class(Trajectory), intent(in) :: this
        integer, intent(in) :: frame

        call trajectory_check_frame(this, frame)
        trajectory_get_box = this%frameArray(frame)%box

    end function trajectory_get_box

    subroutine trajectory_check_frame(this, frame)

        implicit none
        class(Trajectory), intent(in) :: this
        integer, intent(in) :: frame
        character (len=1024) :: message

        if (frame > this%frames_read .or. frame < 1) then
            write(message, "(a,i0,a,i0,a)") "Tried to access frame number ", frame, " when there are ", &
                this%frames_read, " read into memory. Note that Fortran uses one-based indexing."
            call error_stop_program(trim(message))
        end if

    end subroutine trajectory_check_frame

    subroutine trajectory_close(this)

        implicit none
        class(Trajectory), intent(inout) :: this

        call this%dcd%close()

    end subroutine trajectory_close

end module dcdfort_trajectory
