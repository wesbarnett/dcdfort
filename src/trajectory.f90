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
!> @brief Module that contains Trajectory type

module dcdfort_trajectory

    use dcdfort_common
    use dcdfort_index, only: IndexFile
    use dcdfort_reader, only: dcdfile
    use dcdfort_utils, only: vol
    implicit none
    private

    type :: Frame
        real(kind=real32), allocatable :: xyz(:,:)
        real(kind=real64) :: box(6)
    end type

    !> @brief Trajectory class
    type, public :: Trajectory
        !> number of trajectory frames (snapshots) in Trajectory object
        integer(kind=int32) :: NFRAMES
        !> timestep of first frame in trajectory file
        integer(kind=int32) :: ISTART         
        !> timestep of last frame in trajectory file
        integer(kind=int32) :: IEND  
        !> frequency trajectory was saved (in time steps)
        integer(kind=int32) :: NEVERY 
        !> simulation time step
        real(kind=real32) :: timestep
        type(IndexFile), private :: ndx
        type(dcdfile), private :: dcd
        type(Frame), allocatable, private :: frameArray(:)
        integer(kind=int32), private :: NUMATOMS
        integer(kind=int32), private :: N
        integer(kind=int32), private :: frames_read
        integer(kind=int32), private :: FRAMES_REMAINING
        logical(kind=int32), private :: read_only_index_group
    contains
        !> Trajectory class method which opens DCD file and optionally index file.
        procedure :: open => trajectory_open
        !> Trajectory class method which reads a specified number of frames into memory after using the open() method
        procedure :: read_next => trajectory_read_next
        !> Trajectory class method which skips a specified number of frames
        procedure :: skip_next => trajectory_skip_next
        !> Trajectory class method which closes a DCD file which was opened with open()
        procedure :: close => trajectory_close
        !> Trajectory class method which opens, reads, and closes a trajectory file
        procedure :: read => trajectory_read
        !> Trajectory class method which returns the coordinates of a particle
        procedure :: x => trajectory_get_xyz
        !> Trajectory class method which returns the box from a simulation frame
        procedure :: box => trajectory_get_box
        !> Trajectory class method which gets the number of atoms in the system or an index group
        procedure :: natoms => trajectory_get_natoms
        procedure :: vol => trajectory_vol
    end type

contains

    !> @brief Trajectory class method which opens DCD file and optionally index file.
    !
    !> @param[inout] this the Trajectory object
    !> @param[in] filename name of DCD file
    !> @param[in] ndxfile name of GROMACS style index file
    subroutine trajectory_open(this, filename, ndxfile)

        implicit none
        class(Trajectory), intent(inout) :: this
        character(len=*), intent(in) :: filename
        character(len=*), intent(in), optional :: ndxfile

        call this%dcd%open(trim(filename))

        write(error_unit,'(a)') prompt//"Opened "//trim(filename)//" for reading."
        call this%dcd%read_header(this%nframes, this%istart, this%nevery, this%iend, this%timestep, this%NUMATOMS)

        write(error_unit,'(a,i0,a)') prompt, this%NUMATOMS, " atoms present in system."
        write(error_unit,'(a,i0,a)') prompt, this%NFRAMES, " frames present in trajectory file."
        write(error_unit,'(a,i0)') prompt//"First timestep in trajectory file: ", this%ISTART
        write(error_unit,'(a,i0)') prompt//"Last timestep in strajectory file: ", this%IEND
        write(error_unit,'(a,f12.6)') prompt//"Simulation timestep: ", this%timestep
        write(error_unit,'(a,i0,a)') prompt//"Trajectory written every ", this%nevery, " timesteps."

        this%FRAMES_REMAINING = this%NFRAMES

        this%N = this%NUMATOMS ! Save for use when user selects just one group
        if (present(ndxfile)) call this%ndx%read(ndxfile, this%NUMATOMS)

    end subroutine trajectory_open

    !> @brief Trajectory class method which gets the number of atoms in the system or an index group
    !
    !> @param[inout] this the Trajectory object
    !> @param[in] group name of index group
    !> @return number of atoms in index group (if specified) or in system
    pure function trajectory_get_natoms(this, group)

        implicit none
        integer(kind=int32) :: trajectory_get_natoms
        class(Trajectory), intent(in) :: this
        character(len=*), intent(in), optional :: group

        if (this%read_only_index_group .and. present(group)) then
            call error_stop_program("Do not specify an index group in natoms() when already specifying an & 
                &index group with read() or read_next().")
        end if

        trajectory_get_natoms = merge(this%ndx%get_natoms(group), this%NUMATOMS, present(group))

    end function trajectory_get_natoms

    !> @brief Trajectory class method which skips a specified number of frames
    !
    !> @param[inout] this Trajectory class
    !> @param[in] F number of frames to skip; if not specified, 1 frame is skipped
    !> @return number of frames skipped
    function trajectory_skip_next(this, F)

        implicit none
        class(Trajectory), intent(inout) :: this
        integer(kind=int32), intent(in), optional :: F
        integer(kind=int32) :: trajectory_skip_next, N

        ! If the user specified how many frames to read and it is greater than one, use it
        N = merge(F, 1, present(F))

        ! Are we near the end of the file?
        N = min(this%FRAMES_REMAINING, N)
        this%FRAMES_REMAINING = this%FRAMES_REMAINING - N

        call this%dcd%skip_next(N)
        write(error_unit,'(a,i0,a)') prompt//"Skipped ", N, " frames."
        trajectory_skip_next = N

    end function trajectory_skip_next

    !> @brief Trajectory class method which reads a specified number of frames into memory after using the open() method
    !
    !> @param[inout] this Trajectory class
    !> @param[in] F number of frames to read in; if not specified, 1 frame is read
    !> @param[in] ndxgrp read only this index group into memory
    !> @param[in] every Read in every this many frames; default is to read in every frame
    !> @return number of frames read in
    function trajectory_read_next(this, F, ndxgrp, every)

        implicit none
        integer(kind=int32) :: trajectory_read_next
        class(Trajectory), intent(inout) :: this
        integer(kind=int32), intent(in), optional :: F, every
        character(len=*), optional :: ndxgrp
        integer(kind=int32) :: I, J, N, S, K
        real(kind=real32), allocatable :: xyz(:,:)

        ! If the user specified how many frames to read and it is greater than one, use it
        N = merge(F, 1, present(F))
        S = merge(every, 1, present(every)) - 1

        ! Are we near the end of the file?
        N = min(this%FRAMES_REMAINING, N)
        this%FRAMES_REMAINING = this%FRAMES_REMAINING - N

        if (present(every)) N = N / every

        if (allocated(this%frameArray)) deallocate(this%frameArray)
        allocate(this%frameArray(N))

        write(error_unit,*)

        this%read_only_index_group = .false.

        if (present(ndxgrp)) then

            allocate(xyz(this%N,3))
            this%NUMATOMS = this%natoms(trim(ndxgrp))
            do I = 1, N

                if (modulo(I*(S+1), 1000) .eq. 0) call print_frames_saved(I, ndxgrp)

                allocate(this%frameArray(I)%xyz(this%NUMATOMS,3))
                call this%dcd%skip_next(S)
                call this%dcd%read_next(xyz, this%frameArray(I)%box)

                do J = 1, size(this%ndx%group)
                    if (trim(this%ndx%group(J)%title) .eq. trim(ndxgrp)) then
                        ! Tends to be faster to list this out individually
                        do K = 1, this%NUMATOMS
                            this%frameArray(I)%xyz(K,1) = xyz(this%ndx%group(J)%LOC(K),1)
                            this%frameArray(I)%xyz(K,2) = xyz(this%ndx%group(J)%LOC(K),2)
                            this%frameArray(I)%xyz(K,3) = xyz(this%ndx%group(J)%LOC(K),3)
                        end do
                        exit
                    end if
                end do

            end do
            deallocate(xyz)

            this%read_only_index_group = .true.

        else

            do I = 1, N

                if (modulo(I*(S+1), 1000) .eq. 0) call print_frames_saved(I)

                allocate(this%frameArray(I)%xyz(this%NUMATOMS,3))
                call this%dcd%skip_next(S)
                call this%dcd%read_next(this%frameArray(I)%xyz, this%frameArray(I)%box)

            end do

        end if

        this%frames_read = N
        trajectory_read_next = N
        call print_frames_saved(N, ndxgrp)

    end function trajectory_read_next

    subroutine print_frames_saved(I, ndxgrp)

        implicit none
        integer(kind=int32), intent(in) :: I
        character(len=*), intent(in), optional :: ndxgrp

        if (present(ndxgrp)) then
            write(error_unit,'(a,i0)') achar(27)//"[1A"//achar(27)//"[K"//prompt//"Frames saved &
                &from index group '"//trim(ndxgrp)//"': ", I
        else
            write(error_unit,'(a,i0)') achar(27)//"[1A"//achar(27)//"[K"//prompt//"Frames saved: ", I
        end if

    end subroutine print_frames_saved

    !> @brief Trajectory class method which returns the coordinates of a particle
    !> @details Gets the coordinates of a particle from the read in trajectory. Returns a real array with 3 elements (x, y, and z).
    !! An optional index group can be specified; if so, the atomid is in relationship to the group.
    !> @param[inout] this Trajectory class
    !> @param[in] frame
    !> @param[in] atom atomid of particle to get
    !> @param[in] group optional index group
    !> @return coordinate of the particle
    pure function trajectory_get_xyz(this, frame, atom, group)

        implicit none
        real(kind=real32) :: trajectory_get_xyz(3)
        integer(kind=int32), intent(in) :: frame, atom
        integer(kind=int32) :: atom_tmp, natoms
        class(Trajectory), intent(in) :: this
        character(len=1024) :: message
        character(len=*), intent(in), optional :: group

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

        ! Faster to list individually
        trajectory_get_xyz(1) = this%frameArray(frame)%xyz(atom_tmp,1)
        trajectory_get_xyz(2) = this%frameArray(frame)%xyz(atom_tmp,2)
        trajectory_get_xyz(3) = this%frameArray(frame)%xyz(atom_tmp,3)

    end function trajectory_get_xyz

    !> @brief Trajectory class method which opens, reads, and closes a trajectory file
    !> @param[inout] this Trajectory class
    !> @param[in] dcdfile Name of DCD trajectory file
    !> @param[in] ndxfile Name of optional index file
    !> @param[in] ndxgrp Name of optional group. If specified, only that group will be read into memory; otherwise, all particles
    !> read into memory.
    !> @param[in] every Read in every this many frames; default is 1
    !> @param[in] skip Skip this many frames at the beginning of the trajectory file; default is 0
    subroutine trajectory_read(this, dcdfile, ndxfile, ndxgrp, every, skip)

        implicit none
        class(Trajectory), intent(inout) :: this
        character(len=*), optional :: ndxfile, ndxgrp
        character(len=*) :: dcdfile
        integer(kind=int32) :: N
        integer(kind=int32), intent(in), optional :: every, skip

        call this%open(dcdfile, ndxfile)

        if (present(every)) then
            write(error_unit,'(a,i0,a)') prompt//"Saving ", every, " snapshots into memory."
        end if
        if (present(skip)) then
            write(error_unit,'(a,i0,a)') prompt//"Skipping first ", skip, " snapshots."
        end if

        if (present(skip)) then
            N = this%skip_next(skip)
            N = this%read_next(this%NFRAMES-N, ndxgrp, every)
        else
            N = this%read_next(this%NFRAMES, ndxgrp, every)
        end if

        call this%close()

    end subroutine trajectory_read

    !> @brief Trajectory class method which returns the box from a simulation frame
    !> @param[inout] this Trajectory class
    !> @param[in] frame Which frame to return the box dimensions
    !> @return  box dimensions of specified snapshot. Array containing 6 elements, ordered as A, B, C, alpha, beta, gamma.
    !! A = length of unit cell vector along x-axis;
    !! B = length of unit cell vector in xy-plane;
    !! C = length of unit cell vector in yz-plane;
    !! alpha = cosine of angle between B and C;
    !! beta = cosine of angle between A and C;
    !! gamma = cosine of angle between A and B;
    pure function trajectory_get_box(this, frame)

        implicit none
        real(kind=real64) :: trajectory_get_box(6)
        class(Trajectory), intent(in) :: this
        integer(kind=int32), intent(in) :: frame
        integer(kind=int32) :: i

        call trajectory_check_frame(this, frame)
        do i = 1, 6
            trajectory_get_box(i) = this%frameArray(frame)%box(i)
        end do

    end function trajectory_get_box

    pure subroutine trajectory_check_frame(this, frame)

        implicit none
        class(Trajectory), intent(in) :: this
        integer(kind=int32), intent(in) :: frame
        character(len=1024) :: message

        if (frame > this%frames_read .or. frame < 1) then
            write(message, "(a,i0,a,i0,a)") "Tried to access frame number ", frame, " when there are ", &
                this%frames_read, " read into memory. Note that Fortran uses one-based indexing."
            call error_stop_program(trim(message))
        end if

    end subroutine trajectory_check_frame

    !> @brief Trajectory class method which closes a DCD file which was opened with open()
    !> @param[inout] this Trajectory class
    subroutine trajectory_close(this)

        implicit none
        class(Trajectory), intent(inout) :: this

        call this%dcd%close()

    end subroutine trajectory_close

    !> @brief Trajectory class method which returns the volume of the box
    !> @param[inout] this Trajectory class
    !> @param[in] frame snapshot to get box volume of
    !> @return the volume of the box of the frame specified
    pure function trajectory_vol(this, frame)

        implicit none
        class(Trajectory), intent(in) :: this
        integer(kind=int32), intent(in) :: frame
        real(kind=real64) :: trajectory_vol

        trajectory_vol = vol(this%box(frame))

    end function trajectory_vol

end module dcdfort_trajectory
