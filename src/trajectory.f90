! TODO: box, groups/types, lammps file
module dcdfort_trajectory

    use, intrinsic :: iso_c_binding
    use, intrinsic :: iso_fortran_env
    implicit none
    private

    type :: Frame
        real(C_FLOAT), allocatable :: xyz(:,:)
        real(C_FLOAT) :: box(6)
    end type

    type, public :: Trajectory
        type(C_PTR), pointer :: v
        type(Frame), allocatable :: frameArray(:)
        integer :: NFRAMES
        integer :: NUMATOMS, N
        integer :: FRAMES_REMAINING
    contains
        procedure :: open => trajectory_open
        procedure :: natoms => trajectory_get_natoms
        procedure :: read_next => trajectory_read_next
        procedure :: x => trajectory_get_xyz
        procedure :: box => trajectory_get_box
        procedure :: read => trajectory_read
        procedure :: close => trajectory_close
    end type

    interface 

        type(C_PTR) function open_dcd_read(filename, filetype, natoms) bind(C, name='open_dcd_read')
            import
            character(kind=C_CHAR), intent(in) :: filename(*), filetype(*)
            integer(C_INT), intent(out) :: natoms
        end function

        integer(C_INT) function get_nframes(v) bind(C, name='get_nframes')
            import
            type(C_PTR), intent(in) :: v
        end function

        integer(C_INT) function read_next_wrapper(v, natoms, coords, box) bind(C, name='read_next_wrapper')
            import
            type(C_PTR) :: v
            integer(C_INT) :: natoms
            real(C_FLOAT) :: coords(*), box(*)
        end function

        subroutine close_file_read(v) bind(C, name='close_file_read')
            import
            type(C_PTR) :: v
        end subroutine

    end interface

contains

    subroutine error_stop_program(message)

        implicit none
        character (len=*), intent(in) :: message

        write(error_unit,*)
        write(error_unit,'(a, a)') "LIBDCDFORT ERROR: ", message
        write(error_unit,*)
        call abort()

    end subroutine error_stop_program

    subroutine trajectory_open(this, filename_in)

        use, intrinsic :: iso_c_binding, only: C_NULL_CHAR, c_f_pointer

        implicit none
        class(Trajectory), intent(inout) :: this
        type(C_PTR) :: v_c
        character (len=*), intent(in) :: filename_in
        character (len=206) :: filename, filetype
        logical :: ex

        inquire(file=trim(filename_in), exist=ex)

        if (ex .eqv. .false.) then
            call error_stop_program(trim(filename_in)//" does not exist.")
        end if

        ! Set the file name to be read in for C.
        filename = trim(filename_in)//C_NULL_CHAR
        filetype = "dcd"//C_NULL_CHAR

        v_c = open_dcd_read(filename, filetype, this%NUMATOMS)
        call c_f_pointer(v_c, this % v)
        this%NFRAMES = get_nframes(this % v)
        this%FRAMES_REMAINING = this%NFRAMES

        write(error_unit,'(a)') "Opened "//trim(filename)//" for reading."
        write(error_unit,'(i0,a)') this%NUMATOMS, " atoms present in system."
        write(error_unit,'(i0,a)') this%NFRAMES, " frames present in trajectory file."

    end subroutine trajectory_open

    ! TODO: Get number of atoms in a group
    function trajectory_get_natoms(this)

        implicit none
        integer :: trajectory_get_natoms
        class(Trajectory), intent(in) :: this

        trajectory_get_natoms = this%NUMATOMS

    end function trajectory_get_natoms

    function trajectory_read_next(this, F)

        implicit none
        integer :: trajectory_read_next
        class(Trajectory), intent(inout) :: this
        integer, intent(in), optional :: F
        integer :: I, N, STAT

        ! If the user specified how many frames to read and it is greater than one, use it
        N = merge(F, 1, present(F))

        ! Are we near the end of the file?
        N = min(this%FRAMES_REMAINING, N)
        this%FRAMES_REMAINING = this%FRAMES_REMAINING - N

        if (allocated(this%frameArray)) deallocate(this%frameArray)

        allocate(this%frameArray(N))

        write(error_unit,*)

        do I = 1, N

            if (modulo(I, 1000) .eq. 0) call print_frames_saved(I)

            allocate(this%frameArray(I)%xyz(3,this%NUMATOMS))
            STAT = read_next_wrapper(this%v, this%NUMATOMS, this%frameArray(I)%xyz, this%frameArray(I)%box)

        end do

    end function trajectory_read_next

    subroutine print_frames_saved(I)

        implicit none
        integer, intent(in) :: I
        write(error_unit,'(a,i0)') achar(27)//"[1A"//achar(27)//"[K"//"Frames saved: ", I

    end subroutine print_frames_saved

    ! TODO: groups
    function trajectory_get_xyz(this, frame, atom)

        implicit none
        real :: trajectory_get_xyz(3)
        integer, intent(in) :: frame, atom
        integer :: atom_tmp, natoms
        class(Trajectory), intent(inout) :: this
        character (len=1024) :: message

        call trajectory_check_frame(this, frame)

        atom_tmp = atom
        natoms = this%natoms()

        if (atom > natoms .or. atom < 1) then
            write(message, "(a,i0,a,i0,a)") "Tried to access atom number ", atom_tmp, " when there are ", &
                natoms, ". Note that Fortran uses one-based indexing."
            call error_stop_program(trim(message))
        end if

        trajectory_get_xyz = this%frameArray(frame)%xyz(:,atom)

    end function trajectory_get_xyz

    ! TODO: implement lammps file and groups/types
    subroutine trajectory_read(this, dcdfile)

        implicit none
        class(Trajectory), intent(inout) :: this
        character (len=*) :: dcdfile
        integer :: N

        call this%open(dcdfile)

        N = this%read_next(this%NFRAMES)

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

        if (frame > this%NFRAMES .or. frame < 1) then
            write(message, "(a,i0,a,i0,a)") "Tried to access frame number ", frame, " when there are ", &
                this%NFRAMES, ". Note that Fortran uses one-based indexing."
            call error_stop_program(trim(message))
        end if

    end subroutine trajectory_check_frame

    subroutine trajectory_close(this)

        implicit none
        class(Trajectory), intent(inout) :: this

        call close_file_read(this % v) 

    end subroutine trajectory_close

end module dcdfort_trajectory
