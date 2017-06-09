!
module dcdfort_trajectory

    use, intrinsic :: iso_c_binding, only: C_PTR, C_CHAR, C_FLOAT, C_INT
    use, intrinsic :: iso_fortran_env
    implicit none
    private

    type :: Frame
        real(C_FLOAT), allocatable :: xyz(:,:)
        integer(C_INT) :: STEP
        real(C_FLOAT) :: box(3,3), prec, time
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

    end interface

contains

    subroutine error_stop_program(message)

        implicit none
        character (len=*), intent(in) :: message

        write(error_unit,*)
        write(error_unit,'(a, a)') "LIBGMXFORT ERROR: ", message
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

        write(error_unit,'(a)') "Opened "//trim(filename)//" for reading."
        write(error_unit,'(i0,a)') this%NUMATOMS, " atoms present in system."
        !write(error_unit,'(i0,a)') this%NFRAMES, " frames present in trajectory file."

    end subroutine trajectory_open

    ! TODO: Get number of atoms in a group
    function trajectory_get_natoms(this)

        implicit none
        integer :: trajectory_get_natoms
        class(Trajectory), intent(in) :: this

        trajectory_get_natoms = this%NUMATOMS

    end function trajectory_get_natoms

end module dcdfort_trajectory
