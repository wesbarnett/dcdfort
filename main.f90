program main

    use dcdfort_trajectory

    implicit none

    type(Trajectory) :: trj
    integer :: n

    call trj%open("dump.dcd")

    write(*,*) "atoms: ", trj%natoms()
    write(*,*) "frames: ", trj%nframes

    n = trj%read_next(5)

end program main
