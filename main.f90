program main

    use dcdfort_trajectory

    implicit none

    type(Trajectory) :: trj

    call trj%open("dump.dcd")

    write(*,*) "atoms: ", trj%NUMATOMS
    write(*,*) "frames: ", trj%NFRAMES

end program main
