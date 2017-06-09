program main

    use dcdfort_trajectory

    implicit none

    type(Trajectory) :: trj
    integer :: n, i, j, f, k

    call trj%read("dump.dcd")

    do i = 1, trj%nframes

        do j = 1, trj%natoms()

            write(*,*) i, j, trj%x(i,j)

        end do

    end do

end program main
