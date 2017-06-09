program main

    use dcdfort_trajectory

    implicit none

    type(Trajectory) :: trj
    integer :: n, i, j, f, k

    call trj%open("dump.dcd")

    n = trj%natoms()
    f = trj%nframes

    do i = 1, 1

        k = trj%read_next()

        do j = 1, n

            write(*,*) i, j, trj%x(1,j)

        end do

    end do

end program main
