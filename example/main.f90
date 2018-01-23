program main

    use dcdfort_trajectory

    implicit none

    type(Trajectory) :: trj
    real(8) :: box(6)
    real :: x(3)
    integer :: i, j

    call trj%read("dump.dcd")

    ! Can parallelize this loop easily
    do i = 1, trj%nframes

        box = trj%box(i)

        ! Do calculations here
        do j = 1, trj%natoms()

            x = trj%x(i,j)

            ! Just for testing
            write(*,*) x

        end do

    end do

end program main
