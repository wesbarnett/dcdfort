program main

    use dcdfort_trajectory

    implicit none

    type(Trajectory) :: trj
    integer :: i, j
    real :: box(6)

    call trj%read("dump.dcd")

    do i = 1, 1

        do j = 1, trj%natoms()

            write(*,*) trj%x(i,j)

        end do

        box = trj%box(i)
        write(*,*) box(1), box(2), box(3), box(4), box(5), box(6)

    end do

end program main
