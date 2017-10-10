program main

    use dcdfort_trajectory

    implicit none

    type(Trajectory) :: trj
    integer :: i, j
    real(8) :: box(6)

    ! Opens file, reads all atoms in, closes file
    call trj%read("dump.dcd")

    ! Cycles through frames read in
    ! Easy to parallelize this loop
    do i = 1, trj%nframes

        ! Box for this frame
        box = trj%box(i)
        write(*,*) box(1), box(2), box(3), box(4), box(5), box(6)

        ! Cycles through all atoms
        do j = 1, trj%natoms()

            ! Do some calculations here
            write(*,*) trj%x(i,j)

        end do

    end do

end program main
