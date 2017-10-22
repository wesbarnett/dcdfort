program speed

    use dcdfort_reader

    implicit none
    type(dcdfile) :: dcd
    integer :: nframes, istart, nevery, iend, natoms, i, j, k
    character (len=1024) :: filename
    real :: timestep, start, finish
    real, allocatable :: xyz(:,:), xyz2(:)
    real(8) :: box(6)

    call get_command_argument(1,filename)

    call cpu_time(start)

    call dcd%open(trim(filename))
    call dcd%read_header(nframes, istart, nevery, iend, timestep, natoms)

    allocate(xyz(natoms,3))
    allocate(xyz2(3))

    do i = 1, nframes
        call dcd%read_next(xyz, box)
    end do

    call dcd%close()

    call cpu_time(finish)

    print *, real(dcd%filesize)/(finish-start)/(1024*1024), " MiB/sec"

    call cpu_time(start)

    do i = 1, nframes
        do j = 1, natoms
            xyz2(1) = xyz(j,1)
            xyz2(2) = xyz(j,2)
            xyz2(3) = xyz(j,3)
        end do
    end do

    call cpu_time(finish)

    print *, real(natoms*12*nframes)/(finish-start)/(1024*1024), " MiB/sec"

end program speed
