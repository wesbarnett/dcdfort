program speed

    use dcdfort_reader

    implicit none
    type(dcdfile) :: dcd
    integer :: nframes, istart, nevery, iend, natoms, i
    character (len=1024) :: filename
    real :: timestep, start, finish
    real, allocatable :: xyz(:,:)
    real(8) :: box(6)

    call get_command_argument(1,filename)

    call cpu_time(start)

    call dcd%open(trim(filename))
    call dcd%read_header(nframes, istart, nevery, iend, timestep, natoms)

    allocate(xyz(natoms,3))

    do i = 1, nframes
        call dcd%read_next(xyz, box)
    end do

    call dcd%close()

    call cpu_time(finish)

    print *, real(dcd%filesize)/(finish-start)/(1024*1024), " MiB/sec"

end program speed
