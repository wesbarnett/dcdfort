# libdcdfort

Fortran library for natively reading in DCD trajectory files generated
from [LAMMPS](http://lammps.sandia.gov) simulations for analysis. Uses
an object-oriented style.

This is similar to my other project
[libgmxfort](https://github.com/wesbarnett/libgmxfort).

**Note:** DCD files generated from simulation packages other than
LAMMPS may not work correctly with this library.

## Compilation

After cloning the repository, or extracting the release tarball, cd into the
repository. Then:

```bash
mkdir build
cd build
cmake .. -DCMAKE_INSTALL_PREFIX=/usr/local
make
```

## Testing

To test your build, do:

```bash
make test
```

If any tests do not pass, please file an issue.

## Installation

The following will install the library to the location specified by
the cmake flag `-DCMAKE_INSTALL_PREFIX`, which is `/usr/local` by
default.

```bash
make install
```

## Usage

Compile your program with `-ldcdfort`. You may also need to use `-I` to point to
where the modules files are even with all of the right environment variables set
(by default at `/usr/local/include`). 

### Linking other cmake projects

A file is included to easily link other cmake projects to the dcdfort
installation. Use `find_package ( dcdfort )` and the variables
`dcdfort_INCLUDE_DIRS` and `dcdfort_LIBRARIES` in your `CMakeLists.txt`.

### pkg-config

A pkg-config file is included, so that it can
be used in your program compilations. You may need to set the
`PKG_CONFIG_PATH` environment variable to find the file (by default in
the directory `/usr/local/lib/pkgconfig`). See `man 1 pkg-config` for
more information.

### API

Add `use dcdfort_trajectory` to your Fortran program in order to use
the `Trajectory` class and `use dcdfort_utils` in order to use any of
the other utilities.  There is an example in the `example` folder on
how to do this.

#### Reading in trajectory and index files

Typically you will open a trajectory file (and optionally a
corresponding GROMACS-style index file). Then you will read in the
entire trajectory file at once, or you can read it in in chunks. Then
you should close the trajectory file when done.

The simplest way to use this library is to construct a `Trajectory` object and
then use the `read()` method:

```fortran
use dcdfort_trajectory
implicit none
type(Trajectory) :: trj
call trj%read("traj.dcd")
```

If you have a corresponding index file, you can add a second argument to open:

```fortran
call trj%read("traj.dcd", "index.ndx")
```

Now information regarding the index groups is stored in memory and can be used
in some of the following methods.

The `read()` method opens the dcd file, reads in all information, and then
closes it. The `trj` object in this example now stores all of the coordinates and
information from the .dcd file.

If you want to read in the trajectory file in frame-by-frame use `read_next()`
instead of `read()`. To use this, you must additionally open and close the dcd
file on your own. By default it reads in one frame:

```fortran
integer :: n
call trj%open("traj.dcd", "index.ndx")
n = trj%read_next()
call trj%close()
```

To read in more than one, specify an argument. The following reads in 10 frames:

```fortran
n = trj%read_next(10)
```

`read_next()` returns the number of frames actually read in. It is a function,
and not a subroutine. This is useful for using it with a `do while` loop. For
example:

```fortran
use dcdfort_trajectory

implicit none

type(Trajectory) :: trj
integer :: i, n

call trj%open("traj.dcd", "index.ndx")

n = trj%read_next(10)
do while (n > 0)
    do i = 1, n
        ! do some things with the frames read in
    end do
    n = trj%read_next(10)
end do

call trj%close()
```

To skip a frame without reading it into memory use `skip_next()`. You can also
pass an integer argument to indicate how many frames to skip. The function
returns the actual number of frames skipped (you might be near the end of the
file and not able to skip all you specified).

#### Getting simulation information

After calling `read()` or `read_next()` every atom's coordinates are accessible
via the `x()` method. For example, to get the coordinates of the first atom in
the first frame you would do the following. The frame is the first argument and
the atom number is the second argument. 

```fortran
real :: myatom(3)
! ...
myatom = trj%x(1, 1)
```

**Note**: Fortran uses one-based indexing, and that convention is retained here.

If you read in an index file, you can get atom coordinates in relationship to
that. The following gets the fifth atom in index group C in the 10th frame:

```fortran
myatom = trj%x(10, 5, "C")
```

If the index group does not exist, then an error will be thrown, causing the
program to stop.

**Note:** If you have more than one group in your index file with the same name,
this will simply use the first group with that name. It's best not to repeat
group names in your index file. The library will give you a warning if it finds
that an index name is duplicated, but the program will continue.

Note that when you use `x()` you will still have to give it the frame number as
the first argument even if you only read in one frame with `read_next()`.  You
can always get the *total* number of frames in a trajectory file object with the
`nframes` member:

If you want direct access to the object storing a coordinate, do the
following use `trj%frameArray(i)%xyz(j,k)` where `i` is the frame
number, `j` are the x, y, and z coordinates (so `1`, `2`, and `3`),
and `k` is the atom number. The `x()` method is just a convenient way
to get this object.

```fortran
integer :: n
! ...
n = trj%nframes
```

This is distinct from the number of frames read in using `read_next()`. The
frame number passed to the `x()` method, and other methods here, is always in
relationship to the number of frames read in, not the total number of frames in
the file.

To get the timestep corresponding with the first saved frame in the
trajectory file do:

```fortran
integer :: istart
! ...
istart = trj%istart
```

To get the timestep corresponding with the last saved frame in the
trajectory file do:

```fortran
integer :: iend
! ...
iend = trj%iend
```

To get how often frames were saved in your simulation to this
trajectory file use the `nevery` object. This corresponds with the
fifth column in a LAMMPS `dump dcd` line where you indicated to dump
every this many timesteps. It is the column labeled `N` in the LAMMPS
[dump manual page](http://lammps.sandia.gov/doc/dump.html).

```fortran
real(8) :: nevery
! ...
nsavc = trj%nevery
```

To get the simulation timestep, use the `timestep` object. This
corresponds to the `timestep` setting in LAMMPS.

```fortran
real(8) :: timestep
! ...
delta = trj%timestep
```

**Warning:** Some programs such as *catdcd* overwrite time step
information. *dcdfort* outputs this information whenever it opens a
file. If you intend on using this information in your analysis
program, double check that it is correct. If you are only using LAMMPS
output, you shouldn't have to worry about this.

You can also get the number of atoms with the `natoms()` method:

```fortran
integer :: n
! ...
n = trj%natoms()
```

If you want to know how many atoms are in an index group include the group name
as an argument. In this example the group name is "C":

```fortran
n = trj%natoms("C")
```

If that index group does not exist, then the method will simply return 0.

To get the box coordinates, use `box`. The following gets the box of the `2`nd
frame:

```fortran
real(8) :: mybox(6)
! ...
mybox = trj%box(2)
```

The first three elements of the array are the x, y, and z dimensions. The last
three elements are the alpha, beta, and gamma angles.

#### Reading in specific groups only

As shown above, the most common use of this library is to use `read()` or
`read_next()` to save all atom locations and then use getters like `x()` and
`natoms()` to get information about them by specifying an index group as an
argument.

To save memory, you can save just a specific index group with read():

```fortran
trj%read(xtcfile, ndxfile, "C")
```

If you do this, you only have access to the group above, and you should never
pass an index group name to getters like x(), since only one group is available.
If you do specify a group in a getter after already specifying it in `read()` or
`read_next()`, you will get an error, and the program will stop.

#### Utilities

There are several functions and subroutines in the `dcdfort_utils` module,
including periodic boundary and distance calculations. Check out the source file
for what is available.

## License

This project is release under the following license.

```
libdcdfort

Copyright (C) 2017 James W. Barnett

This program is free software; you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free Software
Foundation; either version 2 of the License, or (at your option) any later
version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with
this program; if not, write to the Free Software Foundation, Inc., 51 Franklin
Street, Fifth Floor, Boston, MA 02110-1301 USA.
```

See the file `LICENSE` for full details.
