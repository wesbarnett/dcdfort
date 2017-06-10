# libdcdfort

Personal code for reading in DCD trajectory files for analysis. Uses and
object-oriented style, reading in the DCD format using the VMD plugin (provided
with this source) and DCD.

This is similar to my other project
[libgmxfort](https://github.com/wesbarnett/libgmxfort), and there is the remote
possibility that I could merge these two projects together, but really I just
needed something quick and simple to work with a new project I have using the
DCD file format.

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

If the tests do not pass, please file an issue.

## Installation

The following will install the library to the location specified by
`DCMAKE_INSTALL_PREFIX`.

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
`dcdfort_INCLUDE_DIRS` and `dcdfort_LIBRARIES`.

### pkg-config

A pkg-config file is included, so that it can
be used in your program compilations. You may need to set `PKG_CONFIG_PATH` to
find the file (by default in the directory `/usr/local/lib/pkgconfig`)

### API

To use the library always put `use dcdfort_trajectory` for in order to use the
`Trajectory` class and `use dcdfort_utils` in order to use any of the other
utilities.  There is an example in the `example` folder on how to do this.

Typically you will open a trajectory file (and optionally a corresponding index
file). Then you will read in the entire trajectory file at once, or you can read
it in in chunks. Then you should close the trajectory file when done.

The simplest way to use this library is to construct a `Trajectory` object and
then use the `read()` method:

```fortran
use dcdfort_trajectory
implicit none
type(Trajectory) :: trj
call trj%read("traj.dcd")
```

The `read()` method opens the dcd file, reads in all information, and then
closes it. The `trj` object in this example now stores all of the coordinates and
information from the .dcd file.

If you want to read in the trajectory file in frame-by-frame use `read_next()`
instead of `read()`. To use this, you must additionally open and close the dcd
file on your own. By default it reads in one frame:

```fortran
integer :: n
call trj%open("traj.dcd")
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

call trj%open("traj.dcd")

n = trj%read_next(10)
do while (n > 0)
    do i = 1, n
        ! do some things with the frames read in
    end do
    n = trj%read_next(10)
end do

call trj%close()
```

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

Note that when you use `x()` you will still have to give it the frame number as
the first argument even if you only read in one frame with `read_next()`.  You
can always get the number of frames in a trajectory file object with the
`nframes` member:

```fortran
integer :: n
! ...
n = trj%nframes
```

You can also get the number of atoms with the `natoms()` method:

```fortran
integer :: n
! ...
n = trj%natoms()
```

To get the box coordinates, use `box`. The following gets the box of the `2`nd
frame:

```fortran
real :: mybox(6)
! ...
mybox = trj%box(2)
```

The first three elements of the array are the x, y, and z dimensions. The last
three elements are the alpha, beta, and gamma angles.

As shown above, the most common use of this library is to use `read()` or
`read_next()` to save all atom locations and then use getters like `x()` and
`natoms()` to get information about them by specifying an index group as an
argument.

There are several functions and subroutines in the `dcdfort_utils` module,
including periodic boundary and distance calculations. Check out the source file
for what is available.

## License

Other than the DCD plugin which is released under its own license, this project
is release under the following license.

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

See LICENSE for full details.
