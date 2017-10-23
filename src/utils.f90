!
! This file is part of libdcdfort
! https://github.com/wesbarnett/dcdfort
!
! Copyright (c) 2017 by James W. Barnett
!
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 2 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License along
! with this program; if not, write to the Free Software Foundation, Inc.,
! 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

!> @file
!> @author James W. Barnett, Columbia University
!
!> @brief Module that contains some useful utilities

module dcdfort_utils

    implicit none
    public
    real(8), private, parameter :: pi = 2.0d0*acos(0.0d0)
    real(8), private, parameter :: degreesToRadians = pi/180.0d0
 
contains

    !> @brief Corrects for the periodic boundary condition
    !> @details Moves particle (or vector) the distance of half the box if it is more than half the distance of the box
    !> @param[in] a original coordinates
    !> @param[in] box simulation box
    !> @return the shifted coordinates
    function pbc(a, box)

        implicit none
        real(8), intent(in) :: a(3), box(6)
        real(8) :: pbc(3), tbox(3,3) = 0.0d0
        integer :: I, shift

        ! A = box(1)
        ! B = box(2)
        ! C = box(3)
        ! alpha = box(4)
        ! beta = box(5)
        ! gamma = box(6)

        !  ax  bx  cx
        !   0  by  cy
        !   0   0  cz

        ! convert angles to box vectors
        ! ax = A
        tbox(1,1) = box(1)

        ! bx = B*cos(gamma)
        tbox(1,2) = box(2)*dcos(box(6)*degreesToRadians)

        ! cx = C*cos(beta)
        tbox(1,3) = box(3)*dcos(box(5)*degreesToRadians)

        ! by = B*sin(gamma) = sqrt(B**2 - bx**2)
        tbox(2,2) = dsqrt(box(2)**2-tbox(1,2)**2)

        ! cy = (dot(B,C) - bx*cx) / by = (B*C*cos(alpha) - bx*cx) / by
        tbox(2,3) = (box(2)*box(3)*dcos(box(4)*degreesToRadians) - tbox(1,2)*tbox(1,3))/tbox(2,2)

        ! cz = sqrt(C**2 - cx**2 - cy**2)
        tbox(3,3) = dsqrt(box(3)**2-tbox(1,3)**2-tbox(2,3)**2)

        pbc = a

        shift = nint(pbc(3) / tbox(3,3))
        if (shift .ne. 0) then
            pbc(3) = pbc(3) - tbox(3,3) * shift
            pbc(2) = pbc(2) - tbox(2,3) * shift
            pbc(1) = pbc(1) - tbox(1,3) * shift
        end if

        shift = nint(pbc(2) / tbox(2,2))
        if (shift .ne. 0) then
            pbc(2) = pbc(2) - tbox(2,2) * shift
            pbc(1) = pbc(1) - tbox(1,2) * shift
        end if

        shift = nint(pbc(1) / tbox(1,1))
        if (shift .ne. 0) then
            pbc(1) = pbc(1) - tbox(1,1) * shift
        end if

    end function pbc

    !> @brief Performs cross product between two vectors
    !> @param[in] a first vector
    !> @param[in] b second vector
    !> @result resulting cross product
    function cross(a, b)

        implicit none
        real(8) :: cross(3)
        real(8), intent(in), dimension(3) :: a, b

        cross(1) = a(2) * b(3) - a(3) * b(2)
        cross(2) = a(3) * b(1) - a(1) * b(3)
        cross(3) = a(1) * b(2) - a(2) * b(1)

    end function cross

    !> @brief Calculates the distance squared between two points
    !> @param[in] a first point
    !> @param[in] b second point
    !> @param[in] box box, if pbc to be accounted for
    !> @result distance squared
    function distance2(a, b, box)

        implicit none
        real(8) :: distance2
        real(8), intent(in), dimension(3) :: a, b
        real(8) :: c(3)
        real(8), intent(in), optional :: box(6)

        if (present(box)) then
            c = pbc(a - b, box)
        else
            c = a - b
        end if
        distance2 = dot_product(c, c)

    end function distance2

    !> @brief Calculates the distance between two points
    !> @param[in] a first point
    !> @param[in] b second point
    !> @param[in] box box, if pbc to be accounted for
    !> @result distance 
    function distance(a, b, box)

        implicit none
        real(8) :: distance
        real(8), intent(in), dimension(3) :: a, b
        real(8), intent(in), optional :: box(6)

        if (present(box)) then
            distance = dsqrt(distance2(a, b, box))
        else
            distance = dsqrt(distance2(a, b))
        end if

    end function distance

    !> @brief Calculates the bond vector between two points
    !> @param[in] a first point
    !> @param[in] b second point
    !> @param[in] box box, if pbc to be accounted for
    !> @result bond vector pointing from a to b
    function bond_vector(a, b, box)

        implicit none
        real(8) :: bond_vector(3)
        real(8), intent(in), dimension(3) :: a, b
        real(8), intent(in) :: box(6)

        bond_vector = pbc(a-b, box)

    end function bond_vector

    !> @brief Calculates the magnitude of a vector
    !> @param[in] a vector
    !> @result magnitude of a
    function magnitude(a)

        implicit none
        real(8) :: magnitude
        real(8), intent(in) :: a(3)

        magnitude = dsqrt(dot_product(a, a))

    end function magnitude

    !> @brief Calculates the bond angle between two vectors
    !> @details Calculates the angle between the vector formed by a-b and b-c.
    !> @param[in] a first point
    !> @param[in] b middle point
    !> @param[in] c third point
    !> @param[in] box box, if pbc to be accounted for
    !> @result bond angle between a-b-c
    function bond_angle(a, b, c, box)

        implicit none
        real(8) :: bond_angle
        real(8), intent(in), dimension(3) :: a, b, c
        real(8), intent(in) :: box(6)
        real(8), dimension(3) :: bond1, bond2

        bond1 = bond_vector(b, a, box)
        bond2 = bond_vector(b, c, box)

        bond_angle = dacos(dot_product(bond1, bond2)/(magnitude(bond1)*magnitude(bond2)))

    end function bond_angle

    !> @brief Calculates the dihedral angle between two planes formed by four atoms
    !> @details Calculates the dihedral angle between the vectors formed by i-j, j-k, k-l
    !> @param[in] i first point
    !> @param[in] j middle point
    !> @param[in] k third point
    !> @param[in] l fourth point
    !> @param[in] box box, if pbc to be accounted for
    !> @result dihedral angle forms by i-j-k-l
    function dihedral_angle(i, j, k, l, box)

        implicit none
        real(8) :: dihedral_angle
        real(8), intent(in), dimension(3) :: i, j, k, l
        real(8), intent(in), dimension(6) :: box
        real(8) :: A_mag, B_mag, G_mag
        real(8), dimension(3) :: H, G, F, A, B, cross_BA

        H = bond_vector(k, l, box)
        G = bond_vector(k, j, box)
        F = bond_vector(j, i, box)
        A = cross(F, G)
        B = cross(H, G)
        cross_BA = cross(B, A)
        A_mag = magnitude(A)
        B_mag = magnitude(B)
        G_mag = magnitude(G)

        dihedral_angle = atan2(dot_product(cross_BA, G) / (A_mag*B_mag*G_mag), dot_product(A, B) / (A_mag*B_mag))

    end function dihedral_angle

end module dcdfort_utils
