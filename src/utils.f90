module dcdfort_utils

    implicit none
    public
    real(8), private, parameter :: pi = 2.0d0*acos(0.0d0)
    real(8), private, parameter :: degreesToRadians = pi/180.0d0
 
contains

    function pbc(a, box)

        implicit none
        real(8), intent(in) :: a(3), box(6)
        real(8) :: pbc(3), tbox(3,3) = 0.0
        integer :: I

        ! A = box(1), B = box(2), C = box(3)
        ! alpha = box(4), beta = box(5), gamma = box(6)

        ! convert angles to box vectors
        tbox(1,1) = box(1)

        tbox(1,2) = box(2)*cos(box(6)*degreesToRadians)
        tbox(2,2) = sqrt(box(2)**2-tbox(1,2)**2)

        tbox(1,3) = box(3)*cos(box(5)*degreesToRadians)
        tbox(2,3) = (box(2)*box(3)*cos(box(4)*degreesToRadians) - tbox(1,2)*tbox(1,3))/tbox(2,2)
        tbox(3,3) = sqrt(box(2)**2-tbox(1,3)**2-tbox(2,3)**2)

        pbc = a
        do I = 3, 1, -1
            pbc(1:I) = pbc(1:I) - tbox(1:I,I) * nint(pbc(I) / tbox(I,I))
        end do

    end function pbc


    function cross(a, b)

        implicit none
        real(8) :: cross(3)
        real(8), intent(in), dimension(3) :: a, b

        cross(1) = a(2) * b(3) - a(3) * b(2)
        cross(2) = a(3) * b(1) - a(1) * b(3)
        cross(3) = a(1) * b(2) - a(2) * b(1)

    end function cross


    function distance2(a, b, box)

        implicit none
        real(8) :: distance2
        real(8), intent(in), dimension(3) :: a, b
        real(8) :: c(3)
        real(8), intent(in) :: box(6)

        c = pbc(a - b, box)
        distance2 = dot_product(c, c)

    end function distance2


    function distance(a, b, box)

        implicit none
        real(8) :: distance
        real(8), intent(in), dimension(3) :: a, b
        real(8), intent(in) :: box(6)
        distance = dsqrt(distance2(a, b, box))

    end function distance

    function bond_vector(a, b, box)

        implicit none
        real(8) :: bond_vector(3)
        real(8), intent(in), dimension(3) :: a, b
        real(8), intent(in) :: box(6)

        bond_vector = pbc(a-b, box)

    end function bond_vector


    function magnitude(a)

        implicit none
        real(8) :: magnitude
        real(8), intent(in) :: a(3)

        magnitude = dsqrt(dot_product(a, a))

    end function magnitude


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
