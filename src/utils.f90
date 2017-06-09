module dcdfort_utils

    implicit none
    public
 
contains

    ! TODO: FIX!! Only works for CUBIC boxes
    function pbc(a, box)

        implicit none
        real(8), intent(in) :: a(3), box(6)
        real(8) :: pbc(3)
        integer :: I
        integer :: shift

        pbc = a
        pbc(3) = pbc(3) - box(3) * nint(pbc(3)/box(3))
        pbc(2) = pbc(2) - box(2) * nint(pbc(2)/box(2))
        pbc(1) = pbc(1) - box(1) * nint(pbc(1)/box(1))

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
        real(8), intent(in), dimension(3,3) :: box
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
