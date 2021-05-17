module tfem
    use glb
    implicit none
    type quad4_g
        integer :: nod = 4
        integer :: dim = 2
        integer :: dof = 2
        integer :: nst = 2
        integer :: nip
        integer :: nim
        real(rkd) :: me(8,8)
        real(rkd) :: ce(8,8)
        real(rkd) :: ke(8,8)
        real(rkd) :: dee(2,2)
        real(rkd) :: bee(3,8)
        real(rkd) :: fun(4)
        real(rkd) :: der(2,4)
        real(rkd) :: jac(2,2)
        real(rkd) :: det
        real(rkd) :: derv(2,4)
        real(rkd), allocatable :: gpt(:,:)
        real(rkd), allocatable :: gwe(:)
    end type quad4_g

    type quad4
        integer :: g(8)
        integer :: num(8)
        integer :: typ
    end type quad4

end module tfem