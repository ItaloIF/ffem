subroutine beemat(bee,deriv)
    ! this subroutine forms the bee matrix in 2-d (ih=3 or 4) or 3-d (ih=6)
    use glb
    implicit none
    real(rkd), intent(in) :: deriv(:,:)
    real(rkd), intent(out) :: bee(:,:)
    integer :: k
    integer :: l
    integer :: m
    integer :: n
    integer :: ih
    integer :: nod
    real(rkd) :: x
    real(rkd) :: y
    real(rkd) :: z
    bee = 0.0_rkd
    ih = ubound(bee,1)
    nod = ubound(deriv,2)
    select case(ih)
        case (3,4)
            do m = 1,nod
                k = 2*m
                l = k-1
                x = deriv(1,m)
                y = deriv(2,m)
                bee(1,l) = x
                bee(3,k) = x
                bee(2,k) = y
                bee(3,l) = y
            end do
        case (6)
            do m = 1,nod
                n = 3*m
                k = n-1
                l = k-1
                x = deriv(1,m)
                y = deriv(2,m)
                z = deriv(3,m)
                bee(1,l) = x
                bee(4,k) = x
                bee(6,n) = x
                bee(2,k) = y
                bee(4,l) = y
                bee(5,n) = y
                bee(3,n) = z
                bee(5,k) = z
                bee(6,l) = z
            end do
        case default
            write(*,*) 'wrong dimension for nst in bee matrix'
        end select
        return
end subroutine beemat