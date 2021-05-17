subroutine invert(matrix)
    ! this subroutine inverts a small square matrix onto itself.
    use glb
    implicit none
    real(rkd), intent(in out) :: matrix(:,:)
    real(rkd) :: det
    real(rkd) :: j11
    real(rkd) :: j12
    real(rkd) :: j13
    real(rkd) :: j21
    real(rkd) :: j22
    real(rkd) :: j23
    real(rkd) :: j31
    real(rkd) :: j32
    real(rkd) :: j33
    real(rkd) :: con
    integer :: ndim
    integer :: i
    integer :: k
    ndim = ubound(matrix,1)
    if (ndim.eq.2) then
        det = matrix(1,1)*matrix(2,2) - matrix(1,2)*matrix(2,1)
        j11 = matrix(1,1)
        matrix(1,1) = matrix(2,2)
        matrix(2,2) = j11
        matrix(1,2) = -matrix(1,2)
        matrix(2,1) = -matrix(2,1)
        matrix = matrix/det
    else if (ndim.eq.3) then
        det = matrix(1,1)*(matrix(2,2)*matrix(3,3)-matrix(3,2)*matrix(2,3))
        det = det-matrix(1,2)*(matrix(2,1)*matrix(3,3)-matrix(3,1)*matrix(2,3))
        det = det+matrix(1,3)*(matrix(2,1)*matrix(3,2)-matrix(3,1)*matrix(2,2))
        j11 = matrix(2,2)*matrix(3,3)-matrix(3,2)*matrix(2,3)
        j21 = -matrix(2,1)*matrix(3,3)+matrix(3,1)*matrix(2,3)
        j31 = matrix(2,1)*matrix(3,2)-matrix(3,1)*matrix(2,2)
        j12 = -matrix(1,2)*matrix(3,3)+matrix(3,2)*matrix(1,3)
        j22 = matrix(1,1)*matrix(3,3)-matrix(3,1)*matrix(1,3)
        j32 = -matrix(1,1)*matrix(3,2)+matrix(3,1)*matrix(1,2)
        j13 = matrix(1,2)*matrix(2,3)-matrix(2,2)*matrix(1,3)
        j23 = -matrix(1,1)*matrix(2,3)+matrix(2,1)*matrix(1,3)
        j33 = matrix(1,1)*matrix(2,2)-matrix(2,1)*matrix(1,2)
        matrix(1,1) = j11
        matrix(1,2) = j12
        matrix(1,3) = j13
        matrix(2,1) = j21
        matrix(2,2) = j22
        matrix(2,3) = j23
        matrix(3,1) = j31
        matrix(3,2) = j32
        matrix(3,3) = j33
        matrix = matrix/det
    else
        do k = 1,ndim
            con = matrix(k,k)
            matrix(k,k) = 1.0_rkd
            matrix(k,:) = matrix(k,:)/con
            do i = 1,ndim
                if (i.ne.k) then
                    con = matrix(i,k)
                    matrix(i,k) = 0.0_rkd
                    matrix(i,:) = matrix(i,:) - matrix(k,:)*con
                end if
            end do
        end do
    end if
    return
end subroutine invert