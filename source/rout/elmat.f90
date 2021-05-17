subroutine elmat(area,rho,emm)
    use glb
    implicit none
    real(rkd), intent(in) :: area
    real(rkd), intent(in) :: rho
    real(rkd), intent(out) :: emm(:,:)
    real(rkd) :: zero = 0.0_rkd
    real(rkd) :: pt2 = 0.2_rkd
    real(rkd) :: pt25 = 0.25_rkd
    integer :: i
    integer :: ndof

    ndof = ubound(emm,1)
    emm = zero
    select case(ndof)
        case (8)
            do i = 1,8
                emm(i,i) = pt25*area*rho
            end do
        case (16)
            do i = 1,16
                emm(i,i) = pt2*area*rho
            end do
            do i = 1,13,4
                emm(i,i) = pt25*emm(3,3)
            end do
            do i = 2,14,4
                emm(i,i) = pt25*emm(3,3)
            end do
        case default
            write(*,*) 'Wrong number of nodes for quad element'
    end select
    return
end subroutine elmat