function determinant(jac) result(det)
    ! this function returns the determinant of a 1x1, 2x2 or 3x3 Jacobian matrix.
    use glb
    implicit none
    real(rkd), intent(in) :: jac(:,:)
    real(rkd) ::  det
    integer :: it
    it = ubound(jac,1)
    select case(it)
        case(1)
            det = 1.0_rkd
        case(2)
            det = jac(1,1)*jac(2,2) - jac(1,2)*jac(2,1)
        case(3)
            det = jac(1,1)*(jac(2,2)*jac(3,3)-jac(3,2)*jac(2,3))
            det = det - jac(1,2)*(jac(2,1)*jac(3,3)-jac(3,1)*jac(2,3))
            det = det + jac(1,3)*(jac(2,1)*jac(3,2)-jac(3,1)*jac(2,2))
        case default
            write(*,*) ' wrong dimension for Jacobian matrix'
    end select
    return
end function determinant