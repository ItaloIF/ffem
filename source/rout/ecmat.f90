subroutine ecmat(ecm,fun)
    use glb
    implicit none
    real(rkd), intent(in) :: fun(:)
    real(rkd), intent(out) :: ecm(:,:)
    integer :: nodof
    integer :: ndof
    integer :: nod
    integer :: i
    integer :: j
    real(rkd), allocatable :: nt(:,:)
    real(rkd), allocatable :: tn(:,:)
    real(rkd) :: zero = 0.0_rkd

    nod = size(fun)
    ndof = ubound(ecm,1)
    nodof = ndof/nod
    allocate(nt(ndof,nodof))
    allocate(tn(nodof,ndof))
    nt = zero
    tn = zero
    do i = 1,nod
        do j = 1,nodof
            nt((i-1)*nodof+j,j) = fun(i)
            tn(j,(i-1)*nodof+j) = fun(i)
        end do
    end do
    ecm = matmul(nt,tn)
    return
end subroutine ecmat