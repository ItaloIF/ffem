subroutine num_to_g(num, nf, g)
    implicit none
    integer, intent(in) :: num(:)
    integer, intent(in) :: nf(:,:)
    integer, intent(out) :: g(:)
    integer :: i
    integer :: k
    integer :: nod
    integer :: nodof

    nod = ubound(num,1)
    nodof = ubound(nf,1)
    do i = 1,nod
        k = i*nodof
        g(k-nodof+1:k) = nf(:,num(i))
    end do
    return
end subroutine num_to_g