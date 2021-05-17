subroutine amux(a,ia,ja,x,y)
    use glb
    implicit none
        real(rkd), intent(in) :: a(:)
    integer, intent(in out) :: ia(:)
    integer, intent(in out) :: ja(:)
    real(rkd), intent(in) :: x(:)
    real(rkd), intent(out) :: y(:)
    integer :: n
    integer :: i
    integer :: j
    real(rkd) :: t
    n = size(x)
    y = 0.0_rkd
    do i = 1,n
        do j = ia(i),ia(i+1)-1
            y(i) = y(i) + a(j)*x(ja(j))
        end do
    end do
    return
end subroutine amux