subroutine fkdiag(kdiag, g)
    implicit none
    integer, intent(in) :: g(:)
    integer, intent(out) :: kdiag(:)
    integer :: idof
    integer :: i
    integer :: iwp1
    integer :: j
    integer :: im
    integer :: k
    
    idof = size(g)
    do i = 1,idof
        iwp1 = 1
        if (g(i) /= 0) then
            do j = 1,idof
                if (g(j) /= 0) then
                    im = g(i) - g(j) + 1
                    if (im > iwp1) iwp1 = im
                end if
            end do
            k = g(i)
            if (iwp1 > kdiag(k)) kdiag(k) = iwp1
        end if
    end do
    return
    end subroutine fkdiag