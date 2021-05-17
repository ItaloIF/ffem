subroutine formnf(nf)
    implicit none
    integer, intent(in out) :: nf(:,:)
    integer :: i
    integer :: j
    integer :: m
    
    m = 0
    do  j = 1,UBOUND(nf,2)
        do i = 1,UBOUND(nf,1)
            if (nf(i,j) /= 0) then
                m = m+1
                nf(i,j) = m
            end if
        end do
    end do
    return
end subroutine formnf