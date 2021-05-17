subroutine asemspr(g,ks,ke,ia,ja,isys)
    use glb
    implicit none
    integer, intent(in) :: g(:)
    real(rkd), intent(in out) :: ks(:)
    real(rkd), intent(in) :: ke(:,:)
    integer, intent(in) :: ia(:)
    integer, intent(in) :: ja(:)
    integer, intent(in) :: isys

    integer :: n
    integer :: neq
    integer :: nnz
    integer :: i
    integer :: j
    integer :: k
    integer :: i_glob
    integer :: j_glob

    n = size(ke,1)
    nnz = size(ks,1)
    neq = size(ia,1) - 1

    do i = 1,n
        i_glob = g(i)
        if (i_glob.ne.0) then
            do j = 1,n
                j_glob = g(j)
                if (j_glob.ne.0) then
                    if (isys.eq.1) then
                        ! symmetric matix
                        if (j_glob.ge.i_glob) then
                            write(*,*) 'isys = 1, no code'
                        end if
                    else if (isys.eq.0) then
                        ! nonsymetric matrix
                        do k = ia(i_glob), ia(i_glob+1)-1
                            if (ja(k).eq.j_glob) then
                                ks(k) = ks(k) + ke(i,j)
                                exit
                            end if
                        end do
                    else
                        write(*,*) 'wrong value of isys in asemspr'
                        stop
                    end if
                end if
            end do
        end if
    end do
    return
end subroutine asemspr