subroutine asemcoo1(g,ks,ke,nnz,ia,ja,isys)
    use glb
    implicit none
    integer, intent(in) :: g(:)
    real(rkd), intent(in out) :: ks(:)
    real(rkd), intent(in) :: ke(:,:)
    integer, intent(in out) :: nnz
    integer, intent(in out) :: ia(:)
    integer, intent(in out) :: ja(:)
    integer, intent(in) :: isys

    integer :: n
    integer :: i
    integer :: j
    integer :: i_glob
    integer :: j_glob

    n = size(ke,1)
    do i = 1,n
        i_glob = g(i)
        if (i_glob.ne.0) then
            do j = 1,n
                j_glob = g(j)
                if (ke(i,j).ne.0.0_rkd) then
                    if (j_glob.ne.0) then
                        if (isys.eq.1) then
                            ! symmetric matix
                            if (j_glob.ge.i_glob) then
                                nnz = nnz + 1
                                ks(nnz) = ke(i,j)
                                ia(nnz) = i_glob
                                ja(nnz) = j_glob
                            end if
                        else if (isys.eq.0) then
                            ! nonsymetric matrix
                            nnz = nnz + 1
                            ks(nnz) = ke(i,j)
                            ia(nnz) = i_glob
                            ja(nnz) = j_glob
                        else
                            write(*,*) 'wrong value of isys in asemcoo1'
                            stop
                        end if
                    end if
                end if
            end do
        end if
    end do
    return
end subroutine asemcoo1