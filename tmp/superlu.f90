module mod_superlu

    implicit none
    contains

    subroutine superlu_solver(n,nnz,a,ia,ja,b)
    use glbvar
    implicit none
    integer(ikd), intent(in) :: n
    integer(ikd), intent(in) :: nnz
    real(rkd), intent(in out) :: a(:)
    integer(ikd), intent(in) :: ia(:)
    integer(ikd), intent(in) :: ja(:)
    real(rkd), intent(in out) :: b(:)
    
    integer(4) :: nrhs
    integer(4) :: ldb
    integer :: info
    integer :: iopt
    integer(4) :: factors

    nrhs = 1
    ldb = n
    info = 1
    iopt = 1
    call c_fortran_dgssv( iopt, n, nnz, nrhs, a(1:nnz), ja(1:nnz), ia(1:n+1), b, ldb, factors, info )

    iopt = 2
    call c_fortran_dgssv( iopt, n, nnz, nrhs, a(1:nnz), ja(1:nnz), ia(1:n+1), b, ldb, factors, info )

    end subroutine superlu_solver

end module mod_superlu