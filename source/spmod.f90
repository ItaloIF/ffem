module spmod

    implicit none
    contains

    subroutine sp_solver(idsol,n,nnz,a,ia,ja,b)
        use glb
        implicit none
        integer, intent(in) :: idsol
        integer, intent(in) :: n
        integer, intent(in) :: nnz
        real(rkd), intent(in out) :: a(:)
        integer, intent(in) :: ia(:)
        integer, intent(in) :: ja(:)
        real(rkd), intent(in out) :: b(:)
        
        select case (idsol)
            case (1)
                ! SuperLU solver
                call sp_lu(n,nnz,a,ia,ja,b)
            case (2)
                ! Umfpack Solver
                write(*,*) 'Umfpack no coded'
            case default
                write(*,*) 'Wrong solver id value'
            end select

    end subroutine sp_solver

    subroutine sp_lu(n,nnz,a,ia,ja,b)
        use glb
        implicit none
        integer, intent(in) :: n
        integer, intent(in) :: nnz
        real(rkd), intent(in out) :: a(:)
        integer, intent(in) :: ia(:)
        integer, intent(in) :: ja(:)
        real(rkd), intent(in out) :: b(:)
        
        integer:: nrhs
        integer :: ldb
        integer :: info
        integer :: iopt
        integer :: factors

        nrhs = 1
        ldb = n
        iopt = 1
        info = 1
        call c_fortran_dgssv( iopt, n, nnz, nrhs, a(1:nnz), ja(1:nnz), ia(1:n+1), b, ldb, factors, info )
        iopt = 2
        call c_fortran_dgssv( iopt, n, nnz, nrhs, a(1:nnz), ja(1:nnz), ia(1:n+1), b, ldb, factors, info )
    end subroutine sp_lu

end module spmod