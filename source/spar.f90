module spar
    interface
    subroutine amux(a,ia,ja,x,y)
        use glb
        implicit none
        real(rkd), intent(in) :: x(:)
        real(rkd), intent(out) :: y(:)
        real(rkd), intent(in) :: a(:)
        integer, intent(in out) :: ia(:)
        integer, intent(in out) :: ja(:)
    end subroutine amux

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
    end subroutine asemcoo1

    subroutine asemcoo2(g,ks,ke,nnz,ia,ja,isys)
        use glb
        implicit none
        integer, intent(in) :: g(:)
        real(rkd), intent(in out) :: ks(:)
        real(rkd), intent(in) :: ke(:,:)
        integer, intent(in out) :: nnz
        integer, intent(in out) :: ia(:)
        integer, intent(in out) :: ja(:)
        integer, intent(in) :: isys
    end subroutine asemcoo2

    subroutine asemspr(g,ks,ke,ia,ja,isys)
        use glb
        implicit none
        integer, intent(in) :: g(:)
        real(rkd), intent(in out) :: ks(:)
        real(rkd), intent(in) :: ke(:,:)
        integer, intent(in) :: ia(:)
        integer, intent(in) :: ja(:)
        integer, intent(in) :: isys
    end subroutine asemspr

    subroutine clncsr(job,value,n,a,ia,ja,indu)
        use glb
        implicit none
        integer, intent(in) :: job
        integer, intent(in) :: value
        integer, intent(in) :: n
        real(rkd), intent(in out) :: a(:)
        integer, intent(in out) :: ia(:)
        integer, intent(in out) :: ja(:)
        integer, intent(out) :: indu(n)
    end subroutine clncsr

    subroutine clncsr2(job,value,n,a,ia,ja,indu)
        use glb
        implicit none
        integer, intent(in) :: job
        integer, intent(in) :: value
        integer, intent(in) :: n
        real(rkd), intent(in out) :: a(:)
        integer, intent(in out) :: ia(:)
        integer, intent(in out) :: ja(:)
        integer, intent(out) :: indu(n)
    end subroutine clncsr2

    subroutine coicsr(job,n,a,ia,ja)
        use glb
        implicit none
        integer, intent(in) :: job
        integer, intent(in) :: n
        real(rkd), intent(in out) :: a(:)
        integer, intent(in out) :: ia(:)
        integer, intent(in out) :: ja(:)
    end subroutine coicsr

    end interface
end module spar
include 'spar/amux.f90'
include 'spar/asemcoo1.f90'
include 'spar/asemcoo2.f90'
include 'spar/asemspr.f90'
include 'spar/clncsr.f90'
include 'spar/clncsr2.f90'
include 'spar/coicsr.f90'