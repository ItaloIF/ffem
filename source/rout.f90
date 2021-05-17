module rout
    interface
    subroutine beemat(bee, deriv)
        use glb
        implicit none
        real(rkd), intent(in) :: deriv(:,:)
        real(rkd), intent(out) :: bee(:,:)
    end subroutine beemat

    subroutine deemat(dee,e,v)
        use glb
        implicit none
        real(rkd), intent(in) :: e
        real(rkd), intent(in) :: v
        real(rkd), intent(out) :: dee(:,:)
    end subroutine deemat

    function determinant(jac) result(det)
        use glb
        implicit none
        real(rkd), intent(in) :: jac(:,:)
        real(rkd) :: det
    end function determinant

    subroutine ecmat(ecm,fun)
        use glb
        implicit none
        real(rkd), intent(in) :: fun(:)
        real(rkd), intent(out) :: ecm(:,:)
    end subroutine ecmat

    subroutine elmat(area,rho,emm)
        use glb
        implicit none
        real(rkd), intent(in) :: area
        real(rkd), intent(in) :: rho
        real(rkd), intent(out) :: emm(:,:)
    end subroutine elmat

    subroutine fkdiag(kdiag, g)
        implicit none
        integer, intent(in) :: g(:)
        integer, intent(out) :: kdiag(:)
    end subroutine fkdiag

    subroutine formnf(nf)
        implicit none
        integer, intent(in out) :: nf(:,:)
    end subroutine formnf

    subroutine getname(argv, nlen)
        implicit none
        integer :: narg
        integer, intent(out) :: nlen
        character(len=:), allocatable, intent(out) ::  argv
        integer :: lnblnk
        integer :: iargc
    end subroutine getname

    subroutine invert(matrix)
        use glb
        implicit none
        real(rkd), intent(in out) :: matrix(:,:)
    end subroutine invert

    subroutine num_to_g(num, nf, g)
        implicit none
        integer, intent(in) :: num(:)
        integer, intent(in) :: nf(:,:)
        integer, intent(out) :: g(:)
    end subroutine num_to_g

    subroutine sample(element,s,wt)
        use glb
        implicit none
        character(len=4), intent(in) :: element
        real(rkd), intent(out) :: s(:,:)
        real(rkd), intent(out), optional :: wt(:)
    end subroutine sample

    subroutine shape_der(der,points,i)
        use glb
        implicit none
        integer, intent(in) :: i
        real(rkd), intent(in) :: points(:,:)
        real(rkd), intent(out) :: der(:,:)
    end subroutine shape_der

    subroutine shape_fun(fun,points,i)
        use glb
        implicit none
        integer, intent(in) :: i
        real(rkd), intent(in) :: points(:,:)
        real(rkd), intent(out) :: fun(:)
    end subroutine shape_fun

    end interface
end module rout
include 'rout/beemat.f90'
include 'rout/deemat.f90'
include 'rout/determinant.f90'
include 'rout/ecmat.f90'
include 'rout/elmat.f90'
include 'rout/fkdiag.f90'
include 'rout/formnf.f90'
include 'rout/getname.f90'
include 'rout/invert.f90'
include 'rout/num_to_g.f90'
include 'rout/sample.f90'
include 'rout/shape_der.f90'
include 'rout/shape_fun.f90'