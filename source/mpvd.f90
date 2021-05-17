module mpvd
    interface

    subroutine main_pvd(unit, model, tv)
        use glb
        implicit none
        integer, intent(in) :: unit
        character(*), intent(in) :: model
        real(rkd), intent(in) :: tv(:)
    end subroutine main_pvd

    subroutine step_pvd(unit, coord, cone, disp, nf)
        use glb
        implicit none
        integer, intent(in) :: unit
        real(rkd), intent(in) :: coord(:,:)
        integer, intent(in) ::  cone(:,:)
        real(rkd), intent(in) :: disp(0:)
        integer, intent(in) :: nf(:,:)
    end subroutine step_pvd

    end interface

end module mpvd
include 'mpvd/main_pvd.f90'
include 'mpvd/step_pvd.f90'