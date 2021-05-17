module mgid
    interface

    subroutine mesh_gid(coord,cone,matno,element,unit,iact,ctype,nn)
        use glb
        implicit none
        real(rkd), intent(in) :: coord(:,:)
        integer, intent(in) ::  cone(:,:) ! connectivity
        integer, intent(in) ::  matno(:)
        character(*), intent(in) :: element
        integer, intent(in) ::  unit
        integer, intent(in out) ::  iact
        character(*), intent(in) :: ctype
        integer, intent(in), optional ::  nn
    end subroutine mesh_gid

    subroutine resp_gid(tte,tdisp,points,nf,istep,noutp,unit,element,dd)
        use glb
        implicit none
        real(rkd), intent(in) :: tte(:,:,:)
        real(rkd), intent(in) :: tdisp(0:)
        real(rkd), intent(in) :: points(:,:)
        integer, intent(in) :: nf(:,:)
        integer, intent(in) :: istep
        integer, intent(in) :: noutp
        integer, intent(in) :: unit
        character(*), intent(in) :: element
        real(rkd), intent(in), optional :: dd(:,:,:)
    end subroutine resp_gid

    end interface

end module mgid
include 'mgid/mesh_gid.f90'
include 'mgid/resp_gid.f90'