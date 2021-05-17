subroutine elela10(prop,dee)
    use glb
    use rout
    implicit none
    real(rkd), intent(in) :: prop(:)
    real(rkd), intent(in out) :: dee(:,:)

    real(rkd) :: e
    real(rkd) :: v
    e = prop(1)
    v = prop(2)
    call deemat(dee,e,v)
    return
end subroutine elela10