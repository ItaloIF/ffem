subroutine elstr10(prop,sti,esi,stf)
    use glb
    use rout
    implicit none
    real(rkd), intent(in) :: prop(:)
    real(rkd), intent(in) :: sti(:)
    real(rkd), intent(in) :: esi(:)
    real(rkd), intent(out) :: stf(:)

    real(rkd) :: e
    real(rkd) :: v
    real(rkd), allocatable :: dee(:,:)
    integer ::  nst
    nst = size(sti)
    allocate(dee(nst,nst))
    e = prop(1)
    v = prop(2)
    call deemat(dee,e,v)
    stf = sti + matmul(dee,esi)
    return
end subroutine elstr10