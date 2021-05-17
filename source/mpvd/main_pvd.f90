subroutine main_pvd(unit, model, tv)
    use glb
    implicit none
    integer, intent(in) :: unit
    character(*), intent(in) :: model
    real(rkd), intent(in) :: tv(:)

    integer :: nt
    integer :: i

    nt = size(tv)
    write(unit,101)
    write(unit,102)
    write(unit,103)
    do i = 1,nt
        write(unit,104) tv(i), model, tv(i)
    end do
    write(unit,105)
    write(unit,106)

    101  format('<?xml version="1.0"?>')
    102  format('<VTKFile type="Collection" compressor="vtkZLibDataCompressor">')
    103  format('  <Collection>')
    104  format('    <DataSet timestep="',es16.10,'" group="" part="0" file="',a,'_T',es16.10,'_P0.vtu"/>')
    105  format('  </Collection>')
    106  format('</VTKFile>')
end subroutine main_pvd