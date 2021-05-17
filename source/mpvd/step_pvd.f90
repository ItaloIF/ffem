subroutine step_pvd(unit, coord, cone, disp, nf)
    use glb
    implicit none
    integer, intent(in) :: unit
    real(rkd), intent(in) :: coord(:,:)
    integer, intent(in) ::  cone(:,:)
    real(rkd), intent(in) :: disp(0:)
    integer, intent(in) :: nf(:,:)

    integer :: npoin
    integer :: ndim
    integer :: nelem
    integer :: nnode
    integer :: i
    real(rkd) :: zero = 0.0_rkd

    npoin = ubound(coord,2)
    ndim = ubound(coord,1)
    nelem = ubound(cone,2)
    nnode = ubound(cone,1)

    write(unit,101)
    write(unit,102)
    write(unit,104)
    write(unit,106) npoin, nelem
   
    write(unit,120)
    write(unit,122)
    do i = 1,npoin
        write(unit,124) i
    end do
    write(unit,111)
    write(unit,123)
    do i = 1,npoin
        if (ndim.eq.3) then
            write(unit,110) disp(nf(:,i))
        else
            write(unit,110) disp(nf(:,i)), zero
        end if
    end do
    write(unit,111)
    write(unit,121)
    write(unit,125)
    write(unit,127)
    do i = 1,nelem
        write(unit,119) i
    end do
    write(unit,111)
    write(unit,126)

    write(unit,108)
    write(unit,109)
    do i = 1,npoin
        if (ndim.eq.3) then
            write(unit,110) coord(:,i)
        else
            write(unit,110) coord(:,i), zero
        end if
    end do
    write(unit,111)
    write(unit,112)
    write(unit,113)
    write(unit,115)
    do i = 1,nelem
        write(unit,119) cone(:,i) - 1
    end do
    write(unit,111)
    write(unit,116)
    do i = 1,nelem
        write(unit,119) i*nnode
    end do
    write(unit,111)
    write(unit,117)
    do i = 1,nelem
        write(unit,119) 9
    end do
    write(unit,111)
    write(unit,114)

    write(unit,107)
    write(unit,105)
    write(unit,103)


    101  format('<?xml version="1.0"?>')
    102  format('<VTKFile type="UnstructuredGrid" version="0.1" byte_order="BigEndian">')
    103  format('</VTKFile>')

    104  format('  <UnstructuredGrid>')
    105  format('  </UnstructuredGrid>')

    106  format('    <Piece NumberOfPoints="',i0,'" NumberOfCells="',i0,'">')
    107  format('    </Piece>')

    108  format('      <Points>')
    109  format('        <DataArray type="Float32" NumberOfComponents="3" Format="ascii">')
    110  format('         ',*(f18.8))
    111  format('        </DataArray>')
    112  format('      </Points>')

    113  format('      <Cells>')
    114  format('      </Cells>')
    115  format('        <DataArray type="Int32" Name="connectivity" Format="ascii">')
    116  format('        <DataArray type="Int32" Name="offsets" Format="ascii">')
    117  format('        <DataArray type="Int32" Name="types" Format="ascii">')
    118  format('      </Points>')
    119  format('         ',*(' ',i0))

    120  format('      <PointData Scalars="scalars">')
    121  format('      </PointData>')
    122  format('        <DataArray type="Int32" Name="NodeTag" format="ascii">')
    123  format('        <DataArray type="Float32" Name="Displacement" NumberOfComponents="3" format="ascii">')
    124  format('          ',i0)

    125  format('      <CellData Scalars="scalars">')
    126  format('      </CellData>')
    127  format('        <DataArray type="Int32" Name="EleTag" Format="ascii">')

end subroutine step_pvd