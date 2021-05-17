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

    integer :: npoin
    integer :: ndim
    integer :: idim
    integer :: ipoin
    integer :: ielem
    integer :: nelem
    integer :: nnode

    npoin = ubound(coord,2)
    ndim = ubound(coord,1)
    nelem = ubound(cone,2)
    nnode = ubound(cone,1)

    select case(element)
        case ('line')
            ! agregar
        case ('tri')
            if (nnode.eq.3) write(unit,*) 'MESH'// ctype //'dimension 2 ElemType Triangle Nnode 3'
            if (nnode.eq.6) write(unit,*) 'MESH dimension 2 ElemType Triangle Nnode 6'
            if (iact.eq.0) then
                write(unit,*) 'Coordinates'
                do ipoin = 1,npoin
                    write(unit,400) ipoin, (coord(idim,ipoin), idim = 1,ndim)
                end do
                write(unit,*) 'End Coordinates'
            end if
            write(unit,*) 
            write(unit,*) 'Elements'
            if (nnode.eq.3) then
                do ielem = 1,nelem
                    write(unit,100) ielem,cone(1,ielem),cone(3,ielem),cone(2,ielem),matno(ielem)
                end do
            else if (nnode.eq.6) then
                do ielem = 1,nelem
                    write(unit,100) ielem,cone(3,ielem),cone(1,ielem),cone(5,ielem),cone(2,ielem),cone(6,ielem),cone(4,ielem), &
                                    matno(ielem)
                end do
            end if
        case ('quad')
            if (nnode.eq.4) write(unit,*) 'MESH'// ctype //'dimension 2 ElemType Quadrilateral Nnode 4'
            if (nnode.eq.8) write(unit,*) 'MESH'// ctype //'dimension 2 ElemType Quadrilateral Nnode 8'
            if (nnode.eq.9) write(unit,*) 'MESH'// ctype //'dimension 2 ElemType Quadrilateral Nnode 9'
            if (iact.eq.0) then
                write(unit,*) 'Coordinates'
                do ipoin = 1,npoin
                    write(unit,400) ipoin, (coord(idim,ipoin), idim = 1,ndim)
                end do
                write(unit,*) 'End Coordinates'
            end if
            write(unit,*)
            write(unit,*) 'Elements'
            if (nnode.eq.4) then
                do ielem = 1,nelem
                    write(unit,100) ielem,cone(1,ielem),cone(4,ielem),cone(3,ielem),cone(2,ielem),matno(ielem)
                end do
            elseif (nnode.eq.8) then
                do ielem = 1,nelem
                    write(unit,100) ielem,cone(1,ielem),cone(7,ielem),cone(5,ielem), &
                                    cone(3,ielem),cone(8,ielem),cone(6,ielem), &
                                    cone(4,ielem),cone(2,ielem),matno(ielem)
                end do
            elseif (nnode.eq.9) then
                do ielem = 1,nelem
                    write(unit,100) ielem,cone(1,ielem),cone(7,ielem),cone(5,ielem), &
                                    cone(3,ielem),cone(8,ielem),cone(6,ielem),cone(4,ielem), &
                                    cone(2,ielem),cone(9,ielem),matno(ielem)
                end do
            end if
        case ('hexa')
            if (nnode.eq.8) write(unit,*) 'MESH'// ctype //'dimension 3 ElemType Hexahedra Nnode 8'
            if (nnode.eq.20) write(unit,*) 'MESH'// ctype //'dimension 3 ElemType Hexahedra Nnode 20'
            if (iact.eq.0) then
                write(unit,*) 'Coordinates'
                do ipoin = 1,npoin
                    write(unit,400) ipoin, (coord(idim,ipoin), idim = 1,ndim)
                end do
                write(unit,*) 'End Coordinates'
            end if
            write(unit,*)
            write(unit,*) 'Elements'
            if (nnode.eq.8) then
                do ielem = 1,nelem
                    write(unit,100) ielem,cone(1,ielem),cone(2,ielem),cone(3,ielem), &
                                    cone(4,ielem),cone(5,ielem),cone(6,ielem), &
                                    cone(7,ielem),cone(8,ielem),matno(ielem)
                end do
            else if (nnode.eq.20) then
                ! revisar sentido convention de GID con numeration of Smith
                do ielem = 1,nelem
                    write(unit,100) ielem,&
                                    cone(1,ielem),cone(7,ielem),cone(19,ielem), cone(13,ielem), &
                                    cone(3,ielem),cone(5,ielem),cone(17,ielem), cone(15,ielem), &
                                    cone(8,ielem),cone(12,ielem),cone(20,ielem), cone(9,ielem), &
                                    cone(2,ielem),cone(6,ielem),cone(18,ielem), cone(14,ielem), &
                                    cone(4,ielem),cone(11,ielem),cone(16,ielem), cone(10,ielem), &
                                    matno(ielem)
                end do
            end if
        case default
            write(*,*) 'wrong element in mesh_gid'
    end select
    write(unit,*) 'End Elements'
    write(unit,*)
    iact = iact + 1
    400 format(i10,*(e16.8))
    100 format(i10,*(i7))
    return
end subroutine mesh_gid