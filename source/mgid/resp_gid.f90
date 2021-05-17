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

    integer :: nels
    integer :: iel
    integer :: kelgs
    integer :: kel
    integer :: npoin
    integer :: ipoin
    integer :: nip
    integer :: ndim
    integer :: i
    integer :: idim
    integer :: nst
    real(rkd) :: zero = 0.0_rkd
    real(rkd) :: j2

    ndim = ubound(points,2)
    nip = ubound(points,1)
    npoin = ubound(nf,2)
    nst = ubound(tte,1)
    nels = ubound(tte,3)

    if (noutp.eq.1) then
        write(unit,*) 'GiD Post Results File 1.0'
        write(unit,*)
        if (element.eq.'tri')  write(unit,3000)
        if (element.eq.'quad') write(unit,3001)
        if (element.eq.'hexa') write(unit,3002)
        write(unit,*)
        write(unit,901) nip
        write(unit,*)
        write(unit,902)
        write(unit,*)
        do i = 1,nip
            write(unit,907) (points(i,idim), idim = 1,ndim)
        end do
        write(unit,*)
        write(unit,903)

        901  format('Number of Gauss Points:',i10)
        902  format('Natural Coordinates: Given')
        903  format('End gausspoints')
        907  format(5x,*(e20.8))
        3000 format('GaussPoints "GaussSoil" Elemtype Triangle')
        3001 format('GaussPoints "GaussSoil" Elemtype Quadrilateral')
        3002 format('GaussPoints "GaussSoil" Elemtype Hexahedra')
    end if

    write(unit,*)
    write(unit,700) istep
    write(unit,*)
    if(ndim.eq.2) write(unit,*) 'ComponentNames "Desl-X", "Desl-Y"'
    if(ndim.eq.3) write(unit,*) 'ComponentNames "Desl-X", "Desl-Y", "Desl-Z"'
    write(unit,*)
    write(unit,*) 'Values'
    do ipoin = 1,npoin
        write(unit,910) ipoin, tdisp(nf(:,ipoin))
    end do
    write(unit,*)
    write(unit,*) 'End Values'

    700  format('Result "Displacements" "Load Analysis" ',i10,' Vector OnNodes')
    910  format(i10,*(e14.6))

    if (nels.ne.0) then
        ! stresses
        write(unit,*)
        write(unit,705) istep
        write(unit,*)
        write(unit,*) 'Values'
        write(unit,*)

        do  iel = 1,nels
            kelgs = 0
            do i = 1,nip
                kelgs = kelgs + 1
                kel = iel
                if (kelgs.eq.1) then
                    if (nst.eq.3) write(unit,706) kel,tte(1,i,iel),tte(2,i,iel),tte(3,i,iel),zero,zero,zero
                    if (nst.eq.4) write(unit,706) kel,tte(1,i,iel),tte(2,i,iel),tte(3,i,iel),tte(4,i,iel),zero,zero
                    if (nst.eq.6) write(unit,706) kel,tte(1,i,iel),tte(2,i,iel),tte(3,i,iel),tte(4,i,iel),tte(5,i,iel), &
                                                  tte(6,i,iel)
                else
                    if (nst.eq.3) write(unit,707) tte(1,i,iel),tte(2,i,iel),tte(3,i,iel),zero,zero,zero
                    if (nst.eq.4) write(unit,707) tte(1,i,iel),tte(2,i,iel),tte(3,i,iel),tte(4,i,iel),zero,zero
                    if (nst.eq.6) write(unit,707) tte(1,i,iel),tte(2,i,iel),tte(3,i,iel),tte(4,i,iel),tte(5,i,iel), &
                                                  tte(6,i,iel)
                end if
            end do
        end do
        write(unit,*)
        write(unit,*) 'End Values'
        write(unit,*)

        ! von misses stresses
        if (present(dd)) then
            write(unit,*)
            write(unit,805) istep
            write(unit,*)
            write(unit,806)
            write(unit,*)
            write(unit,*) 'Values'
            write(unit,*)
            do  iel = 1,nels
                kelgs = 0
                do i = 1,nip
                    kelgs = kelgs + 1
                    kel = iel
                    if (kelgs.eq.1) then
                        write(unit,706) kel, dd(1,i,iel)
                    else
                        write(unit,707) dd(1,i,iel)
                    end if
                end do
            end do
            write(unit,*)
            write(unit,*) 'End Values'
            write(unit,*)
        end if

        705  format('Result "IntStresses" "Load Analysis" ',i10,' Matrix OnGaussPoints "GaussSoil"')
        805  format('Result "VonMises" "Load Analysis" ',i10,' Scalar OnGaussPoints "GaussSoil"')
        806  format('  ComponentNames "Svm "')
        706  format(i10,1x,*(e16.8))
        707  format(11x,*(e16.8))
    end if
end subroutine resp_gid