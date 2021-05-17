module tool
    use glb
    implicit none

    contains

    subroutine readmsh(fname,mfix,nodof,ndim,nod,nn,nels,g_coord,nf,g_num)
        use glb
        implicit none
        character(len=60), intent(in):: fname
        integer, intent(in) :: mfix(:,:)
        integer, intent(in) :: nodof
        integer, intent(in) :: ndim
        integer, intent(in) :: nod
        integer, intent(out) :: nn
        integer, intent(out) :: nels
        real(rkd), allocatable, intent(out) :: g_coord(:,:)
        integer, allocatable, intent(out) :: nf(:,:)
        integer, allocatable, intent(out) :: g_num(:,:)

        character(len=500) :: commt
        integer :: nent
        integer :: nnent
        integer :: i
        integer :: j
        integer :: k
        integer :: inode
        integer :: idim
        integer :: ient
        integer :: dumm
        integer :: nelet
        integer :: iels
        integer, allocatable :: g_tot(:,:)

        real(rkd) :: c(2,3)
        real(rkd) :: d1, d2, d3
        integer :: tric
        integer :: nfix
        integer :: tfix(nodof)

        nfix = ubound(mfix,2)
        open(100, file = trim(fname), status='unknown')
        do
            read(100,'(a)') commt
            if (trim(commt).eq.'$Nodes') exit
        end do
        read(100,*) nent, nn
        allocate(g_coord(ndim,nn))
        allocate(nf(nodof,nn))
        nf = 1
        inode = 1
        do i = 1,nent
            tfix = 1
            read(100,*) idim, ient, dumm, nnent
            do k = 1,nfix
                if (mfix(1,k).eq.idim .and. mfix(2,k).eq.ient) then
                    tfix = mfix(3:nodof+2,k)
                    exit
                end if
            end do
            do j = 1,nnent
                read(100,*)
            end do
            do j = 1,nnent
                read(100,*) g_coord(:,inode)
                nf(:,inode) = tfix
                inode = inode + 1
            end do
        end do
        do
            read(100,'(a)') commt
            if (trim(commt).eq.'$Elements') exit
        end do
        read(100,*) nent, nelet
        allocate(g_tot(nod,nelet))
        iels = 0
        do i = 1,nent
            read(100,*) idim, ient, dumm, nnent
            if (idim.eq.ndim) then
                do j = 1,nnent
                    iels = iels + 1
                    read(100,*) dumm, g_tot(:,iels)
                    if (nod.eq.6) then
                        c = g_coord(:,g_tot(1:3,iels))
                        d1 = (c(1,1)-c(1,2))**2 + (c(2,1)-c(2,2))**2
                        d2 = (c(1,2)-c(1,3))**2 + (c(2,2)-c(2,3))**2
                        d3 = (c(1,3)-c(1,1))**2 + (c(2,3)-c(2,1))**2
                        if (max(d1,d2,d3).eq.d1) then
                            tric = 3
                        else if (max(d1,d2,d3).eq.d2) then
                            tric = 1
                        else
                            tric = 2
                        end if
                        call gmsh_to_fem(ndim,nod,g_tot(:,iels),tric)
                    else
                        call gmsh_to_fem(ndim,nod,g_tot(:,iels))
                    end if
                end do
            else
                do j = 1,nnent
                    read(100,*)
                end do
            end if
        end do
        nels = iels
        allocate(g_num(nod,nels))
        g_num = g_tot(:,1:nels)
        close(100)
    end subroutine readmsh

    subroutine gmsh_to_fem(ndim,nod,elepts,tric)
        use glb
        implicit none
        integer, intent(in) :: ndim
        integer, intent(in) :: nod
        integer, intent(in out) :: elepts(nod)
        integer, intent(in), optional :: tric
        integer :: dum(nod)

        select case(ndim)
            case (2)
                select case(nod)
                    case (3)
                        dum(1) = elepts(1)
                        dum(2) = elepts(3)
                        dum(3) = elepts(2)
                        elepts = dum
                    case (6)
                        if (tric.eq.3) then
                            dum(1) = elepts(1)
                            dum(2) = elepts(6)
                            dum(3) = elepts(3)
                            dum(4) = elepts(5)
                            dum(5) = elepts(2)
                            dum(6) = elepts(4)
                        else if (tric.eq.1) then 
                            dum(1) = elepts(2)
                            dum(2) = elepts(4)
                            dum(3) = elepts(1)
                            dum(4) = elepts(6)
                            dum(5) = elepts(3)
                            dum(6) = elepts(5)
                        else
                            dum(1) = elepts(3)
                            dum(2) = elepts(5)
                            dum(3) = elepts(2)
                            dum(4) = elepts(4)
                            dum(5) = elepts(1)
                            dum(6) = elepts(6)
                        end if
                        elepts = dum
                    case (10)
                        dum(1) = elepts(1)
                        dum(2) = elepts(9)
                        dum(3) = elepts(8)
                        dum(4) = elepts(3)
                        dum(5) = elepts(7)
                        dum(6) = elepts(6)
                        dum(7) = elepts(2)
                        dum(8) = elepts(5)
                        dum(9) = elepts(4)
                        dum(10) = elepts(10)
                        elepts = dum
                    case (15)
                        dum(1) = elepts(1)
                        dum(2) = elepts(12)
                        dum(3) = elepts(11)
                        dum(4) = elepts(10)
                        dum(5) = elepts(3)
                        dum(6) = elepts(9)
                        dum(7) = elepts(8)
                        dum(8) = elepts(7)
                        dum(9) = elepts(2)
                        dum(10) = elepts(6)
                        dum(11) = elepts(5)
                        dum(12) = elepts(4)
                        dum(13) = elepts(13)
                        dum(14) = elepts(15)
                        dum(15) = elepts(14)
                        elepts = dum
                    case (4)
                        dum(1) = elepts(1)
                        dum(2) = elepts(4)
                        dum(3) = elepts(3)
                        dum(4) = elepts(2)
                        elepts = dum
                    case (5)
                        write(*,*) 'not code of 5 nodes'
                        stop
                    case (8)
                        dum(1) = elepts(1)
                        dum(2) = elepts(8)
                        dum(3) = elepts(4)
                        dum(4) = elepts(7)
                        dum(5) = elepts(3)
                        dum(6) = elepts(6)
                        dum(7) = elepts(2)
                        dum(8) = elepts(5)
                        elepts = dum
                    case (9)
                        dum(1) = elepts(1)
                        dum(2) = elepts(8)
                        dum(3) = elepts(4)
                        dum(4) = elepts(7)
                        dum(5) = elepts(3)
                        dum(6) = elepts(6)
                        dum(7) = elepts(2)
                        dum(8) = elepts(5)
                        dum(9) = elepts(9)
                        elepts = dum
                    case default
                        write(*,*) 'wrong number of nodes'
                        stop
                end select
            case (3)
                select case(nod)
                    case (4)
                        dum(1) = elepts(1)
                        dum(2) = elepts(2)
                        dum(3) = elepts(4)
                        dum(4) = elepts(3)
                        elepts = dum
                    case (8)
                        dum(1) = elepts(1)
                        dum(2) = elepts(4)
                        dum(3) = elepts(8)
                        dum(4) = elepts(5)
                        dum(5) = elepts(2)
                        dum(6) = elepts(3)
                        dum(7) = elepts(7)
                        dum(8) = elepts(6)
                        elepts = dum
                    case (14)
                        write(*,*) 'not code of 14 nodes'
                        stop
                    case (20)
                        dum(1) = elepts(1)
                        dum(2) = elepts(10)
                        dum(3) = elepts(4)
                        dum(4) = elepts(16)
                        dum(5) = elepts(8)
                        dum(6) = elepts(18)
                        dum(7) = elepts(5)
                        dum(8) = elepts(11)
                        dum(9) = elepts(9)
                        dum(10) = elepts(14)
                        dum(11) = elepts(20)
                        dum(12) = elepts(17)
                        dum(13) = elepts(2)
                        dum(14) = elepts(12)
                        dum(15) = elepts(3)
                        dum(16) = elepts(15)
                        dum(17) = elepts(7)
                        dum(18) = elepts(19)
                        dum(19) = elepts(6)
                        dum(20) = elepts(13)
                        elepts = dum
                    case default
                        write(*,*) 'wrong number of nodes'
                        stop
                end select  
            case default
                write(*,*) 'wrong number of dimensions'
                stop
        end select
    end subroutine


end module tool