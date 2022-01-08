program pc1
    use glbvar
    use main
    use geom
    use tool
    use gidrou
    implicit none
    integer(ikd) :: fixed_freedoms
    integer(ikd) :: i
    integer(ikd) :: iel
    integer(ikd) :: k
    integer(ikd) :: loaded_nodes
    integer(ikd) :: ndim
    integer(ikd) :: ndof
    integer(ikd) :: nels
    integer(ikd) :: neq
    integer(ikd) :: nip
    integer(ikd) :: nsp
    integer(ikd) :: nlen
    integer(ikd) :: nn
    integer(ikd) :: nod
    integer(ikd) :: nodof
    integer(ikd) :: nprops = 3
    integer(ikd) :: np_types
    integer(ikd) :: nst
    integer(ikd) :: unit
    integer(ikd) :: istep
    real(rkd) :: det
    real(rkd) :: penalty = 1.0e20_rkd
    real(rkd) :: zero = 0.0_rkd
    character(len=:), allocatable :: argv
    character(len=15) :: element
    character(len=20) :: fname
    character(len=20) :: ctype
    real(rkd), allocatable :: strsg(:,:,:)
    real(rkd), allocatable :: dd(:,:,:)
    integer(ikd), allocatable :: etype(:)
    integer(ikd), allocatable :: g(:)
    integer(ikd), allocatable :: g_g(:,:)
    integer(ikd), allocatable :: g_num(:,:)
    integer(ikd), allocatable :: kdiag(:)
    integer(ikd), allocatable :: nf(:,:)
    integer(ikd), allocatable :: no(:)
    integer(ikd), allocatable :: node(:)
    integer(ikd), allocatable :: num(:)
    integer(ikd), allocatable :: sense(:)
    real(rkd), allocatable :: bee(:,:)
    real(rkd), allocatable :: coord(:,:)
    real(rkd), allocatable :: dee(:,:)
    real(rkd), allocatable :: der(:,:)
    real(rkd), allocatable :: deriv(:,:)
    real(rkd), allocatable :: eld(:)
    real(rkd), allocatable :: fun(:)
    real(rkd), allocatable :: gc(:)
    real(rkd), allocatable :: gravlo(:)
    real(rkd), allocatable :: g_coord(:,:)
    real(rkd), allocatable :: jac(:,:)
    real(rkd), allocatable :: km(:,:)
    real(rkd), allocatable :: kv(:)
    real(rkd), allocatable :: loads(:)
    real(rkd), allocatable :: points(:,:)
    real(rkd), allocatable :: prop(:,:)
    real(rkd), allocatable :: sigma(:)
    real(rkd), allocatable :: value(:)
    real(rkd), allocatable :: weights(:)
    integer(ikd) :: nent
    integer(ikd), allocatable :: mfix(:,:)

    call getname(argv,nlen)
    open(10, file=argv(1:nlen)//'.dat')
    open(11, file=argv(1:nlen)//'.res')
    read(10,*) element, nod, nip, nodof, nst, ndim, np_types
    ndof = nod*nodof
    read(10,*) fname
    read(10,*) nent
    allocate(mfix(nodof+2,nent))
    read(10,*) mfix

    call readmsh(fname,mfix,nodof,ndim,nod,nn,nels,g_coord,nf,g_num)
    
    allocate(points(nip,ndim))
    allocate(dee(nst,nst))
    allocate(coord(nod,ndim))
    allocate(jac(ndim,ndim))
    allocate(weights(nip))
    allocate(num(nod))
    allocate(der(ndim,nod))
    allocate(deriv(ndim,nod))
    allocate(bee(nst,ndof))
    allocate(km(ndof,ndof))
    allocate(eld(ndof))
    allocate(sigma(nst))
    allocate(g(ndof))
    allocate(g_g(ndof,nels))
    allocate(gc(ndim))
    allocate(fun(nod))
    allocate(etype(nels))
    allocate(prop(nprops,np_types))
    
    read(10,*) prop
    etype = 1
    if (np_types.gt.1) read(10,*) etype

    call formnf(nf)
    neq = maxval(nf)
    allocate(kdiag(neq))
    allocate(loads(0:neq))
    allocate(gravlo(0:neq))
    kdiag = 0

    do iel = 1,nels
        num = g_num(:,iel)
        call num_to_g(num,nf,g)
        g_g(:,iel) = g
        call fkdiag(kdiag,g)
    end do
    do i = 2,neq
        kdiag(i) = kdiag(i) + kdiag(i-1)
    end do
    allocate(kv(kdiag(neq)))
    write(11,'(2(a,i5))') 'There are',neq,' equations and the skyline storage is',kdiag(neq)
    call sample(element,points,weights)
    kv = zero
    gravlo = zero
    do iel = 1,nels
        ! call fmdsig(dee,prop(1,etype(iel)),prop(2,etype(iel)))
        call deemat(dee,prop(1,etype(iel)),prop(2,etype(iel)))
        num = g_num(:,iel)
        coord = transpose(g_coord(:,num))
        g = g_g(:,iel)
        km = zero
        eld = zero
        do i = 1,nip
            call shape_fun(fun,points,i)
            call shape_der(der,points,i)
            jac = matmul(der,coord)
            det = determinant(jac)
            call invert(jac)
            deriv = matmul(jac,der)
            call beemat(bee,deriv)
            km = km + matmul(matmul(transpose(bee),dee),bee)*det*weights(i)
            eld(nodof:ndof:nodof) = eld(nodof:ndof:nodof) + fun(:)*det*weights(i)
        end do
        call fsparv(kv,km,g,kdiag)
        gravlo(g) = gravlo(g) - eld*prop(3,etype(iel))
    end do
    loads = zero
    read(10,*) loaded_nodes, (k, loads(nf(:,k)), i = 1,loaded_nodes)
    loads = loads + gravlo

    read(10,*) fixed_freedoms
    if (fixed_freedoms.ne.0) then
        allocate(node(fixed_freedoms))
        allocate(sense(fixed_freedoms))
        allocate(value(fixed_freedoms))
        allocate(no(fixed_freedoms))
        read(10,*) (node(i), sense(i), value(i), i = 1,fixed_freedoms)
        do i = 1,fixed_freedoms
            no(i) = nf(sense(i), node(i))
        end do
        kv(kdiag(no)) = kv(kdiag(no)) + penalty
        loads(no) = kv(kdiag(no))*value
    end if
    close(10)
    call sparin(kv,kdiag)
    call spabac(kv,loads,kdiag)
    loads(0) = zero
    if (ndim.eq.3) then
        write(11,'(/a)') '  Node   x-disp      y-disp      z-disp'
    else
        write(11,'(/a)') '  Node   x-disp      y-disp'
    end if
    do k = 1,nn
        write(11,'(i5,3e12.4)') k, loads(nf(:,k))
    end do

    ! stresses
    ! nip = 1
    allocate(strsg(nst,nip,nels))
    ! deallocate(points)
    ! deallocate(weights)
    ! allocate(points(nip,ndim))
    ! allocate(weights(nip))
    ! call sample(element,points,weights)

    ! vm streses
    nsp = 1
    allocate(dd(nsp,nip,nels))

    write(11,'(/a,i2,a)') 'The integration point (nip=',nip,') stresses are:'
    if (ndim.eq.3) then
        write(11,'(a,a)') '    Element     x-coord     y-coord     z-coord', &
                        '    sig_x       sig_y       sig_z       tau_xy      tau_yz      tau_zx'
    else
        write(11,'(a,a)') '    Element x-coord     y-coord', &
                          '          sig_x       sig_y       tau_xy'
    end if
    
    do iel = 1,nels
        ! call fmdsig(dee,prop(1,etype(iel)),prop(2,etype(iel)))
        call deemat(dee,prop(1,etype(iel)),prop(2,etype(iel)))
        num = g_num(:,iel)
        coord = transpose(g_coord(:,num))
        g = g_g(:,iel)
        eld = loads(g)
        do i = 1,nip
            call shape_der(der,points,i)
            call shape_fun(fun,points,i)
            gc = matmul(fun,coord)
            jac = matmul(der,coord)
            call invert(jac)
            deriv = matmul(jac,der)
            call beemat(bee,deriv)
            sigma = matmul(dee,matmul(bee,eld))
            if (ndim.eq.3) then
                write(11,'(i8,4x,3e12.4)',advance='no') iel, gc
                write(11,'(6e12.4)') sigma
            else
                write(11,'(i8,2e12.4,5x,3e12.4)') iel,gc,sigma
            end if
            strsg(:,i,iel) = sigma
            ! vm streses
            call invar(sigma, dd(1,i,iel))
            ! call invar(sigma, dd(1,i,iel), sigm, theta)
        end do
    end do
    close(11)

    ! GiD
    open(12,file = argv(1:nlen)//'.post.msh')
    open(13,file = argv(1:nlen)//'.post.res')

    if (element.eq.'quad') ctype = ' "Quadrilateral" '
    if (element.eq.'hexa') ctype = ' "Hexahedron" '
    if (element.eq.'tri') ctype = ' "Triangle" '
    if (element.eq.'tetra') ctype = ' "Tetrahedron" '
    i = 0
    unit = 12
    call mesh_gid(g_coord,g_num,etype,element,unit,i,ctype,nn)
    unit = 13
    call output_gid(strsg,loads,points,nf,1,1,unit,element,dd)
    istep = 1
    call smooth(g_coord,g_num,strsg,istep,nip,element,unit,etype,prop,dd)
    close(12)
    close(13)
    
end program pc1