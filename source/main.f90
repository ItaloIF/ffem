program main
    use glb
    use rout
    use tool
    use spar
    use mater
    use spmod
    use mgid
    use tsmod
    implicit none
! Variables
    ! Scalar integers
    integer :: nlen
    integer :: vers
    integer :: nod
    integer :: nip
    integer :: nodof
    integer :: nst
    integer :: ndim
    integer :: nptps
    integer :: nent
    integer :: nprops = 4
    integer :: nentp
    integer :: nstep
    integer :: miter
    integer :: npri
    integer :: nres
    integer :: nload
    integer :: i
    integer :: j
    integer :: it
    integer :: iter
    integer :: nn
    integer :: nels
    integer :: neq
    integer :: ndof
    integer :: nim
    integer :: isys
    integer :: nnz
    integer :: ntot
    integer :: ns
    integer :: idmod
    integer :: idsol
    integer :: unit

    ! Scalar reals
    real(rkd) :: toler
    real(rkd) :: dtim
    real(rkd) :: alfa
    real(rkd) :: beta
    real(rkd) :: gamma
    real(rkd) :: delta
    real(rkd) :: zero = 0.0_rkd
    real(rkd) :: one = 1.0_rkd
    real(rkd) :: pt5 = 0.5_rkd
    real(rkd) :: area
    real(rkd) :: rhow
    real(rkd) :: det
    real(rkd) :: dvol
    real(rkd) :: sy
    real(rkd) :: em
    real(rkd) :: pv
    real(rkd) :: time
    real(rkd) :: a0, a1, a2, a3, a4, a5, a6, a7
    real(rkd) :: resid
    real(rkd) :: retot
    real(rkd) :: ratio
    real(rkd) :: ft

    ! Dynamic integer arrays
    integer, allocatable :: mfix(:,:)
    integer, allocatable :: mprop(:,:)
    integer, allocatable :: tldid(:)
    integer, allocatable :: nf(:,:)
    integer, allocatable :: g_num(:,:)
    integer, allocatable :: etype(:)
    integer, allocatable :: g_g(:,:)
    integer, allocatable :: g(:)
    integer, allocatable :: num(:)
    integer, allocatable :: ia(:)
    integer, allocatable :: ja(:)
    integer, allocatable :: indu(:)
    ! Dynamic real arrays
    real(rkd), allocatable :: prop(:,:)
    real(rkd), allocatable :: tload(:,:)
    real(rkd), allocatable :: g_coord(:,:)
    real(rkd), allocatable :: loads(:)
    real(rkd), allocatable :: kv(:)
    real(rkd), allocatable :: km(:,:)
    real(rkd), allocatable :: mv(:)
    real(rkd), allocatable :: mm(:,:)
    real(rkd), allocatable :: cv(:)
    real(rkd), allocatable :: kp(:)
    real(rkd), allocatable :: weights(:)
    real(rkd), allocatable :: points(:,:)
    real(rkd), allocatable :: coord(:,:)
    real(rkd), allocatable :: der(:,:)
    real(rkd), allocatable :: deriv(:,:)
    real(rkd), allocatable :: dee(:,:)
    real(rkd), allocatable :: bee(:,:)
    real(rkd), allocatable :: fun(:)
    real(rkd), allocatable :: ecm(:,:)
    real(rkd), allocatable :: jac(:,:)
    real(rkd), allocatable :: effst(:,:)
    real(rkd), allocatable :: strsg(:,:,:)
    real(rkd), allocatable :: estrt(:,:,:)
    real(rkd), allocatable :: estri(:)
    real(rkd), allocatable :: stres(:)
    real(rkd), allocatable :: desig(:)
    real(rkd), allocatable :: osv(:,:,:)
    real(rkd), allocatable :: pl(:,:)
    real(rkd), allocatable :: eld(:)
    real(rkd), allocatable :: eload(:)
    real(rkd), allocatable :: bload(:)

    real(rkd), allocatable :: d(:)
    real(rkd), allocatable :: v(:)
    real(rkd), allocatable :: a(:)
    real(rkd), allocatable :: vc(:)
    real(rkd), allocatable :: kd(:)
    real(rkd), allocatable :: rt(:)
    real(rkd), allocatable :: u0(:)
    real(rkd), allocatable :: u10(:)
    real(rkd), allocatable :: u20(:)
    real(rkd), allocatable :: out(:)
    ! Static characters
    character(len=10) :: dumm
    character(len=10) :: fetp
    character(len=60) :: fmesh
    character(len=20) :: ctype
    ! Dynamic characters
    character(len=:), allocatable :: fdat               ! file input
    ! Logical
    logical :: mcons = .true.
    
! Input Data
    call getname(fdat,nlen)
    open(1, file=fdat(1:nlen)//'.dat')
    vers = 1
    idmod = 1
    idsol = 1
    if (vers.eq.1) then
        read(1,'(32x)',advance='no'); read(1,*) fetp
        read(1,'(32x)',advance='no'); read(1,*) nod
        read(1,'(32x)',advance='no'); read(1,*) nip
        read(1,'(32x)',advance='no'); read(1,*) nodof
        read(1,'(32x)',advance='no'); read(1,*) nst
        read(1,'(32x)',advance='no'); read(1,*) ndim
        read(1,'(32x)',advance='no'); read(1,*) fmesh
        read(1,'(32x)',advance='no'); read(1,*) nent
        allocate(mfix(nodof+2,nent))
        read(1,*)
        read(1,*) mfix
        read(1,'(32x)',advance='no'); read(1,*) nptps
        allocate(prop(nprops,nptps))
        read(1,*)
        read(1,*) prop
        read(1,'(32x)',advance='no'); read(1,*) nentp
        allocate(mprop(nodof+1,nentp))
        read(1,*)
        read(1,*) mprop
        read(1,'(32x)',advance='no'); read(1,*) toler
        read(1,'(32x)',advance='no'); read(1,*) nstep
        read(1,'(32x)',advance='no'); read(1,*) miter
        read(1,'(32x)',advance='no'); read(1,*) dtim
        read(1,'(32x)',advance='no'); read(1,*) alfa
        read(1,'(32x)',advance='no'); read(1,*) beta
        read(1,'(32x)',advance='no'); read(1,*) gamma
        read(1,'(32x)',advance='no'); read(1,*) delta
        read(1,'(32x)',advance='no'); read(1,*) npri
        read(1,'(32x)',advance='no'); read(1,*) nres
        read(1,'(32x)',advance='no'); read(1,*) nload
        allocate(tldid(nload))
        allocate(tload(nodof,nload))
        read(1,*)
        read(1,*) (tldid(i),tload(:,i),i=1,nload)
    end if
    close(1)
! Read msh file and initialisation
    call readmsh(fmesh,mfix,nodof,ndim,nod,nn,nels,g_coord,nf,g_num)
    allocate(etype(nels))  ! temporal
    etype = 1
    call formnf(nf)
    neq = maxval(nf)
    
    ndof = nod*nodof
    allocate(g_g(ndof,nels))

    allocate(loads(0:neq))
    loads = zero
    do i = 1,nload
        loads(nf(:,tldid(i))) = tload(:,i)
    end do
    
! Find global array sizes and loop elements
    allocate(num(nod))
    allocate(g(ndof))
    do i = 1,nels
        num = g_num(:,i)
        call num_to_g(num,nf,g)
        g_g(:,i) = g
    end do
    ntot = nod*ndim
    ns = nels*(ntot*ntot)
! Allocate arrays
    allocate(ia(ns))
    allocate(ja(ns))
    allocate(mv(ns))
    allocate(indu(neq))
    allocate(mm(ndof,ndof))
    allocate(km(ndof,ndof))
    allocate(coord(nod,ndim))
    allocate(der(ndim,nod))
    allocate(jac(ndim,ndim))
    allocate(deriv(ndim,nod))
    allocate(dee(nst,nst))
    allocate(bee(nst,ndof))
    allocate (fun(nod))
    allocate(ecm(ndof,ndof))
    allocate(effst(nip,nels))
    allocate(strsg(nst,nip,nels))
    allocate(estrt(nst,nip,nels))
    allocate(estri(nst))
    allocate(desig(nst))
    allocate(stres(nst))
    allocate(osv(nst,nip,nels))
    allocate(pl(nst,nst))
    allocate(eld(ndof))
    allocate(eload(neq))
    allocate(bload(neq))

    strsg = zero
    estrt = zero
    osv = zero

    allocate(d(neq))
    allocate(v(neq))
    allocate(a(neq))
    allocate(vc(neq))
    allocate(kd(neq))
    allocate(rt(neq))
    allocate(u0(neq))
    allocate(u10(neq))
    allocate(u20(neq))

    allocate(out(0:neq))

! Mass matrix assembly (with another number of integration points)
    nim = nip  ! modify for read this value
    nnz = 0
    isys = 0
    allocate(points(nim,ndim))
    allocate(weights(nim))
    call sample(fetp,points,weights)
    mv = zero
    do i  = 1,nels
        num = g_num(:,i)
        coord = transpose(g_coord(:,num))
        g = g_g(:,i)
        mm = zero
        area = zero
        rhow = prop(3,etype(i))
        do j = 1,nim
            call shape_fun(fun,points,j)
            call shape_der(der,points,j)
            jac = matmul(der,coord)
            det = determinant(jac)
            dvol = det*weights(j)
            area = area + det*weights(j)
            if (mcons) then
                call ecmat(ecm,fun)
                mm = mm + ecm*rhow*dvol
            end if
        end do
        if (.not.mcons) call elmat(area,rhow,mm)
        call asemcoo2(g,mv,mm,nnz,ia,ja,isys)
    end do
! Sparse form
    call coicsr(1,neq,mv(1:nnz),ia(1:nnz),ja(1:nnz))
    call clncsr2(3,1,neq,mv,ia,ja,indu)
    nnz = ia(neq+1) - 1
    allocate(kv(nnz))
    allocate(cv(nnz))
    allocate(kp(nnz))
! Initial stiffness matrix assembly
    deallocate(points)
    deallocate(weights)
    allocate(points(nip,ndim))
    allocate(weights(nip))
    call sample(fetp,points,weights)
    kv = zero
    do i = 1,nels
        em = prop(1,etype(i))
        pv = prop(2,etype(i))
        call deemat(dee,em,pv)
        num = g_num(:,i)
        coord = transpose(g_coord(:,num))
        g = g_g(:,i)
        km = zero
        do j = 1,nip
            call shape_fun(fun,points,j)
            call shape_der(der,points,j)
            jac = matmul(der,coord)
            det = determinant(jac)
            call invert(jac)
            deriv = matmul(jac,der)
            call beemat(bee,deriv)
            km = km + matmul(matmul(transpose(bee),dee),bee)*det*weights(j)
        end do
        call asemspr(g,kv,km,ia(1:neq+1),ja(1:nnz),isys)
    end do

! Damping matrix
    cv = zero
    cv = alfa*mv + beta*kv

! Initial conditions: d, v, a (displacement, velocity, aceleration)
    d = zero
    v = zero
    a = zero
    rt = zero
    vc = zero
    kd = zero
    u0 = zero
    u10 = zero
    u20 = zero

    a0 = dtim*dtim*(pt5 - delta)
    a1 = dtim*(one - gamma)
    a2 = dtim*dtim*delta
    a3 = dtim*gamma
    a4 = one/a2
    a5 = dtim*beta*gamma
    a6 = dtim*alfa*gamma
    a7 = one + a6

    call amux(cv,ia,ja,v,vc)
    call amux(kv,ia,ja,d,kd)
    call tseries(3,zero,1.0_rkd,ft,0.0_rkd,1.0_rkd,0.1_rkd)
    rt = loads(1:neq)*ft
    a = rt - vc - kd
    call sp_solver(idsol,neq,nnz,mv,ia,ja,a)

! Time stepping loop
    time = zero
    do it = 1,nstep
        time =  time + dtim
        ! Stiffness matrix asesembly
        kv = zero
        do i = 1,nels
            em = prop(1,etype(i))
            pv = prop(2,etype(i))
            sy = prop(4,etype(i)) 
            num = g_num(:,i)
            coord = transpose(g_coord(:,num))
            g = g_g(:,i)
            km = zero
            do j = 1,nip
                call shape_der(der,points,j)
                jac = matmul(der,coord)
                det = determinant(jac)
                call invert(jac)
                deriv = matmul(jac,der)
                call beemat(bee,deriv)
                call models(5,idmod,prop(:,etype(i)),osv(:,j,i),strsg(:,j,i),dee)
                km = km + matmul(matmul(transpose(bee),dee),bee)*det*weights(j)
            end do
            call asemspr(g,kv,km,ia(1:neq+1),ja(1:nnz),isys)
        end do
        ! Predict displacement, velocity and acceleration
        u0 = d
        u10 = v
        d = d + dtim * v + a0 * a
        v = v + a1 * a
        u20 = a4 * (u0 - d)

        rt = zero
        call tseries(3,time,1.0_rkd,ft,0.0_rkd,1.0_rkd,0.1_rkd)
        rt = loads * ft
        cv = alfa*mv + beta*kv
        kp = a7 * a4 * mv + kv * (a5 * a4 + one)
        ! Iteration loop
        do iter = 1,miter
            ! Update stresses
            eload = zero
            do i = 1,nels
                em = prop(1,etype(i))
                pv = prop(2,etype(i))
                sy = prop(4,etype(i)) 
                num = g_num(:,i)
                coord = transpose(g_coord(:,num))
                g = g_g(:,i)
                eld = u0(g)
                bload = zero
                do j = 1,nip
                    call shape_der(der,points,j)
                    jac = matmul(der,coord)
                    det = determinant(jac)
                    call invert(jac)
                    deriv = matmul(jac,der)
                    call beemat(bee,deriv)
                    estri = matmul(bee,eld)
                    estri = estri - estrt(:,j,i)
                    estrt(:,j,i) = estrt(:,j,i) + estri
                    desig = matmul(dee,estri)
                    stres = strsg(:,j,i) + desig
                    call models(3,idmod,prop(:,etype(i)),osv(:,j,i),strsg(:,j,i),dee, &
                                stres,estri,estrt(:,j,i))
                    strsg(:,j,i) = stres
                    bload  = bload + matmul(transpose(bee),stres)*det*weights(j)
                end do
                eload(g) = eload(g) + bload
            end do
            ! Update displacement and check convergence
            call amux(cv,ia,ja,u10,vc)
            call amux(mv,ia,ja,u20,kd)
            kd = rt - vc - kd - eload
            call sp_solver(idsol,neq,nnz,kp,ia,ja,kd)
            u0 = u0 + kd
            u20 = (u0 - d)*a4
            u10 = v + a3*u20
            resid = dot_product(kd,kd)
            resid = dsqrt(resid)
            retot = dot_product(u0,u0)
            retot = dsqrt(retot)
            ratio = resid/retot
            write(*,*) iter,ratio
            if (ratio.le.toler) exit
            if (iter.eq.miter) stop
        end do
        d = u0
        v = u10
        a = u20
    end do

    ! GiD
    open(12,file = fdat(1:nlen)//'.post.msh')
    open(13,file = fdat(1:nlen)//'.post.res')

    if (fetp.eq.'quad') ctype = ' "Quadrilateral" '
    if (fetp.eq.'hexa') ctype = ' "Hexahedron" '
    if (fetp.eq.'tri') ctype = ' "Triangle" '
    if (fetp.eq.'tetra') ctype = ' "Tetrahedron" '
    i = 0
    unit = 12
    call mesh_gid(g_coord,g_num,etype,fetp,unit,i,ctype,nn)
    unit = 13
    out = zero
    out(1:) = d
    call resp_gid(strsg,out,points,nf,1,1,unit,fetp)
    close(12)
    close(13)

end program main

