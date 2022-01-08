program main
    use glb
    use rout
    use tool
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
    integer :: nprops = 3
    integer :: nentp
    integer :: nstep
    integer :: miter
    integer :: npri
    integer :: nres
    integer :: nload
    integer :: i
    integer :: nn
    integer :: nels
    integer :: neq
    integer :: ndof
    integer :: nim
    ! Scalar reals
    real(rkd) :: toler
    real(rkd) :: dtim
    real(rkd) :: alfa
    real(rkd) :: beta
    real(rkd) :: gamma
    real(rkd) :: delta
    real(rkd) :: zero = 0.0_rkd
    ! Dynamic integer arrays
    integer, allocatable :: mfix(:,:)
    integer, allocatable :: mprop(:,:)
    integer, allocatable :: tldid(:)
    integer, allocatable :: nf(:,:)
    integer, allocatable :: g_num(:,:)
    integer, allocatable :: etype(:)
    integer, allocatable :: kdiag(:)
    integer, allocatable :: g_g(:,:)
    integer, allocatable :: g(:)
    integer, allocatable :: num(:)
    ! Dynamic real arrays
    real(rkd), allocatable :: prop(:,:)
    real(rkd), allocatable :: tload(:,:)
    real(rkd), allocatable :: g_coord(:,:)
    real(rkd), allocatable :: loads(:)
    real(rkd), allocatable :: kv(:)
    real(rkd), allocatable :: weights(:)
    real(rkd), allocatable :: points(:,:)
    ! Static characters
    character(len=10) :: dumm
    character(len=10) :: fetp
    character(len=20) :: fmesh
    ! Dynamic characters
    character(len=:), allocatable :: fdat               ! file input
    
! Input Data
    call getname(fdat,nlen)
    open(1, file=fdat(1:nlen)//'.dat')
    vers = 1
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
    
! Find global array sizes
    allocate(kdiag(neq))
    allocate(num(nod))
    allocate(g(ndof))
    kdiag = 0
    do i = 1,nels
        num = g_num(:,i)
        call num_to_g(num,nf,g)
        g_g(:,i) = g
        call fkdiag(kdiag,g)
    end do
    do i = 2,neq
        kdiag(i) = kdiag(i) + kdiag(i-1)
    end do
! Allocate arrays
    allocate(kv(kdiag(neq)))

! Mass matrix assembly with another number of integration points
    nim = nip  ! modify for read this value
    allocate(points(nim,ndim))
    allocate(weights(nim))
    do i  = 1,nels
      num = g_num(:,i)
    end do

end program main