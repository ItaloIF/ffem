subroutine geela10(prop,dee)
    use glb
    implicit none
    real(rkd), intent(in) :: prop(:)
    real(rkd), intent(out) :: dee(:,:)

    real(rkd) :: e
    real(rkd) :: v

    e = prop(1)
    v = prop(2)
    call deemat(dee,e,v)
    return
end subroutine geela10

subroutine gestr10(prop,tau,deps,ddd,st)
    use glb
    implicit none
    real(rkd), intent(in) :: prop(:)
    real(rkd), intent(in out) :: tau(:)
    real(rkd), intent(in) :: deps(:)
    real(rkd), intent(in out) :: ddd(:)
    real(rkd), intent(out) :: st(:)

    real(rkd) :: e
    real(rkd) :: v
    real(rkd) :: uniax
    real(rkd) :: zero = 0.0_rkd
    real(rkd) :: one = 1.0_rkd
    real(rkd) :: d3 = 3.0_rkd
    real(rkd) :: d6 = 6.0_rkd
    real(rkd) :: root3
    real(rkd), allocatable :: dee(:,:)
    real(rkd), allocatable :: desig(:,:)
    real(rkd), allocatable :: devia(:,:)
    real(rkd) :: steff
    real(rkd) :: theta
    real(rkd) :: varj2
    real(rkd) :: yield
    real(rkd) :: preys
    real(rkd) :: hards
    real(rkd) :: epstc
    real(rkd) :: escur
    integer ::  nst
    integer ::  ncrit
    integer ::  iasoc
    logical :: flag = .false.

    ih = size(tau)
    root3 = dsqrt(d3)
    allocate(dee(ih,ih))
    allocate(desig(ih))
    allocate(devia(ih))
    
    ! material parameters
    e = prop(1)
    v = prop(2)
    uniax = prop(5)
    hards = prop(6)
    frict = prop(7)
    ncrit = nint(prop(8))
    iasoc = nint(prop(9))
    ! state variables
    effsc = ddd(1)
    epstc = ddd(2)
    nplas = nint(ddd(8))

    if (ncrit.eq.3) uniax = uniax*cos(frict*0.017453292_rkd)
    if (ncrit.eq.4) uniax = d6*uniax*cos(frict*0.017453292_rkd) / &
                            (root3*(d3-sin(frict*0.017453292_rkd)))
    ! Calculate effective stress and trial values
    call deemat(dee,e,v)
    desig = matmul(dee,deps)
    st = tau + desig
    call invar10(prop,st,ncrit,devia,steff,theta,varj2,yield)
    preys = uniax + hards*abs(epstc)

    ! if the gauss point(g.p.) has yield in the previous interation
    if (nplas.eq.1) then
        ! ask if the current g.p. is unloading in this iteration
        escur = yield - effsc
        if (escur.le.zero) then
            ! elastic behavior - unloading
            nplas = 0
            effsc = yield
            if (escur.eq.zero) nplas = 1
            flag = .false.
        else
            ! plastic behaviour - loading
            nplas = 1
            efact = one
            flag = .true.
        end if
    else
        ! point has not previously yield
        escur = yield - preys
        if (escur.le.zero) then
            ! elastic behaviour - loading or unloading
            nplas = 0
            st = tau + desig
            effsc = yield
            if (escur.eq.zero) nplas = 1
            flag = .false.
        else
            nplas = 1
            rfact = escur/(yield - effsc)
            fy0 = effsc - preys
            fy1 = yield - preys
            alfa0 = zero
            alfa1 = one
            call pegasus()
            rfact = one - alpha
            flag = .true.
        end if
    end if

    return
end subroutine gestr10

subroutine invar10(prop,stemp,ncrit,devia,steff,theta,varj2,yield)
    use glb
    implicit none
    real(rkd), intent(in) :: prop(:)
    real(rkd), intent(in) :: stemp(:)
    integer, intent(in) :: ncrit
    real(rkd), intent(out) :: devia(:)
    real(rkd), intent(out) :: steff
    real(rkd), intent(out) :: theta
    real(rkd), intent(out) :: varj2
    real(rkd), intent(out) :: yield

    integer :: ih
    real(rkd) :: root3
    real(rkd) :: zero = 0.0_rkd
    real(rkd) :: one = 1.0_rkd
    real(rkd) :: d2 = 2.0_rkd
    real(rkd) :: d3 = 3.0_rkd
    real(rkd) :: d6 = 6.0_rkd
    real(rkd) :: pt5 = 0.5_rkd
    real(rkd) :: smean
    real(rkd) :: varj3
    real(rkd) :: phira
    real(rkd) :: snphi
    real(rkd) :: sint3

    ih = size(devia)
    root3 = dsqrt(d3)
    select case (ih)
        case (4)
            smean = (stemp(1) + stemp(2) + stemp(4))/d3
            devia(1) = stemp(1) - smean
            devia(2) = stemp(2) - smean
            devia(3) = stemp(3)
            devia(4) = stemp(4) - smean
            varj2 = devia(3)*devia(3) + pt5*(devia(1)*devia(1)+devia(2)*devia(2)+devia(4)*devia(4))
            varj3 = devia(4)*(devia(4)*devia(4) - varj2)
        case (6)
            smean = (stemp(1) + stemp(2) + stemp(3))/d3
            devia(1) = stemp(1) - smean
            devia(2) = stemp(2) - smean
            devia(3) = stemp(3) - smean
            devia(4) = stemp(4)
            devia(5) = stemp(5)
            devia(6) = stemp(6)
            varj2 = devia(4)*devia(4) + devia(5)*devia(5) + devia(6)*devia(6) + &
                    pt5*(devia(1)*devia(1)+devia(2)*devia(2)+devia(3)*devia(3))
            varj3 = devia(1)*devia(2)*devia(3) + d2*devia(4)*devia(5)*devia(6) - &
                    devia(1)*devia(5)**2 - devia(2)*devia(6)**2 - devia(3)*devia(4)**2
    end select
    steff = dsqrt(varj2)
    if (steff.eq.zero) then
        sint3 = zero
    else
        sint3 = -d3*dsqrt(d3)*varj3/(d2*varj2*steff)
    end if
    if (sint3.lt.-one) sint3 = -one
    if (sint3.gt.one) sint3 = one
    theta = asin(sint3)/d3
    select case (ncrit)
        case (1)
            ! Tresca
            yield = d2*cos(theta)*steff
        case (2)
            ! Von-Mises
            yield = root3*steff
        case (3)
            ! Mohr Coulomb
            phira = prop(7)*0.017453292_rkd
            snphi = sin(phira)
            yield = smean*snphi + steff*(cos(theta) - sin(theta)*snphi/root3)
        case (4)
            ! Drucker-Prager
            phira = prop(7)*0.017453292_rkd 
            snphi = sin(phira)
            yield = d6*smean*snphi/(root3*(d3-snphi)) + steff
    end select
    return
end subroutine invar10

subroutine pegasus()
    use glb

    integer :: i
    do while (i.le.maxit)
        i = i + 1
        
    end do
end subroutine pegasus