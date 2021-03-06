subroutine shape_der(der,points,i)
    ! this subroutine produces derivatives of shape functions withe respect to local coordinates.
    use glb
    implicit none
    integer, intent(in) :: i
    real(rkd), intent(in) :: points(:,:)
    real(rkd), intent(out) :: der(:,:)
    real(rkd) :: eta
    real(rkd) :: xi
    real(rkd) :: zeta
    real(rkd) :: xi0
    real(rkd) :: eta0
    real(rkd) :: zeta0
    real(rkd) :: etam
    real(rkd) :: etap
    real(rkd) :: xim
    real(rkd) :: xip
    real(rkd) :: c1, c2, c3
    real(rkd) :: t1, t2, t3, t4, t5, t6, t7, t8, t9
    real(rkd) :: x2p1
    real(rkd) :: x2m1
    real(rkd) :: e2p1
    real(rkd) :: e2m1
    real(rkd) :: zetam
    real(rkd) :: zetap
    real(rkd), parameter :: zero = 0.0_rkd
    real(rkd), parameter :: pt125 = 0.125_rkd
    real(rkd), parameter :: pt25 = 0.25_rkd
    real(rkd), parameter :: pt5 = 0.5_rkd
    real(rkd), parameter :: pt75 = 0.75_rkd
    real(rkd), parameter :: one = 1.0_rkd
    real(rkd), parameter :: two = 2.0_rkd
    real(rkd), parameter :: d3 = 3.0_rkd
    real(rkd), parameter :: d4 = 4.0_rkd
    real(rkd), parameter :: d5 = 5.0_rkd
    real(rkd), parameter :: d6 = 6.0_rkd
    real(rkd), parameter :: d8 = 8.0_rkd
    real(rkd), parameter :: d9 = 9.0_rkd
    real(rkd), parameter :: d10 = 10.0_rkd
    real(rkd), parameter :: d11 = 11.0_rkd
    real(rkd), parameter :: d12 = 12.0_rkd
    real(rkd), parameter :: d16 = 16.0_rkd
    real(rkd), parameter :: d18 = 18.0_rkd
    real(rkd), parameter :: d27 = 27.0_rkd
    real(rkd), parameter :: d32 = 32.0_rkd
    real(rkd), parameter :: d36 = 36.0_rkd
    real(rkd), parameter :: d54 = 54.0_rkd
    real(rkd), parameter :: d64 = 64.0_rkd
    real(rkd), parameter :: d128 = 128_rkd
    integer :: xii(20)
    integer :: etai(20)
    integer :: zetai(20)
    integer :: l
    integer :: ndim
    integer :: nod

    ndim = ubound(der,1)
    nod = ubound(der,2)
    select case(ndim)
        case (1)
            ! one dimensional elements
            xi = points(i,1)
            select case(nod)
                case (2)
                    der(1,1) = -pt5
                    der(1,2) = pt5
                case (3)
                    t1 = -one-xi
                    t2 = -xi
                    t3 = one-xi
                    der(1,1) = -(t3+t2)/two
                    der(1,2) = (t3+t1)
                    der(1,3) = -(t2+t1)/two
                case (4)
                    t1 = -one-xi
                    t2 = -one/d3-xi
                    t3 = one/d3-xi
                    t4 = one-xi
                    der(1,1) = -(t3*t4+t2*t4+t2*t3)*d9/d16
                    der(1,2) = (t3*t4+t1*t4+t1*t3)*d27/d16
                    der(1,3) = -(t2*t4+t1*t4+t1*t2)*d27/d16
                    der(1,4) = (t2*t3+t1*t3+t1*t2)*d9/d16
                case (5)
                    t1 = -one-xi 
                    t2 = -pt5-xi
                    t3 = -xi
                    t4 = pt5-xi
                    t5 = one-xi
                    der(1,1) = -(t3*t4*t5+t2*t4*t5+t2*t3*t5+t2*t3*t4)*two/d3
                    der(1,2) = (t3*t4*t5+t1*t4*t5+t1*t3*t5+t1*t3*t4)*d8/d3
                    der(1,3) = -(t2*t4*t5+t1*t4*t5+t1*t2*t5+t1*t2*t4)*d4
                    der(1,4) = (t2*t3*t5+t1*t3*t5+t1*t2*t5+t1*t2*t3)*d8/d3
                    der(1,5) = -(t2*t3*t4+t1*t3*t4+t1*t2*t4+t1*t2*t3)*two/d3
                case default
                    write(*,*) 'wrong number of nodes in shape_der'
                end select
        case (2)
            ! two dimensional elements
            xi=points(i,1)
            eta = points(i,2)
            c1 = xi 
            c2 = eta 
            c3 = one - c1 - c2
            etam = pt25*(one-eta)
            etap = pt25*(one+eta)
            xim = pt25*(one-xi)
            xip = pt25*(one+xi)
            x2p1 = two*xi + one
            x2m1 = two*xi - one
            e2p1 = two*eta + one
            e2m1 = two*eta - one
            select case(nod)
                case (3)
                    der(1,1) = one
                    der(1,3) = zero
                    der(1,2) = -one
                    der(2,1) = zero
                    der(2,3) = one
                    der(2,2) = -one
                case (6)
                    der(1,1) = d4*c1-one
                    der(1,6) = d4*c2
                    der(1,5) = zero
                    der(1,4) = -d4*c2
                    der(1,3) = -(d4*c3-one)
                    der(1,2) = d4*(c3-c1)
                    der(2,1) = zero
                    der(2,6) = d4*c1
                    der(2,5) = d4*c2-one
                    der(2,4) = d4*(c3-c2)
                    der(2,3) = -(d4*c3-one)
                    der(2,2) = -d4*c1
                case (10)
                    der(1,1) = (d27*c1**2-d18*c1+two)/two
                    der(1,9) = (d9*(d6*c1-one)*c2)/two
                    der(1,8) = (d9*(d3*c2-one)*c2)/two
                    der(1,7) = zero
                    der(1,6) = -(d9*(d3*c2-one)*c2)/two
                    der(1,5) = (d9*(d6*c1+d6*c2-d5)*c2)/two
                    der(1,4) = -(d27*c1**2+d54*c1*c2-d36*c1+d27*c2**2-d36*c2+d11)/two
                    der(1,3) = (d9*(d9*c1**2+d12*c1*c2-d10*c1+d3*c2**2-d5*c2+two))/two
                    der(1,2) = -(d9*(d9*c1**2+d6*c1*c2-d8*c1-c2+one))/two
                    der(1,10) = -d27*(((c2-one)+c1)+c1)*c2
                    der(2,1) = zero
                    der(2,9) = (d9*(d3*c1-one)*c1)/two
                    der(2,8) = (d9*(d6*c2-one)*c1)/two
                    der(2,7) = (d27*c2**2-d18*c2+two)/two
                    der(2,6) = -(d9*((c1+c2-one)*(d6*c2-one)+(d3*c2-one)*c2))/two
                    der(2,5) = (d9*(d3*c1**2+d12*c1*c2-d5*c1+d9*c2**2-d10*c2+two))/two
                    der(2,4) = -(d27*c1**2+d54*c1*c2-d36*c1+d27*c2**2-d36*c2+d11)/two
                    der(2,3) = (d9*(d6*c1+d6*c2-d5)*c1)/two
                    der(2,2) = -(d9*(d3*c1-one)*c1)/two
                    der(2,10) = -d27*(((c2-one)+c1)+c2)*c1
                case (15)
                    t1 = c1-pt25  
                    t2 = c1-pt5 
                    t3 = c1-pt75
                    t4 = c2-pt25
                    t5 = c2-pt5
                    t6 = c2-pt75
                    t7 = c3-pt25
                    t8 = c3-pt5
                    t9 = c3-pt75
                    der(1,1) = d32/d3*(t2*t3*(t1+c1)+c1*t1*(t3+t2))
                    der(1,12) = d128/d3*c2*(t2*(t1+c1)+c1*t1)
                    der(1,11) = d64*c2*t4*(t1+c1)
                    der(1,10) = d128/d3*c2*t4*t5
                    der(1,9) = zero
                    der(1,8) = -d128/d3*c2*t4*t5
                    der(1,7) = -d64*c2*t4*(t7+c3)
                    der(1,6) = -d128/d3*c2*(t8*(t7+c3)+c3*t7)
                    der(1,5) = -d32/d3*(t8*t9*(t7+c3)+c3*t7*(t8+t9))
                    der(1,4) = d128/d3*(c3*t7*t8-c1*(t8*(t7+c3)+c3*t7))
                    der(1,3) = d64*(c3*t7*(t1+c1)-c1*t1*(t7+c3))
                    der(1,2) = d128/d3*(c3*(t2*(t1+c1)+c1*t1)-c1*t1*t2)
                    der(1,13) = d128*c2*(c3*(t1+c1)-c1*t1)
                    der(1,15) = d128*c2*t4*(c3-c1)
                    der(1,14) = d128*c2*(c3*t7-c1*(t7+c3))
                    der(2,1) = zero
                    der(2,12) = d128/d3*c1*t1*t2
                    der(2,11) = d64*c1*t1*(t4+c2)
                    der(2,10) = d128/d3*c1*(t5*(t4+c2)+c2*t4)
                    der(2,9) = d32/d3*(t5*t6*(t4+c2)+c2*t4*(t6+t5))
                    der(2,8) = d128/d3*((c3*(t5*(t4+c2)+c2*t4))-c2*t4*t5)
                    der(2,7) = d64*(c3*t7*(t4+c2)-c2*t4*(t7+c3))
                    der(2,6) = d128/d3*(c3*t7*t8-c2*(t8*(t7+c3)+c3*t7))
                    der(2,5) = -d32/d3*(t8*t9*(t7+c3)+c3*t7*(t8+t9))
                    der(2,4) = -d128/d3*c1*(t8*(t7+c3)+c3*t7)
                    der(2,3) = -d64*c1*t1*(t7+c3)
                    der(2,2) = -d128/d3*c1*t1*t2
                    der(2,13) = d128*c1*t1*(c3-c2)
                    der(2,15) = d128*c1*(c3*(t4+c2)-c2*t4)
                    der(2,14) = d128*c1*(c3*t7-c2*(c3+t7))
                case (4)
                    der(1,1) = -etam
                    der(1,2) = -etap
                    der(1,3) = etap
                    der(1,4) = etam
                    der(2,1) = -xim
                    der(2,2) = xim
                    der(2,3) = xip
                    der(2,4) = -xip
                case (5)
                    der(1,1) = -etam+pt5*xi*(one-eta**2)
                    der(1,2) = -etap+pt5*xi*(one-eta**2)
                    der(1,3) = etap+pt5*xi*(one-eta**2)
                    der(1,4) = etam+pt5*xi*(one-eta**2)
                    der(1,5) = -two*xi*(one-eta**2)
                    der(2,1) = -xim+pt5*eta*(one-xi**2)
                    der(2,2) = xim+pt5*eta*(one-xi**2)
                    der(2,3) = xip+pt5*eta*(one-xi**2)
                    der(2,4) = -xip+pt5*eta*(one-xi**2)
                    der(2,5) = -two*eta*(one-xi**2)
                case (8)
                    der(1,1) = etam*(two*xi+eta)
                    der(1,2) = -d8*etam*etap
                    der(1,3) = etap*(two*xi-eta)
                    der(1,4) = -d4*etap*xi
                    der(1,5) = etap*(two*xi+eta)
                    der(1,6) = d8*etap*etam
                    der(1,7) = etam*(two*xi-eta)
                    der(1,8) = -d4*etam*xi
                    der(2,1) = xim*(xi+two*eta)
                    der(2,2) = -d4*xim*eta
                    der(2,3) = xim*(two*eta-xi)
                    der(2,4) = d8*xim*xip
                    der(2,5) = xip*(xi+two*eta)
                    der(2,6) = -d4*xip*eta
                    der(2,7) = xip*(two*eta-xi)
                    der(2,8) = -d8*xim*xip 
                case (9)
                    etam = eta - one
                    etap = eta + one
                    xim = xi - one
                    xip = xi + one
                    der(1,1) = pt25*x2m1*eta*etam
                    der(1,2) = -pt5*x2m1*etap*etam
                    der(1,3) = pt25*x2m1*eta*etap
                    der(1,4) = -xi*eta*etap
                    der(1,5) = pt25*x2p1*eta*etap
                    der(1,6) = -pt5*x2p1*etap*etam
                    der(1,7) = pt25*x2p1*eta*etam
                    der(1,8) = -xi*eta*etam
                    der(1,9) = two*xi*etap*etam
                    der(2,1) = pt25*xi*xim*e2m1
                    der(2,2) = -xi*xim*eta
                    der(2,3) = pt25*xi*xim*e2p1
                    der(2,4) = -pt5*xip*xim*e2p1
                    der(2,5) = pt25*xi*xip*e2p1
                    der(2,6) = -xi*xip*eta
                    der(2,7) = pt25*xi*xip*e2m1
                    der(2,8) = -pt5*xip*xim*e2m1
                    der(2,9) = two*xip*xim*eta
                case default
                    write(*,*) 'wrong number of nodes in shape_der'
            end select
        case (3)
            ! three dimensional elements
            xi = points(i,1)
            eta = points(i,2)
            zeta = points(i,3)
            etam = one - eta
            xim = one - xi
            zetam = one - zeta
            etap = eta + one
            xip = xi + one
            zetap = zeta + one
            select case(nod)
                case (4)
                    der(1:3,1:4) = zero
                    der(1,1) = one
                    der(2,2) = one
                    der(3,3) = one
                    der(1,4) = -one
                    der(2,4) = -one
                    der(3,4) = -one
                case (8)
                    der(1,1) = -pt125*etam*zetam
                    der(1,2) = -pt125*etam*zetap
                    der(1,3) = pt125*etam*zetap
                    der(1,4) = pt125*etam*zetam
                    der(1,5) = -pt125*etap*zetam
                    der(1,6) = -pt125*etap*zetap
                    der(1,7) = pt125*etap*zetap
                    der(1,8) = pt125*etap*zetam
                    der(2,1) = -pt125*xim*zetam
                    der(2,2) = -pt125*xim*zetap
                    der(2,3) = -pt125*xip*zetap
                    der(2,4) = -pt125*xip*zetam
                    der(2,5) = pt125*xim*zetam
                    der(2,6) = pt125*xim*zetap
                    der(2,7) = pt125*xip*zetap
                    der(2,8) = pt125*xip*zetam
                    der(3,1) = -pt125*xim*etam
                    der(3,2) = pt125*xim*etam
                    der(3,3) = pt125*xip*etam
                    der(3,4) = -pt125*xip*etam
                    der(3,5) = -pt125*xim*etap
                    der(3,6) = pt125*xim*etap
                    der(3,7) = pt125*xip*etap
                    der(3,8) = -pt125*xip*etap
                case (14)
                    der(1,1) = (two*xi*eta+two*xi*zeta+d4*xi+eta*zeta+eta+zeta)*(eta-one)*(zeta-one)/d8
                    der(1,2) = -(two*xi*eta-two*xi*zeta+d4*xi-eta*zeta+eta-zeta)*(eta-one)*(zeta+one)/d8
                    der(1,3) = -(two*xi*eta-two*xi*zeta+d4*xi+eta*zeta-eta+zeta)*(eta-one)*(zeta+one)/d8
                    der(1,4) = (two*xi*eta+two*xi*zeta+d4*xi-eta*zeta-eta-zeta)*(eta-one)*(zeta-one)/d8
                    der(1,5) = -(eta-one)*(zeta+one)*(zeta-one)*xi
                    der(1,6) = -(eta+one)*(eta-one)*(zeta+one)*(zeta-one)/two
                    der(1,7) =  (eta+one)*(eta-one)*(zeta+one)*xi
                    der(1,8) = (eta+one)*(eta-one)*(zeta+one)*(zeta-one)/two
                    der(1,9) = -(eta+one)*(eta-one)*(zeta-one)*xi
                    der(1,10) = (two*xi*eta-two*xi*zeta-d4*xi+eta*zeta+eta-zeta)*(eta+one)*(zeta-one)/d8
                    der(1,11) = -(two*xi*eta+two*xi*zeta-d4*xi-eta*zeta+eta+zeta)*(eta+one)*(zeta+one)/d8
                    der(1,12) = -(two*xi*eta+two*xi*zeta-d4*xi+eta*zeta-eta-zeta)*(eta+one)*(zeta+one)/d8
                    der(1,13) = (two*xi*eta-two*xi*zeta-d4*xi-eta*zeta-eta+zeta)*(eta+one)*(zeta-one)/d8
                    der(1,14) =  (eta+one)*(zeta+one)*(zeta-one)*xi
                    der(2,1) = (two*xi*eta+xi*zeta+xi+two*eta*zeta+d4*eta+zeta)*(xi-one)*(zeta-one)/d8
                    der(2,2) = -(two*xi*eta-xi*zeta+xi-two*eta*zeta+d4*eta-zeta)*(xi-one)*(zeta+one)/d8
                    der(2,3) = -(two*xi*eta-xi*zeta+xi+two*eta*zeta-d4*eta+zeta)*(xi+one)*(zeta+one)/d8
                    der(2,4) = (two*xi*eta+xi*zeta+xi-two*eta*zeta-d4*eta-zeta)*(xi+one)*(zeta-one)/d8
                    der(2,5) = -(xi+one)*(xi-one)*(zeta+one)*(zeta-one)/two
                    der(2,6) = -(xi-one)*(zeta+one)*(zeta-one)*eta
                    der(2,7) =  (xi+one)*(xi-one)*(zeta+one)*eta
                    der(2,8) =  (xi+one)*(zeta+one)*(zeta-one)*eta
                    der(2,9) = -(xi+one)*(xi-one)*(zeta-one)*eta
                    der(2,10) = (two*xi*eta-xi*zeta-xi+two*eta*zeta+d4*eta-zeta)*(xi-one)*(zeta-one)/d8
                    der(2,11) = -(two*xi*eta+xi*zeta-xi-two*eta*zeta+d4*eta+zeta)*(xi-one)*(zeta+one)/d8
                    der(2,12) = -(two*xi*eta+xi*zeta-xi+two*eta*zeta-d4*eta-zeta)*(xi+one)*(zeta+one)/d8
                    der(2,13) = (two*xi*eta-xi*zeta-xi-two*eta*zeta-d4*eta+zeta)*(xi+one)*(zeta-one)/d8
                    der(2,14) = (xi+one)*(xi-one)*(zeta+one)*(zeta-one)/two
                    der(3,1) = (xi*eta+two*xi*zeta+xi+two*eta*zeta+eta+d4*zeta)*(xi-one)*(eta-one)/d8
                    der(3,2) = -(xi*eta-two*xi*zeta+xi-two*eta*zeta+eta-d4*zeta)*(xi-one)*(eta-one)/d8
                    der(3,3) = -(xi*eta-two*xi*zeta+xi+two*eta*zeta-eta+d4*zeta)*(xi+one)*(eta-one)/d8
                    der(3,4) = (xi*eta+two*xi*zeta+xi-two*eta*zeta-eta-d4*zeta)*(xi+one)*(eta-one)/d8
                    der(3,5) = -(xi+one)*(xi-one)*(eta-one)*zeta  
                    der(3,6) = -(xi-one)*(eta+one)*(eta-one)*zeta  
                    der(3,7) = (xi+one)*(xi-one)*(eta+one)*(eta-one)/two
                    der(3,8) =  (xi+one)*(eta+one)*(eta-one)*zeta
                    der(3,9) = -(xi+one)*(xi-one)*(eta+one)*(eta-one)/two
                    der(3,10) = (xi*eta-two*xi*zeta-xi+two*eta*zeta+eta-d4*zeta)*(xi-one)*(eta+one)/d8
                    der(3,11) = -(xi*eta+two*xi*zeta-xi-two*eta*zeta+eta+d4*zeta)*(xi-one)*(eta+one)/d8
                    der(3,12) = -(xi*eta+two*xi*zeta-xi+two*eta*zeta-eta-d4*zeta)*(xi+one)*(eta+one)/d8
                    der(3,13) = (xi*eta-two*xi*zeta-xi-two*eta*zeta-eta+d4*zeta)*(xi+one)*(eta+one)/d8
                    der(3,14) = (xi+one)*(xi-one)*(eta+one)*zeta
                case (20)
                    xii = (/-1,-1,-1,0,1,1,1,0,-1,-1,1,1,-1,-1,-1,0,1,1,1,0/)
                    etai = (/-1,-1,-1,-1,-1,-1,-1,-1,0,0,0,0,1,1,1,1,1,1,1,1 /)
                    zetai = (/-1,0,1,1,1,0,-1,-1,-1,1,1,-1,-1,0,1,1,1,0,-1,-1/)
                    do l = 1,20
                        xi0 = xi*xii(l)
                        eta0 = eta*etai(l)
                        zeta0 = zeta*zetai(l)
                        if (l.eq.4 .or. l.eq.8 .or. l.eq.16 .or. l.eq.20) then
                            der(1,l) = -pt5*xi*(one+eta0)*(one+zeta0)
                            der(2,l) = pt25*etai(l)*(one-xi*xi)*(one+zeta0)
                            der(3,l) = pt25*zetai(l)*(one-xi*xi)*(one+eta0)
                        else if (l.ge.9 .and. l.le.12) then
                            der(1,l) = pt25*xii(l)*(one-eta*eta)*(one+zeta0)
                            der(2,l) = -pt5*eta*(one+xi0)*(one+zeta0)
                            der(3,l) = pt25*zetai(l)*(one+xi0)*(one-eta*eta)
                        else if (l.eq.2 .or. l.eq.6 .or. l.eq.14 .or. l.eq.18) then
                            der(1,l) = pt25*xii(l)*(one+eta0)*(one-zeta*zeta)
                            der(2,l) = pt25*etai(l)*(one+xi0)*(one-zeta*zeta)
                            der(3,l) = -pt5*zeta*(one+xi0)*(one+eta0)
                        else
                            der(1,l) = pt125*xii(l)*(one+eta0)*(one+zeta0)*(two*xi0+eta0+zeta0-one)
                            der(2,l) = pt125*etai(l)*(one+xi0)*(one+zeta0)*(xi0+two*eta0+zeta0-one)
                            der(3,l) = pt125*zetai(l)*(one+xi0)*(one+eta0)*(xi0+eta0+two*zeta0-one)
                        end if
                    end do
                case default
                    write(*,*) 'wrong number of nodes in shape_der'
                end select
        case default
            write(*,*) 'wrong number of dimensions in shape_der'
    end select
    return
end subroutine shape_der