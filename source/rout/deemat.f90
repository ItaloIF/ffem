subroutine deemat(dee,e,v)
    ! this subroutine returns the elastic dee matrix for ih=3 (plane strain), ih=4 (axisymmetry or plane strain elastoplasticity) or ih=6 (three dimensions).
    use glb
    implicit none
    real(rkd), intent(in) :: e
    real(rkd), intent(in) :: v
    real(rkd), intent(out) :: dee(:,:)
    real(rkd) :: v1
    real(rkd) :: v2
    real(rkd) :: c
    real(rkd) :: vv
    real(rkd) :: zero = 0.0_rkd
    real(rkd) :: pt5 = 0.5_rkd
    real(rkd) :: one = 1.0_rkd
    real(rkd) :: two = 2.0_rkd
    integer :: i
    integer :: ih
    dee = zero
    ih = ubound(dee,1)
    v1 = one - v
    c = e/((one+v)*(one-two*v))
    select case(ih)
        case (3)
            dee(1,1) = v1*c
            dee(2,2) = v1*c
            dee(1,2) = v*c
            dee(2,1) = v*c
            dee(3,3) = pt5*c*(one-two*v)
        case (4)
            dee(1,1) = v1*c
            dee(2,2) = v1*c
            dee(4,4) = v1*c
            dee(3,3) = pt5*c*(one-two*v) 
            dee(1,2) = v*c
            dee(2,1) = v*c
            dee(1,4) = v*c
            dee(4,1) = v*c
            dee(2,4) = v*c
            dee(4,2) = v*c
        case (6)
            v2 = v/(one-v)
            vv = (one-two*v)/(one-v)*pt5
            do i = 1,3
                dee(i,i) = one
            end do
            do i = 4,6
                dee(i,i) = vv
            end do
            dee(1,2) = v2
            dee(2,1) = v2
            dee(1,3) = v2
            dee(3,1) = v2
            dee(2,3) = v2
            dee(3,2) = v2
            dee = dee*e/(two*(one+v)*vv)
        case default
            write(*,*) 'wrong size for dee matrix'
            stop
        end select
        return    
end subroutine deemat