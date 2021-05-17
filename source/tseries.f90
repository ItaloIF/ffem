module tsmod
    contains
    subroutine tseries(idts,t,fac,f,p1,p2,p3)
        use glb
        implicit none
        integer, intent(in) :: idts
        real(rkd), intent(in) :: t
        real(rkd), intent(in) :: fac
        real(rkd), intent(out) :: f
        real(rkd), intent(in), optional :: p1
        real(rkd), intent(in), optional :: p2
        real(rkd), intent(in), optional :: p3

        real(rkd) :: pi
        real(rkd) :: ts
        real(rkd) :: te
        real(rkd) :: per        

        f = 0.0_rkd
        pi = 4*atan(1.0_rkd)
        select case (idts)
            case (1)
                ! Constant
                f = fac
            case (2)
                ! Linear
                ts = p1
                f = fac*(t - ts)
            case (3)
                ! Sinusoidal
                ts = p1
                te = p2
                per = p3
                if (t.ge.ts.and.t.le.te) f = fac*sin(2*pi*(t-ts)/(per))
            case default
                write(*,*) 'Incorrect time series option'
        end select

    end subroutine tseries

end module tsmod



