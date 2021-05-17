module mater
    interface

    subroutine elela10(prop,dee)
        use glb
        implicit none
        real(rkd), intent(in) :: prop(:)
        real(rkd), intent(out) :: dee(:,:)
    end subroutine elela10

    subroutine elstr10(prop,sti,esi,stf)
        use glb
        implicit none
        real(rkd), intent(in) :: prop(:)
        real(rkd), intent(in) :: sti(:)
        real(rkd), intent(in) :: esi(:)
        real(rkd), intent(out) :: stf(:)
    end subroutine elstr10

    end interface

    contains
    subroutine models(job,model,prop,osv,sti,dee,stf,esi,est)
        use glb
        implicit none
        integer, intent(in) :: job
        integer, intent(in) :: model
        real(rkd), intent(in) :: prop(:)
        real(rkd), intent(in out) :: osv(:)
        real(rkd), intent(in) :: sti(:)
        real(rkd), intent(in out) :: dee(:,:)
        real(rkd), intent(in out), optional :: stf(:)
        real(rkd), intent(in out), optional :: esi(:)
        real(rkd), intent(in out), optional :: est(:)

        ! (1) Elastic matrix
        ! (2) Initial conditions
        ! (3) Evaluate stresses
        ! (5) Elastoplastic matrix
        select case (model)
            case (1)
                ! Elastic model
                select case (job)
                    case (1,5)
                        call elela10(prop,dee)
                    case (3)
                        call elstr10(prop,sti,esi,stf)
                    case default
                        write(*,*) 'wrong job in mater model'
                        stop
                end select
            ! case (2)
            !     ! General model Hinton & Owen
            !     select case (job)
            !         case (1)
            !             call geela10(prop,dee)
            !         case (3)
            !             call gestr10(prop,tau,deps,dee,st)
            !         case (5)
            !             call geelp10(prop,dee)
            !         case default
            !             write(*,*) 'wrong job in mater model'
            !             stop
            !     end select
            case default
                write(*,*) 'wrong model id'
                stop
        end select
    end subroutine models

end module mater

include 'mater/elela10.f90'
include 'mater/elstr10.f90'