program femg
    use glb
    type pts
        integer :: id
        real(rkd) :: coor(3)
        integer :: ndof
        integer, allocatable :: cons(:)
    end type pts
    implicit none

! Variables
    ! Scalar integers
    type(pts) :: 
    integer :: nn


! Input Data
    call getname(fdat,nlen)
    open(1, file=fdat(1:nlen)//'.dat')
    ! read points
    call read_points(1, nn, pts)
    ! read constrains
    call read


end program