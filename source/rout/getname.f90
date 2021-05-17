subroutine getname(argv, nlen)
    implicit none
    integer :: narg
    integer, intent(out) :: nlen
    integer :: lnblnk
    integer :: iargc
    character(len=:), allocatable, intent(out) :: argv
    character(len=500) :: buff
    logical :: found

    narg = IARGC()
    if (narg < 1) then
        write(*,*) 'Please enter the base name of data file: '
        read(*,*) buff
        nlen = len(trim(buff))
        allocate(character(nlen) :: argv)
        argv = buff(1:nlen)
    else
        call get_command_argument(number=1, length=nlen)
        allocate(character(nlen) :: argv)
        call get_command_argument(number=1, value=argv)
    end if
    nlen = LNBLNK(argv)
    inquire(file=argv(1:nlen)//'.dat', exist=found)
    if (.NOT.found) then
        write(*,*) 'Data file not found: ',argv(1:nlen)//'.dat'
        write(*,*) 'Please create or check spelling.'
        stop
    end if
    return
end subroutine getname