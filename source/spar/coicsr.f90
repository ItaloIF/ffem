subroutine coicsr(job,n,a,ia,ja)
    !------------------------------------------------------------------------
    ! IN-PLACE coo-csr conversion routine.
    !------------------------------------------------------------------------
    ! this subroutine converts a matrix stored in coordinate format into 
    ! the csr format. The conversion is done in place in that the arrays 
    ! a,ja,ia of the result are overwritten onto the original arrays.
    !------------------------------------------------------------------------
    ! on entry:
    !--------- 
    ! n	= integer. row dimension of A.
    ! nnz	= integer. number of nonzero elements in A.
    ! job   = integer. Job indicator. when job=1, the real values in a are
    !         filled. Otherwise a is not touched and the structure of the
    !         array only (i.e. ja, ia)  is obtained.
    ! a	= real array of size nnz (number of nonzero elements in A)
    !         containing the nonzero elements 
    ! ja	= integer array of length nnz containing the column positions
    ! 	  of the corresponding elements in a.
    ! ia	= integer array of length nnz containing the row positions
    ! 	  of the corresponding elements in a.
    ! iwk	= integer work array of length n+1 
    ! on return:
    !----------
    ! a
    ! ja 
    ! ia	= contains the compressed sparse row data structure for the 
    !         resulting matrix.
    ! Note: 
    !-------
    !         the entries of the output matrix are not sorted (the column
    !         indices in each are not in increasing order) use coocsr
    !         if you want them sorted.
    !----------------------------------------------------------------------!
    !  Coded by Y. Saad, Sep. 26 1989
    !  Modify by II, Jan 31 2021                                          
    !----------------------------------------------------------------------!
    use glb
    implicit none
    integer, intent(in) :: job
    integer, intent(in) :: n
    real(rkd), intent(in out) :: a(:)
    integer, intent(in out) :: ia(:)
    integer, intent(in out) :: ja(:)

    logical val     
    integer :: nnz
    integer :: init
    integer :: iwk(n+1)
    integer :: k
    integer :: i
    integer :: j
    integer :: inext
    integer :: jnext
    integer :: ipos
    real(rkd) :: t
    real(rkd) :: tnext
    nnz = ubound(a,1)
    val = job.eq.1
    ! find pointer array for resulting matrix
    iwk = 0
    do k = 1,nnz
        i = ia(k)
        iwk(i+1) = iwk(i+1) + 1
    end do
    iwk(1) = 1
    do i = 2,n
        iwk(i) = iwk(i-1) + iwk(i)
    end do
    ! loop for a cycle in chasing process
    init = 1
    k = 0
    do
        if (val) t = a(init)
        i = ia(init)
        j = ja(init)
        ia(init) = -1
        do 
            k = k + 1
            ! current row number is i.  determine  where to go
            ipos = iwk(i)
            ! save the chased element
            if (val) tnext = a(ipos)
            inext = ia(ipos)
            jnext = ja(ipos)
            ! then occupy its location
            if (val) a(ipos)  = t
            ja(ipos) = j
            ! update pointer information for next element to come in row i
            iwk(i) = ipos + 1
            ! determine  next element to be chased
            if (ia(ipos).lt.0) exit
            t = tnext
            i = inext
            j = jnext 
            ia(ipos) = -1
            if (k.ge.nnz) exit
        end do
        if (k.ge.nnz) exit
        do 
            init = init + 1
            if (init.gt.nnz) exit
            if (ia(init).ge.0) exit
        end do
        if (init.gt.nnz) exit
    end do
    ia(2:n+1) = iwk(1:n)
    ia(1) = 1
    return
end subroutine coicsr