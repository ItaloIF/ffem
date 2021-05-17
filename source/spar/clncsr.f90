subroutine clncsr(job,value,n,a,ia,ja,indu)
    ! this routine performs two tasks to clean up a CSR matrix
    ! -- remove duplicate/zero entries,
    ! -- perform a partial ordering, new order lower triangular part,
    !    main diagonal, upper triangular part.
    ! on entry :
    !   job :
    !       0 : nothing is done
    !       1 : eliminate duplicate entries, zero entries
    !       2 : eliminate duplicate entries and perform partial ordering
    !       3 : eliminate duplicate entries, sort the entries in the
    !           increasing order of column indices
    !   value :
    !       0 : 
    !       1 :
    !   n : row dimension of the matrix
    !   a,ja,ia : input matrix in CSR format
    ! on return :
    !   a,ja,ia : cleaned matrix
    !----------------------------------------------------------------------!
    !  Modify by II, Feb 01 2021                                          
    !----------------------------------------------------------------------!

    use glb
    implicit none
    integer, intent(in) :: job
    integer, intent(in) :: value
    integer, intent(in) :: n
    real(rkd), intent(in out) :: a(:)
    integer, intent(in out) :: ia(:)
    integer, intent(in out) :: ja(:)
    integer, intent(out) :: indu(n)

    integer :: i
    integer :: j
    integer :: k
    integer :: ipos
    integer :: kfirst
    integer :: klast
    integer :: iwk(n+1)
    real(rkd) :: tmp
    integer :: ko
    if (job.le.0) return
    indu = 0
    iwk = ia(1:n+1)
    k = 1
    do i = 1,n
        ia(i) = k
        ipos = iwk(i)
        klast = iwk(i+1)
        do
            if (ipos.ge.klast) exit
            j = ja(ipos)
            if (indu(j).eq.0) then
                ! new entry
                if (value.ne.0) then
                    if (a(ipos).ne.0.0D0) then
                        indu(j) = k
                        ja(k) = ja(ipos)
                        a(k) = a(ipos)
                        k = k + 1
                     endif
                else
                    indu(j) = k
                    ja(k) = ja(ipos)
                    k = k + 1
                end if
            else if (value.ne.0) then
                ! duplicate entry
                a(indu(j)) = a(indu(j)) + a(ipos)
            end if
            ipos = ipos + 1
        end do
        ! remove marks before working on the next row
        do ipos = ia(i),k-1
            indu(ja(ipos)) = 0
        end do
    end do
    ia(n+1) = k
    if (job.le.1) return

    ! partial ordering
    ! split the matrix into strict upper/lower triangular
    ! parts, INDU points to the the beginning of the upper part
    do i = 1,n
        klast = ia(i+1) - 1
        kfirst = ia(i)
        do
            if (klast.le.kfirst) exit
            if (ja(klast).lt.i.and.ja(kfirst).ge.i) then
                ! swap klast with kfirst
                j = ja(klast)
                ja(klast) = ja(kfirst)
                ja(kfirst) = j
                if (value.ne.0) then
                    tmp = a(klast)
                    a(klast) = a(kfirst)
                    a(kfirst) = tmp
                end if 
            end if
            if (ja(klast).ge.i) klast = klast - 1
            if (ja(kfirst).lt.i) kfirst = kfirst + 1
        end do
        if (ja(klast).lt.i) then
            indu(i) = klast + 1
         else
            indu(i) = klast
         endif
    end do
    if (job.le.2) return

    ! order the entries according to column indices
    ! burble-sort is used
    do i = 1,n
        do ipos = ia(i),indu(i)-1
            do j = indu(i)-1,ipos+1,-1
                k = j - 1
                if (ja(k).gt.ja(j)) then
                    ko = ja(k)
                    ja(k) = ja(j)
                    ja(j) = ko
                    if (value.ne.0) then
                        tmp = a(k)
                        a(k) = a(j)
                        a(j) = tmp
                    end if
                end if
            end do
        end do
        do ipos = indu(i),ia(i+1)-1
            do j = ia(i+1)-1,ipos+1,-1
                k = j - 1
                if (ja(k).gt.ja(j)) then
                    ko = ja(k)
                    ja(k) = ja(j)
                    ja(j) = ko
                    if (value.ne.0) then
                       tmp = a(k)
                       a(k) = a(j)
                       a(j) = tmp
                    end if
                 end if
            end do
        end do
    end do
    return
end subroutine clncsr