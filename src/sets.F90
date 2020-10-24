
module sets

  public :: intersection, union

contains


  !===================================================
  subroutine sort_lists(lists)
    !-------------------------------------------------
    use healpix_types, only: I4B
    use healpix_modules, only: assert_alloc, iindexx
    implicit none
    integer(i4b), intent(inout), dimension(0:,0:) :: lists
    !
    integer(i4b) :: i,j,nd, nl, status
    integer(i4b), dimension(:), allocatable :: idx, buffer
    character(len=*), parameter :: code = 'sort_lists'
    !-------------------------------------------------

    nd = size(lists,1)
    nl = size(lists,2)

    allocate(idx(1:nd), stat=status)
    call assert_alloc(status, code, 'idx')
    allocate(buffer(1:nd), stat=status)
    call assert_alloc(status, code, 'buffer')

    call iindexx(nd, lists(0:nd-1,0), idx) ! idx is 1-based
    do j=0,nl-1
       do i=1,nd
          buffer(i) = lists(idx(i)-1, j)
       enddo
       lists(0:nd-1, j) = buffer(1:nd)
    enddo
    deallocate(idx)
    deallocate(buffer)
    return
  end subroutine sort_lists

  !===================================================
  function binary_search(list, value, istart) result (ipos)
    !-------------------------------------------------
    ! if list is a (ascending) sorted listed
    ! ipos = binary_search(list, value, istart)
    ! is such that list(ipos) = value
    ! ipos=-1 in case of failure, ipos is 0-based,
    ! istart: lowest value of ipos allowed (0 is recommended)
    !
    use healpix_types, only: I4B
    implicit none
    integer(i4b), dimension(0:), intent(in) :: list
    integer(i4b), intent(in) :: value, istart
    integer(i4b) :: ipos
    integer(i4b) :: ilow, i, ihi, n
    !-------------------------------------------------

    n = size(list)
    ilow = istart
    ihi  = n - 1
    ipos = -1 ! failure
    do 
       if (ilow == ihi) exit
       i = (ilow+ihi+1)/2  ! ceil( (low+hi)/2 )
       if (list(i) > value) then
          ihi = i - 1
       else
          ilow = i
       end if
    enddo
    if (list(ilow) == value) ipos = ilow
    return
    
  end function binary_search

  !===================================================
  function union(list1, list2, sort) result (lout)
    !-------------------------------------------------
    ! listany => union(list1, list2 [,sort=])
    ! for 2 lists of integer numbers (the 2nd must be sorted)
    ! returns a threefold list
    !  - the union list numbers
    !  - their (0-based) location in the 1st list, with -1 for those coming from the second list
    !  - their (0-based) location in the 2nd list, with -1 for those coming from the 1st list
    ! if sort is set to .true., the output list is sorted 
    ! according to its 1st field
    !------------------------------------------------
    use healpix_types, only: I4B, LGT
    use healpix_modules, only: iindexx, assert_alloc
    implicit none
    
    integer(i4b), dimension(0:),  intent(inout) :: list1,list2
    logical(LGT), optional,       intent(in)    :: sort
    integer(i4b), dimension(:,:), pointer    :: lout
    !
    integer(i4b), dimension(:,:),   pointer :: linter
    integer(i4b) :: n1,n2,ni,nu, i1,i2,i,j, status, badval
    character(len=*), parameter :: code='union'
    logical(LGT) :: sort_out

    sort_out = .false.
    if (present(sort)) sort_out = sort

    linter => intersection(list1, list2)
    n1 = size(list1)
    n2 = size(list2)
    ni = size(linter,1)
    nu = n1 + n2 - ni
    allocate(lout(0:nu-1,0:2), stat=status)
    call assert_alloc(status, code, 'lout')


    ! fill union list with intersection
    do i=0,ni-1
       lout(i, 0) = linter(i, 0)
       lout(i, 1) = linter(i, 1)
       lout(i, 2) = linter(i, 2)
    enddo

    if (nu > ni) then 
       ! edit temporarily lists to flag out elements in intersection
       badval = min(0, minval(list1), minval(list2)) - 1
       do i=0,ni-1
          list1(linter(i,1)) = badval
          list2(linter(i,2)) = badval
       enddo
       ! fill union list with part of list1 not in intersection
       i = ni
       do i1=0,n1-1
          if (list1(i1) > badval) then
             lout(i,0) = list1(i1)
             lout(i,1) = i1
             lout(i,2) = -1
             i = i+1
          endif
       enddo
       ! fill union list with part of list2 not in intersection
       do i2=0,n2-1
          if (list2(i2) > badval) then
             lout(i,0) = list2(i2)
             lout(i,1) = -1
             lout(i,2) = i2
             i = i+1
          endif
       enddo
       ! restore lists
       do i=0,ni-1
          list1(linter(i,1)) = linter(i,0)
          list2(linter(i,2)) = linter(i,0)
       enddo

    endif
    deallocate(linter)
    ! sort union list
    if (sort_out) then
       call sort_lists(lout)
    endif

    return
  end function union
  !===================================================
  function intersection(list1, list2, sort) result (lout)
    !-------------------------------------------------
    ! listcom => intersection(list1, list2 [,sort=])
    ! for 2 lists of integer numbers (the 2nd must be sorted)
    ! returns a threefold list containing
    !  - the list of common numbers
    !  - their (0-based) location in the 1st list
    !  - their (0-based) location in the 2nd list
    ! if sort is set to .true., the output list is sorted 
    ! according to its 1st field
    !------------------------------------------------
    use healpix_types, only: I4B, LGT
    use healpix_modules, only: fatal_error, assert_alloc
    implicit none
    integer(i4b), dimension(0:),  intent(in) :: list1, list2
    logical(LGT), optional,       intent(in)    :: sort
    integer(i4b), dimension(:,:), pointer :: lout
    !
    integer(i4b), dimension(:), allocatable :: diff
    integer(i4b) :: n1,n2,n,i1,i2,i2t,k, min2, max2, status
    logical(LGT) :: done, sort_out
    character(len=*), parameter :: code='intersection'

    done     = .false.
    sort_out = .false.
    if (present(sort)) sort_out = sort
    
    n1 = size(list1)
    n2 = size(list2)
    ! use short-cut if lists are identical
    if (            n1 == n2        .and. &
         & list1(0)    == list2(0)  .and. &
         & list1(n1-1) == list2(n2-1)) then
       allocate(diff(0:n1-1),stat=status)
       call assert_alloc(status, code, 'diff')
       diff = abs(list1-list2)
       if (maxval(diff) == 0) then
          allocate(lout(0:n1-1,0:2),stat=status)
          call assert_alloc(status, code, 'lout')
          do i1=0,n1-1
             lout(i1,0) = list1(i1)
             lout(i1,1) = i1
             lout(i1,2) = i1
          enddo
          done=.true.
       endif
       deallocate(diff)
    endif
    if (done) then
       if (sort_out) call sort_lists(lout)
       return
    endif

    ! quick-check that 2nd list is sorted
    min2 = minval(list2)
    max2 = maxval(list2)
    if (min2 /= list2(0) .or. max2 /= list2(n2-1)) then
       print*,min2,max2
       print*,list2(0),list2(n2-1)
       call fatal_error('List2 is not sorted in '//code)
    endif

    ! general case
    do k=0,1 ! k=0: count matches, k=1: fill arrays
       if (k==1) then
          allocate(lout(0:n-1,0:2),stat=status)
          call assert_alloc(status, code, 'lout')
       endif
       n  = 0
       i2 = 0
       do i1=0, n1-1
          i2t = binary_search(list2, list1(i1), i2)
          if (i2t >= 0) then
             if (k == 1) then
                lout(n,0) = list1(i1)!common values
                lout(n,1) = i1       !0 based location in 1st list
                lout(n,2) = i2t      !0 based location in 2nd list
             endif
             n = n + 1
             !i2 = i2t
          endif
       enddo
    enddo

    if (sort_out) call sort_lists(lout)

    return 
  end function intersection
  
end module sets


