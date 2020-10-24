module kernels_tools
  use healpix_types
  implicit none
contains


  !===============================================================================================
  subroutine wl2kernels(wl, kernels)
    !===============================================================================================
    ! This routine computes the window functions of the decoupled estimator
    ! as specified in Chon et al. Eq 62
    ! faster than previous Spice's compute windows
    !       - use kernel symmetry to reduce number of calculations
    !       - OPenMP: use dynamic rather than static scheduling of loops
    ! ---------------------------------------------------------------------------------------------
    use apodize_mod,   only: apodizefunction
    use misc_utils,    only: assert_alloc, wall_clock_time
    use rec3jj_module, only: rec3jj
    !-------------------
    real(DP), dimension(0:,0:,1:), intent(out) :: kernels
    real(DP), dimension(0:),       intent(in)  :: wl

    real(DP), parameter :: zero =  0.00000000000_dp
    real(DP), parameter :: one  =  1.00000000000_dp
    real(DP), parameter :: two  =  2.00000000000_dp
    real(DP), parameter :: mtwo = -2.00000000000_dp
    real(DP), parameter :: eps = 1.e-12_dp
    integer                                   :: l1max, l2max, l3max, ncl, nkernels
    real(DP),     dimension(:)  , allocatable :: wigner00, wigner22, work
    real(DP),     dimension(:,:), allocatable :: clw
    real(DP)     :: dl1, dl2, dl3min, dl3max, tmp00, tmp02, tmp22e, tmp22o
    integer(i4b) :: l1, l2, l, l3eff, ierror1, ierror2, l3e, l3o, off_even, off_odd, j, j1
    integer(i4b) :: nl3min, nl3max, l2start, l2end, lmin, i
    integer(i4b) :: status
    integer(I4B) :: incl,inclpercent,nlog
    character(len=*), parameter :: code = 'wl2kernels'
    real(SP) :: wt0, wt1, wtt, wttold

    real(DP), allocatable, dimension(:) :: fl
    real(DP) :: kcross
    logical(LGT) :: fullkernels, megaverbose=.true.
    !===============================================================================================
    call wall_clock_time(wt0)
    wttold = wt0


    l1max    = size(kernels,1)-1 ! first dim of kernel
    l2max    = size(kernels,2)-1 ! 2nd dim of kernel
    nkernels = size(kernels,3)   ! 3rd dim of kernel
    l3max    = size(wl,     1)-1 ! dimension of apodization window
    l3eff = l1max + l2max ! largest ell dealt with

    allocate(fl(0:l3max), stat=status)
    call assert_alloc(status, code, 'fl vector')
    ! f_L <- (2L+1)w_L
    do l=0,l3max
       fl(l) = wl(l)*( TWO * l + ONE) / FOURPI
    enddo

    ! fl(l) = delta(l) then kernel = Identity
    if (abs(fl(0) - ONE) < EPS .and. maxval(abs(fl(1:l3max))) <  EPS) then
       kernels = ZERO
       do l=0, min(l1max, l2max)
          kernels(l,l,1:nkernels) = ONE
       enddo
       if (megaverbose) print*,'shortcut kernel calculations'
       return
    endif

    ! matrix is symmetric but can be rectangular
    ! compute upper half if #columns  > #rows
    ! compute lower half if #columns  < #rows
    ! start loops

    nlog = max(nint( 2.*(l1max/600.)**3), 1 ) ! number of logs propto number of calculations
    nlog = min(nlog, 25) ! no more than 25 logs
    inclpercent=max((l1max+1)/nlog,1)
    incl=0
!$OMP PARALLEL DEFAULT(NONE) &
!$OMP SHARED(l1max, l2max, l3max, l3eff, fl, kernels, kcross, fullkernels, incl, inclpercent, megaverbose, wtt, wttold) &
!$OMP PRIVATE(l1, dl1, l2, dl2, l2start, l2end, l3e, l3o, dl3min, dl3max, nl3min, nl3max, wigner00, wigner22, work, &
!$OMP         ierror1, ierror2, tmp00, tmp22o, tmp22e, tmp02, off_even, off_odd, status)

    allocate(wigner00(0:l3eff), stat=status)
    call assert_alloc(status,code,'wigner00')

    allocate(wigner22(0:l3eff), stat=status)
    call assert_alloc(status,code,'wigner22')

    allocate(work(1:l3eff+1), stat=status)
    call assert_alloc(status,code,'work')

    ierror1 = 0
    ierror2 = 0
!$OMP DO schedule(dynamic, 2)
!     do l1 = l1max, 0, -1
    do l1 = 0, l1max
       dl1 = real(l1, kind=DP)
       if (l1max >= l2max) then
          l2start = 0
          l2end = min(l1, l2max)
       else
          l2start = l1
          l2end = l2max
       endif

       do l2 = l2start, l2end
          dl2 = real(l2, kind=DP)
          
          if (abs(l1-l2) <= l3max) then
             work(1:l1+l2+1)     = ZERO
             call rec3jj(work, dl1, dl2, ZERO, ZERO, dl3min, dl3max, l3eff+1, ierror1)
             nl3min = nint(dl3min)
             nl3max = min(nint(dl3max), l3max) ! truncate to ell's where window is non-zero
             wigner00(nl3min:nl3max) = work(1:nl3max-nl3min+1)
             
             work(1:l1+l2+1)     = ZERO
             if (l1 > 1 .and. l2 > 1) then
                call rec3jj(work, dl1, dl2, TWO,  MTWO, dl3min, dl3max, l3eff+1, ierror2)
             endif
             wigner22(nl3min:nl3max) = work(1:nl3max-nl3min+1)

             tmp00  = ZERO
             tmp22o = ZERO
             tmp22e = ZERO
             tmp02  = ZERO
             
             off_even = 1
             if ( mod(l1+l2,2) == mod(nl3min,2) ) off_even = 0
             off_odd = 1 - off_even
             
             do l3e = nl3min + off_even, nl3max, 2 ! all l2 such that l1+l2+l3 is even
                tmp00  = tmp00 +  wigner00(l3e)*wigner00(l3e) * fl(l3e) ! TT
                tmp22e = tmp22e + wigner22(l3e)*wigner22(l3e) * fl(l3e) ! EE or BB
                tmp02  = tmp02 +  wigner00(l3e)*wigner22(l3e) * fl(l3e) ! TE or TB
             enddo
             do l3o = nl3min + off_odd, nl3max, 2  ! all l3 such that l1+l2+l3 is odd
                tmp22o = tmp22o + wigner22(l3o)*wigner22(l3o) * fl(l3o) ! EE or BB
             enddo
             kernels(l1, l2, 1) = kernels(l1, l2, 1) + tmp00   ! TT
             kernels(l1, l2, 2) = kernels(l1, l2, 2) + tmp22e + tmp22o ! EE, BB, EB
             kernels(l1, l2, 3) = kernels(l1, l2, 3) + tmp22e - tmp22o ! EE -> BB, BB -> EE
             kernels(l1, l2, 4) = kernels(l1, l2, 4) + tmp02   ! TE, TB
          endif

       enddo ! loop on l2
!$OMP ATOMIC
      incl=incl+1
      if (megaverbose .and. mod(incl,inclpercent) == 0) then
         call wall_clock_time(wtt)
         write(*,*) '% kernel done =', &
               &           nint( incl * 100.0_dp/real(l1max+1,kind=DP)), wtt-wttold
         wttold = wtt
      endif
    enddo ! loop on l1
!$OMP END DO
    deallocate(work)
    deallocate(wigner00, wigner22)
!$OMP END PARALLEL

    ! fill other half of matrix and apply 2*l + 1 factor in one direction
    lmin = min(l1max, l2max)
!$OMP PARALLEL DEFAULT(NONE) &
!$OMP SHARED(l1max, l2max, lmin, kernels, kcross, fullkernels) &
!$OMP PRIVATE(l1, l2, l2start, l2end)
!$OMP DO schedule(dynamic, 8)
    do l1 = 0, lmin
       if (l1max >= l2max) then
          l2start = 0
          l2end = min(l1-1, lmin)
       else
          l2start = l1+1
          l2end = lmin
       endif
       do l2 = l2start, l2end
          kernels(l2, l1, 1:4) = kernels(l1, l2, 1:4)
       enddo
    enddo
!$OMP END DO
!$OMP END PARALLEL

!$OMP PARALLEL DEFAULT(NONE) SHARED(l2max, kernels) PRIVATE(l2)
!$OMP DO schedule(dynamic, 8)
    do l2 = 0, l2max
       kernels(:,l2,:) = kernels(:,l2,:) * (TWO * l2 + ONE)
    enddo
!$OMP END DO
!$OMP END PARALLEL

    deallocate(fl)
    call wall_clock_time(wt1)
    if (megaverbose) then
       print*,'Kernel calculation (wall time) = ',wt1-wt0
       print*,l1max,l2max,nkernels
    endif

  end subroutine wl2kernels
end module kernels_tools

!============================================================================
program test
  use healpix_types
  use healpix_modules
  use kernels_tools
  implicit none
  integer(i4b), parameter :: nlmax = 3000
  integer(i4b) :: i, status
  real(dp), dimension(:,:,:), allocatable :: kernels
  real(dp), dimension(0:2*nlmax) :: wl


  allocate(kernels(0:nlmax, 0:nlmax, 1:4), stat=status)
  call assert_alloc(status, 'main', 'kernels')

  do i=0, nlmax
     wl(i) = exp(-(i*0.1d0)**2)
  enddo
  call wl2kernels(wl, kernels)
!   do i=1,4
!      print*,kernels(0:4,0:4,i)
!   enddo

end program test
