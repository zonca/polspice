module windows
  implicit none
contains


  !===============================================================================================
  subroutine compute_kernels
    !===============================================================================================
    ! This routine computes the window functions of the decoupled estimator
    ! as specified in Chon et al. Eq 62
    ! faster than previous Spice's compute windows
    !       - use kernel symmetry to reduce number of calculations
    !       - OPenMP: use dynamic rather than static scheduling of loops
    ! ---------------------------------------------------------------------------------------------
    use spice_common, only: sp, dp, i4b, PI, TWOPI, FOURPI, &
         &    nlmax, ncor, &
         &    fullkernels, kernels, nkernels, kcross, TEnorm, tenormoutput, &
         &    verbose, megaverbose, apodizetype, apodizesigma, thetamax, &
         &    polarization
    use apodize_mod,   only: apodizefunction
    use misc_utils,    only: assert_alloc, wall_clock_time
    use do_legendre_mod, only: do_legendre
#ifdef OLDREC3JJ
    use rec3jj_module, only: rec3jj
#else
    use rec3jjcmb_module, only: w3j_220, init_w3j
#endif
!     real(DP),     intent(in), dimension(0:l3max, 1:ncl)         :: cl_mask
!     real(DP),     intent(out), dimension(0:l1max, 0:l2max, 1:4) :: kernels
    !-------------------
    real(DP), parameter :: zero =  0.00000000000_dp
    real(DP), parameter :: one  =  1.00000000000_dp
    real(DP), parameter :: two  =  2.00000000000_dp
    real(DP), parameter :: mtwo = -2.00000000000_dp
    real(DP), parameter :: eps = 1.e-12_dp
    integer                                   :: l1max, l2max, l3max
    real(DP),     dimension(:),     allocatable :: wigner00, wigner22, work
    real(DP),     dimension(:,:)  , allocatable :: work2
    !real(DP),     dimension(:,:), allocatable :: clw
    real(DP)     :: dl1, dl2, dl3min, dl3max, tmp00, tmp02, tmp22e, tmp22o
    integer(i4b) :: l1, l2, l, l3eff, ierror1, ierror2, l3e, l3o, off_even, off_odd
    integer(i4b) :: nl3min, nl3max, l2start, l2end, lmin, i
    integer(i4b) :: status
    integer(I4B) :: incl,inclpercent,nlog
    character(len=*), parameter :: code = 'cl2kernels'
    real(SP) :: ct0,ct1,wt0, wt1, wti
    real(DP), allocatable, dimension(:)     :: mu_,w_
    real(DP), allocatable, dimension(:,:,:) :: Pl_
    real(DP) :: fcompletion

    real(DP) :: mumin, wxi
    real(DP), allocatable, dimension(:) :: apfunction, fl
    !===============================================================================================
    call cpu_time(ct0)
    call wall_clock_time(wt0)


    l1max =   nlmax ! first dim of kernel
    l2max = 2*nlmax ! 2nd dim of kernel
    l3max =   nlmax ! dimension of apodization window
    l3eff = l1max + l2max ! largest ell dealt with

    ! generate Legendre polynomials (Pl_) and Gauss-Legendre nodes (mu_) and weights (w_)
    ! should be done before calling apodizefunction
    allocate(mu_(1:nlmax+1), stat=status)
    call assert_alloc(status, code, 'mu_')
    allocate(w_(1:nlmax+1), stat=status)
    call assert_alloc(status, code, 'w_')
    ! allocate(Pl_(0:nlmax,0:nlmax,0:ncor-1), stat=status) ! non-local variables (spice_common)
    allocate(Pl_(0:nlmax,0:nlmax,0:0), stat=status) ! 2019-12-05
    call assert_alloc(status, code, 'Pl_')
    mumin = cos(thetamax*PI/180.0_DP)
!    call do_legendre(nside, nlmax, mumin, mu_, w_, Pl_, megaverbose, ncor) ! 2010-03-12
!    call do_legendre(nlmax, mumin, mu_, w_, Pl_, megaverbose, ncor) ! 2019-12-02
    call do_legendre(nlmax, mumin, mu_, w_, Pl_, megaverbose, 1) 

    ! Compute the apodization function in real space
    allocate(apfunction(0:nlmax),stat=status)
    apfunction = ZERO
    do i=0,nlmax
       call apodizefunction(mu_(i+1), apodizesigma, thetamax, apfunction(i), type=apodizetype)
    enddo

    ! Get the Legendre transform of the apodization function
    allocate(fl(0:nlmax), stat=status)
    fl = ZERO
    do i=0,nlmax
       wxi = w_(i+1)*apfunction(i)*TWOPI
       do l=0,nlmax
          fl(l) = fl(l) + Pl_(i,l,0)*wxi
       enddo
    enddo

    ! f_L <- (2L+1)f_L
    do l=0,nlmax 
       fl(l) = fl(l)*( TWO * l + ONE) / FOURPI
    enddo

    deallocate(apfunction)
    deallocate(mu_, w_, Pl_)

    ! fullkernels, kernels, kcross and nkernels are defined in spice_common
    if (fullkernels) then
       if (polarization) then
          nkernels = 4
          if (megaverbose) write(*,*) 'Compute full kernel'
       else
          nkernels = 1
          if (megaverbose) write(*,*) 'Compute TT kernel'
       endif
       allocate(kernels(0:nlmax, 0:2*nlmax, 1:nkernels), stat=status)
       call assert_alloc(status, code, 'can not allocate kernels')
       kernels = ZERO
    else
       if (megaverbose .or. verbose) then
          write(*,*) 'Compute TE kernel for apodization'
          write(*,*) '(to skip this calculation, see the options'
          write(*,*) ' -tenormfileout and -tenormfilein'
          write(*,*) ' in code documentation)'
       endif
       allocate(kcross(0:nlmax, 0:2*nlmax), stat=status)
       call assert_alloc(status, code, 'can not allocate kcross')
       kcross = ZERO
       nkernels = 1
    endif

    ! fl(l) = delta(l) then kernel = Identity
    if (abs(fl(0) - ONE) < EPS .and. maxval(abs(fl(1:nlmax))) <  EPS) then
       if (fullkernels) then
          do l=0, nlmax
             kernels(l,l,1:nkernels) = ONE
          enddo
       else
          do l=0, nlmax
             kcross(l,l) = ONE
          enddo
       endif
       if (megaverbose) print*,'shortcut kernel calculations'
!        deallocate(apfunction, fl)
!        deallocate(mu_, w_, Pl_)
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
!$OMP SHARED(l1max, l2max, l3max, l3eff, fl, kernels, kcross, fullkernels, incl, inclpercent, megaverbose, nkernels, wt0) &
!$OMP PRIVATE(l1, dl1, l2, dl2, l2start, l2end, l3e, l3o, dl3min, dl3max, nl3min, nl3max, wigner00, wigner22, work, work2, &
!$OMP         ierror1, ierror2, tmp00, tmp22o, tmp22e, tmp02, off_even, off_odd, status, wti, fcompletion)

    allocate(wigner00(0:l3eff), stat=status)
    call assert_alloc(status,code,'wigner00')

    allocate(wigner22(0:l3eff), stat=status)
    call assert_alloc(status,code,'wigner22')

#ifdef OLDREC3JJ
    allocate(work(1:l3eff+1), stat=status)
    call assert_alloc(status,code,'work')
#else
    call init_w3j(l3max)
    !call init_w3j()
    allocate(work2(0:l3eff,0:2), stat=status)
    call assert_alloc(status,code,'work2')
#endif
    ierror1 = 0
    ierror2 = 0
!$OMP DO schedule(dynamic, 2)
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
#ifdef OLDREC3JJ
             work(1:2*min(l1,l2)+1)     = ZERO
             call rec3jj(work, dl1, dl2, ZERO, ZERO, dl3min, dl3max, l3eff+1, ierror1)
             nl3min = nint(dl3min)
             nl3max = min(nint(dl3max), l3max) ! truncate to ell's where window is non-zero
             wigner00(nl3min:nl3max) = work(1:nl3max-nl3min+1)
             
             work(1:2*min(l1,l2)+1)     = ZERO
             if (l1 > 1 .and. l2 > 1) then
                call rec3jj(work, dl1, dl2, TWO,  MTWO, dl3min, dl3max, l3eff+1, ierror2)
             endif
             wigner22(nl3min:nl3max) = work(1:nl3max-nl3min+1)
#else
             call w3j_220(work2, l1, l2, nl3min, nl3max)
             nl3max = min(nl3max, l3max) ! truncate to ell's where window is non-zero
             wigner00(nl3min:nl3max) = work2(nl3min:nl3max, 0)
             wigner22(nl3min:nl3max) = work2(nl3min:nl3max, 2) 
             ! returns (l1,l2,l3;-2,2,0) instead of expected (l1,l2,l3;2,-2,0), 
             ! differing on the sign of odd l1+l2+l3.
             ! Does not matter since below those terms are either squarred or dropped
#endif

             tmp00  = ZERO
             tmp22o = ZERO
             tmp22e = ZERO
             tmp02  = ZERO
             
             off_even = 1
             if ( mod(l1+l2,2) == mod(nl3min,2) ) off_even = 0
             off_odd = 1 - off_even
             
             if (fullkernels) then
                do l3e = nl3min + off_even, nl3max, 2 ! all l3 such that l1+l2+l3 is even
                   tmp00  = tmp00 +  wigner00(l3e)*wigner00(l3e) * fl(l3e) ! TT
                   tmp22e = tmp22e + wigner22(l3e)*wigner22(l3e) * fl(l3e) ! EE or BB
                   tmp02  = tmp02 +  wigner00(l3e)*wigner22(l3e) * fl(l3e) ! TE or TB
                enddo
                do l3o = nl3min + off_odd, nl3max, 2  ! all l3 such that l1+l2+l3 is odd
                   tmp22o = tmp22o + wigner22(l3o)*wigner22(l3o) * fl(l3o) ! EE or BB
                enddo
                kernels(l1, l2, 1) = kernels(l1, l2, 1) + tmp00   ! TT
                if (nkernels == 4) then
                   kernels(l1, l2, 2) = kernels(l1, l2, 2) + tmp22e + tmp22o ! EE, BB, EB
                   kernels(l1, l2, 3) = kernels(l1, l2, 3) + tmp22e - tmp22o ! EE -> BB, BB -> EE
                   kernels(l1, l2, 4) = kernels(l1, l2, 4) + tmp02   ! TE, TB
                endif
             else
                do l3e = nl3min + off_even, nl3max, 2 ! all l3 such that l1+l2+l3 is even
                   tmp02  = tmp02 +  wigner00(l3e)*wigner22(l3e) * fl(l3e) ! TE or TB
                enddo
                kcross(l1, l2) = kcross(l1, l2) + tmp02
             endif
          endif

       enddo ! loop on l2
!$OMP ATOMIC
      incl=incl+1
      if (megaverbose .and. mod(incl,inclpercent) == 0) then
         call wall_clock_time(wti)
#ifdef OLDREC3JJ
         fcompletion = (incl /real(l1max+1,kind=DP))**2
#else
         fcompletion = (incl /real(l1max+1,kind=DP))
#endif
         if (fcompletion > 1.e-3) write(*,'(a,f6.1,a,f8.1,a)') 'kernel: ',100.0_dp*fcompletion,'% in ',wti-wt0,'s'
!          write(*,*) '% kernel done =', &
!                &           nint( incl * 100.0_dp/real(l1max+1,kind=DP))
      endif
    enddo ! loop on l1
!$OMP END DO
    if (allocated(work))  deallocate(work)
    if (allocated(work2)) deallocate(work2)
    deallocate(wigner00, wigner22)
!$OMP END PARALLEL

    ! fill other half of matrix and apply 2*l + 1 factor in one direction
    lmin = min(l1max, l2max)
!$OMP PARALLEL DEFAULT(NONE) &
!$OMP SHARED(l1max, l2max, lmin, kernels, kcross, fullkernels, nkernels) &
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
       if (fullkernels) then
          do l2 = l2start, l2end
             kernels(l2, l1, 1:nkernels) = kernels(l1, l2, 1:nkernels)
          enddo
       else
          do l2 = l2start, l2end
             kcross(l2, l1) = kcross(l1, l2)
          enddo
       endif
    enddo
!$OMP END DO
!$OMP END PARALLEL

    if (fullkernels) then
!$OMP PARALLEL DEFAULT(NONE) SHARED(l2max, kernels) PRIVATE(l2)
!$OMP DO schedule(dynamic, 8)
       do l2 = 0, l2max
          kernels(:,l2,:) = kernels(:,l2,:) * (TWO * l2 + ONE)
       enddo
!$OMP END DO
!$OMP END PARALLEL
    else
!$OMP PARALLEL DEFAULT(NONE) SHARED(l2max, kcross) PRIVATE(l2)
!$OMP DO schedule(dynamic, 8)
       do l2 = 0, l2max
          kcross(:,l2) = kcross(:,l2) * (TWO * l2 + ONE)
       enddo
!$OMP END DO
!$OMP END PARALLEL
    endif

    if (tenormoutput) then
       allocate(TEnorm(0:nlmax,1:1), stat=status)
       call assert_alloc(status, code, 'TEnorm')
       if (fullkernels) then
          do l1=0,nlmax
             TEnorm(l1,1) = sum(kernels(l1,:,4))
          enddo
       else
          do l1=0,nlmax
             TEnorm(l1,1) = sum(kcross (l1,:))
          enddo
       endif
    endif

    deallocate(fl)
!     deallocate(apfunction)
!     deallocate(mu_, w_, Pl_)
    call cpu_time(ct1)
    call wall_clock_time(wt1)
    if (verbose) print*,'Kernel calculation (CPU,wall) time [s]: ',ct1-ct0,wt1-wt0

  end subroutine compute_kernels

!====================================================================================================
  subroutine get_TE_factor(Crossl)
!====================================================================================================
    ! 2010-07-07: automatically triggers allocation and calculation of kernels or kcross
    !-----------------------------------------------------------------------------------------------
    use spice_common, only: dp, i4b, nlmax, kernels, kcross, TEnorm, megaverbose
    use spice_parameters, only: fullkernels
    use misc_utils,       only: assert
    REAL(DP) :: Crossl(0:nlmax)
    
    integer(i4b) :: i

    if (allocated(TEnorm)) then
       Crossl(0:nlmax) = TEnorm(0:nlmax,1)
    else
       if (fullkernels) then
          if (.not. allocated(kernels)) then
             if (megaverbose) print*,'Compute kernels required for TE factor'
             call compute_kernels
          endif
          do i=0,nlmax
             Crossl(i) = sum(kernels(i,:,4))
          enddo
       else
          if (.not. allocated(kcross)) then
             if (megaverbose) print*,'Compute kcross required for TE factor'
             call compute_kernels
          endif
          do i=0,nlmax
             Crossl(i) = sum(kcross(i,:))
          enddo
       endif
    endif
  end subroutine get_TE_factor
    

end module windows
