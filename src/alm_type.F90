module alm_type
  public
  private::cheap_isqrt ! already public in Healpix 3.60

  interface findzrange
     module procedure findzrange_1, findzrange_n
  end interface

  interface alm2cl
     module procedure alm2cl_spice, alm2cl_type
  end interface

  interface sub_alm2cl
     module procedure sub_alm2cl, sub_alm2cl_type
  end interface
contains
  !=======================================================================
  function cheap_isqrt(lin) result (lout)
    !=======================================================================
    use healpix_types
    implicit none
    integer(i4b), intent(in) :: lin
    integer(i4b) :: lout
    lout = floor(sqrt(dble(lin)), kind=I4B)
    return
  end function cheap_isqrt
    
  !=======================================================================
  subroutine pix2ringinfo_ring(nside, pixel, iz, kshift, pmin, pmax, nr, z)
    !=======================================================================
    use healpix_types
    implicit none
    integer(i4b), parameter :: MKD = I4B
    integer(kind=i4b), intent(in)            :: nside 
    integer(kind=MKD), intent(in)            :: pixel  ! ring ordered pixel index
    integer(kind=i4b), intent(out), optional :: iz     ! ring in {1 .. 4nside-1}
    integer(kind=i4b), intent(out), optional :: kshift !  shift in {0,1}
    integer(kind=MKD), intent(out), optional :: pmin   ! lowest pixel in ring
    integer(kind=MKD), intent(out), optional :: pmax   ! highest pixel in ring
    integer(kind=MKD), intent(out), optional :: nr     ! number of pixels in ring
    real(kind=DP),     intent(out), optional :: z      ! cos(theta) in ]-1,1[
    
    integer(MKD)    :: npix, ncap, nl2, nl4, ip, ir
    real(DP)        :: fact1, fact2
    integer(i4b)    :: iz_in, kshift_in
    integer(MKD)    :: pmin_in, pmax_in, nr_in
    real(DP)        :: z_in
    
    !=======================================================================
    
    npix  = (12_MKD * nside) * nside
    ncap  = 2*nside*(nside-1_MKD) ! number of pixels in the north polar cap
    nl4   = 4_MKD*nside
    nl2   = 2_MKD*nside
    fact2 = (3.000000000000_dp*nside)*nside
    fact1 =  1.500000000000_dp*nside
    !     ------------ identifies ring number --------------
    if (pixel < ncap) then ! North Polar cap -------------
       iz_in = (cheap_isqrt(2*pixel+2) + 1)/2
       z_in  =  1.0_dp - (iz_in / fact2) * iz_in
    elseif (pixel < npix-ncap) then ! Equatorial region ------
       ip    = pixel - ncap
       iz_in = INT( ip / nl4 ) + nside ! counted from North pole
       z_in  = (nl2 - iz_in) / fact1
    else ! South Polar cap -----------------------------------
       ip    = npix - pixel
       iz_in = (cheap_isqrt(2*ip) + 1) / 2     ! counted from South pole
       z_in  =  - 1.0_dp + (iz_in / fact2) * iz_in
       iz_in = nl4 - iz_in ! counted from North Pole
    endif
    
    !     ------------ compute pixel range in ring --------------
    if (iz_in >= nside .and. iz_in <= 3*nside) then ! equatorial region
       ir = iz_in - nside + 1  ! in {1, 2*nside + 1}
       pmin_in = ncap  + nl4*(ir-1_MKD) !  lowest pixel number in the ring
       pmax_in = pmin_in + nl4 - 1_MKD    ! highest pixel number in the ring
       kshift_in = MODULO(ir,2)
       nr_in = nl4
    else
       if (iz_in < nside) then       !    north pole
          ir = iz_in
          pmin_in = 2*ir*(ir-1_MKD)        !  lowest pixel number in the ring
          pmax_in = pmin_in + 4*ir - 1_MKD   ! highest pixel number in the ring
       else                          !    south pole
          ir = nl4 - iz_in
          pmin_in = npix - 2*ir*(ir+1_MKD) !  lowest pixel number in the ring
          pmax_in = pmin_in + 4*ir - 1_MKD   ! highest pixel number in the ring
       endif
       nr_in = ir*4
       kshift_in = 1
    endif

    if  (present(iz))     iz     = iz_in
    if  (present(kshift)) kshift = kshift_in
    if  (present(pmin))   pmin   = pmin_in
    if  (present(pmax))   pmax   = pmax_in
    if  (present(nr))     nr     = nr_in
    if  (present(z))      z      = z_in

  end subroutine pix2ringinfo_ring
  
  !=======================================================================
  !         FINDZRANGE
  !=======================================================================
  function findzrange_n(nside, map) result (zrange)
    !=======================================================================
    ! find out the z-range on which the pixels are non-zero
    ! assuming the map to be RING ordered
    !=======================================================================
    use healpix_types
    use spice_parameters, only: KMAP
    implicit none
    integer(i4b),             intent(in)  :: nside
    real(KMAP),               intent(in)  :: map(0:,1:)
    real(DP), dimension(1:2)              :: zrange
    integer(i4b)                          :: i,nd
    real(DP), dimension(:,:), allocatable :: zr
    
    nd = size(map,2)
    allocate(zr(1:2,1:nd))
    do i=1,nd
       zr(1:2,i) = findzrange_1(nside, map(0:,i))
    enddo
    zrange(1) = minval(zr(1,:))
    zrange(2) = maxval(zr(2,:))
    
    return
  end function findzrange_n
  !=======================================================================
  function findzrange_1(nside, map) result (zrange)
    use healpix_types
    use healpix_modules, only: nside2npix, pix2ang_ring, fatal_error
    use spice_parameters, only: KMAP
    implicit none
    integer(i4b),             intent(in)  :: nside
    real(KMAP),               intent(in)  :: map(0:)
    real(DP), dimension(1:2)              :: zrange
    integer(i8b)                 :: p, npix, pmin, pmax, np
    real(DP)                     :: thnorth,thsouth, phi
    real(DP),         parameter  :: eps=1.d-6
    character(len=*), parameter  :: code='findzrange'
    !logical(LGT),     parameter :: optimize_zrange = .false.
    logical(LGT),     parameter :: optimize_zrange = .true.
    
    zrange = (/ -1.d0, 1.d0 /)
    if (optimize_zrange) then
       npix = nside2npix(nside)
       np   = size(map)
       if (np /= npix) then
          write(*,*) code//'> inconsistent map size and Nside'
          print*,nside,npix,np
          call fatal_error
       endif
       
       do pmin=0,npix-1 ! find first non-zero pixel, starting from North Pole
          if (map(pmin) /= 0.0_KMAP) exit
       enddo
       if (pmin >= npix) then
          write(*,*) code//'> No non-zero pixel found'
          call fatal_error
       endif
       
       do pmax=npix-1,pmin,-1 ! find last non-zero pixel
          if (map(pmax) /= 0.0_KMAP) exit
       enddo
       if (pmin > pmax) then
          write(*,*) code//'> Inconsistent pixel indexes'
          call fatal_error
       endif
       
       ! pmin < pmax -> thnorth  < thsouth
       call pix2ang_ring(nside, pmin, thnorth, phi)
       call pix2ang_ring(nside, pmax, thsouth, phi)
       ! zrange(1) < zrange(2): analyze on the range [zrange(1),zrange(2)]
       zrange = (/ max(cos(thsouth)-eps,-1.d0), min(cos(thnorth)+eps,1.d0) /)
       
!!!print*,'Fsky = ',(zrange(2)-zrange(1))/2.d0
    endif
    
    return
  end function findzrange_1


  !=======================================================================
  !         SUB_ALM2CL
  !========================================================
  subroutine sub_alm2cl(alm1, i1, alm2, i2, cl, i3)
    !========================================================
    use healpix_types
    use spice_parameters, only: KMAPC
    use misc_utils, only: fatal_error
    implicit none
    integer(I4B),                        intent(in) :: i1, i2, i3
    complex(KMAPC), dimension(1:,0:,0:), intent(in) :: alm1, alm2
    real(DP),       dimension(0:,1:),    intent(out):: cl

    integer(I4B) :: nlmax, nmmax, l, mm, m, j1, j2, j3
    complex(DPC) :: dc
    real(DP), parameter :: two = 2.000000000000000000_dp
    real(DP), parameter :: one = 1.000000000000000000_dp

    nlmax = min( size(alm1, 2), size(alm2, 2), size(cl, 1)) - 1
    nmmax = min( size(alm1, 3), size(alm2, 3)) - 1
    j1 = size(alm1, 1)
    j2 = size(alm2, 1)
    j3 = size(cl,   2)

    if (i1 > j1 .or. i2 > j2 .or. i3 > j3) then
       call fatal_error('invalid index in alm -> C(l)')
    endif
    do l = 0, nlmax
       mm = min(l, nmmax)
       dc = cmplx(0.0_dp, 0.0_dp, kind=DP)
       do m=1,mm
          dc =   dc + cmplx(      alm1(i1,l,m) , kind=DP) &
               &     *cmplx(conjg(alm2(i2,l,m)), kind=DP)
       enddo
       dc = (dc + conjg(dc)) + cmplx(      alm1(i1,l,0) , kind=DP) &
            &                 *cmplx(conjg(alm2(i2,l,0)), kind=DP)
       cl(l,i3) = real(dc, kind=DP) / (two*l + one)
    enddo

    return
  end subroutine sub_alm2cl
  !========================================================
  subroutine sub_alm2cl_type(talm1, i1, talm2, i2, cl, i3)
    !========================================================
    use healpix_types,    only: I4B, DP, DPC, LGT
    use healpix_modules,  only: fatal_error
    use spice_parameters, only: almtype
    implicit none
    integer(I4B),                        intent(in) :: i1, i2, i3
    type(almtype),                       intent(in) :: talm1, talm2
    real(DP),       dimension(0:,1:),    intent(out):: cl

    integer(I4B) :: nlmax, nmmax, l, mm, m, j1, j2, j3
    integer(I4B) :: nl1max, nm1max, nl2max, nm2max
    integer(I4B) :: nlm1, nlm2
    complex(DPC) :: dc
    complex(DPC), allocatable, dimension(:) :: x1, x2
    real(DP), parameter :: two = 2.000000000000000000_dp
    real(DP), parameter :: one = 1.000000000000000000_dp
    logical(LGT) :: tab1, tab2, vec1, vec2

    nl1max = talm1%lmax
    nl2max = talm2%lmax
    nm1max = talm1%mmax
    nm2max = talm2%mmax
    nlm1   = 2*nl1max + 1
    nlm2   = 2*nl2max + 1
    nlmax = min( nl1max, nl2max, size(cl, 1) - 1)
    nmmax = min( nm1max, nm2max                 )
    j1 = talm1%npol
    j2 = talm2%npol
    j3 = size(cl,   2)
    vec1 = (talm1%type == talm1%VECT)
    tab1 = (talm1%type == talm1%TAB)
    vec2 = (talm2%type == talm2%VECT)
    tab2 = (talm2%type == talm2%TAB)

    if (vec1 .eqv. tab1 .or. vec2 .eqv. tab2) then
       call fatal_error('Must be either TAB or VECT alm in alm2cl')
    endif

    if (i1 > j1 .or. i2 > j2 .or. i3 > j3) then
       call fatal_error('invalid index in alm -> C(l)')
    endif

    allocate(x1(0:nlmax),x2(0:nlmax))
    do l = 0, nlmax
       mm = min(l, nmmax)
       ! prepare alm's
       if (tab1) then
          x1(0:mm) = cmplx(talm1%almtab(i1, l, 0:mm), kind=DP)
       else
          do m=0,mm
             x1(m) = cmplx(talm1%almvec(i1, l + (m * (nlm1 - m))/2), kind=DP)
          enddo
       endif
       if (tab2) then
          x2(0:mm) = cmplx(talm2%almtab(i2, l, 0:mm), kind=DP)
       else
          do m=0,mm
             x2(m) = cmplx(talm2%almvec(i2, l + (m * (nlm2 - m))/2), kind=DP)
          enddo
       endif
       ! sum the (small) terms m>0
       dc = cmplx(0.0_dp, 0.0_dp, kind=DP)
       do m=1,mm
          dc =   dc + x1(m) * conjg(x2(m))
       enddo
       ! add the (possibly) leading m=0 term
       dc = (dc + conjg(dc)) + x1(0) * conjg(x2(0))
       ! normalize
       cl(l,i3) = real(dc, kind=DP) / (two*l + one)
    enddo
    deallocate(x1, x2)

    return
  end subroutine sub_alm2cl_type

  !=======================================================================
  !         ALM2CL
  !========================================================
  subroutine alm2cl_type(talm1, talm2, cl, symmetric)
  !========================================================
    ! computes C(l) from a_lm, in the order
    ! TT, [EE, TE, BB, [TB, EB, [ET, BT, BE]]]
    !
    ! TE= alm1_T * alm2_E 
    ! unless symmetric is set : TE = (alm1_T * alm2_E + alm1_E * alm2_T)/2
    !=======================================================
    use healpix_types
    use spice_parameters, only: almtype
    implicit none
    type(almtype),                       intent(in) :: talm1, talm2
    real(DP)    ,   dimension(0:, 1: ),  intent(out):: cl
    logical(LGT),   optional,            intent(in):: symmetric
    ! 
    integer(I4B) :: nlmax, nmmax
    integer(I4B) :: ncl, na1, na2, k1, k2
    real(DP), parameter :: half = 0.500000000000000000_dp
    logical(LGT) :: polarisation, bcoupling, do_sym, asympol

    real(DP), allocatable, dimension(:,:) :: cl_work
    
    !========================================================

    ncl   = size(cl, 2)
    nlmax = size(cl, 1) - 1
    na1 = talm1%npol
    na2 = talm2%npol
    polarisation = (na1 >= 3 .and. na2 >= 3 .and. ncl >=4)
    bcoupling    = (ncl >=6) .and. polarisation
    asympol      = (ncl >=9) .and. polarisation
    do_sym = .false.
    cl = 0.0_DP
    if (present(symmetric)) do_sym = symmetric
    if (polarisation .and. do_sym) then
       print*,'Symmetric TE C(l)'
    endif
    if (polarisation) then
       allocate(cl_work(0:nlmax, 1:2))
    endif

     ! TT power spectrum
    call sub_alm2cl(talm1, 1, talm2, 1, cl, 1)
    if (polarisation) then
       ! GG or EE power spectrum
       call sub_alm2cl(talm1, 2, talm2, 2, cl, 2)
       ! CC or BB power spectrum
       call sub_alm2cl(talm1, 3, talm2, 3, cl, 3)

       ! TG or TE power spectrum
       call sub_alm2cl(talm1, 1, talm2, 2, cl_work, 1)
       call sub_alm2cl(talm1, 2, talm2, 1, cl_work, 2)
       k1 = 4  ;   k2 = k1 +3
                    cl(0:,k1) =  cl_work(0:,1)
       if (asympol) cl(0:,k2) =                cl_work(0:,2)
       if (do_sym ) cl(0:,k1) = (cl_work(0:,1)+cl_work(0:,2)) * half

       if (bcoupling) then
          ! TC or TB power spectrum
          call sub_alm2cl(talm1, 1, talm2, 3, cl_work, 1)
          call sub_alm2cl(talm1, 3, talm2, 1, cl_work, 2)
          k1 = 5  ;   k2 = k1 +3
                       cl(0:,k1) =  cl_work(0:,1)
          if (asympol) cl(0:,k2) =                cl_work(0:,2)
          if (do_sym ) cl(0:,k1) = (cl_work(0:,1)+cl_work(0:,2)) * half


          ! GC or EB power spectrum
          call sub_alm2cl(talm1, 2, talm2, 3, cl_work, 1)
          call sub_alm2cl(talm1, 3, talm2, 2, cl_work, 2)
          k1 = 6  ;   k2 = k1 +3
                       cl(0:,k1) =  cl_work(0:,1)
          if (asympol) cl(0:,k2) =                cl_work(0:,2)
          if (do_sym ) cl(0:,k1) = (cl_work(0:,1)+cl_work(0:,2)) * half

       endif
    endif

    if (allocated(cl_work)) deallocate(cl_work)

    return
  end subroutine alm2cl_type
  !========================================================
  subroutine alm2cl_spice(nlmax, nmmax, alm1, alm2, cl, symmetric)
  !========================================================
    ! adapted from Healpix alm2cl with SPC alm and DP cl
    !
    ! computes C(l) from a_lm, in the order
    ! TT, [EE, TE, BB, [TB, EB, [ET, BT, BE]]]
    !
    ! TE= alm1_T * alm2_E 
    ! unless symmetric is set : TE = (alm1_T * alm2_E + alm1_E * alm2_T)/2
    !=======================================================
    use healpix_types
    use spice_parameters, only: KMAPC
    implicit none
    integer(I4B),                        intent(in) :: nlmax, nmmax
    complex(KMAPC), dimension(1:,0:,0:), intent(in) :: alm1, alm2
    real(DP)    ,   dimension(0:, 1: ),  intent(out):: cl
    logical(LGT),   optional,            intent(in):: symmetric
    ! 
    integer(I4B) :: ncl, na1, na2, k1, k2
    !complex(DPC)     :: dc, dcs
    !real(DP), parameter :: two = 2.000000000000000000_dp
    !real(DP), parameter :: one = 1.000000000000000000_dp
    real(DP), parameter :: half = 0.500000000000000000_dp
    logical(LGT) :: polarisation, bcoupling, do_sym, asympol

    real(DP), allocatable, dimension(:,:) :: cl_work
    
    !========================================================

    ncl = size(cl, 2)
    na1 = size(alm1, 1)
    na2 = size(alm2, 1)
    polarisation = (na1 >= 3 .and. na2 >= 3 .and. ncl >=4)
    bcoupling    = (ncl >=6) .and. polarisation
    asympol      = (ncl >=9) .and. polarisation
    do_sym = .false.
    cl = 0.0_DP
    if (present(symmetric)) do_sym = symmetric
    if (polarisation .and. do_sym) then
       print*,'Symmetric TE C(l)'
    endif
    if (polarisation) then
       allocate(cl_work(0:nlmax, 1:2))
    endif

     ! TT power spectrum
    call sub_alm2cl(alm1, 1, alm2, 1, cl, 1)
    if (polarisation) then
       ! GG or EE power spectrum
       call sub_alm2cl(alm1, 2, alm2, 2, cl, 2)
       ! CC or BB power spectrum
       call sub_alm2cl(alm1, 3, alm2, 3, cl, 3)

       ! TG or TE power spectrum
       call sub_alm2cl(alm1, 1, alm2, 2, cl_work, 1)
       call sub_alm2cl(alm1, 2, alm2, 1, cl_work, 2)
       k1 = 4  ;   k2 = k1 +3
                    cl(0:,k1) =  cl_work(0:,1)
       if (asympol) cl(0:,k2) =                cl_work(0:,2)
       if (do_sym ) cl(0:,k1) = (cl_work(0:,1)+cl_work(0:,2)) * half

       if (bcoupling) then
          ! TC or TB power spectrum
          call sub_alm2cl(alm1, 1, alm2, 3, cl_work, 1)
          call sub_alm2cl(alm1, 3, alm2, 1, cl_work, 2)
          k1 = 5  ;   k2 = k1 +3
                       cl(0:,k1) =  cl_work(0:,1)
          if (asympol) cl(0:,k2) =                cl_work(0:,2)
          if (do_sym ) cl(0:,k1) = (cl_work(0:,1)+cl_work(0:,2)) * half


          ! GC or EB power spectrum
          call sub_alm2cl(alm1, 2, alm2, 3, cl_work, 1)
          call sub_alm2cl(alm1, 3, alm2, 2, cl_work, 2)
          k1 = 6  ;   k2 = k1 +3
                       cl(0:,k1) =  cl_work(0:,1)
          if (asympol) cl(0:,k2) =                cl_work(0:,2)
          if (do_sym ) cl(0:,k1) = (cl_work(0:,1)+cl_work(0:,2)) * half

       endif
    endif

    if (allocated(cl_work)) deallocate(cl_work)

    return
  end subroutine alm2cl_spice

  !================================================================
  !================================================================
  subroutine alloc_alm_type(talm, npol, lmax, mmax, vector)
    !--------------------------------------------------
    use healpix_types,    only: I4B, LGT
    use healpix_modules,  only: assert_alloc, fatal_error
    use spice_parameters, only: almtype, KMAPC
    implicit none
    type(almtype), intent(out)          :: talm
    integer(I4B),  intent(in)           :: npol, lmax, mmax
    logical(LGT),  intent(in), optional :: vector
    
    logical(LGT) :: myvector
    integer(I4B) :: status

    myvector = .false.
    if (present(vector)) myvector = vector

    talm%npol = npol
    talm%lmax = lmax
    talm%mmax = mmax
    if (myvector) then
       talm%type = talm%VECT
       allocate(talm%almvec(1:npol, 0:lmax+(mmax*(2*lmax-mmax+1))/2 ), stat=status)
       talm%almvec = 0.0_KMAPC
    else 
       talm%type = talm%TAB
       allocate(talm%almtab(1:npol, 0:lmax, 0:mmax), stat=status)
       talm%almtab = 0.0_KMAPC
    endif
    !print*,'alloc_alm_type',myvector,talm%type,talm%VECT,talm%TAB,&
    !     & talm%lmax,talm%mmax

    return
  end subroutine alloc_alm_type
  !================================================================
  subroutine dealloc_alm_type(talm)
    !--------------------------------------------------
    use spice_parameters, only: almtype
    implicit none
    type(almtype), intent(inout)          :: talm

    if (allocated(talm%almvec)) deallocate(talm%almvec)
    if (allocated(talm%almtab)) deallocate(talm%almtab)
    talm%type = talm%VOID
    return
  end subroutine dealloc_alm_type

  !================================================================
  subroutine map2alm(tmap, talm, ringweights, index)
    !================================================================
    use spice_parameters, only: maptype, almtype, KMAP, KMAPC
    use healpix_types,    only: I4B, I8B, SP, DP, LGT
    use healpix_modules,  only: fatal_error
    use alm_tools,        only: m2a_hpx => map2alm ! rename Healpix's map2alm into m2a_hpx
    use healpix_sharp_f90_cut
    implicit none
    type(maptype), intent(inout),target         :: tmap
    type(almtype), intent(inout)                :: talm
    real(DP),      intent(in), dimension(1:,1:) :: ringweights
    integer(I4B),  intent(in), optional         :: index
    !
    logical(LGT)                                :: slice
    real(DP), dimension(1:2)                    :: zrange
    character(len=*),        parameter          :: code="map2alm"
    real(KMAP),   dimension(:,:), allocatable   :: mapring
    integer(i4b), dimension(:),   allocatable   :: rings
    integer(i4b) :: pmin, pmax, ring, rmin, rmax, nrings, nmaps 
    integer(i8b) :: np
    real(DP)     :: znorth, zsouth
    real(KMAP), pointer, dimension(:) :: pmap1

    slice = present(index)

    if (tmap%type == tmap%FULL) then
       if(talm%type == talm%TAB) then
          if (slice) then
             pmap1  => tmap%map(0:,index)
             zrange = findzrange(tmap%nside, pmap1)
             call m2a_hpx(tmap%nside, talm%lmax, talm%mmax, pmap1, talm%almtab, zrange, ringweights)
             nullify(pmap1)
          else
             zrange = findzrange(tmap%nside, tmap%map)
             call m2a_hpx(tmap%nside, talm%lmax, talm%mmax, tmap%map, talm%almtab, zrange, ringweights)
          endif
          tmap%zbounds = zrange
       else
          call fatal_error("configuration (full map, vec alm) not implemented in "//code)
       endif
    else
       call pix2ringinfo_ring(tmap%nside, minval(tmap%pixel), iz=rmin, pmin=pmin, z=znorth)
       call pix2ringinfo_ring(tmap%nside, maxval(tmap%pixel), iz=rmax, pmax=pmax, z=zsouth)
       tmap%zbounds = (/ zsouth, znorth /)
       np     = pmax - pmin + 1
       nrings = rmax - rmin + 1
       allocate(rings(1:nrings))
       rings  = (/ (ring, ring=rmin, rmax) /) ! list of useful rings
       nmaps  = tmap%nmaps
       if (slice) nmaps = 1
       allocate(mapring(0:np-1,1:nmaps)) ! put input map in rings
       mapring = 0.
       if(talm%type == talm%TAB) then
          if (slice) then
             mapring(tmap%pixel - pmin,1:1) = tmap%map(:,index:index)
#if DOUBLE
             call sharp_hpcut_map2alm_x_d(&
#else
             call sharp_hpcut_map2alm_x_s(&
#endif
             & tmap%nside, talm%lmax, talm%mmax, talm%type, nrings, rings, np, &
             & mapring, talm%almtab, ringweights)
          else
             mapring(tmap%pixel - pmin,:) = tmap%map(:,:)
#if DOUBLE
             call sharp_hpcut_map2alm_pol_x_d(&
#else
             call sharp_hpcut_map2alm_pol_x_s(&
#endif
             & tmap%nside, talm%lmax, talm%mmax, talm%type, nrings, rings, np, &
             & mapring, talm%almtab, ringweights)
          endif
       else
          if (slice) then
             mapring(tmap%pixel - pmin,1:1) = tmap%map(:,index:index)
#if DOUBLE
             call sharp_hpcut_map2alm_x_d(&
#else
             call sharp_hpcut_map2alm_x_s(&
#endif
             & tmap%nside, talm%lmax, talm%mmax, talm%type, nrings, rings, np, &
             & mapring, talm%almvec, ringweights)
          else
             mapring(tmap%pixel - pmin,:) = tmap%map(:,:)
#if DOUBLE
             call sharp_hpcut_map2alm_pol_x_d(&
#else
             call sharp_hpcut_map2alm_pol_x_s(&
#endif
             & tmap%nside, talm%lmax, talm%mmax, talm%type, nrings, rings, np, &
             & mapring, talm%almvec, ringweights)
          endif          
       endif
       deallocate(rings)
       deallocate(mapring)
    endif
    return
  end subroutine map2alm
  !================================================================
  subroutine dump_alm(talm, almfile, verbose)
  !================================================================
    use spice_parameters, only: almtype, KMAP, KMAPC, version
    use healpix_types,    only: FILENAMELEN, I4B, I8B, LGT
    use healpix_modules,  only: add_card, assert_alloc, dump_alms, &
         &                      fatal_error, write_minimal_header
    implicit none
    type(almtype),    intent(in)           :: talm
    character(len=*), intent(in)           :: almfile
    logical(LGT),     intent(in), optional :: verbose

    character(len=FILENAMELEN)   :: fn_bang
    integer(I4B)       :: i, l, nlheader, nlmax, nmmax, npol, status
    integer(I8B)       :: m, j0, j1
    character(len=80), dimension(1:120)  :: header
    logical(LGT)                         :: polar, myverbose
    complex(KMAPC), allocatable, dimension(:,:,:) :: almbuff
    !-----------------------------------------------------------------

    if (trim(almfile) == '' .or. trim(almfile) == '!') return
    myverbose = .false.
    if (present(verbose)) myverbose = verbose

    fn_bang = trim(adjustl(almfile))
    if (fn_bang(1:1) /= '!') fn_bang = '!'//trim(fn_bang)! add leading '!'
    if (myverbose) write(*,*) 'Dumping alm into '//trim(almfile)

    nlmax = talm%lmax
    nmmax = talm%mmax
    npol  = talm%npol
    polar = (npol == 3)
    do i = 1, npol
       header = ''
       call write_minimal_header(header,'alm', &
            creator = 'Spice', version = VERSION, &
            nlmax = nlmax, nmmax = nmmax, polar = polar) 
       if (i == 1) then
          call add_card(header,"EXTNAME","'ANALYSED a_lms (TEMPERATURE)'")
       elseif (i == 2) then
          call add_card(header,"EXTNAME","'ANALYSED a_lms (GRAD / ELECTRIC component)'")
       elseif (i == 3) then
          call add_card(header,"EXTNAME","'ANALYSED a_lms (CURL / MAGNETIC component)'")
       else
          call fatal_error('Too many spins in dumped alms to '//trim(almfile))
       endif
       nlheader = SIZE(header)

       if (talm%type == talm%TAB) then
          call dump_alms(fn_bang, talm%almtab(i,0:nlmax,0:nmmax), &
               &                  nlmax, header, nlheader, i-1_i4b)

       elseif (talm%type == talm%VECT) then
          if (i == 1) then
             allocate(almbuff(1:1, 0:nlmax, 0:nmmax), stat=status)
             call assert_alloc(status, 'dump_alm', 'almbuff')
             almbuff = 0.0_KMAPC
          endif
          do m=0,nmmax
             j0 = (m * (2*nlmax+1 - m))/2
             j1 = j0 + nlmax - m - 1
             almbuff(1,m:nlmax,m) = talm%almvec(i,j0:j1)
          enddo
          call dump_alms(fn_bang, almbuff(1,0:nlmax,0:nmmax), &
               &                  nlmax, header, nlheader, i-1_i4b)
          if (i == npol) deallocate(almbuff)

       else
          call fatal_error('Wrong type of dumped alms to '//trim(almfile))
       endif
    enddo ! i: extension
    return
  end subroutine dump_alm
  !================================================================
  subroutine edit_alm(talm, almfile, verbose)
  !================================================================
    use spice_parameters, only: almtype, KMAP, KMAPC, version
    use healpix_types,    only: FILENAMELEN, I4B, I8B, LGT
    use healpix_modules,  only: assert_alloc, fatal_error, &
         &                      number_of_alms, read_conbintab
    implicit none
    type(almtype),    intent(inout)        :: talm
    character(len=*), intent(in)           :: almfile
    logical(LGT),     intent(in), optional :: verbose

    integer(I4B) :: iext, next, status, l, m
    integer(I4B) :: nalms, lm, j
    integer(I4B) :: nlmax, nmmax, npol
    real(KMAP), dimension(:,:), allocatable :: alm_read
    complex(KMAPC) :: z
    logical(LGT)                         :: myverbose
    !-----------------------------------------------------------------
    if (trim(almfile) == '') return

    myverbose = .false.
    if (present(verbose)) myverbose = verbose
    if (myverbose) write(*,*) 'Modifing alm with weights read from '//trim(almfile)

    nlmax = talm%lmax
    nmmax = talm%mmax
    npol  = talm%npol

    nalms = number_of_alms(almfile, next)
    allocate(alm_read(0:nalms, 1:6),stat=status)
    call assert_alloc(status, 'edit_alm', 'alm_read')
    do iext=1, min(next,npol)
       call read_conbintab(almfile, alm_read, nalms, extno=iext-1)
       do lm=0, nalms-1
          l = nint( alm_read(lm, 1) )
          m = nint( alm_read(lm, 2) )
          if (l <= nlmax .and. m <= nmmax .and. abs(m) <= l .and. m >= 0) then
             z = cmplx(alm_read(lm, 3), alm_read(lm, 4), kind=KMAP)
             if (talm%type == talm%TAB) then
                talm%almtab(iext, l, m) = talm%almtab(iext, l, m) * z
             else
                j = l + (m * (2*nlmax+1 - m))/2
                talm%almvec(iext, j) = talm%almvec(iext, j) * z
             endif
          endif
       enddo
    enddo
    deallocate(alm_read)

    return
  end subroutine edit_alm
end module alm_type
