module map_type

  use healpix_types,    only: DP, HPX_DBADVAL
  use spice_parameters, only: KMAP
  type tcombine
     character(len=4) :: operator = 'junk'
     real(DP)         :: weight   = 1.0_DP
     real(KMAP)       :: fmissval = HPX_DBADVAL
  end type tcombine

  private :: tcombine
  public

contains
!=======================================================================
subroutine myfitsreadmap(mapfile, map, npix, nmaps, fmissval, code)
  !----------------------------------------------------------
  use healpix_types
  use spice_parameters, only: KMAP
  use fits_spice, only: read_fits_cutx, fix_polcconv
  use healpix_modules, only: input_map, getsize_fits, &
       & getnumext_fits, fatal_error
  implicit none

  character(len=*),                        intent(in)  :: mapfile
  integer(I4B),                            intent(in)  :: npix, nmaps
  real(KMAP), dimension(0:npix-1,1:nmaps), intent(out) :: map
  real(KMAP),       optional,              intent(in)  :: fmissval
  character(len=*), optional,              intent(in)  :: code
  integer(I4B) :: njunk, map1type, polcconv, polarisation, ncols, n_ext, extno, i

  integer(I4B), dimension(:),   allocatable :: pixel
  real(KMAP),   dimension(:,:), allocatable :: buffer

  n_ext = getnumext_fits(mapfile)
  njunk = getsize_fits(mapfile, type=map1type, polcconv=polcconv, polarisation=polarisation, nmaps=ncols)
  if (map1type == 3 .and. n_ext == 1) then
     ! cut sky (polarized) map, with I (Q and U) in same single extension
     allocate(pixel(0:njunk-1))
     allocate(buffer(0:njunk-1, 1:nmaps))
     call read_fits_cutx(mapfile, pixel, buffer, extno=0)
     if (nmaps == 3) then
        call fix_polcconv(mapfile, buffer, polcconv)
     endif
     map(:,:) = fmissval
     do i=1, nmaps
        map(pixel,i) = buffer(:,i)
     enddo
     deallocate(pixel)
     deallocate(buffer)

  elseif (map1type == 3 .and. n_ext >= 3 .and. nmaps == 3 &
       ! cut sky polarized map, with I, Q and U in separate extensions
       & .or. map1type == 2 &
       ! full sky map (polarized or not)
       & ) then
     call input_map(mapfile, map, npix, nmaps, fmissval=fmissval) ! 2016-08-09
  else
     call fatal_error(code//" can not read "//trim(mapfile))
  endif
return
end subroutine myfitsreadmap

!=======================================================================
subroutine myfitsreadcutmap(file, tmap, index, fmissval, code)
  !----------------------------------------------------------
  use healpix_types
  use spice_parameters, only: KMAP, maptype
  use fits_spice, only: read_fits_cutx, fix_polcconv
  use healpix_modules, only: fatal_error, getsize_fits, read_fits_cut4, getnumext_fits
  implicit none

  character(len=*),             intent(in)    :: file
  type(maptype),                intent(inout) :: tmap
  integer(I4B),     optional,   intent(in)    :: index
  real(KMAP),       optional,   intent(in)    :: fmissval
  character(len=*), optional,   intent(in)    :: code
  integer(I4B) :: njunk, map1type, istokes

  integer(i8b), dimension(:), allocatable :: obs_npix
  integer(i4b), dimension(:), allocatable :: pixel
  real(SP),     dimension(:), allocatable :: signal
  integer(i4b)  :: extno, imap, imapout, nmaps, status, polcconv, polarisation, n_ext
  logical(LGT)  :: increment
  !character(len=), parameter :: code = "myfitsreadcutmap"


  increment = present(index)
  nmaps = tmap%nmaps
  if (increment) nmaps = 1
  

  n_ext = getnumext_fits(file)
  allocate(obs_npix(1:nmaps))
  obs_npix(1) = getsize_fits(file, extno=0, polcconv=polcconv, polarisation=polarisation)
  if (n_ext == 1) then
     obs_npix(1:nmaps) = obs_npix(1)
  else
     do imap=1,nmaps
        extno = imap - 1
        obs_npix(imap) = getsize_fits(file, extno=extno)
     enddo
  endif


  if (increment .and. obs_npix(1) /= tmap%npix &
       & .or. minval(obs_npix) /= maxval(obs_npix)) then
     print*,nmaps,tmap%nmaps,tmap%npix,increment
     print*,obs_npix
     call fatal_error("Mismatched T,Q,U map size in "&
          & //trim(file)//" in "//code)
  endif

  if ((n_ext == 1  .and. polarisation == 1)  .or.  nmaps == 1) then
     ! read FITS file containing PIXEL, INTENSITY [, Q_POLARISATION, U_POLARISATION]
     if (increment) then
        call read_fits_cutx(file, tmap%pixel, tmap%map(0:, index:index), extno=0)
     else
        call read_fits_cutx(file, tmap%pixel, tmap%map                 , extno=0)
     endif
  else
     ! read FITS file containing PIXEL, SIGNAL [, N_OBS, SERROR]
     ! with 1 extension per Stokes (,Q,U) parameter
     allocate(pixel (0:obs_npix(1)-1), stat = status)
     allocate(signal(0:obs_npix(1)-1), stat = status)
     do imap = 1,nmaps
        extno = imap - 1
        imapout = imap
        if (increment) imapout = index
        call read_fits_cut4(file, int(obs_npix(imap),kind=i4b), &
             &              pixel, signal, extno=extno)
        if (imapout == 1) then
           tmap%pixel = pixel
           tmap%npix  = obs_npix(imap)
           tmap%map(0:,imapout) = signal
        else
           pixel = abs(pixel - tmap%pixel)
           if (maxval(pixel) /= 0) then
              call fatal_error("Mismatched pixels in T,Q,U maps in "&
                   & //trim(file)//" in "//code)
           endif
           tmap%map(0:,imapout) = signal
        endif
     enddo
     deallocate(pixel )
     deallocate(signal)
  endif

  ! ------------ Deal with polcconv -----------
  if (nmaps == 3) then
     call fix_polcconv(file, tmap%map, polcconv)
  endif


  deallocate(obs_npix)
  

  return
end subroutine myfitsreadcutmap
!=======================================================================
subroutine alloc_map_type(tmap, nside, nmaps, npixtot, obsnpix, ordering, coordsys, value)
  !----------------------------------------------------------
  use healpix_types
  use healpix_modules,  only: assert_alloc
  use spice_parameters, only: KMAP, maptype
  implicit none
  type(maptype),    intent(inout) :: tmap
  integer(i4b),     intent(in)    :: nside, nmaps, npixtot, obsnpix, ordering
  character(len=*), intent(in)    :: coordsys
  real(KMAP), optional, intent(in):: value
  !
  integer(i4b)                    :: status
  character(len=*), parameter     :: code='alloc_map_type'

  tmap%nside    = nside
  tmap%nmaps    = nmaps
  tmap%ordering = ordering
  tmap%fmissval = HPX_DBADVAL
  if (obsnpix < int(tmap%THRESHOLD * npixtot)) then
     tmap%type = tmap%CUT4B
     tmap%npix = obsnpix
     allocate(tmap%map(0:tmap%npix-1,1:tmap%nmaps),stat=status)
     call assert_alloc(status, code, 'map')
     allocate(tmap%pixel(0:tmap%npix-1),stat=status)
     call assert_alloc(status, code, 'pixel')
     tmap%pixel = -1
  else
     tmap%type = tmap%FULL
     tmap%npix = npixtot
     allocate(tmap%map(0:tmap%npix-1,1:tmap%nmaps),stat=status)
     call assert_alloc(status, code, 'map')
     tmap%zbounds = (/ -1.d0, 1.d0 /)
  endif

  if (present(value)) then
     tmap%map = value
  endif

  return
end subroutine alloc_map_type
!=======================================================================
subroutine copy_map_type(tmap1, tmap2)
  !----------------------------------------------------------
  use healpix_types
  use healpix_modules,  only: assert_alloc
  use spice_parameters, only: maptype
  implicit none
  type(maptype), intent(out) :: tmap1
  type(maptype), intent(in)  :: tmap2
  integer(I4B) :: status
  character(len=*), parameter :: code = 'copy_map_type'

  tmap1%nside    = tmap2%nside
  tmap1%npix     = tmap2%npix
  tmap1%nmaps    = tmap2%nmaps
  tmap1%ordering = tmap2%ordering
  tmap1%type     = tmap2%type
  tmap1%coordsys = tmap2%coordsys
  tmap1%zbounds  = tmap2%zbounds

  if (allocated(tmap2%pixel)) then
     allocate(tmap1%pixel(0:tmap1%npix-1),stat=status)
     call assert_alloc(status, code, 'pixel')
     tmap1%pixel = tmap2%pixel
  endif
  if (allocated(tmap2%pixel8)) then
     allocate(tmap1%pixel8(0:tmap1%npix-1),stat=status)
     call assert_alloc(status, code, 'pixel8')
     tmap1%pixel8 = tmap2%pixel8
  endif
  if (allocated(tmap2%map)) then
     allocate(tmap1%map(0:tmap1%npix-1,1:tmap1%nmaps),stat=status)
     call assert_alloc(status, code, 'map')
     tmap1%map = tmap2%map
  endif
end subroutine copy_map_type

!=======================================================================
subroutine dealloc_map_type(tmap)
  !----------------------------------------------------------
  use spice_parameters, only: maptype
  implicit none
  type(maptype),intent(inout) :: tmap  
  if (allocated(tmap%pixel)) deallocate(tmap%pixel)
  if (allocated(tmap%map  )) deallocate(tmap%map  )
  tmap%type = tmap%VOID
  return
end subroutine dealloc_map_type
!=======================================================================
subroutine fill_map_type(file, tmap, fmissval, index)
  !----------------------------------------------------------
  use healpix_types
  use healpix_modules,  only: fatal_error
  use spice_parameters, only: KMAP, maptype
  implicit none
  character(len=*), intent(in)           :: file
  type(maptype),    intent(inout)        :: tmap
  real(KMAP),       intent(in), optional :: fmissval
  integer(I4b),     intent(in), optional :: index
  character(len=*), parameter            :: code='fill_map_type'
  integer(i4b)                           :: npix

  if (present(fmissval)) then
     tmap%fmissval = fmissval
  endif

  if (tmap%type == tmap%FULL) then
     npix = tmap%npix
     if (present(index)) then
        call myfitsreadmap(file, tmap%map(:,index:index), npix, 1,          fmissval=fmissval)
     else
        call myfitsreadmap(file, tmap%map,                npix, tmap%nmaps, fmissval=fmissval)
     endif
  elseif (tmap%type == tmap%CUT4B) then
     call myfitsreadcutmap(file, tmap, index, fmissval, code)
  else
     call fatal_error('Unknown type of map in '//code)
  endif
  return
end subroutine fill_map_type
!=======================================================================
subroutine add_map_type(tmapsum, tmap, weight) ! tmapsum = tmapsum + weight * tmap
  !----------------------------------------------------------
  use healpix_types
  use healpix_modules,  only: fatal_error
  use spice_parameters, only: KMAP, maptype
  implicit none
  type(maptype),  intent(inout)    :: tmapsum
  type(maptype),  intent(in)       :: tmap
  real(DP),       intent(in)       :: weight
  character(len=*), parameter      :: code='add_map_type'
  type(tcombine)                   :: tc
  !----------------------------------------------------------
  tc%operator = 'WSUM'
  tc%weight  = weight
  call combine_map_type(tmapsum, tmap, tc)
  return
end subroutine add_map_type
!=======================================================================
subroutine mult_map_type(tmapmult, tmap, fmissval) ! tmapmult = tmapmult * tmap
  !----------------------------------------------------------
  use healpix_types
  use healpix_modules,  only: fatal_error
  use spice_parameters, only: KMAP, maptype
  implicit none
  type(maptype),  intent(inout)        :: tmapmult
  type(maptype),  intent(in)           :: tmap
  real(KMAP),     intent(in), optional :: fmissval
  character(len=*),     parameter      :: code='mult_map_type'
  type(tcombine)                       :: tc
  !----------------------------------------------------------
  tc%operator = 'MULT'
  if (present(fmissval)) then
     tc%fmissval = fmissval
  else
     tc%fmissval = tmap%fmissval
  endif
  call combine_map_type(tmapmult, tmap, tc)
  return
end subroutine mult_map_type
!=======================================================================
subroutine combine_map_type(tmapout, tmap, tc)
  !----------------------------------------------------------
  use healpix_types
  use healpix_modules,  only: fatal_error, strupcase
  use spice_parameters, only: KMAP, maptype
  use sets,             only: intersection
  implicit none
  type(maptype),  intent(inout)    :: tmapout
  type(maptype),  intent(in)       :: tmap
  type(tcombine), intent(in)       :: tc
  character(len=*), parameter      :: code='combine_map_type'
  integer(i8b)                     :: p, p1, p2, ninter
  integer(i4b)                     :: j1, j2, lb1, lb2, q
  logical(LGT)                     :: do_wsum, do_mult, &
       &                              full1, full2, cut1, cut2, defined1
  integer(i4b), dimension(:,:), pointer     :: linter
  real(KMAP),   allocatable, dimension(:,:) :: buffer
  real(KMAP) :: thr

  do_wsum = (strupcase(tc%operator) == 'WSUM')
  do_mult = (strupcase(tc%operator) == 'MULT')
  if (do_wsum .eqv. do_mult) then
     call fatal_error('Choose either WSUM or MULT in '//code)
  endif
  if (tmapout%nside /= tmap%nside) then
     call fatal_error('Inconsistent Nsides in '//code)
  endif
  if (tmapout%ordering /= tmap%ordering) then
     call fatal_error('Inconsistent Orderings in '//code)
  endif
  if (do_wsum .and. (tmapout%nmaps /= tmap%nmaps)) then
     call fatal_error('Added maps must have same dimensions in '//code)
  endif
  thr = abs(1.e-5 * tc%fmissval)

  full1    = (tmapout%type == tmapout%FULL)
  full2    = (   tmap%type ==    tmap%FULL)
  cut1     = (tmapout%type == tmapout%CUT4B)
  cut2     = (   tmap%type ==    tmap%CUT4B)
  defined1 = full1
  if (.not. defined1) defined1 = (maxval(tmapout%pixel) .ge. 0)
  !print*,'full1, cut1, full2, cut2',full1, cut1, full2, cut2
  !print*,'defined1',defined1

  if (full1 .and. full2) then
     do j1=1,tmapout%nmaps
        j2 = min(j1, tmap%nmaps)
        do p=0,tmapout%npix-1
           if (   abs(   tmap%map(p,j2) - tc%fmissval) > thr .and. &
                & abs(tmapout%map(p,j1) - tc%fmissval) > thr) then
              if (do_mult) then ! multiplication
                 tmapout%map(p,j1) = tmapout%map(p,j1) * tmap%map(p,j2)
              else ! weighted sum
                 tmapout%map(p,j1) = tmapout%map(p,j1) + tc%weight * tmap%map(p,j2)
              endif
           else ! invalid pixel
              tmapout%map(p,j1) = tc%fmissval
           endif
        enddo
     enddo
  elseif ((cut1 .and. cut2) .or. (cut1 .and. full2) .or. (full1 .and. cut2)) then
     if (defined1) then
        if (cut1 .and. cut2) then
           linter => intersection(tmapout%pixel, tmap%pixel)
           ninter = size(linter,1)
        elseif (cut1 .and. full2) then
           ninter = tmapout%npix
           allocate(linter(0:ninter-1, 0:2))
           linter(:,0) = tmapout%pixel
           linter(:,1) = (/ (q, q=0,ninter-1) /)
           linter(:,2) = linter(:,0)
        else  ! (full1 .and. cut2) then
           ninter = tmap%npix
           allocate(linter(0:ninter-1, 0:2))
           linter(:,0) = tmap%pixel
           linter(:,1) = linter(:,0)
           linter(:,2) = (/ (q, q=0,ninter-1) /)
        endif
        lb1    = lbound(linter,1)
        lb2    = lbound(linter,2)
        ! pixel
        if (cut2) then
           if (cut1) deallocate(tmapout%pixel)
           allocate(tmapout%pixel(0:ninter-1))
           tmapout%pixel = linter(lb1:ninter-1+lb1,lb2)
        endif
        ! map
        allocate(buffer(0:ninter-1,1:tmapout%nmaps))
        do j1=1,tmapout%nmaps
           j2 = min(j1, tmap%nmaps)
           do p=0,ninter-1
              p1 = linter(p+lb1,1+lb2)
              p2 = linter(p+lb1,2+lb2)
              if (    abs(   tmap%map(p2,j2) - tc%fmissval) > thr .and. &
                   &  abs(tmapout%map(p1,j1) - tc%fmissval) > thr) then ! corrected 2020-04-08
                 if (do_mult) then ! multiplication
                    buffer(p,j1) = tmapout%map(p1,j1) * tmap%map(p2,j2)
                 else ! weighted sum
                    buffer(p,j1) = tmapout%map(p1,j1) + tc%weight * tmap%map(p2,j2)
                 endif
              else
                 buffer(p,j1) = tc%fmissval
              endif
           enddo
        enddo
        deallocate(tmapout%map)
        allocate(tmapout%map(0:ninter-1,1:tmapout%nmaps))
        tmapout%map   = buffer
        tmapout%type  = tmapout%CUT4B
        tmapout%npix  = ninter
        deallocate(buffer)
     else ! undefined (cut) tmapout -> copy weight*tmap or tmap into tmapout
        ! pixel
        deallocate(tmapout%pixel)
        allocate(tmapout%pixel(0:tmap%npix-1))
        tmapout%pixel =          tmap%pixel
        tmapout%npix  =          tmap%npix
        ! map
        deallocate(tmapout%map)
        allocate(tmapout%map(0:tmap%npix-1,1:tmapout%nmaps))
        do j1=1,tmapout%nmaps
           j2 = min(j1, tmap%nmaps)
           do p=0,tmapout%npix-1
              if (abs(   tmap%map(p,j2) - tc%fmissval) > thr) then
                 if (do_mult) then ! multiplication
                    tmapout%map(p,j1) = tmap%map(p,j2)
                 else ! weighted sum
                    tmapout%map(p,j1) = tc%weight * tmap%map(p,j2)
                 endif
              else ! invalid pixel
                 tmapout%map(p,j1) = tc%fmissval
              endif
           enddo
        enddo
!         if (do_mult) then ! multiplication
!            do j1=1,tmapout%nmaps
!               j2 = min(j1, tmap%nmaps)
!               tmapout%map(:,j1)   = tmap%map(:,j2)
!            enddo
!         else ! weighted sum
!            tmapout%map   = tc%weight * tmap%map
!         endif
     endif
  else
     call fatal_error('Unknown type of map in '//code)
  endif
  return
end subroutine combine_map_type
!=======================================================================
subroutine intersect_map_type(tmap1, tmap2)
  !----------------------------------------------------------
  use healpix_types
  use healpix_modules,  only: fatal_error
  use spice_parameters, only: KMAP, maptype
  use sets,             only: intersection
  implicit none
  type(maptype),  intent(inout)    :: tmap1,tmap2
  !
  character(len=*), parameter      :: code='intersect_map_type'
  integer(i4b), dimension(:,:), pointer     :: linter
  integer(i8b) :: p, p1, p2, ninter
  integer(i4b) :: j1, j2, lb1, lb2, q
  real(KMAP),   allocatable, dimension(:,:) :: buffer
  logical(LGT)                     :: full1, full2, cut1, cut2, defined1

!   if (tmap1%type /= tmap2%type) then
!      call fatal_error('Inconsistent map types in '//code)
!   endif
  if (tmap1%nside /= tmap2%nside) then
     call fatal_error('Inconsistent Nsides in '//code)
  endif
  if (tmap1%ordering /= tmap2%ordering) then
     call fatal_error('Inconsistent Orderings in '//code)
  endif
  full1    = (tmap1%type == tmap1%FULL)
  full2    = (tmap2%type == tmap2%FULL)
  cut1     = (tmap1%type == tmap1%CUT4B)
  cut2     = (tmap2%type == tmap2%CUT4B)
  defined1 = full1
  if (.not. defined1) defined1 = (maxval(tmap1%pixel) .ge. 0)

  if (full1 .and. full2) then
     return
  elseif ((cut1 .and. cut2) .or. (cut1 .and. full2) .or. (full1 .and. cut2)) then
     if (defined1) then
        if (cut1 .and. cut2) then
           if (tmap1%npix == tmap2%npix) then
              if (maxval(abs(tmap1%pixel-tmap2%pixel)) == 0) then ! identical maps
                 return ! already matching maps
              endif
           endif
           linter => intersection(tmap1%pixel, tmap2%pixel)
           ninter = size(linter,1)
        elseif (cut1 .and. full2) then
           ninter = tmap1%npix
           allocate(linter(0:ninter-1, 0:2))
           linter(:,0) = tmap1%pixel
           linter(:,1) = (/ (q, q=0,ninter-1) /)
           linter(:,2) = linter(:,0)
        else  ! (full1 .and. cut2) then
           ninter = tmap2%npix
           allocate(linter(0:ninter-1, 0:2))
           linter(:,0) = tmap2%pixel
           linter(:,1) = linter(:,0)
           linter(:,2) = (/ (q, q=0,ninter-1) /)
        endif
        lb1    = lbound(linter,1)
        lb2    = lbound(linter,2)
        ! pixels 1 & 2
        deallocate(tmap1%pixel, tmap2%pixel)
        allocate(tmap1%pixel(0:ninter-1))
        allocate(tmap2%pixel(0:ninter-1))
        tmap1%pixel = linter(lb1:ninter-1+lb1,lb2)
        tmap2%pixel = tmap1%pixel
        ! map 1
        allocate(buffer(0:ninter-1,1:tmap1%nmaps))
        do j1=1,tmap1%nmaps
           do p=0,ninter-1
              p1 = linter(p+lb1,1+lb2)
              buffer(p,j1) = tmap1%map(p1,j1)
           enddo
        enddo
        deallocate(tmap1%map)
        allocate(tmap1%map(0:ninter-1,1:tmap1%nmaps))
        tmap1%map = buffer
        ! map 2
        if (tmap2%nmaps /= tmap1%nmaps) then
           deallocate(buffer)
           allocate(buffer(0:ninter-1,1:tmap2%nmaps))
        endif
        do j2=1,tmap2%nmaps
           do p=0,ninter-1
              p2 = linter(p+lb1,2+lb2)
              buffer(p,j2) = tmap2%map(p2,j2)
           enddo
        enddo
        deallocate(tmap2%map)
        allocate(tmap2%map(0:ninter-1,1:tmap2%nmaps))
        tmap2%map = buffer
        deallocate(buffer)
        !
        tmap1%npix  = ninter
        tmap2%npix  = ninter
        tmap1%type  = tmap1%CUT4B
        tmap2%type  = tmap2%CUT4B
        !print*,'      ',tmap1%npix,tmap2%npix
     else ! tmap1 not filled
        call fatal_error('map1 must be defined in '//code)
!         tmap1%npix = tmap2%npix
!         deallocate(tmap1%map)
!         allocate(tmap1%map(0:tmap1%npix-1,1:tmap1%nmaps),stat=status)
!         call assert_alloc(status, code, 'map')
!         deallocate(tmap1%pixel)
!         allocate(tmap1%pixel(0:tmap1%npix-1),stat=status)
!         call assert_alloc(status, code, 'pixel')
!         tmap1%pixel = tmap2%pixel
     endif
  else
     call fatal_error('Unknown type of map in '//code)
  endif
  return
end subroutine intersect_map_type
!=======================================================================
subroutine process_map_type(tmap, process, verbose)
  !----------------------------------------------------------
  use healpix_types
  use healpix_modules,  only: fatal_error, strlowcase, assert_alloc, &
       & convert_nest2ring, convert_ring2nest, nest2ring, ring2nest, &
       & iindexx, string
  use spice_parameters, only: maptype
  implicit none
  type(maptype),    intent(inout)        :: tmap
  character(len=*), intent(in)           :: process
  logical(LGT),     intent(in), optional :: verbose
  character(len=*), parameter            :: code='process_map_type'
  character(len=FILENAMELEN)             :: proc

  logical(LGT) :: myverbose, do_nest2ring, do_ring2nest, do_boolean
  logical(LGT) :: do_ring2ring, do_nest2nest
  integer(i8b) :: i, ipixbad, pold
  integer(i4b) :: j, p, status
  integer(i4b), allocatable, dimension(:) :: pout

  myverbose=.false.
  if (present(verbose)) myverbose = verbose
  proc = strlowcase(process)
  do_nest2ring = (proc(1:9) == 'nest2ring' .and. tmap%ordering == tmap%NESTED)
  do_ring2ring = (proc(1:9) == 'nest2ring' .and. tmap%ordering == tmap%RING)
  do_ring2nest = (proc(1:9) == 'ring2nest' .and. tmap%ordering == tmap%RING)
  do_nest2nest = (proc(1:9) == 'ring2nest' .and. tmap%ordering == tmap%NESTED)
  do_boolean   = (proc(1:7) == 'boolean')


  if (do_boolean) then
     do j=1,tmap%nmaps
        ipixbad = 0
        do i=0,tmap%npix-1
           if (tmap%map(i,j) /= 0. .and. tmap%map(i,j) /= 1.) then
              tmap%map(i,j) = 1.
              ipixbad       = ipixbad + 1
           endif
        enddo
        if (ipixbad > 0 .and. myverbose) then
           write(*,*)
           write(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
           write(*,*) 'WARNING in input_files :'
           write(*,*) 'The mask file has pixels different from 0 or 1 '
           write(*,*) 'Pixels different from 0 are set equal to 1'
           write(*,*) 'Fraction of pixels changed :',real(ipixbad)/real(tmap%npix)
           write(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
           write(*,*)
        endif
     enddo
     return
  endif

  if (tmap%type == tmap%FULL) then
     if (do_nest2ring) then
        if (myverbose) write(*,*) '... Nested to Ring conversion'
        call convert_nest2ring(tmap%nside, tmap%map)
        tmap%ordering = tmap%RING
        return
     elseif (do_ring2nest) then
        if (myverbose) write(*,*) '... Ring to Nested conversion'
        call convert_ring2nest(tmap%nside, tmap%map)
        tmap%ordering = tmap%NESTED
        return
     elseif (do_ring2ring .or. do_nest2nest) then
        return
     else
        call fatal_error('unknown processing '//trim(process)//' for full sky map in '//code)
        return
     endif
  endif

  if (tmap%type == tmap%CUT4B) then
     if (do_nest2ring .or. do_ring2nest .or. do_ring2ring .or. do_nest2nest) then
        allocate(pout(0:tmap%npix-1), stat=status)
        call assert_alloc(status,code,"pout")
        if (do_nest2ring .or. do_ring2nest) then
           ! convert pixel index
           if (do_nest2ring) then
              if (myverbose) write(*,*) '... Nested to Ring conversion (and pixel sorting)'
              do p=0,tmap%npix-1
                 call nest2ring(tmap%nside, tmap%pixel(p), pout(p))
              enddo
              tmap%ordering = tmap%RING
           else
              if (myverbose) write(*,*) '... Ring to Nested conversion (and pixel sorting)'
              do p=0,tmap%npix-1
                 call ring2nest(tmap%nside, tmap%pixel(p), pout(p))
              enddo
              tmap%ordering = tmap%NESTED
           endif
           tmap%pixel = pout
        else
           if (myverbose) write(*,*) '... pixel sorting (and duplicates detection)'
        endif
        ! sort pixels and maps on (new) pixel index
        call iindexx(int(tmap%npix,I4B),tmap%pixel,pout) ! pout is 1-based
        tmap%pixel = tmap%pixel(pout-1)
        do j=1,tmap%nmaps
           tmap%map(:,j) = tmap%map(pout-1,j)
        enddo
        deallocate(pout)
        ! check for duplicates
        pold = tmap%pixel(0)
        do i=1,tmap%npix-1
           p = tmap%pixel(i)
           if (p == pold) then
              call fatal_error('multiple pixels with same index ('&
                   & //trim(string(p,format='(i0)'))//') detected by '//code)
           endif
           pold = p
        enddo
        return
     else
        call fatal_error('unknown processing '//trim(process)//' for cut sky map in '//code)
     endif
  endif

  call fatal_error('Unknown type of map in '//code)
  return
end subroutine process_map_type
!=======================================================================
subroutine remove_dipole_map_type(tmap, index, degree, multipoles, &
     & fmissval, mask, weights, silent)
  !----------------------------------------------------------
  use healpix_types
  use healpix_modules,   only: fatal_error
  use spice_parameters,  only: maptype, KMAP
  use remove_dipole_mod, only: remove_dipole
  implicit none
  type(maptype), intent(inout)            :: tmap 
  integer(i4b),  intent(in)               :: index, degree
  real(DP),      intent(out),dimension(0:):: multipoles
  real(KMAP),    intent(in),   optional   :: fmissval
  type(maptype), intent(inout),optional   :: mask, weights
  logical(LGT),  intent(in),  optional    :: silent
!!!  real(KMAP),    parameter                :: fbad_value = -1.6375e30_KMAP

  multipoles = 0.0_dp
  if (tmap%type == tmap%FULL) then
     call remove_dipole(tmap%nside, tmap%map(:,index), tmap%ordering, degree, &
          & multipoles, tmap%zbounds, &
          & fmissval, mask%map(:,index), weights%map(:,index), silent)
  else
     call sub_remove_dipole_map_type(tmap, index, degree, multipoles, &
          & fmissval, mask, weights, silent)
     !call fatal_error('Dipole removal not implemented for partial map')
  endif
  return
end subroutine remove_dipole_map_type
!=======================================================================
subroutine sub_remove_dipole_map_type(tmap, index, degree, multipoles, &
     & fmissval, tmask, tweights, silent)
  !----------------------------------------------------------
  use healpix_types
  use healpix_modules,   only: fatal_error, pix2vec_ring, pix2vec_nest, &
       & dsvdcmp, dsvbksb, nside2npix
  ! nside2npix, npix2nside, query_strip
  use spice_parameters,  only: maptype, KMAP
  use remove_dipole_mod, only: remove_dipole
  implicit none
  type(maptype), intent(inout)            :: tmap 
  integer(i4b),  intent(in)               :: index, degree
  real(DP),      intent(out),dimension(0:):: multipoles
  real(KMAP),    intent(in), optional     :: fmissval
  type(maptype), intent(inout), optional     :: tmask, tweights
  logical(LGT),  intent(in), optional     :: silent

  ! local
  real   (kind=KMAP)                :: fmiss_effct
  integer(kind=i4b)                 :: ipix, npix, n_mult, n_mult1, nside
  integer(kind=i4b)                 :: i, ordering
  logical(lgt)                      :: dodipole, do_mask, do_weights
  real(kind=dp)                     :: flag, temp, wmin
  real(kind=dp), dimension(1:3)     :: vec
  real(kind=dp), dimension(0:3)     :: b, wdiag
  real(kind=dp), dimension(0:3,0:3) :: mat, umat, vmat, imat
  real(kind=dp), dimension(0:3)     :: dmultipoles, tmpvec
  real(kind=dp)                     :: theta1, theta2
  integer(kind=i4b)                 :: ncpix, ncp, nbad, npfull
  !integer(kind=i4b), dimension(:), allocatable :: cut_pixel
  logical :: be_verbose
  character(len=30)                 :: legend

  character(len=*),      parameter :: code = "SUB_REMOVE_DIPOLE_MAP_TYPE"
  real   (kind=KMAP),    parameter :: fbad_value = -1.6375e30_KMAP

  !============================================================

  !!npix = nside2npix(nside)
  multipoles = 0.0_dp
  fmiss_effct = fbad_value
  if (present(fmissval)) fmiss_effct = fmissval
  be_verbose = .true.
  if (present(silent)) be_verbose = .not.silent

  do_mask = .false.
  if (present(tmask)) then
     if (tmask%npix > 100) then
        do_mask = .true.
     else
        if (tmask%npix > 1) print*,'WARNING: '//code//' mask ignored'
     endif
  endif

  do_weights = .false.
  if (present(tweights)) then
     if (tweights%npix > 100) then
        do_weights = .true.
     else
        if (tweights%npix > 1) print*,'WARNING: '//code//' weights ignored'
     endif
  endif

  nside    = tmap%nside
  ordering = tmap%ordering
  npfull   = nside2npix(nside)
  if (ordering < 1 .or. ordering > 2) then
     print*,code//': ordering should be one of 1 (RING) or 2 (NESTED)'
     print*,'       it is: ', ordering
     print*,code//"> ABORT ! "
     call fatal_error       
  endif

  if (   do_mask    .and. tmask%ordering    /= ordering  .or. &
       & do_weights .and. tweights%ordering /= ordering) then
     print*,code//': Mismatched ordering between map, mask and weights'
     call fatal_error
  endif
  if (   do_mask    .and. tmask%nside    /= nside  .or. &
       & do_weights .and. tweights%nside /= nside) then
     print*,code//': Mismatched Nside between map, mask and weights'
     call fatal_error
  endif

  if (degree == 0) then
     if (be_verbose) print*,code//"> No monopole nor dipole removal"
     return
  elseif (degree == 1) then
     dodipole = .false.
     legend = 'monopole'
  else if (degree == 2) then
     dodipole = .true.
     legend = 'monopole and dipole'
  else
     print*,code//"> degree can only be "
     print*,"      0: nothing done, or"
     print*,"      1: monopole (l=0) removal or "
     print*,"      2: monopole and dipole (l=0,1) removal"
     print*,code//"> ABORT ! "
     call fatal_error
  endif

  n_mult = (degree)**2
  n_mult1 = size(multipoles)
  if (n_mult1 < n_mult) then 
     print*,'WARNING: '//code//': multipoles array is not large enough to accomodate results:'
     print*, 'size=',n_mult1, ';   needs:',n_mult
     n_mult = n_mult1
  endif


  if (do_mask .and. do_weights) call intersect_map_type(tmask, tweights)
  if (do_mask                 ) call intersect_map_type(tmap,  tmask)
  if (              do_weights) call intersect_map_type(tmap,  tweights)


! !  !!!!! TO BE DONE !!!!
!   !----------------------------------------------
!   ! flag out pixels excluded by cut
!   !----------------------------------------------
!   if (zbounds(2) <= zbounds(1)) then ! remove one strip
!      ncpix = npix * (zbounds(1)-zbounds(2))/2. * 1.1 + 10* nside
!   else ! remove 2 strips
!      ncpix = npix * (2.0 + zbounds(1)-zbounds(2))/2 * 1.1 + 10*nside
!   endif
!   allocate(cut_pixel(0:ncpix))
!   theta1 = acos(zbounds(1)) ; theta2 = acos(zbounds(2))
!   call query_strip(nside, theta1, theta2, cut_pixel, ncp, nest=ordering-1, inclusive=0)
!   if (ncp > 0) then
!      if (be_verbose) print*,code,"> "//trim(string(ncp,format='(i0)'))//&
!           & " pixels flagged out in galactic cut"
!      call alloc_map_type(tcut, nside, 1, ncp, ordering, '')
!      tcut%pixel = cut_pixel(0:ncp-1)
!      tcut%map   = fmiss_effct
!   endif
!   deallocate(cut_pixel)

  !----------------------------------------------
  ! generate least square linear system
  !----------------------------------------------
  mat = 0.0_dp
  b   = 0.0_dp
  nbad = npfull - tmap%npix
  do ipix = 0, tmap%npix-1

     ! flag = 1 for good values
     !      = 0 for bad values
     flag = 1.0_dp
     if (do_weights) flag = flag * tweights%map(ipix,index)

     if ( abs(tmap%map(ipix,index) - fmiss_effct) <= abs(1.e-5*fmiss_effct) ) then
        nbad = nbad + 1
        goto 20
     endif
     if (do_mask) then
        if (abs(tmask%map(ipix,index)) <= 1.e-10) then
           nbad = nbad + 1
           goto 20
        endif
     endif

     ! pure monopole related quantities
     temp = tmap%map(ipix,index) * flag
     b(0) = b(0) + temp
     mat(0,0) = mat(0,0) + flag
     !!if (temp /= temp) print*,ipix,temp,'  <- NaN'

     ! dipole related quantities
     if (dodipole) then
        ! computes dipole basis functions
        ! pixel -> vector
        if (ordering == 1) then 
           call pix2vec_ring( nside, tmap%pixel(ipix), vec)
        else
           call pix2vec_nest( nside, tmap%pixel(ipix), vec)
        endif

        ! construct vector T*(1,x,y,z)
        ! computes dipole basis functions
        b(1:3) = b(1:3) + temp * vec(1:3)

        ! construct (half of) matrix (1,x,y,z)#(1,x,y,z)
        mat(1:3,0) = mat(1:3,0) + vec(1:3)*flag
        do i = 1, 3
           mat(i:3,i) = mat(i:3,i) + vec(i:3) * (vec(i) * flag)
        enddo
     endif

20   continue
  enddo

  if (be_verbose) then
9001 format(a,i0,a,g12.5,a)
     write(*,9001) code//' Excluding ',nbad, &
          & ' pixels when computing '//trim(legend)//' (', &
          &  100.d0*nbad/dble(npfull),'%)'
  endif

  if (mat(0,0) == 0.d0) then
     print*,code//'> Sum of weights is zero'
     call fatal_error
  endif

  if (dodipole) then
     ! symmetric matrix
     do i=0,2
        mat(i,i+1:3) = mat(i+1:3,i)
     enddo
     
     if ( abs( (mat(1,1) + mat(2,2) + mat(3,3))/mat(0,0) - 1.d0) > 1.d-10) then
        print*,code//'> Error in matrix construction.'
        call fatal_error              
     endif
  endif

  !----------------------------------------------
  ! solve system    mat . (mono, dip_x, dip_y, dip_z) = b
  !----------------------------------------------

  if (dodipole) then
     ! SVD decomposition
     umat = mat
     call dsvdcmp(umat, 4, 4, 4, 4, wdiag, vmat)
     ! thresholding
     wmin = maxval(wdiag)* 1.e-6_dp
     where (wdiag < wmin)
        wdiag = 0.0_dp
     end where
     ! back substitution
     call dsvbksb(umat, wdiag, vmat, 4, 4, 4, 4, b, dmultipoles)

  else
     dmultipoles(0) = b(0) / mat(0,0) ! average temperature
  endif

  !----------------------------------------------
  ! remove monopole and dipole
  !----------------------------------------------
  do ipix = 0, tmap%npix-1
     if ( abs(tmap%map(ipix,index) - fmiss_effct) < abs(1.e-5*fmiss_effct) ) goto 10
     temp = dmultipoles(0)
     if (dodipole) then
        ! computes dipole basis functions
        ! pixel -> vector
        if (ordering == 1) call pix2vec_ring( nside, tmap%pixel(ipix), vec)
        if (ordering == 2) call pix2vec_nest( nside, tmap%pixel(ipix), vec)
        temp = temp + SUM( dmultipoles(1:3) * vec(1:3))
     endif
     tmap%map(ipix,index) = tmap%map(ipix,index) - temp
10   continue
  enddo

  multipoles(0:n_mult-1) = dmultipoles(0:n_mult-1)

  return
end subroutine sub_remove_dipole_map_type
!=======================================================================

end module map_type
