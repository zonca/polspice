! edited to run with Healpix 3.*
module deal_with_files
! subroutine check_headers
! subroutine assert_read_piofile
! subroutine set_bad_pixels_to_0
! subroutine input_files
! subroutine output_results
! subroutine check_header
! subroutine set_nlmax
! subroutine check_header_asc
! subroutine check_open_output_files
! subroutine write_corcl_file
! subroutine read_corcl_file
! subroutine check_mylread
! subroutine generate_gaussian_beam
! subroutine my_generate_beam
! function my_isfits
! subroutine get_pixwin_filename
! subroutine read_twod_fits_DP
! subroutine write_twod_fits_DP
! subroutine write_nd_fits_DP
! subroutine check_twod_fits_header
! subroutine check_nd_fits_header
! subroutine deletefile
!!!!!!!! function my_getnumext_fits
! subroutine write_spice_header

  interface set_nlmax
     module procedure set_nlmax_scal, set_nlmax_vect
  end interface
contains
!=======================================================================
subroutine check_headers
!=======================================================================
  use spice_common
  use misc_utils, only: fatal_error, string
#ifdef PIO
  use piolib
#endif
  implicit none
  integer(I4B) :: my_npixtot
  logical(LGT) :: head_defined !must be .false. to get nside, npixtot, nlmax from file
  ! then subsequent calls with .true. will compare previsou nside and npixtot to file value
  integer(I8B) :: myerr,my_group,my_nx8,firstline,lastline
  character(len=FILENAMELEN) :: grpname='' !SP
  integer(i4b) :: imap, jstokes, k

  do imap=1,2
     head_defined=.false.
     do jstokes=1,3
        do k=1,nlmaps(jstokes,imap)
           call check_header(listmapfile(k,jstokes,imap), verbose, megaverbose, .true., &
                &            nsides(imap), npixtots(imap), obsnpixs(imap), &
                &            nlmax, ordering_tqu(k,jstokes,imap), dry,&
                &            polarization, head_defined) 
        enddo
     enddo
  enddo


  if (clmapinput) call check_header_asc(cl_inmap_file, nsides, nlmax, ncor,&
 &                                      head_defined, verbose, megaverbose, npixtots)
  if (masks_present) then
     call check_header(maskfile, verbose, megaverbose, .false.,&
          &            nsides(1), npixtots(1), obsnpixs_m(1), nlmax, ordering_m1, dry,&
          &            .false., head_defined)
  endif
  if (masks2_present) then
     call check_header(maskfile2, verbose, megaverbose, .false.,&
 &                     nsides(2), npixtots(2), obsnpixs_m(2), nlmax, ordering_m2, dry, &
 &                     .false., head_defined)
  endif
  if (clmaskinput) then
     ncmask = -1 ! read from file actual number of columns
     call check_header_asc(cl_inmask_file, nsides, nlmax, ncmask,&
 &                         head_defined, verbose, megaverbose, npixtots)
  endif
  if (weights_present) then
     call check_header(weightfile, verbose, megaverbose, .false.,&
 &                     nsides(1), npixtots(1), obsnpixs_w(1), nlmax, ordering_w1, dry, &
 &                     .false., head_defined)
  endif
  if (weights2_present) then
     call check_header(weightfile2, verbose, megaverbose, .false.,&
 &                     nsides(2), npixtots(2), obsnpixs_w(2), nlmax, ordering_w2, dry, &
 &                     .false., head_defined)
  endif

  if (tenorminput) then
     call check_te_header(tenorm_in_file, verbose, megaverbose, nlmax, &
          thetamax, apodize, apodizetype, apodizesigma)
  endif

  if (windowinput) then
     if (megaverbose) write(*,*) 'Check header for file '//trim(window_in_file)
#ifdef PIO
     myerr=PIOgetgrpname(grpname,window_in_file)
     my_group=PIOopentab2dgrp(grpname,"r")
     myerr=PIOgettab2dlines(firstline,lastline,window_in_file,my_group)
     my_nx8=PIOgettab2dcolumngrp(my_group)
     if (myerr /= 0) then
        write(*,*) 'ERROR in check_headers'
        write(*,*) 'Could not check header for file'
        write(*,*) trim(window_in_file)
        write(*,*) PIOerrmess(myerr)
     else
        windowfile_ny=lastline-firstline ! +1 !SP
        windowfile_nx=my_nx8
     endif
#else 
     call check_twod_fits_header(window_in_file,my_npixtot,windowfile_nx,windowfile_ny)
#endif
     if (windowfile_nx < nlmax+1 .or. windowfile_ny < 2*nlmax+1) then
        write(*,*) 'ERROR in check_headers'
        write(*,*) 'The input window function dimensions are too small for current study'
        write(*,*) 'nx found =',windowfile_nx,',   nx needed =',nlmax+1
        write(*,*) 'ny found =',windowfile_ny,',   ny needed =',2*nlmax+1
        CALL FATAL_ERROR
     endif
     if ((windowfile_nx > nlmax+1 .or. windowfile_ny > 2*nlmax+1) .and. verbose) then
        write(*,*) 'WARNING in check_headers'
        write(*,*) 'The input window function dimensions are too large for current study'
        write(*,*) 'only a subset will be used'
        write(*,*) 'nx found =',windowfile_nx,',   nx needed =',nlmax+1
        write(*,*) 'ny found =',windowfile_ny,',   ny needed =',2*nlmax+1
     endif
  endif

  ! This can be put only here because the default name of pixwin_file depends
  ! on the value of nside

  if (.not.map2_present) then
     nsides(2)   = nsides(1)
     npixtots(2) = npixtots(1)
  endif

  ! print*,map2_present,pwf2_set
  if (map2_present) then ! 2 maps
     if (.not.pwf2_set) then ! 2nd pixel, not set: copy 1st one
        pixwin_file_name_input2 = pixwin_file_name_input
        default_pix_file2       = default_pix_file
        correct_pix2            = correct_pix
     endif
  else ! a single map, same pixels
     pixwin_file_name_input2 = pixwin_file_name_input
     default_pix_file2       = default_pix_file
     correct_pix2            = correct_pix
  endif
  ! print*,trim(pixwin_file_name_default),trim(pixwin_file_name_input), default_pix_file, correct_pix
  ! print*,trim(pixwin_file_name_default2),trim(pixwin_file_name_input2),default_pix_file2,correct_pix2

  if (correct_pix) call get_pixwin_filename(pixwin_file_name_default, &
 &                                          pixwin_file_name_input, &
 &                                          pixwin_file_name, &
 &                                          default_pix_file,nsides(1),healpix_data,stringlength)

  if (correct_pix2) call get_pixwin_filename(pixwin_file_name_default2, &
 &                                          pixwin_file_name_input2, &
 &                                          pixwin_file_name2, &
 &                                          default_pix_file2,nsides(2),healpix_data,stringlength)

end subroutine check_headers

!--------------------------------------------------------------------
subroutine assert_read_piofile(nread, nexpected, file, code)
  use healpix_types
  use misc_utils, only: fatal_error
  implicit none
  integer(I8B) :: nread
  integer(I4B) :: nexpected
  character(len=*) :: file, code

  if (nread /= nexpected) then
     write(*,*) 'ERROR in '//code
     write(*,*) 'Inconsistent number of data for file'
     write(*,*) trim(file)
     write(*,*) 'Number of pixels detected :',nread
     write(*,*) 'Npixtot should be :',nexpected
     CALL FATAL_ERROR
  endif
  return
end subroutine assert_read_piofile
!--------------------------------------------------------------------
subroutine set_bad_pixels_to_0(map, npix, nmaps, bad_value, verbose, name)
  use spice_parameters, only: KMAP
  use healpix_types
  implicit none

  !integer, parameter :: KMAP = SP
  integer(i4b),               intent(in)    :: npix, nmaps
  real(KMAP), dimension(0:npix-1,1:nmaps), intent(inout) :: map
  real(KMAP),                 intent(in)    :: bad_value
  logical(lgt),               intent(in)    :: verbose
  character(len=*),           intent(in)    :: name

  integer(i4b) :: n_nan, i, j

  n_nan = 0
  if (verbose) write(*,*) 'Setting bad pixels to 0 ('//trim(name)//')'
  do j=1,nmaps
     do i=0,npix-1
        if (abs(map(i,j)/bad_value-1.0) < 1.e-6) map(i,j) = 0.
        if (map(i,j) /= map(i,j)) then ! to detect NaN
           n_nan = n_nan + 1
           map(i,j) = 0.
        endif
     enddo
  enddo
  if (n_nan > 0) then
     write(*,*)          '=============================================================================='
     write(*,'(a,i0,a)') ' **WARNING**: Detected ',n_nan,' NaN-valued pixels in input '//trim(name)//', replaced with 0.'
     write(*,*)          '             Note that the mask/weight were NOT updated.'
     write(*,*)          '=============================================================================='
  endif

  return
end subroutine set_bad_pixels_to_0
!=======================================================================
subroutine input_files
!=======================================================================
  use spice_common
  use map_type, only: alloc_map_type, fill_map_type,  &
       & add_map_type, process_map_type, dealloc_map_type
  use misc_utils, only: assert_alloc, fatal_error, string, wall_clock_time
  use pix_tools,  only: convert_nest2ring
#ifdef PIO
  use piolib
  use my_pio_routines
  use fitstools, only: read_dbintab
#else
  use fitstools, only: read_dbintab, getsize_fits, fits2cl
!!  use fitstools, only: read_dbintab,input_map, getsize_fits, fits2cl
#endif
  implicit none

  integer(I4B) :: ipixbad,i,j,my_npixtot, status,n_nan
  integer(I4B) :: istokes, map1type, map2type, njunk, kem
  real(DP), allocatable, dimension(:) :: kc_buff
  real(DP) :: dnullval
  logical(LGT) :: anynull

  integer(I8B) :: myerr,mygroup,nbdata, p
  character(len=FILENAMELEN) :: grpname
  !real(SP),pointer,dimension(:) :: tmpmap
  real(DP),pointer,dimension(:) :: tmpmapDP
  character(len=*), parameter :: code = 'input_files'
  !real(KMAP), allocatable, dimension(:,:) :: tmp_buffer
  type(maptype) :: tmp_buffer
  integer(I4B) :: nmap_buffer
  integer(I4B) :: jmap, kstokes, k
  character(len=80), allocatable, dimension(:) :: header
  real(SP) :: ctime0, ctime1, wtime0, wtime1
  character(len=8) :: coordsys = '' ! place holder
  real(KMAP), parameter :: ZERO = 0.0_KMAP

  call wall_clock_time(wtime0)
  call cpu_time(ctime0)

! -------------------------------------------------------------------
! read all the maps (if applicable), re-order them and combine them
! -------------------------------------------------------------------
  nmap_buffer = nmap
  if (maxval(nlmaps(1:3,1)) > 0) then
     call alloc_map_type(tmap1_in, nsides(1), nmap, npixtots(1), obsnpixs(1), 1, coordsys, value=ZERO) ! 2019-10
  endif
  if (maxval(nlmaps(1:3,2)) > 0) then
     call alloc_map_type(tmap2_in, nsides(2), nmap, npixtots(2), obsnpixs(2), 1, coordsys, value=ZERO) ! 2019-10
  endif

  if (allocated(tmap1_in%map)) then
     do jmap=1,2 ! 1st (and optional second) map(s)
        call alloc_map_type(tmp_buffer, nsides(jmap), nmap, npixtots(jmap), obsnpixs(jmap), 1, coordsys)

        do i=1,nlmaps(1, jmap)
           tmp_buffer%ordering = ordering_tqu(i,1,jmap)
           ! read (polarized) map, and set bad pixels to 0
#ifdef PIO
           do kstokes=1,nmap ! T,Q,U
              if (megaverbose) write(*,'(a,i1,a,a)') &
                   & ' Importing (',jmap,')' ,trim(listmapfile(i,kstokes,jmap))
              call mypioreadmap (listmapfile(i,kstokes,jmap), tmp_buffer%map(0:,kstokes), npixtots(jmap), code)
           enddo
           call set_bad_pixels_to_0(tmp_buffer%map, npixtots(1), nmap, HPX_SBADVAL, .false., 'map')
#else
           kstokes = 1
           if (megaverbose) write(*,'(a,i1,a,a)') &
                & ' Importing (',jmap,') ',trim(listmapfile(i,kstokes,jmap))
           !call myfitsreadmap(listmapfile(i,kstokes,jmap), tmp_buffer, npixtots(jmap), nmap, fmissval=0._KMAP)
           call fill_map_type(listmapfile(i,kstokes,jmap), tmp_buffer, fmissval=ZERO)
#endif
           ! reorder to RING if necessary
!            if (ordering_tqu(i,1,jmap) == 2) then
!               if (megaverbose) print*,'... and reordering ...'
!               call process_map_type(tmp_buffer,'nest2ring')
!            endif
           call process_map_type(tmp_buffer,'nest2ring',verbose=verbose) ! 2020-05-26
           ! combine
           if (megaverbose) print*,'... and combining '//trim(string(listmapw8(i,1,jmap)))
!            if (jmap == 1) map_in  = map_in  + listmapw8(i,1,jmap) * tmp_buffer
!            if (jmap == 2) map2_in = map2_in + listmapw8(i,1,jmap) * tmp_buffer
           if (jmap == 1) call add_map_type(tmap1_in, tmp_buffer, listmapw8(i,1,jmap))
           if (jmap == 2) call add_map_type(tmap2_in, tmp_buffer, listmapw8(i,1,jmap))
        enddo
        call dealloc_map_type(tmp_buffer)
     enddo
  endif

! -------------------------------------------------------------------
! read all the binary masks (if applicable) and re-order them
! -------------------------------------------------------------------
  if (masks_present) then
     if (megaverbose) write(*,*) 'Read input mask file'
     call alloc_map_type(tmask1_map, nsides(1), nmask1, npixtots(1), obsnpixs_m(1), ordering_m1, coordsys) ! 2019-10
#ifdef PIO
     call mypioreadmap(maskfile, tmask1_map%map(0:,1), npixtots(1), code)
     if (maskfilep_toread) then
        call mypioreadmap(maskfilep, tmask1_map%map(0:,2), npixtots(1), code)
     endif
#else
     call fill_map_type(maskfile, tmask1_map, fmissval=ZERO, index=1)
     if (maskfilep_toread) then
        call fill_map_type(maskfilep, tmask1_map, fmissval=ZERO, index=2)
     endif
#endif	
     ! guaranty that the masks are zero or one
     call process_map_type(tmask1_map, 'boolean', verbose=verbose)
     ! convert mask1 to RING if necessary
     ! if (tmask1_map%ordering == 2) call process_map_type(tmask1_map,'nest2ring')
     call process_map_type(tmask1_map,'nest2ring',verbose=verbose) ! 2020-05-26
  endif

  if (masks2_present) then
     if (megaverbose) write(*,*) 'Read input mask file (2)'
     ! allocate(mask2_map(0:npixtots(2)-1,1:nmask2))
     call alloc_map_type(tmask2_map, nsides(2), nmask2, npixtots(2), obsnpixs_m(2), ordering_m2, coordsys) ! 2019-10
#ifdef PIO
     call mypioreadmap(maskfile2, tmask2_map%map(0:,1), npixtots(2), code)
     if (maskfilep2_toread) then
        call mypioreadmap(maskfilep2, tmask2_map%map(0:,2), npixtots(2), code)
     endif
#else
     call fill_map_type(maskfile2, tmask2_map, fmissval=ZERO, index=1)
     if (maskfilep2_toread) then
        call fill_map_type(maskfilep2, tmask2_map, fmissval=ZERO, index=2)
     endif
#endif	
     ! guaranty that the masks are zero or one
     call process_map_type(tmask2_map, 'boolean', verbose=verbose)
     ! convert mask2 to RING if necessary
     ! if (tmask2_map%ordering == 2) call process_map_type(tmask2_map,'nest2ring')
     call process_map_type(tmask2_map,'nest2ring', verbose=verbose) ! 2020-05-26
  endif

! -------------------------------------------------------------------
! read all the weighted masks (if applicable) and re-order them
! -------------------------------------------------------------------
  if (weights_present) then
     if (megaverbose) write(*,*) 'Read input weight file'
     ! allocate(weight_map(0:npixtots(1)-1,1:nweight1))
     call alloc_map_type(tweight1_map, nsides(1), nweight1, npixtots(1), obsnpixs_w(1), ordering_w1, coordsys) ! 2019-10
#ifdef PIO
     call mypioreadmap(weightfile, tweight1_map%map(0:,1), npixtots(1), code)
     if (weightfilep_toread) then
        call mypioreadmap(weightfilep, tweight1_map%map(0:,2), npixtots(1), code)
     endif     
#else
     call fill_map_type(weightfile, tweight1_map, fmissval=ZERO, index=1)
     if (weightfilep_toread) then
        !call fill_map_type(weightfile, tweight1_map, fmissval=ZERO, index=2) ! 2019-01-20
        call fill_map_type(weightfilep, tweight1_map, fmissval=ZERO, index=2)
     endif
#endif
     if (minval(tweight1_map%map) < 0.) then
        write(*,*) 'ERROR in input_files :'
        write(*,*) 'The weight map has negative values.'
        CALL FATAL_ERROR
     endif
     ! convert weight1 to RING if necessary
     ! if (tweight1_map%ordering == 2) call process_map_type(tweight1_map,'nest2ring')
     call process_map_type(tweight1_map,'nest2ring', verbose=verbose) ! 2020-05-26
  endif

  if (weights2_present) then
     if (megaverbose) write(*,*) 'Read input weight file (2)'
     call alloc_map_type(tweight2_map, nsides(2), nweight2, npixtots(2), obsnpixs_w(2), ordering_w2, coordsys) ! 2019-10
#ifdef PIO
     call mypioreadmap(weightfile2, tweight2_map%map(0:,1), npixtots(2), code)
     if (weightfilep2_toread) then
        call mypioreadmap(weightfilep2, tweight2_map%map(0:,2), npixtots(2), code)
     endif
#else
     call fill_map_type(weightfile2, tweight2_map, fmissval=ZERO, index=1)
     if (weightfilep2_toread) then
        !call fill_map_type(weightfile2, tweight2_map, fmissval=ZERO, index=2) ! 2019-01-20
        call fill_map_type(weightfilep2, tweight2_map, fmissval=ZERO, index=2)
     endif
#endif
     if (minval(tweight2_map%map) < 0.) then
        write(*,*) 'ERROR in input_files :'
        write(*,*) 'The weight map (2) has negative values.'
        CALL FATAL_ERROR
     endif
     ! convert weight2 to RING if necessary
     ! if (tweight2_map%ordering == 2) call process_map_type(tweight2_map,'nest2ring')
     call process_map_type(tweight2_map,'nest2ring', verbose=verbose) ! 2020-05-26
  endif

! -------------------------------------------------------------------
! read window matrix (if applicable)
! -------------------------------------------------------------------
  if (windowinput) then
     if (megaverbose) write(*,*) 'Read input window function'
#ifdef PIO
     nbdata=PIOreadtab2dobject(tmpmapDP,window_in_file,' ')
     if (nbdata < 0_I8B) then
        write(*,*) 'ERROR in input_files :'
        write(*,*) 'Something went wrong while reading file'
        write(*,*) trim(window_in_file)
        write(*,*) PIOerrmess(nbdata)
     endif
     allocate(kcross(0:nlmax,0:2*nlmax))
     do j=0,2*nlmax
        do i=0,nlmax
           kcross(i,j)=tmpmapDP(1 + i + j*windowfile_nx)
        enddo
     enddo
     myerr=PIOdeletetab2dtable(tmpmapDP)
#else
     my_npixtot= windowfile_nx * windowfile_ny
     
     allocate(kc_buff(0:my_npixtot-1))
     call read_twod_fits_DP(window_in_file,my_npixtot,windowfile_nx,windowfile_ny,kc_buff)
     allocate(kcross(0:nlmax,0:2*nlmax))
     do j=0,2*nlmax
        do i=0,nlmax
           kcross(i,j)=kc_buff(i + j*windowfile_nx)
        enddo
     enddo
     deallocate(kc_buff)
#endif
  endif

! -------------------------------------------------------------------
! read smaller l-dependent files
! -------------------------------------------------------------------
  if (clmapinput) then 
     if (megaverbose) write(*,*) 'Read input precomputed Cls for the map'
     allocate(cl_map_precomp (0:nlmax,1:ncor))
     call read_corcl_file(cl_map_precomp,nlmax,cl_inmap_file,.false.,ncor)
  endif
  
  if (clmaskinput) then
     if (megaverbose) write(*,*) 'Read input precomputed Cls for the masks/weights'
     allocate(cl_mask_precomp(0:nlmax,1:ncmask))
     call read_corcl_file(cl_mask_precomp,nlmax,cl_inmask_file,.false.,ncmask)
  endif

  if (tenorminput) then
     if (megaverbose) write(*,*) 'Read previously computed TE normalization function'
     allocate(TEnorm(0:nlmax,1))
     allocate(header(1:200))
     call fits2cl(tenorm_in_file, TEnorm, nlmax, 1, header)
     deallocate(header)
  endif

  if (correct_pix .or. correct_pix2) then
     !lwpix = 4*minval(nsides)
     if (correct_pix) then
        if (megaverbose) write(*,*) 'Read pixel window function ',trim(pixwin_file_name)
        allocate(wlpix1(0:4*nsides(1),1)) 
        call read_dbintab(pixwin_file_name, wlpix1,nsides(1)*4+1,1,dnullval,anynull)
     endif
     if (correct_pix2) then
        if (megaverbose) write(*,*) 'Read pixel window function (2)',trim(pixwin_file_name2)
        allocate(wlpix2(0:4*nsides(2),1)) 
        call read_dbintab(pixwin_file_name2,wlpix2,nsides(2)*4+1,1,dnullval,anynull)
     endif
  endif

  if (subtract_noise_cor.or.subtract_noise_cl) then
     allocate(xi_noise(0:nlmax,1:ncor))
     allocate(cl_noise(0:nlmax,1:ncor))
  endif

  if (subtract_noise_cor) then
     if (megaverbose) write(*,*) 'Read input noise correlation function'
     call read_corcl_file(xi_noise,nlmax,noisecor_file_name,.true.,ncor)
  endif

  if (subtract_noise_cl) then
     if (megaverbose) write(*,*) 'Read input noise Cls'
     call read_corcl_file(cl_noise,nlmax,noisecl_file_name,.false.,ncor)
  endif

  if (correct_transfer_function) then 
     if (megaverbose) write(*,*) 'Read transfer function'
     allocate(transfer_function(0:nlmax,1:ncor))
     call read_corcl_file(transfer_function, nlmax, tf_file, .false., ncor)
  endif

  if (correct_beam .or. correct_beam2) then
     allocate(gb(0:nlmax,1:1)) !<<<<<<<<<<<<<<<<
     gb = 1.0_dp
     allocate(gb2(0:nlmax,1:1)) !<<<<<<<<<<<<<<<<
     gb2 = 1.0_dp
     if (map2_present .and. (fwhm /= fwhm2 .or. trim(beam_file) /= trim(beam_file2))) then
        write(*,*) 'Warning: the 2 beams are different'
        if (fwhm /= fwhm2) print*,fwhm,fwhm2
        if (trim(beam_file) /= trim(beam_file2)) print*,trim(beam_file),trim(beam_file2)
     endif
  endif

  if (correct_beam) then
     if (beam_present) then
        if (megaverbose) write(*,*) 'Read input beam'
        call my_generate_beam(nlmax, gb, megaverbose, beam_file)
     else
        call generate_gaussian_beam(fwhm, nlmax, gb)
     endif
  endif

  if (correct_beam2) then
     if (beam2_present) then
        if (megaverbose) write(*,*) 'Read input beam (2)'
        call my_generate_beam(nlmax, gb2, megaverbose, beam_file2)
     else
        call generate_gaussian_beam(fwhm2, nlmax, gb2)
     endif
  endif

  call wall_clock_time(wtime1)
  call cpu_time(ctime1)

  if (verbose) then
     write(*,*) 'Data input (CPU,Wall) time [s]: ',ctime1-ctime0,wtime1-wtime0
  endif

end subroutine input_files

!=======================================================================
subroutine output_results
!=======================================================================
  use spice_common, only:I4B, I8B, DP, LGT, FILENAMELEN, &
       cloutput, coroutput, clmapoutput, clmaskoutput, megaverbose, &
       cl, cl_final, cl_mask, xi_final, mu, &
       ncor, ncov, nlmax, ncmask, nsides,  &
       cl_file_name, cor_file_name, cov_out_file, &
       cl_outmap_file, cl_outmask_file,  &
       kernels, kernels_out_file, kernelsoutput, nkernels, &
       kcross,  window_out_file,  windowoutput,  &
       TEnorm,  tenorm_out_file,  tenormoutput, &
       do_cov, cov_matrix
  use healpix_modules, only: fatal_error, assert_alloc, add_card, &
       write_minimal_header, write_asctab, del_card
#ifdef PIO
  use piolib
#else
  use spice_parameters, only : decouple, version, thetamax, &
       & apodize, apodizesigma, apodizetype, map_present, map2_present, &
       masks_present, masks2_present, weights_present, weights2_present, &
       weightpower, weightpower2, &
       weightpowerp, weightpowerp2, &
       fwhm, fwhm2, correct_beam, correct_beam2, beam_present, beam2_present,&
       beam_file, beam_file2, &
       correct_transfer_function, subtract_average
#endif
  implicit none
  integer(I4B) :: i,j,k,my_npixtot,my_nx,my_ny,status
  real(DP), allocatable, dimension(:) :: kc_buff

  integer(I8B) :: my_npixtot8,my_nx8,my_ny8, my_nz8, myerr,mygroup, inc !SP
  ! integer(I8B), allocatable, dimension(:) :: index1,index2 !SP
  character(len=FILENAMELEN) :: grpname='' !SP
  character(len=FILENAMELEN) :: command='' !SP
  character(len=10) :: str !SP
  character(len=80), dimension(1:120) :: header
  integer(i4b), dimension(1:3) :: dims
  logical(LGT) :: no_kernels, no_kcross

  if (coroutput) then
     if (megaverbose) write(*,*) 'Output correlation function file'
     call write_corcl_file(xi_final,mu,nlmax,cor_file_name,.true.,ncor,nsides)
  endif

  if (cloutput) then
     if (megaverbose) write(*,*) 'Output Cl file'
!!!!     print*,size(cl_final,1),size(cl_final,2),size(mu),nlmax,ncor,nside
     call write_corcl_file(cl_final,mu,nlmax,cl_file_name,.false.,ncor,nsides)
  endif

  if (clmapoutput) then
     if (megaverbose) write(*,*) 'Output raw Cls for the map'
     call write_corcl_file(cl,mu,nlmax,cl_outmap_file,.false.,ncor,nsides)
  endif

  if (clmaskoutput) then
     if (megaverbose) write(*,*) 'Output raw Cls for the weights/masks'
     call write_corcl_file(cl_mask,mu,nlmax,cl_outmask_file,.false.,ncmask,nsides)
  endif

  !---------------------------------------------------------------------
  if (tenormoutput) then
     if (megaverbose) write(*,*) 'Output the TE normalization function'
     status = 0
     call deletefile(tenorm_out_file,status)
     header = ''
     call write_minimal_header(header,'cl', &
          creator = 'Spice', version = version, nlmax=nlmax)
     call del_card(header, (/ 'POLAR   ','BCROSS  ','ASYMCL  '/) )
     call add_card(header,'COMMENT','=============================================================')
     call add_card(header,'COMMENT','Spice internal product')
     call add_card(header,'COMMENT','Normalisation function of TE spectrum induced by apodization')
     call add_card(header,'COMMENT','of Angular Correlation Function')
     call add_card(header,'COMMENT','==============================================================')
     call add_card(header,'APODIZE',apodize, "Apodization of Xi  (True/False)")
     if (apodize) then
        call add_card(header,'APOTYPE', apodizetype,  "Apodization type: 0=Gaussian, 1=Cosine")
        call add_card(header,'APOSIGMA',apodizesigma, "[Deg] Apodization parameter")
     endif
     call add_card(header,'THETAMAX',thetamax, "[Deg] Largest lag in angul. correl. fct Xi")
     call add_card(header,'EXTNAME',"'TE NORMALISATION'")
     call add_card(header,'TTYPE1',"'TE NORMALISATION'")
     call add_card(header,'TUNIT1','none')
     call write_asctab (TEnorm, nlmax, 1, header, size(header), tenorm_out_file)
     header = ''
  endif
  !---------------------------------------------------------------------
  if (windowoutput) then
     no_kcross  = .not. allocated(kcross)
     no_kernels = .not. allocated(kernels)
     if (no_kcross .and. no_kernels) then
        print*,'Window function not computed. Skipping output'
     else
        if (megaverbose) then
           if (no_kcross) then
              write(*,*) 'Output the window function (extracted from full kernel)'
           else
              write(*,*) 'Output the window function'
           endif
        endif
#ifdef PIO
        my_npixtot=(nlmax+1)*(2*nlmax+1)

        allocate(kc_buff(0:my_npixtot-1),stat=status)
        call assert_alloc(status,'output_results','kc_buff')
        my_npixtot8=my_npixtot
        inc=0
        do j=0,2*nlmax
           do i=0,nlmax
              if (no_kcross) then
                 kc_buff(inc)=kernels(i,j,4)
              else
                 kc_buff(inc)=kcross(i,j)
              endif
              inc=inc+1
           enddo
        enddo
        myerr=PIOgetgrpname(grpname,window_out_file)
        my_nx8 = INT(nlmax+1,kind=I8B) !SP
        myerr=PIOcreatetab2dgrp(grpname,my_nx8) !SP
        myerr=PIOcreatetab2dobject(window_out_file,'PIODOUBLE') !SP
        mygroup=PIOopentab2dgrp(grpname,'w')
!!$     write(str,*) 2*nlmax !SP
!!$     command = 'tab=*,0:'//trim(adjustl(str)) !SP
        write(command,'(''tab=*,'',I6,'':'',I6)') 0,2*nlmax !SP
        myerr=PIOwritetab2dobject(kc_buff(0:my_npixtot-1), & !SP
             & window_out_file,trim(command),mygroup) !SP
        myerr=PIOclosetab2dgrp(mygroup)
        deallocate(kc_buff) !SP
#else
        my_npixtot=(nlmax+1)*(2*nlmax+1)
        my_nx=nlmax+1
        my_ny=2*nlmax+1
        allocate(kc_buff(0:my_npixtot-1),stat=status)
        call assert_alloc(status,'output_results','kc_buff')
        inc=0
        do j=0,2*nlmax
           do i=0,nlmax
              if (no_kcross) then
                 kc_buff(inc)=kernels(i,j,4)
              else
                 kc_buff(inc)=kcross(i,j)
              endif
              inc=inc+1
           enddo
        enddo
        status=0
        call deletefile(window_out_file,status)
        call write_twod_fits_DP(window_out_file,my_npixtot,my_nx,my_ny,kc_buff)
        deallocate(kc_buff)
#endif
     endif
  endif

  !---------------------------------------------------------------------
  if (kernelsoutput) then
     if (megaverbose) write(*,*) 'Output the full kernels into '//trim(kernels_out_file)
     if (nkernels /= 1 .and. nkernels /= 4) then
        print*,nkernels
        print*,'ERROR: Expected either 1 or 4 kernels'
        call fatal_error
     endif
     dims(1:3) = (/ nlmax+1, 2*nlmax+1, nkernels /)
     my_npixtot8 = product(dims(1:3))
     allocate(kc_buff(0:my_npixtot8-1))
     inc=0
     do k=1, nkernels
        do j=0,2*nlmax
           do i=0,nlmax
              kc_buff(inc)=kernels(i,j, k)
              inc=inc+1
           enddo
        enddo
     enddo
#ifdef PIO
     myerr=PIOGetGrpname(grpname, kernels_out_file)
     my_nx8 = INT(  dims(1), kind=I8B)
     my_ny8 = INT(  dims(2), kind=I8B)
     my_nz8 = INT(  dims(3), kind=I8B)
     write(command,'(''tab=*,*,'',I6,'':'',I6)') 0, my_nz8
     myerr=PIOCreateTAB3DGrp(grpname, my_nx8, my_ny8)
     mygroup=PIOOpenTAB3DGrp(grpname,'w')
     myerr=PIOCreateTAB3DObject(kernels_out_file,'PIODOUBLE')
     myerr=PIOWriteTAB3DObject(kc_buff, kernels_out_file, trim(command), mygroup)
     myerr=PIOWriteKeywordObject(my_nx8, '1st_dim', 'nx', kernels_out_file, mygroup)
     myerr=PIOWriteKeywordObject(my_ny8, '2nd_dim', 'ny', kernels_out_file, mygroup)
     myerr=PIOWriteKeywordObject(my_nz8, '3rd_dim', 'nz', kernels_out_file, mygroup)
     myerr=PIOCloseTAB3DGrp(mygroup)
#else
     header = ''
     if (nkernels == 4) then
        call add_card(header,'COMMENT','------------------------------------------------------------')
        call add_card(header,'COMMENT','          Coupling Kernels for Spice                        ')
        call add_card(header,'COMMENT','   These 4 kernels relate the average Spice C(l) estimator  ')
        call add_card(header,'COMMENT','       to the underlying ''true'' CMB power spectra.        ')
        call add_card(header,'COMMENT','                                                            ')
        call add_card(header,'COMMENT',' These kernels are not required internally by Spice         ')
        call add_card(header,'COMMENT','           (except for Kern(*,*,4))                         ')
        call add_card(header,'COMMENT','  and are only provided for convenience                     ')
        call add_card(header,'COMMENT','  (eg, for cosmological interpretation of the Spice C(l))   ')
        call add_card(header,'COMMENT','                                                            ')
        call add_card(header,'COMMENT',' These kernels depend on the apodization scheme, thetamax   ')
        call add_card(header,'COMMENT','  and Lmax, all listed in this header.                      ')
        call add_card(header,'COMMENT',' The way to use them depends on the estimators              ')
        call add_card(header,'COMMENT','    (decoupled or not) being considered                     ')
        call add_card(header,'COMMENT','   (option ''decouple'' of Spice).                          ')
        call add_card(header,'COMMENT','                                                            ')
        call add_card(header,'COMMENT','     E-B Decoupled Estimators                               ')
        call add_card(header,'COMMENT','  <C_TT(l1)>           = Sum_l2 Kern(l1, l2, 1) C_TT(l2)_true    ')
        call add_card(header,'COMMENT','  <C_EE(l1)>           = Sum_l2 Kern(l1, l2, 3) C_EE(l2)_true    ')
        call add_card(header,'COMMENT','  <C_BB(l1)>           = Sum_l2 Kern(l1, l2, 3) C_BB(l2)_true    ')
        call add_card(header,'COMMENT','  <C_TE(l1)>           = Sum_l2 Kern(l1, l2, 4) C_TE(l2)_true    ')
        call add_card(header,'COMMENT','  <C_TB(l1)>           = Sum_l2 Kern(l1, l2, 4) C_TB(l2)_true    ')
        call add_card(header,'COMMENT','  <C_EB(l1)>           = Sum_l2 Kern(l1, l2, 3) C_EB(l2)_true    ')
        call add_card(header,'COMMENT','          with l1 in {0,lmax}, l2 in {0,2*lmax}                  ')
        call add_card(header,'COMMENT','     NB: In the decoupled case, Kern(*,*,2) is *NOT* used   ')
        call add_card(header,'COMMENT','       (see Eq. (91) of Chon et al. 2004)                   ')
        call add_card(header,'COMMENT','                                                            ')
        call add_card(header,'COMMENT','                                                            ')
        call add_card(header,'COMMENT','                                                            ')
        call add_card(header,'COMMENT','     Coupled estimators                                              ')
        call add_card(header,'COMMENT','  <C_TT(l1)>           = Sum_l2 Kern(l1, l2, 1)  C_TT(l2)_true       ')
        call add_card(header,'COMMENT','  <C_EE(l1)+C_BB(l1)>  = Sum_l2 Kern(l1, l2, 2) (C_EE(l2)+C_BB(l2))  ')
        call add_card(header,'COMMENT','  <C_EE(l1)-C_BB(l1)>  = Sum_l2 Kern(l1, l2, 3) (C_EE(l2)-C_BB(l2))  ')
        call add_card(header,'COMMENT','  <C_TE(l1)>           = Sum_l2 Kern(l1, l2, 4)  C_TE(l2)_true       ')
        call add_card(header,'COMMENT','  <C_TB(l1)>           = Sum_l2 Kern(l1, l2, 4)  C_TB(l2)_true       ')
        call add_card(header,'COMMENT','  <C_EB(l1)>           = Sum_l2 Kern(l1, l2, 3)  C_EB(l2)_true       ')
        call add_card(header,'COMMENT','          with l1 in {0,lmax}, l2 in {0,2*lmax}                      ')
        call add_card(header,'COMMENT','       (see Eqs. (56)-(60) of Chon et al. 2004)                      ')
        call add_card(header,'COMMENT','                                                                     ')
        call add_card(header,'COMMENT','                                                             ')
        call add_card(header,'COMMENT','                                                             ')
        call add_card(header,'COMMENT','   Sanity checks:                                            ')
        call add_card(header,'COMMENT','   1) The Kernels have dimension (Lmax+1)*(2*Lmax+1)*4       ')
        call add_card(header,'COMMENT','   2) The Kernels each have the form                         ')
        call add_card(header,'COMMENT','   Kern(l1,l2,i) = S(l1,l2,i) * (2*l2+1)                     ')
        call add_card(header,'COMMENT','   where S(l1,l2,i) = S(l2,l1,i) is symmetric for i=1,2,3,4  ')
        call add_card(header,'COMMENT','        (see Eqs. (58) and (62) of Chon et al. 2004)         ')
        call add_card(header,'COMMENT','                                                             ')
     else
        call add_card(header,'COMMENT','------------------------------------------------------------')
        call add_card(header,'COMMENT','          Coupling Kernel for Spice                         ')
        call add_card(header,'COMMENT','   This kernel relates the average Spice C(l) TT estimator  ')
        call add_card(header,'COMMENT','       to the underlying ''true'' CMB TT power spectrum.    ')
        call add_card(header,'COMMENT','   The kernels specific to polarization                     ')
        call add_card(header,'COMMENT','     (different from this one)                              ')
        call add_card(header,'COMMENT','       can also be generated by Spice                       ')
        call add_card(header,'COMMENT','                                                            ')
        call add_card(header,'COMMENT',' This kernel is not required internally by Spice            ')
        call add_card(header,'COMMENT','  and is only provided for convenience                      ')
        call add_card(header,'COMMENT','  (eg, for cosmological interpretation of the Spice C(l))   ')
        call add_card(header,'COMMENT','                                                            ')
        call add_card(header,'COMMENT',' This kernel depend on the apodization scheme, thetamax     ')
        call add_card(header,'COMMENT','  and Lmax, all listed in this header.                      ')
        call add_card(header,'COMMENT','                                                            ')
        call add_card(header,'COMMENT','  <C_TT(l1)>           = Sum_l2 Kern(l1, l2, 1) C_TT(l2)_true    ')
        call add_card(header,'COMMENT','          with l1 in {0,lmax}, l2 in {0,2*lmax}                  ')
        call add_card(header,'COMMENT','                                                            ')
        call add_card(header,'COMMENT','   Sanity checks:                                           ')
        call add_card(header,'COMMENT','   1) The Kernel has dimension (Lmax+1)*(2*Lmax+1)          ')
        call add_card(header,'COMMENT','   2) The Kernel has the form                               ')
        call add_card(header,'COMMENT','   Kern(l1,l2) = S(l1,l2) * (2*l2+1)                        ')
        call add_card(header,'COMMENT','   where S(l1,l2) = S(l2,l1) is symmetric                   ')
        call add_card(header,'COMMENT','        (see Eqs. (58) and (62) of Chon et al. 2004)        ')
        call add_card(header,'COMMENT','                                                            ')
     endif
     call add_card(header,'COMMENT','                        ')
     call add_card(header,'COMMENT','                        ')
     call add_card(header,'COMMENT','                        ')
     call add_card(header,'COMMENT','                        ')
     call add_card(header,'COMMENT','                        ')
     call add_card(header,'COMMENT','  ref: Chon et al, 2004, MNRAS 350, 914                      ')
     
     call write_spice_header(header, size(header), nlmax, ncor, nsides, .false.)
     status=0
     call deletefile(kernels_out_file,status)
     call write_nd_fits_DP(kernels_out_file, my_npixtot8, 3, dims, kc_buff, &
          &  header, size(header))
#endif
     deallocate(kc_buff)
  endif

  !---------------------------------------------------------------------
  if (do_cov) then
     if (megaverbose) write(*,*) 'Output the covariance matrix '//trim(cov_out_file)
     if (ncov /= 1) then
        print*,ncov
        print*,'ERROR: Expected 1'
        call fatal_error
     endif
     dims(1:3) = (/ nlmax+1, nlmax+1, ncov /)
     my_npixtot8 = product(dims(1:3))
     allocate(kc_buff(0:my_npixtot8-1))
     inc=0
     do k=1, dims(3)
        do j=0, dims(2)-1
           do i=0, dims(1)-1
              kc_buff(inc)=cov_matrix(i,j, k)
              inc=inc+1
           enddo
        enddo
     enddo
     print*,'matrix dimensions: ',dims
     print*,'matrix min and max: ', minval(kc_buff), maxval(kc_buff)
#ifdef PIO
     myerr=PIOGetGrpname(grpname, cov_out_file)
     my_nx8 = INT(  dims(1), kind=I8B)
     my_ny8 = INT(  dims(2), kind=I8B)
     my_nz8 = INT(  dims(3), kind=I8B)
     write(command,'(''tab=*,*,'',I6,'':'',I6)') 0, my_nz8-1
     myerr=PIOCreateTAB3DGrp(grpname, my_nx8, my_ny8)
     if (myerr /= 0) write(*,*), myerr, 'creategrp ',PIOerrmess(myerr)
     mygroup=PIOOpenTAB3DGrp(grpname,'w')
     myerr=PIOCreateTAB3DObject(cov_out_file,'PIODOUBLE')
     if (myerr /= 0) write(*,*), myerr, 'createobj ',PIOerrmess(myerr)
     myerr=PIOWriteTAB3DObject(kc_buff, cov_out_file, trim(command), mygroup)
     if (myerr < 0) write(*,*), myerr, 'writeobj ',PIOerrmess(myerr)
     myerr=PIOWriteKeywordObject(my_nx8, '1st_dim', 'nx', cov_out_file, mygroup)
     myerr=PIOWriteKeywordObject(my_ny8, '2nd_dim', 'ny', cov_out_file, mygroup)
     myerr=PIOWriteKeywordObject(my_nz8, '3rd_dim', 'nz', cov_out_file, mygroup)
     if (myerr /= 0) write(*,*), myerr, 'writekw ',PIOerrmess(myerr)
     myerr=PIOCloseTAB3DGrp(mygroup)
     if (myerr /= 0) write(*,*), myerr, 'closegrp ',PIOerrmess(myerr)
#else
     header = ''
     call add_card(header,'COMMENT','------------------------------------------------------------')
     call add_card(header,'COMMENT','          Covariance Matrix for Spice                        ')
     call add_card(header,'COMMENT','                        ')
     call add_card(header,'COMMENT','                        ')
     call add_card(header,'COMMENT','                        ')
     call add_card(header,'COMMENT','                        ')
     call add_card(header,'COMMENT','                        ')
     
     call write_spice_header(header, size(header), nlmax, ncor, nsides, .false.)
     status=0
     call deletefile(cov_out_file,status)
     call write_nd_fits_DP(cov_out_file, my_npixtot8, 3, dims, kc_buff, &
          &  header, size(header))
#endif
     deallocate(kc_buff)
  endif

end subroutine output_results

!=======================================================================
subroutine check_header(file_name,verbose,megaverbose,mapfile, &
 &                      nside,npixtot,obsnpix,&
 &                      nlmax,ordering,dry,polarization, &
 &                      head_defined)
!=======================================================================
  use healpix_types
  use misc_utils, only: fatal_error
#ifdef PIO
  use piolib  
#else
  use fitstools, only: getsize_fits, getnumext_fits
#endif
  implicit none

  character(len=*), intent(IN)    :: file_name
  logical,          intent(IN)    :: verbose,megaverbose,mapfile
  integer(I4B),     intent(INOUT) :: nside,npixtot,obsnpix
  integer(I4B),     intent(OUT)   :: nlmax,ordering
  logical,          intent(IN)    :: dry,polarization
  logical,          intent(INOUT) :: head_defined
  ! if head_defined=.false. get nside, nlmax and npixtot from file, and switch head_defined to .true.
  ! if head_defined=.true. compare provided nside and npixtot to those read in file

  integer(I4B) :: nside_in,nmap_in,npixtot_in,type,obsnpix_in, polar_in, n_ext
  logical :: first_time=.true.
  save first_time 
  
  integer(I8B) :: mygroup,myerr
  character(len=FILENAMELEN) :: grpname='' !SP
  character(len=FILENAMELEN), dimension(:), pointer :: flgnamePIO,maptypePIO,mapnamePIO
  character(len=FILENAMELEN) :: coordsys_in='',ordering_in='' !SP
  character(len=FILENAMELEN), save :: coordsys='' !SP


#ifdef PIO
  myerr   = PIOgetgrpname(grpname,file_name)
  mygroup = PIOopenmapgrp(grpname,"r")
  nside_in= PIOnsidegrp(mygroup)
  myerr   = PIOcoordsysgrp(coordsys_in,mygroup)
  myerr   = PIOorderinggrp(ordering_in,mygroup)
  myerr   = PIOclosemapgrp(mygroup) 
  npixtot_in=12*nside_in**2
  obsnpix_in=npixtot_in ! 2019-10
  nmap_in=1

  if (ordering_in=='RING') then
     ordering=1
  else
     ordering=2
  endif
#else  
  obsnpix_in=getsize_fits(file_name,nmaps=nmap_in,ordering=ordering, &
 &                        nside=nside_in,type=type,polarisation=polar_in) !EH 2003-06, 2019-10
  npixtot_in = 12*nside_in**2 ! cut sky file  !EH 2003-06
  if (type == 3) then
     n_ext = getnumext_fits(file_name)
     if (polar_in == 1) then
        if (n_ext == 1) nmap_in = 3      ! (Pixel,I,Q,U) cut sky map, 2020-01
        if (n_ext >  2) nmap_in = n_ext  ! (Pixel,Signal,N_obs,Serror) * 3 extensions , 2007-04-03
     else
        if (n_ext <  3) nmap_in = 1 ! unpolarised cut sky map
     endif
  endif
#endif


  if (dry) then
     nside_in=nint(sqrt(dble(npixtot_in)/12.d0))
     if (12*nside_in**2 /= npixtot_in) then
        write(*,*) 'ERROR in check header for file'
        write(*,*) trim(file_name)
        write(*,*) 'npixtot is not a multiple of 12*nside**2'
        write(*,*) 'Even with the dry option, I REFUSE to go further.'
        CALL FATAL_ERROR
     endif
     ordering=1
     nmap_in=1
     if (mapfile .and. polarization) nmap_in=3

     if (first_time) then
        if (verbose) then
           write(*,*)
           write(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
           write(*,*) 'WARNING in check_header:'
           write(*,*) 'You wrote the fits file like a pig so you are using'
           write(*,*) 'the dry option at your own risk.'
           write(*,*) 'Only LOUSY consistency check is done.'
           write(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
           write(*,*)
        endif
        first_time=.false.
     endif
  endif

  if (megaverbose) then
     if (dry) then
        write(*,*) 'Check approximately header for file '//trim(file_name)
     else
        write(*,*) 'Check header for file '//trim(file_name)
     endif
  endif

#ifdef PIO
     coordsys=coordsys_in
#endif

  if (.not. head_defined) then
     nside=nside_in
     npixtot=12*nside**2
     call set_nlmax(nlmax,nside)
!!$#ifdef PIO !SP
!!$     coordsys=coordsys_in !SP
!!$#endif !SP
     if (megaverbose) write(*,*) 'nside for input map file =',nside
     head_defined=.true.
  endif

  if (nmap_in > 1 .and. type==2 .and. verbose .and. ((mapfile .and. .not. polarization) &
 &    .or.(.not. mapfile))) then
     write(*,*)
     write(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
     write(*,*) 'WARNING in check_header:'
     write(*,*) 'I detect multiple maps in file'
     write(*,*) trim(file_name)
     write(*,*) 'Only the first map will be read.'
     write(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
     write(*,*)
  endif

#ifndef PIO
!  if (nmap_in /= 3 .and. mapfile .and. polarization) then
  if (mapfile .and. polarization) then
     if (nmap_in < 3) then
          write(*,*) 'ERROR in check_header for file'
          write(*,*) trim(file_name)
          write(*,*) 'You had requested polarization analysis'
          write(*,*) 'while the input file does not seem to contain'
          write(*,*) 'Stokes parameters maps'
          CALL FATAL_ERROR
     endif
     if (nmap_in > 3) then
          write(*,*) 'WARNING in check_header for file'
          write(*,*) trim(file_name)
          write(*,*) 'has more than 3 columns.'
          write(*,*) 'Will assume that first 3 are Stokes parameters maps.'
     endif
  endif
#endif

  if (nside_in /= nside) then
     write(*,*) 'ERROR in check_header for file'
     write(*,*) trim(file_name)
     write(*,*) 'The input file has wrong value of nside :'
     write(*,*) 'nside_in =',nside_in
     write(*,*) 'Other files have Nside = ',nside
     CALL FATAL_ERROR
  endif

  if (npixtot_in /= npixtot) then
     write(*,*) 'ERROR in check_header for file'
     write(*,*) trim(file_name)
     write(*,*) 'mismatch between npixtot and nside'
     write(*,*) 'npixtot_in =',npixtot_in
     write(*,*) 'expected npixtot =',npixtot,' for Nside =',nside_in
     CALL FATAL_ERROR
  endif

  obsnpix=obsnpix_in ! 2019-10
!   if (obsnpix_in /= obsnpix) then
!      write(*,*) 'ERROR in check_header for file'
!      write(*,*) trim(file_name)
!      write(*,*) 'obsnpix_in =',obsnpix_in
!      write(*,*) 'expected obsnpix =',obsnpix
!      CALL FATAL_ERROR
!   endif

  if (ordering < 1 .or. ordering > 2) then
     write(*,*) 'ERROR in check_header for file'
     write(*,*) trim(file_name)
     write(*,*) 'pixel ordering is not RING nor NESTED :'
     write(*,*) 'ordering =',ordering
     CALL FATAL_ERROR
  endif

  if (ordering == 2) then
     write(*,*) 'WARNING in check_header for file'
     write(*,*) trim(file_name)
     write(*,*) 'pixel ordering is NESTED :'
     write(*,*) 'data will be reordered before usage'
  endif

#ifdef PIO
  if (coordsys_in /= coordsys) then
     write(*,*) 'ERROR in check_header for file'
     write(*,*) trim(file_name)
     write(*,*) 'mismatch between coordinate systems'
     CALL FATAL_ERROR
  endif
#endif
end subroutine check_header

!=======================================================================
subroutine check_te_header(tenorm_in_file, verbose, megaverbose, nlmax, &
     thetamax, apodize, apodizetype, apodizesigma)
!=======================================================================
  use healpix_types
  use healpix_modules, only: fatal_error, getsize_fits, get_card, fits2cl
  implicit none

  character(len=*), intent(in) :: tenorm_in_file
  logical(LGT),     intent(in) :: verbose, megaverbose, apodize
  integer(I4B),     intent(in) :: nlmax, apodizetype
  real(DP),         intent(in) :: thetamax, apodizesigma
  integer(I4B)                 :: nlread, nlmax_, apodizetype_
  logical(LGT)                 :: apodize_
  real(DP )                    :: thetamax_, apodizesigma_
  character(len=80), dimension(1:200)   :: header
  real(DP), allocatable, dimension(:,:) :: buffer
  character(len=*), parameter :: code = 'check_te_header'

  header = ''
  nlread = getsize_fits(tenorm_in_file, mlpol=nlmax_)

  if (nlmax_ < nlmax) then
     write(*,*) 'ERROR in '//code//' for '//trim(tenorm_in_file)
     write(*,*) 'The input TE normalization function lmax is too small for current study'
     write(*,*) 'lmax found =',nlmax_,',   lmax needed =',nlmax
     CALL FATAL_ERROR
  endif
  if (nlmax_ > nlmax .and. verbose) then
     write(*,*) 'WARNING in '//code//' for '//trim(tenorm_in_file)
     write(*,*) 'The input TE normalization function lmax is too small for current study'
     write(*,*) 'Only a subset will be used.'
     write(*,*) 'lmax found =',nlmax_,',   lmax needed =',nlmax
  endif

  allocate(buffer(0:nlmax,1:1))
  call fits2cl(tenorm_in_file, buffer, nlmax, 1, header)
  deallocate(buffer)

  call get_card(header, 'APODIZE',  apodize_)
  if (apodize_ .neqv. apodize) then
     write(*,*) 'ERROR in '//code//' for '//trim(tenorm_in_file)
     write(*,*) 'Mismatched APODIZE: '
     write(*,*) 'found =',apodize_,',    expected =',apodize
     call fatal_error
  endif
  if (apodize_) then
     call get_card(header, 'APOTYPE',  apodizetype_)
     if (apodizetype_ .ne. apodizetype_) then 
        write(*,*) 'ERROR in '//code//' for '//trim(tenorm_in_file)
        write(*,*) 'Mismatched APODIZETYPE: '
        write(*,*) 'found =',apodizetype_,',    expected =',apodizetype
        call fatal_error
     endif
     call get_card(header, 'APOSIGMA', apodizesigma_)
     if (abs(apodizesigma_ / apodizesigma - 1.0_DP) > 1.d-4) then
        write(*,*) 'ERROR in '//code//' for '//trim(tenorm_in_file)
        write(*,*) 'Mismatched APODIZESIGMA: '
        write(*,*) 'found =',apodizesigma_,',    expected =',apodizesigma
        call fatal_error
     endif
  endif
  call get_card(header, 'THETAMAX', thetamax_)
  if (abs(thetamax_ / thetamax - 1.0_DP) > 1.d-4) then
     write(*,*) 'ERROR in '//code//' for '//trim(tenorm_in_file)
     write(*,*) 'Mismatched THETAMAX: '
     write(*,*) 'found =',thetamax_,',    expected =',thetamax
     call fatal_error
  endif

  return
end subroutine check_te_header

!=======================================================================
subroutine set_nlmax_vect(nlmax,nsides)
!=======================================================================
  use healpix_types
  implicit none
  integer(I4B), intent(inout) :: nlmax
  integer(I4B), intent(in), dimension(1:2)    :: nsides
  integer(I4b) :: nside

  nside = minval(nsides)
  call set_nlmax_scal(nlmax, nside)

end subroutine set_nlmax_vect
!=======================================================================
subroutine set_nlmax_scal(nlmax,nside)
!=======================================================================
  use healpix_types
  implicit none
  integer(I4B), intent(inout) :: nlmax
  integer(I4B), intent(in)    :: nside

  if (nlmax > 0) then
     nlmax = min(nlmax,3*nside-1)
  else
     nlmax = 3*nside-1
  endif
end subroutine set_nlmax_scal

!=======================================================================
subroutine check_header_asc(file_name, nsides, nlmax, ncor, head_defined, &
 &                          verbose, megaverbose, npixtots)
!=======================================================================
  ! check header of ASCII file  or  FITS ASCII table containing C(l)
  use healpix_types
  use misc_utils, only: fatal_error
  use head_fits,  only: get_card
#ifdef PIO
  use piolib
  use my_pio_routines, only: pio_objecttype
#else
  use fitstools, only: getsize_fits
#endif
  implicit none

  character(len=*), intent(in)    :: file_name
  logical,          intent(in)    :: verbose,megaverbose
  logical,          intent(inout) :: head_defined
  integer(I4B),     intent(inout) :: nlmax,ncor
  integer(I4B),     intent(inout), dimension(1:2) :: nsides
  integer(I4B),     intent(out),   dimension(1:2) :: npixtots
  character(len=*), parameter :: code='check_header_asc'

  integer(I4B) :: lunit=10
  character(len=FILENAMELEN) :: header
  character(len=22) :: stringjunk
  integer(I4B) :: nlmax_in,ncor_in,nside_in1,nside_in2
  character(len=80), dimension(:), allocatable :: headfits

  integer(I4B)     :: junk, ncol, polarisation, type
  character(len=8) :: starter
  logical          :: header_parsed
!  logical :: my_isfits
#ifdef PIO
  character(len=20) :: object_type
  integer(I8B) :: myerr,mygroup,nlmaxD,ncorD, nsideD
  integer(I8B) :: nsideD1,nsideD2,myerr0,myerr1,myerr2
  character(len=FILENAMELEN) :: comment='',grpname='' !SP
#endif

#ifdef PIO
  myerr   = PIOgetgrpname(grpname,file_name)
  object_type = trim(pio_objecttype(file_name))
  if (trim(object_type) == 'VECT') mygroup = PIOopenVECTGrp(grpname,'r')
  if (trim(object_type) == 'CL'  ) mygroup = PIOopenClGrp  (grpname,'r')
  myerr   = PIOreadkeywordgrp(nlmaxD,comment,'nlmax',mygroup)
  myerr   = PIOreadkeywordgrp(ncorD, comment,'ncor', mygroup)
  myerr0  = PIOreadkeywordgrp(nsideD, comment,'nside', mygroup)
  if (myerror0 < 0) then
     myerr1  = PIOreadkeywordgrp(nsideD1,comment,'nside1',mygroup)
     myerr2  = PIOreadkeywordgrp(nsideD2,comment,'nside2',mygroup)
  else
     nsideD1 = nsideD
     nsideD2 = nsideD
  endif
  if (trim(object_type) == 'VECT') myerr   = PIOcloseVECTGrp(mygroup)
  if (trim(object_type) == 'CL'  ) myerr   = PIOcloseClGrp  (mygroup)
  nlmax_in  = nlmaxD
  ncor_in   = ncorD
  nside_in1 = nsideD1
  nside_in2 = nsideD2

#else
  if (my_isfits(file_name)) then
     junk = getsize_fits(file_name, nside=nside_in1, mlpol=nlmax_in, &
          &             nmaps=ncol, polarisation=polarisation, type=type)
     nside_in2 = nside_in1
     if (nside_in1 < 0) then
        allocate(headfits(1:200))
        call get_fits_header(file_name, headfits)
        call get_card(headfits, 'NSIDE1', nside_in1)
        call get_card(headfits, 'NSIDE2', nside_in2)
        deallocate(headfits)
     endif
     if (type /= 1) then
        write(*,*) trim(file_name)//' is not an ASCII FITS file.'
        call fatal_error
     endif
     if (polarisation == 1) then
        ncor_in = ncol
        ! correlation file : one extra column
!        if (ncol == 2) ncor_in = 1 ! unpolarized
        if (ncol == 5) ncor_in = 4  ! polarized
!        if (ncol == 7) ncor_in = 6 ! polarized + B coupling
     else
        ncor_in = 1
     endif
  else
     nside_in2 = -1
     open(file=file_name,unit=lunit,form='formatted',status='old',action='read')

     header_parsed = .false. ! 1st attempt at header parsing
     read(lunit,'(A)',          err=1, end=1)  header
     read(header,'(A22,4(I12))',err=10,end=10) stringjunk,nlmax_in,ncor_in,nside_in1,nside_in2
     header_parsed = .true.
10   if (.not. header_parsed) then ! 2nd attempt at header parsing
        read(header,'(A22,3(I12))',err=1, end=1 ) stringjunk,nlmax_in,ncor_in,nside_in1
        !!!!!nside_in1 = nside_in2 ! corrected 2020-04-08
        nside_in2 = nside_in1
        header_parsed = .true.
     endif
     if (nside_in2 <= 0) nside_in2 = nside_in1 ! added 2020-04-10
1    close(lunit)
     if (.not. header_parsed) then ! header parsing failed
        write(*,*) 'ERROR while reading/parsing '//trim(file_name)
        write(*,*) 'Is it really an ASCII file generated by Spice?'
        call fatal_error
     endif
  endif
#endif

  if (head_defined) then
     if (nside_in1 /= nsides(1)) then
        write(*,*) 'ERROR in '//code
        write(*,*) 'file : '//trim(file_name)
        write(*,*) 'The value of nside corresponding to this file is inconsistent'
        write(*,*) 'with the one of the input map (1).'
        write(*,*) 'nside of the file :',nside_in1
        write(*,*) 'nside of the input map :',nsides(1)
        CALL FATAL_ERROR
     endif
     if (nside_in2 /= nsides(2)) then
        write(*,*) 'ERROR in '//code
        write(*,*) 'file : '//trim(file_name)
        write(*,*) 'The value of nside corresponding to this file is inconsistent'
        write(*,*) 'with the one of the input map (2).'
        write(*,*) 'nside of the file :',nside_in2
        write(*,*) 'nside of the input map :',nsides(2)
        CALL FATAL_ERROR
     endif
     if (nlmax_in < nlmax) then
        write(*,*) 'ERROR in '//code
        write(*,*) 'file : '//trim(file_name)
        write(*,*) 'This file has the following value of nlmax :',nlmax_in
        write(*,*) 'This value is inconsistent with the current value of nlmax.'
        write(*,*) 'We indeed have nlmax =',nlmax
        CALL FATAL_ERROR
     elseif (nlmax_in > nlmax .and. verbose) then
        write(*,*) 'WARNING in '//code
        write(*,*) 'file : '//trim(file_name)
        write(*,*) 'This file has the following value of nlmax :',nlmax_in
        write(*,*) 'This value is larger than the chosen value of nlmax.'
        write(*,*) 'We indeed have nlmax =',nlmax
     endif
     ! if ncor <= 0, adapt it to the one in file
     if (ncor <= 0) ncor = ncor_in
     ! otherwise, be without mercy
     if (ncor_in < ncor) then
        if (ncor_in /= 1) then
           write(*,*) 'ERROR in '//code           
           write(*,*) 'Inconsitent header in '
           write(*,*) 'file : '//trim(file_name)
           write(*,*) 'ncor read :',ncor
           CALL FATAL_ERROR
        endif
        if (ncor == 4) then
           write(*,*) 'ERROR in '//code           
           write(*,*) 'This file corresponds to unpolarized data, while you have'
           write(*,*) 'asked for a calculation including polarization.'
           CALL FATAL_ERROR
        endif
     elseif (ncor_in > ncor) then
        if (ncor_in /= 4) then
           write(*,*) 'ERROR in '//code           
           write(*,*) 'Inconsitent header in '
           write(*,*) 'file : '//trim(file_name)
           write(*,*) 'ncor read :',ncor
           CALL FATAL_ERROR
        endif
        if (verbose) then
           write(*,*) 'WARNING in '//code
           write(*,*) 'file : '//trim(file_name)
           write(*,*) 'This file corresponds to polarized data, while you have'
           write(*,*) 'asked for a calculation without polarization.'
        endif
     endif
  else
     nsides(1)   = nside_in1
     nsides(2)   = nside_in2
     npixtots(:)= 12*nsides(:)**2
     call set_nlmax(nlmax,nsides)
     if (nsides(1) == nsides(2)) then
        if (megaverbose) write(*,*) 'nside for input map file =',nsides(1)
     else
        if (megaverbose) write(*,*) 'nsides for input map file =',nsides(1:2)
     endif
     head_defined=.true.
  endif
  return

end subroutine check_header_asc

!=======================================================================
subroutine check_open_output_files
!=======================================================================
  use spice_common
  use misc_utils, only: fatal_error
  implicit none

  integer(I4B) :: lunit=10
  character(len=filenamelen) :: xx_file_name
  character(len=10) :: status

  if (overwrite) then
     status = 'unknown'
  else
     status='new'
  endif
  if (coroutput) then
     xx_file_name = trim(adjustl(cor_file_name))
     if (xx_file_name(1:1) == '!') then
        xx_file_name=xx_file_name(2:filenamelen)
        status = 'unknown'
     endif
     open(file=xx_file_name,unit=lunit,form='formatted',status=status,err=1)
     close(lunit)
  endif

  if (cloutput) then
     xx_file_name = trim(adjustl(cl_file_name))
     if (xx_file_name(1:1) == '!') then
        xx_file_name=xx_file_name(2:filenamelen)
        status = 'unknown'
     endif
     open(file=xx_file_name,unit=lunit,form='formatted',status=status,err=2)
     close(lunit)
  endif

  if (clmapoutput) then
     xx_file_name = trim(adjustl(cl_outmap_file))
     if (xx_file_name(1:1) == '!') then
        xx_file_name=xx_file_name(2:filenamelen)
        status = 'unknown'
     endif
     open(file=xx_file_name,unit=lunit,form='formatted',status=status,err=3)
     close(lunit)        
  endif

  if (clmaskoutput) then
     xx_file_name = trim(adjustl(cl_outmask_file))
     if (xx_file_name(1:1) == '!') then
        xx_file_name=xx_file_name(2:filenamelen)
        status = 'unknown'
     endif
     open(file=xx_file_name,unit=lunit,form='formatted',status=status,err=4)
     close(lunit)        
  endif

  return

1 write(*,*) 'ERROR in check_open_output_files'
  write(*,*) 'Impossible to create file'
  write(*,*) trim(cor_file_name)
  CALL FATAL_ERROR

2 write(*,*) 'ERROR in check_open_output_files'
  write(*,*) 'Impossible to create file'
  write(*,*) trim(cl_file_name)
  CALL FATAL_ERROR

3 write(*,*) 'ERROR in check_open_output_files'
  write(*,*) 'Impossible to create file'
  write(*,*) trim(cl_outmap_file)
  CALL FATAL_ERROR

4 write(*,*) 'ERROR in check_open_output_files'
  write(*,*) 'Impossible to create file'
  write(*,*) trim(cl_outmask_file)
  CALL FATAL_ERROR
end subroutine check_open_output_files

!=======================================================================
subroutine write_corcl_file(data, mu, nlmax, file_name, cor_file, ncor, nsides)
!=======================================================================
  use healpix_types
#ifdef PIO
  use piolib
  use piolib_tuple
  use spice_parameters, only : decouple, mytype_cl_tuple, mytype_vect_tuple
  use spice_common, only : out_cl_grp ! EH
  use misc_utils, only: string, assert, fatal_error
#else
  use spice_parameters, only : decouple
  use spice_common, only : fits_out !, &
       !mapfile, maskfile, weightfile, maskfilep, weightfilep, &
       !mapfile2, maskfile2, weightfile2, maskfilep2, weightfilep2, &
       !tf_file
  use fitstools, only: write_asctab
  use head_fits, only: add_card
  use misc_utils, only: fatal_error
#endif
  implicit none

  integer(I4B), intent(in) :: nlmax,ncor
  integer(I4b), intent(in), dimension(1:2) :: nsides
  real(DP), intent(in), dimension(0:nlmax,1:ncor), target :: data
!  data(0:nlmax,1:ncor),mu(nlmax+1)
  real(DP), intent(in), dimension(1:nlmax+1) :: mu
  character(len=*), intent(in) :: file_name
  logical, intent(in) :: cor_file

  integer(I4B) :: l,i,nh,iw
  integer(I4B) :: lunit=10
  character(len=FILENAMELEN) :: fn_bang, fn_nobang
!  character(len=10) :: status
  character(len=filenamelen) :: header

  character(len=80), dimension(1:120)            :: headfits
  real(DP), allocatable, dimension(:,:), target :: tmpdata

#ifdef PIO
  character(len=FILENAMELEN) :: command='',grpname='',object='' !SP
  integer(I8B) :: nlmaxD,ncorD,nsideD1,nsideD2  !,nlmaxDp1
  integer(I8B) :: myerr,mygroup
  real(DP),allocatable,dimension(:) :: cos_tab
  character*2 :: theta_exten='TH'
  character*2 :: npairs_exten='NP'
  character*2, dimension(1:9) :: all_exten
  character(len=DMCPIOSTRINGMAXLEN) :: decoupleS

  type(PIOPtrOnArrayDble), dimension(:), pointer        :: InData 
  character(len=512) :: line
  integer(PIOSHORT) :: n2_short 
  integer(I4B) :: mytype
  character(len=8) :: stype
#endif

  !-----------------------------------------------------------------------
  if (ncor > 9) then
     print*,'ncor (number of C(l) or C(theta) to write) is larger than 9: ', ncor
     print*,'Aborting'
     call fatal_error
  endif

  fn_bang = trim(adjustl(file_name))
  if (fn_bang(1:1) == '!') then 
     fn_nobang = fn_bang(2:filenamelen)  !  no leading '!' (for PIOlib/DMC)
  else
     fn_nobang = fn_bang            !  no leading '!' (for PIOlib/DMC)
     fn_bang = '!'//trim(fn_bang)   !  with leading '!' (for FITS)
  endif

#ifdef PIO  
  nlmaxD=nlmax
  ncorD=ncor
  nsideD1=nside(1)
  nsideD2=nside(2)
!  nlmaxDp1=nlmaxD+1_I8B


  if (PIOCheckTemporary(fn_nobang) == 0) then
     ! if temporary file, create default group
     print*,trim(fn_nobang)
     print*,'temporary output file. Use default group'
     grpname = out_cl_grp
  else
     ! otherwise, get group from full path name
     myerr=PIOgetgrpname(grpname,fn_nobang)
  endif
  print*, '<'//trim(grpname)//'>'


  mytype = mytype_cl_tuple
  write(command,'(''begin='',I6,'';end='',I6)') 0, nlmax

  if (mytype == mytype_cl_tuple .or. mytype == mytype_vect_tuple) then
     ! ----------------  write in Tuple C(l) ----------------
     if (mytype == mytype_cl_tuple) then
        myerr   = PIOCreateClGrp(grpname)
        mygroup = PIOOpenClGrp(grpname,'w')
        stype = 'CL'
     else
        myerr   = PIOCreateVectGrp(grpname)
        mygroup = PIOOpenVectGrp(grpname,'w')
        stype = 'VECT'
     endif
     ! group keywords
     myerr=PIOWriteKeywordGrp(nlmaxD,'Max_Multipole',                'nlmax',mygroup)
     myerr=PIOWriteKeywordGrp(ncorD, 'Number_of_Fields',             'ncor', mygroup)
     if (nsideD1 == nsideD2) then
        myerr=PIOWriteKeywordGrp(nsideD1,'Nside_Parameter_of_Input_Map', 'nside',mygroup)
     else
        myerr=PIOWriteKeywordGrp(nsideD1,'Nside_Parameter_of_Input_Map_1', 'nside1',mygroup)
        myerr=PIOWriteKeywordGrp(nsideD2,'Nside_Parameter_of_Input_Map_2', 'nside2',mygroup)
     endif
     ! object keywords
     decoupleS = 'FALSE'
     if (decouple) decoupleS = 'TRUE'
     myerr=PIOWriteKeywordObject(decoupleS  ,'EB_Decoupled','decouple', trim(fn_nobang),mygroup)

     ! prepare data
     if (cor_file) then
        print*,'write cor_file'
        ! write theta and power spectra
        allocate(InData(1:ncor+1))
        allocate(tmpdata(0:nlmax,1:1))
        tmpdata(0:nlmax,1)       = dacos(mu(1:nlmax+1)) ! theta
        InData(1)%IdxSple        => tmpdata(0:nlmax,1)
        do iw = 1, ncor
           InData(iw+1)%IdxSple  => data(0:nlmax,iw)
        enddo
        n2_short = ncor + 1
     else
        print*,'write cl_file'
        ! write power spectra
        allocate(InData(ncor))
        do iw = 1, ncor
           InData(iw)%IdxSple => data(0:nlmax, iw)
        enddo
        n2_short = ncor
     endif

     ! create tuple object
     myerr = PIOCreateTupleObject(trim(fn_nobang), 'PIODOUBLE', n2_short, mygroup, stype, stype) 
     line = 'creating tuple: '//PIOErrMess(myerr)
     call assert(myerr >= 0, line)

     ! write tuple object
     if (mytype == mytype_cl_tuple) then
        myerr = PIOWriteClTupleObject(InData, trim(fn_nobang), command, mygroup)
     else
        myerr = PIOWriteVectTupleObject(InData, trim(fn_nobang), command, mygroup)
     endif
     line = 'writing tuple: '//PIOErrMess(myerr)
     call assert(myerr >= 0, line)

     if (allocated(tmpdata)) deallocate(tmpdata)
     deallocate(InData)
     ! close group
     if (mytype == mytype_cl_tuple) then
        myerr = PIOcloseClGrp(mygroup)
     else
        myerr = PIOcloseVectGrp(mygroup)
     endif

  else
     ! ---------------- write in flat Vect -------------------
     if (decouple) then
        all_exten = (/ '  ','E ','B ','TE','TB','EB', 'ET', 'BT', 'BE' /)
     else
        all_exten = (/ '  ','Q ','U ','TQ','TU','QU', 'QT', 'UT', 'UQ' /)
     endif
     write(command,'(''begin='',I6,'';end='',I6)') 0,nlmax

     myerr=PIOcreateVECTGrp(grpname) 
     mygroup=PIOopenVECTGrp(grpname,'w')

     myerr=PIOwritekeywordgrp(nlmaxD,' ','nlmax',mygroup)
     myerr=PIOwritekeywordgrp(ncorD, ' ','ncor', mygroup)
     if (nsideD1 == nsideD2) then
        myerr=PIOwritekeywordgrp(nsideD1,' ','nside',mygroup)
     else
        myerr=PIOwritekeywordgrp(nsideD1,' ','nside1',mygroup)
        myerr=PIOwritekeywordgrp(nsideD2,' ','nside2',mygroup)
     endif

     if (cor_file) then
        print*,'write cor_file'
        ! write correlation functions
        do iw = 1, ncor
           object=trim(fn_nobang)//trim(all_exten(iw))
           myerr=PIOcreateVECTObject(object,'PIODOUBLE',mygroup)
           myerr=PIOwriteVECTObject(data(0:nlmax,iw),object,command,mygroup)
        enddo
        ! write cos(theta)
        object=trim(fn_nobang)//theta_exten
        myerr=PIOcreateVECTObject(object,'PIODOUBLE',mygroup)
        myerr=PIOwriteVECTObject(mu(1:nlmax+1),object,command,mygroup)
        ! write theta
        object=trim(fn_nobang)//npairs_exten
        myerr=PIOcreateVECTObject(object,'PIODOUBLE',mygroup)
        allocate(cos_tab(0:nlmax)) !SP
        do l=0,nlmax
           cos_tab(l)=dacos(mu(l+1))
        enddo
        myerr=PIOwriteVECTObject(cos_tab(0:nlmax),object,command,mygroup)
        deallocate(cos_tab)
     else
        print*,'write cl_file'
        ! write power spectra
        do iw = 1, ncor
           object=trim(fn_nobang)//trim(all_exten(iw))
           myerr=PIOcreateVECTObject(object,'PIODOUBLE',mygroup)
           myerr=PIOwriteVECTObject(data(0:nlmax,iw),object,command,mygroup)
        enddo
     endif
     ! close group
     myerr = PIOcloseVECTGrp(mygroup)
  endif
#else
  !----------------- Non Piolib ---------------------
  if (fits_out) then 
     ! FITS file
     ! ---- create FITS header --------
     headfits = ''
     if (cor_file) then
        call add_card(headfits,'EXTNAME',"'ANGULAR CORRELATION'")
        call add_card(headfits,'TTYPE1',"ANGLE","Angular Separation")
        call add_card(headfits,'TUNIT1',"RAD")
        call add_card(headfits,'TTYPE2',"TT","Temperature Angular Correlation")
        if (ncor > 1) then
           call add_card(headfits,'TTYPE3',"QQ")
           call add_card(headfits,'TTYPE4',"UU")
           call add_card(headfits,'TTYPE5',"TQ")
        endif
        if (ncor > 4) then
           call add_card(headfits,'TTYPE6',"TU")
           call add_card(headfits,'TTYPE7',"QU")
        endif
        if (ncor > 6) then
           call add_card(headfits,'TTYPE8',"QT")
           call add_card(headfits,'TTYPE9',"UT")
           call add_card(headfits,'TTYPE10',"UQ")
        endif
     else
        call add_card(headfits,'EXTNAME',"'ANGULAR POWER SPECTRUM'")
        call add_card(headfits,'TTYPE1',"TT")
        if (ncor > 1) then
           call add_card(headfits,'TTYPE2',"EE")
           call add_card(headfits,'TTYPE3',"BB")
           call add_card(headfits,'TTYPE4',"TE")
        endif
        if (ncor > 4) then
           call add_card(headfits,'TTYPE5',"TB")
           call add_card(headfits,'TTYPE6',"EB")
        endif
        if (ncor > 6) then
           call add_card(headfits,'TTYPE7',"ET")
           call add_card(headfits,'TTYPE8',"BT")
           call add_card(headfits,'TTYPE9',"BE")
        endif
     endif
     call write_spice_header(headfits, size(headfits), nlmax, ncor, nsides, .true.)

     ! ---- write FITS file --------
     nh = size(headfits)
     if (cor_file) then
        allocate(tmpdata(0:nlmax,1:ncor+1))
        tmpdata(0:nlmax,1)        = dacos(mu(1:nlmax+1)) ! theta
        tmpdata(0:nlmax,2:ncor+1) = data(0:nlmax,1:ncor)
        call write_asctab(tmpdata, nlmax, ncor+1, headfits, nh, fn_bang)
        deallocate(tmpdata)
     else
        call write_asctab(data, nlmax, ncor, headfits, nh, fn_bang)
     endif

  else
     ! ordinary ASCII file
     open(file=fn_nobang,unit=lunit,form='formatted',status='unknown',action='write')

     if (nsides(1) == nsides(2)) then
        write(header,'(A22,3(I12))') '# nlmax, ncor, nside =',&
             &                          nlmax, ncor, nsides(1)
     else
        !write(header,'(A22,4(I12))') '# nlmax, ncor, nside1, nside2 =',&
        write(header,'(A22,4(I12))') '# nlmax, ncor, nsides=',&
             &                          nlmax, ncor, nsides(1),nsides(2)
     endif
     write(lunit,'(A)') trim(header)
     do l=0,nlmax
        if (cor_file) then
           if (ncor == 1) then
              write(lunit, '(3(E24.16))') dacos(mu(l+1)), mu(l+1), data(l,1)
           else
              write(lunit, '(11(E24.16))') dacos(mu(l+1)), mu(l+1), (data(l,i),i=1,ncor)
           endif
        else
           if (ncor == 1) then
              write(lunit,'(I5,1X,E24.16)') l, data(l,1)
           else
              write(lunit,'(I5,1X,9(E24.16))') l, (data(l,i),i=1,ncor)
           endif
        endif
     enddo
  
     close(lunit)
  endif
#endif
end subroutine write_corcl_file

!=======================================================================
subroutine read_corcl_file(data, nlmax, file_name, cor_file, ncor)
!=======================================================================
!   routine to read C(l) files, C(theta) (correlation) files
!   or transfer function files
!
!=======================================================================

  use healpix_types
  use misc_utils, only: fatal_error
#ifdef PIO
  use piolib
  use piolib_tuple
  use my_pio_routines, only: pio_objecttype
  use spice_parameters, only : decouple
#else
  use fitstools, only: fits2cl
#endif
  implicit none

  integer(I4B),     intent(in)  :: nlmax, ncor
  real(DP),         intent(out) :: data(0:nlmax,1:ncor)
  character(len=*), intent(in)  :: file_name
  logical,          intent(in)  :: cor_file

  real(DP), dimension(1:ncor) ::  data_in
  real(DP)     :: xjunk1,xjunk2
  integer(I4B) :: lin,i,nlmax_in,ncor_in,nside_in,iw
  integer(I4B) :: lunit=10
  character(len=FILENAMELEN) :: header
  character(len=15) :: stringjunk

  integer(I8B) :: lread

  character(len=80), dimension(1:60) :: headfits
  real(DP), allocatable, dimension(:,:) :: tmpdata
!  logical :: my_isfits
#ifdef PIO
  character(len=FILENAMELEN) :: object
  character(len=20) :: object_type
  integer(PIOLONG) :: myerr, mygroup
  integer(PIOSHORT) :: ncor_short
  character(len=DMCPIOSTRINGMAXLEN) :: grpname
  character(len=20) :: command
  type(PIOPtrOnArrayDble), dimension(:), pointer  :: OutData 

  real(DP),pointer,dimension(:) :: data_in_arr
  character(len=2) :: theta_exten='TH'
  character(len=2) :: npairs_exten='NP'
  character(len=2), dimension(1:9) :: all_exten
#endif

  !-----------------------------------------------------------------------
  if (ncor > 9) then
     print*,'ncor (number of C(l) or C(theta) to read) is larger than 9: ', ncor
     print*,'Aborting'
     call fatal_error
  endif

  data(0:nlmax,1:ncor)=0._DP

#ifdef PIO
  if (decouple) then
     all_exten = (/ '  ','E ','B ','TE','TB','EB', 'ET', 'BT', 'BE' /)
  else
     all_exten = (/ '  ','Q ','U ','TQ','TU','QU', 'QT', 'UT', 'UQ' /)
  endif
  command = ' '

  ! identify object type (VECT vs CL, TUPLE vs FLAT) and open it
  myerr       = PIOGetGRPName(grpname,file_name)
  object_type = trim(pio_objecttype(file_name))
  if (trim(object_type) == 'VECT') mygroup     = PIOOpenVECTGrp(grpname,'r')
  if (trim(object_type) == 'CL'  ) mygroup     = PIOOpenCLGrp(  grpname,'r')
  myerr = PIOGetNbVecTuple(ncor_short, file_name, mygroup)

  if (ncor_short > 1) then
     ! ---- read tuple CL or VECT object ----
     if (ncor /= ncor_short) then
        print*,'Was expecting ',ncor,' fields in '//trim(file_name)
        print*,'Found ',ncor_short,' in tuple.'
        call fatal_error
     endif
     if (trim(object_type) == 'CL'  ) myerr = PIOReadClTupleObject  (OutData, trim(file_name), command, mygroup)
     if (trim(object_type) == 'VECT') myerr = PIOReadVectTupleObject(OutData, trim(file_name), command, mygroup)
     lread = myerr - 1
     call check_mylread(lread,trim(file_name),nlmax)
     do iw = 1, ncor
        data(0:lread, iw) = OutData(iw)%IdxSple
     enddo
     myerr = PIODeleteTupleTable(OutData, mygroup)
  else
     ! ---- read flat CL or VECT object ----
     do iw = 1, ncor
        object=trim(file_name)//trim(all_exten(iw))
        if (trim(object_type) == 'VECT') myerr=PIOreadVECTObject(data_in_arr, object, command)
        if (trim(object_type) == 'CL')   myerr=PIOreadCLObject  (data_in_arr, object, command)
        lread=myerr-1
        call check_mylread(lread,object,nlmax)
        data(0:lread,iw)=data_in_arr(1:lread+1)
        myerr=PIOdeleteVECTTable(data_in_arr)
        !      if (object_type == 'VECT') myerr=PIOdeleteVECTTable(data_in_arr)
        !      if (object_type == 'CL'  ) myerr=PIOdeleteCLTable(  data_in_arr)
     enddo
  endif

#else
  ! Non PIOlib
  if (my_isfits(file_name)) then 
     ! FITS file
     if (cor_file) then
        allocate(tmpdata(0:nlmax, 1:ncor+1))
        call fits2cl(file_name, tmpdata, nlmax, ncor+1, headfits)
        data(0:nlmax, 1:ncor) = tmpdata(0:nlmax, 2:ncor+1)
        deallocate(tmpdata)
     else
        call fits2cl(file_name, data, nlmax, ncor, headfits)
     endif
  else
     ! ordinary ASCII file
     open(file=file_name,unit=lunit,form='formatted',status='old')
     read(lunit,'(A)',err=1,end=2) header
     read(header,'(A22,3(I12))',err=1,end=2) stringjunk,nlmax_in,ncor_in,nside_in
     
     do lread=0,nlmax
        if (cor_file) then
           read(lunit,*,err=1,end=2) xjunk1,xjunk2,(data_in(i),i=1,ncor)
        else
           read(lunit,*,err=1,end=2) lin,(data_in(i),i=1,ncor)
        endif
        data(lread,1:ncor)=data_in(1:ncor)
     enddo

     close(lunit)
  endif
#endif

  return

1 write(*,*) 'ERROR in read_corcl_file'
  write(*,*) 'file '//trim(file_name)
  write(*,*) 'seems to be corrupted.'
  close(lunit)
  CALL FATAL_ERROR

2 write(*,*) 'ERROR in read_corcl_file'
  write(*,*) 'file '//trim(file_name)
  write(*,*) 'seems to be truncated.'
  close(lunit)
  CALL FATAL_ERROR

end subroutine read_corcl_file

!=======================================================================
subroutine check_mylread(lread,file_name,nlmax, warn_only)
!=======================================================================
  ! checks number of multipoles read from PIOLIB/DMC C(l) file
  !---------------------------------------------------------------------
  use healpix_types
  use misc_utils, only: fatal_error
  implicit none
  integer(I8B) :: lread
  integer(I4B) :: nlmax
  character(len=*) :: file_name
  logical(LGT), optional :: warn_only
  logical(LGT) :: crash

  crash = .true.
  if (present(warn_only)) crash = .not. warn_only

  if (lread < 0) then
     write(*,*) 'ERROR while reading C(l)/correlation/beam file: '
     write(*,*) 'file '//trim(file_name)
     write(*,*) 'This file is not valid'
     CALL FATAL_ERROR
  elseif (lread > nlmax) then
     if (crash) then
        write(*,*) 'ERROR while reading C(l)/correlation/beam file: '
     else
        write(*,*) 'WARNING while reading C(l)/correlation/beam file: '
     endif
     write(*,*) 'file '//trim(file_name)
     write(*,'(i12,a3,i12)') lread, ' > ',nlmax
     write(*,*) 'This file corresponds to a larger nlmax'
     if (crash) CALL FATAL_ERROR
  endif
end subroutine check_mylread

!=======================================================================
subroutine generate_gaussian_beam(fwhm,nlmax,gb)
!=======================================================================
  ! this routine was moved from deal_with_xi_and_cl.f90 to here, for clarity
  ! 2013-03-28: also returns polarized B(l)
  !======================================================================
  use healpix_types
  implicit none
  real(DP),     intent(in) :: fwhm
  integer(I4B), intent(in) :: nlmax
  real(DP), dimension(0:,1:), intent(out) :: gb

  real(DP) :: arcmin2rad,sigma2fwhm,sigma, fact_pol
  integer(I4B) :: l, nd
  !======================================================================

  nd   = size(gb,2)
  arcmin2rad = PI / (180.0_dp * 60.0_dp)
  sigma2fwhm = sqrt(8.0_dp * log(2.0_dp))

  sigma    = fwhm * arcmin2rad / sigma2fwhm ! in radians
  fact_pol = exp(2.0_dp*sigma**2) ! correction for polarised fields

  ! temperature
  do l=0,nlmax
     gb(l,1) = exp(-0.5_dp * l*(l+1.0_dp) * sigma**2)
  enddo    
  ! electric or gradient
  if (nd > 1) gb(0:nlmax,2) = gb(0:nlmax,1) * fact_pol
  ! magnetic or curl
  if (nd > 2) gb(0:nlmax,3) = gb(0:nlmax,1) * fact_pol

end subroutine generate_gaussian_beam

!=======================================================================
subroutine my_generate_beam(lmax, gb, megaverbose, beam_file)
!=======================================================================
! Hacked from generate_beam of Healpix 1.2, modified to take into account verbose mode
! This routine should be cleaned up.
! adapted to Healpix 2.0, can read DMC file
! 2013-03-28: adapted to Healpix 3.0
!    bref is not the same as in generate_beam of Healpix 3.10
! 2017-03-07: can read multi-column plain ASCII file (as in Healpix 3.4)
!==========================================================================
  use healpix_types
  use misc_utils, only: fatal_error, string
#ifdef PIO
  use piolib
  use my_pio_routines, only: pio_objecttype
#else
  use fitstools, only : fits2cl, getsize_fits
#endif
  implicit none
  real(kind=DP), dimension(0:,1:), intent(out) :: gb
  integer(kind=I4B), intent(in) :: lmax
  character(len=*),  intent(in) :: beam_file
  logical(LGT),      intent(in) :: megaverbose
  
  integer(kind=i4b) :: nsize, type, nlheader, nl, nd, lunit, il, i, junk
  integer(i4b) :: l100, iline, nc
  character(len=80), dimension(1:180) :: header
  character(len=1600) :: str
  character(len=80) :: card
!  logical :: my_isfits

  real(DP),pointer,dimension(:) :: tmpvectDP
  integer(I8B) :: myerr
  character(len=20) :: object_type
  real(DP):: bref
  real(dp), dimension(1:3) :: buffer
  character(len=*), parameter :: mycode = 'my_generate_beam'
  !==========================================================================
  ! test if name of external is given and valid

  nl = size(gb, 1)
  nd = size(gb, 2)
  gb = 0.0_dp
  
  if (nl <= lmax.and.megaverbose) then
     print 9000, 'WARNING in '//mycode//':'
     print 9001, 'beam array only available up to ',nl
  endif
  nl = min(nl, lmax+1)
  
#ifdef PIO
 ! ------- piolib case ------------
  myerr = -1
  object_type = trim(pio_objecttype(beam_file))
  if (object_type == 'VECT') myerr=PIOreadVECTObject(tmpvectDP,beam_file,' ')
  if (object_type == 'CL')   myerr=PIOreadCLObject  (tmpvectDP,beam_file,' ')
  call check_mylread(myerr, beam_file, nl, warn_only=.true.)
  gb(0:nl-1,1) = tmpvectDP(1:nl)
  myerr=PIODeleteVECTTable(tmpvectDP)
!   if (object_type == 'VECT') myerr=PIOdeleteVECTTable(tmpvectDP)
!   if (object_type == 'CL'  ) myerr=PIOdeleteCLTable(  tmpvectDP)
#else
 !------------ Non piolib case ----------
  lunit = 15
  
  ! read file according to its type
  if (.not.my_isfits(beam_file)) then 
     ! ordinary ascii file ?
     lunit = 32
     iline = 0 ! line index
     open(unit=lunit,file=beam_file,status='old', &
          &          form='formatted',action='read')
     do
        iline = iline+1
        read(lunit,'(a)', end=100, err=100) str
        if (len_trim(str)>0 .and. str(1:1) /= '#') then
           buffer = 0_dp
           
           ! try 3 columns
3          continue
           read(str,*, err=2, end=2) il, buffer(1:3)
           nc = 3
           goto 99
           
           ! try 2 columns
2          continue
           read(str,*, err=1, end=1) il, buffer(1:2)
           nc = 2
           goto 99
           
           ! try 1 column
1          continue
           read(str,*, err=666, end=666) il, buffer(1)
           nc = 1
           goto 99
           
           ! failure
666        continue
           print 9000, 'ERROR in '//mycode
           print 9000, 'Unable to parse ASCII file '//trim(beam_file)
           print 9000, 'at line '//trim(adjustl(string(iline)))//' containing'
           print 9000, trim(str)
           print 9000, 'Expected format :'
           print 9000, ' Multipole,  column_1 [, column_2, column_3]'
           print 9000, ' with OPTIONAL colum_2 and column_3'
           call fatal_error
           
           ! record data, keep going until end of file
99         continue
           gb(il, 1:nd) = buffer(1:nd)
           if (il == nl-1) exit
        endif
     enddo
100  continue
     close(lunit)
     if (il < (nl-1).and.megaverbose) then
        print 9000, 'WARNING in '//mycode//' :'
        print 9001, 'Beam transfer function only available up to l= ',il
        print 9000, 'The larger multipoles will be set to 0'
     endif
       
  else 
     junk = getsize_fits(beam_file,type=type,nmaps=nc)
     if (type == 1 .or. type == 2) then ! ASCII or BINARY fits table
        ! FITS file with ascii table
        call fits2cl(beam_file, gb, nl-1, nd, header, fmissval=0.0_dp)
     else
        print 9000, 'ERROR in '//mycode
        print 9000, 'the file '//trim(beam_file)//' is of unknown type,'
        print 9000, ' or does not exist.'
        call fatal_error
     endif
  endif

  if (nc < nd) then
     if (megaverbose) then
        print 9000,'WARNING in '//mycode//' :'
        print 9000,'Not enough columns found in '//trim(beam_file)
        print 9000,'Expected '//trim(adjustl(string(nd)))&
             //', found '//trim(adjustl(string(nc)))
     endif
     do i=nc+1, nd
        if (megaverbose) then
           print 9000,' column #'//trim(adjustl(string(i))) &
                &                //' empty, fill in with column #' &
                &                //trim(adjustl(string(i-1)))
        endif
        gb(:,i) = gb(:,i-1)
     enddo
  else
     if (megaverbose) then
        print 9000,'Read '//trim(adjustl(string(nd)))&
             //' columns from '//trim(beam_file)
        print 9000,'(out of '//trim(adjustl(string(nc)))//')'
     endif
  endif

#endif
 !--------- end ----------
9000 format(a)
9001 format(a,i6)
9002 format(a,i3,a,i3)

  return
end subroutine my_generate_beam


!=======================================================================
function my_isfits(filename) result(fits)
  !---------------------------------------------------------------------
  use misc_utils, only: file_present, fatal_error

  character(len=*), intent(in) :: filename
  logical :: fits, found_unix, found_extended
  !
  integer :: lunit
  character(len=8) :: card

  found_extended = file_present(filename) ! test existence of physical or virtual file
  inquire(file=filename, exist=found_unix) ! test existence of physical file

  if (.not.found_extended) then
     call fatal_error(trim(filename)//' not found.')
  endif

  if (found_unix) then
     ! real file: try to read its first line
     lunit = 100
     open(unit=lunit, file=filename, status='old', form='formatted', action='read')
     read(lunit,'(a)') card
     close(lunit)

     card = adjustl(card)
     fits = (card(1:8) == 'SIMPLE  ' .OR. card(1:8) == 'XTENSION')
  else
     ! virtual file: assume it to be FITS
     fits = .true.
  endif

  return
end function my_isfits
!=======================================================================
subroutine get_pixwin_filename(pixwin_file_name_default, &
 &                             pixwin_file_name_input, &
 &                             pixwin_file_name, &
 &                             default,nside,healpix_data,stringlength)
!=======================================================================
! pixwin_file_name_default : default name for the pixel window file
! pixwin_file_name_input  : pixel window file if default is not taken
! default     (T/F) : whether or not default name is taken for pixel
!                     window file.
!=======================================================================
  use healpix_types
  use misc_utils, only: fatal_error
  implicit none

  character(len=*), intent(in)  :: pixwin_file_name_default
  character(len=*), intent(in)  :: pixwin_file_name_input
  character(len=*), intent(out) :: pixwin_file_name
  logical,          intent(in)  :: default
  integer(I4B),     intent(in)  :: nside
  character(len=*), intent(in)  :: healpix_data
  integer(I4B),     intent(in)  :: stringlength

  character(len=8) :: string='00000000'
  integer(I4B) :: leng
  logical :: exist

! Default name for input window file
  if (default) then
     call convert_to_ascii(string,nside)
     leng=len(trim(healpix_data)//trim(pixwin_file_name_default))
     if (leng > stringlength-9) then
        write(*,*) 'ERROR in get_pixwin_filename :'
        write(*,*) 'string pixwin_file_name_default is too long.'
        CALL FATAL_ERROR
     endif
     pixwin_file_name=trim(healpix_data)//trim(pixwin_file_name_default) &
 &                    //string(5:8)//'.fits'
! Customized file name
  else
     leng=len(trim(pixwin_file_name_input))
     if (leng > stringlength) then
        write(*,*) 'ERROR in get_pixwin_filename :'
        write(*,*) 'string pixwin_file_name_input is too long.'
        CALL FATAL_ERROR
     endif
     pixwin_file_name=trim(pixwin_file_name_input)
  endif

  inquire(file=pixwin_file_name,exist=exist)
  if (.not.exist) then
     write(*,*) 'ERROR in get_pixwin_filename'
     write(*,*) 'INPUT file : '//trim(pixwin_file_name)
     write(*,*) 'does not exist.'
     CALL FATAL_ERROR
  endif
  !print*,'get_pixwin_filename:'
  !print*,trim(pixwin_file_name_default),trim(pixwin_file_name_input),&
  !     & trim(pixwin_file_name),default,nside

end subroutine get_pixwin_filename

!===============================================================================
subroutine read_twod_fits_DP(file_name,npixtot,nx,ny,map)
!===============================================================================
  use healpix_types
  use misc_utils, only: fatal_error
  implicit none
  character*(*), intent(in) :: file_name
  integer(i4b),  intent(in) :: npixtot
  integer(i4b),  intent(in) :: nx,ny ! use-less
  real(DP),      intent(out) :: map(0:npixtot-1)

  integer(I4B) :: naxes(2)
  integer(I4B) :: status,unit,readwrite,blocksize,nfound,nbuffer
  integer(I4B) :: group,firstpix,naxis1,naxis2
  logical :: anynull
  real(DP) :: rnullval
  integer(I4B) :: i,j,inc
  character(len=30) :: errtext

  status=0
  call ftgiou(unit,status)
  readwrite=0

  call ftopen(unit,file_name,readwrite,blocksize,status)
  call ftgknj(unit,'NAXIS',1,2,naxes,nfound,status)
  
  naxis1 = naxes(1)
  naxis2 = naxes(2)

  nbuffer=naxis1*naxis2
  if (nbuffer > npixtot) then
     print*,'array too small in read_twod_fits'
     print*,nbuffer,npixtot
  endif
  group=1
  firstpix=1
  rnullval=max_sp

  call ftgpvd(unit,group,firstpix,nbuffer,rnullval,map,anynull,status)

  call ftgerr(status,errtext)
  call ftclos(unit, status)
  call ftfiou(unit, status)
  if (status > 0) then
     write(*,*) 'ERROR in read_twod_fits :',errtext
     CALL FATAL_ERROR
  endif

end subroutine read_twod_fits_DP


!===============================================================================
subroutine write_twod_fits_DP(file_name,npixtot,nx,ny,map)
!===============================================================================
  use healpix_types
  use misc_utils, only: fatal_error
  implicit none
  character*(*), intent(in) :: file_name
  integer(I4B),  intent(in) :: npixtot,nx,ny
  real(DP),      intent(in) :: map(0:npixtot-1)

  integer(I4B) :: status,unit,blocksize,bitpix,naxis,naxes(2)
  integer(I4B) :: i,j,group,fpixel,nelements
  logical :: simple,extend
  character(len=30) :: errtext

  status=0
  call ftgiou(unit,status)
  blocksize=1
  call ftinit(unit,file_name,blocksize,status)
  simple=.true.
  bitpix=-64
  naxis=2
  naxes(1)=nx
  naxes(2)=ny
  extend=.true.
  call ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)
  group=1
  fpixel=1
  nelements=npixtot
  call ftpprd(unit,group,fpixel,nelements,map,status)
  call ftclos(unit, status)
  call ftfiou(unit, status)
  if (status > 0) then
     write(*,*) 'ERROR in write_twod_fits :',errtext
     CALL FATAL_ERROR
  endif

end subroutine write_twod_fits_DP

!===============================================================================
subroutine write_nd_fits_DP(file_name, npixtot, naxes, naxis, map, header, nlheader)
!===============================================================================
  use healpix_types
  use misc_utils, only: fatal_error
  use fitstools, only: printerror, putrec
  implicit none
  character(len=*), intent(in) :: file_name
  integer(I8B),     intent(in) :: npixtot
  integer(I4B),     intent(in) :: naxes
  integer(I4B),     intent(in) :: naxis(1:naxes)
  real(DP),         intent(in) :: map(0:npixtot-1)
  integer(i4b),     intent(in) :: nlheader
  character(len=80),intent(in) :: header(1:nlheader)

  integer(I4B) :: status, unit, blocksize, bitpix
  integer(I4B) :: i, j, group, fpixel, nelements
!  logical      :: simple, extend
!  character*30 :: errtext

  ! open file
  status=0
  call ftgiou(unit,status)
  call printerror(status)
  blocksize=1
  call ftinit(unit, file_name, blocksize, status)
  call printerror(status)

  if (npixtot /= product(naxis(1:naxes)) ) then
     print*,'Error ',naxes, npixtot, naxis(1:naxes)
     call fatal_error
  endif
  ! write minimal header for image
!  simple=.true.
  bitpix=-64
!  extend=.true.
  call ftphps(unit,         bitpix, naxes, naxis,                        status)
!  call ftphpr(unit, simple, bitpix, naxes, naxis(1:naxes), 0, 1, extend, status)
  call printerror(status)

  ! write user provided header
  do i=1, nlheader
     if (trim(header(i)) /= '') call putrec(unit, header(i), status)
  enddo
  call printerror(status)

  ! write data for image
  group=1
  fpixel=1
  nelements=npixtot
  call ftpprd(unit, group, fpixel, nelements, map, status)
  call printerror(status)

  ! close and exit
  call ftclos(unit, status)
  call ftfiou(unit, status)
  call printerror(status)

end subroutine write_nd_fits_DP

!===============================================================================
subroutine check_twod_fits_header(file_name,npixtot,nx,ny)
!===============================================================================
  use healpix_types
  use misc_utils, only: fatal_error
  implicit none
  character*(*), intent(IN)  :: file_name
  integer(I4B),  intent(OUT) :: npixtot,nx,ny
  
  logical(LGT) :: ultraverbose=.false.
  integer(I4B) :: naxes(2)
  integer(I4B) :: status,unit,readwrite,blocksize,nfound
   
  status=0
  call ftgiou(unit,status)
  readwrite=0
  call ftopen(unit,file_name,readwrite,blocksize,status)
  call ftgknj(unit,'NAXIS',1,2,naxes,nfound,status)
  
  if (nfound .ne. 2) then
     write(*,*) 'ERROR in check_twod_fits_header :'
     write(*,*) 'Failure in reading the NAXIS keywords'
     CALL FATAL_ERROR
  endif

  nx=naxes(1)
  ny=naxes(2)
  npixtot=nx*ny

  if (ultraverbose) write(*,*) 'Dimensions of the image detected ',nx,ny
  
  call ftclos(unit, status)
  call ftfiou(unit, status)

end subroutine check_twod_fits_header

!===============================================================================
subroutine check_nd_fits_header(file_name, npixtot, naxes, naxis)
!===============================================================================
  ! determine image dimension, EH, 2009-01-27
  !
  use healpix_types
  use misc_utils, only: fatal_error
  implicit none
  character(len=*), intent(IN)  :: file_name
  integer(I4B),  intent(OUT) :: npixtot ! number of pixels
  integer(I4B),  intent(OUT) :: naxes   ! number of dimensions
  integer(I4B),  intent(OUT), dimension(1: ) :: naxis ! size in each dimension
  
  logical(LGT) :: ultraverbose=.false.
  integer(I4B) :: status, unit, readwrite, blocksize, nfound
  character(len=80) :: comment
   
  status=0
  call ftgiou(unit,status)
  readwrite=0
  call ftopen(unit, file_name, readwrite, blocksize, status)

  ! number of dimension
  call ftgkyj(unit,'NAXIS', naxes, comment, status)
  if (size(naxis) < naxes) then
     write(*,*) 'ERROR in check_nd_fits_header :'
     write(*,*) 'Found ',naxes,' dimension'
     write(*,*) 'Receiving array has size ',size(naxis)
     call fatal_error
  endif

  ! size in each dimension
  naxis(1:) = 1
  call ftgknj(unit,'NAXIS', 1, naxes, naxis, nfound, status)
  if (nfound .ne. naxes) then
     write(*,*) 'ERROR in check_nd_fits_header :'
     write(*,*) 'Failure in reading the NAXIS keywords'
     CALL FATAL_ERROR
  endif

  ! total number of elements
  npixtot = product(naxis(1:naxes))

  if (ultraverbose) write(*,*) 'Dimensions of the image detected ',naxis(1:naxes)
  
  call ftclos(unit, status)
  call ftfiou(unit, status)

end subroutine check_nd_fits_header

!===============================================================================
subroutine deletefile(filename,status)
!===============================================================================
! Extracted from cookbood.f of the cfitsio library
! A simple little routine to delete a FITS file
!===============================================================================
  use healpix_types
  integer(I4B) :: status,unit,blocksize
  character*(*) :: filename

! Simply return if status is greater than zero
  if (status .gt. 0)return

! Get an unused Logical Unit Number to use to open the FITS file
  call ftgiou(unit,status)

! Try to open the file, to see if it exists
  call ftopen(unit,filename,1,blocksize,status)

  if (status .eq. 0)then
!    file was opened;  so now delete it 
     call ftdelt(unit,status)
  else if (status .eq. 103)then
!    file doesn't exist, so just reset status to zero and clear errors
     status=0
     call ftcmsg
  else
!    there was some other error opening the file; delete the file anyway
     status=0
     call ftcmsg
     call ftdelt(unit,status)
  end if

!  Free the unit number for later reuse
  call ftfiou(unit, status)
end subroutine deletefile

! !=======================================================================
! function my_getnumext_fits(filename)
!   !=======================================================================
!   ! patch to bugged Healpix 2.01 getnumext_fits
!   !
!   !  result = my_getnumext_fits(filename)
!   !    returns the number of extensions present in FITS file 'filename'
!   !
!   ! EH, Nov 2004
!   ! April 2007: close file on exit
!     !=======================================================================
!   use healpix_types
!   implicit none
!   character(LEN=*), intent(IN)             :: filename
!   integer(i4b)                             :: my_getnumext_fits
!   !
!   integer(i4b) :: status, unit, readwrite, blocksize, nhdu
!   !-----------------------------------------------------------------------
!   status         = 0
!   unit           = 149
!   my_getnumext_fits = 0
!   readwrite      = 0 ! Read only
!   call ftopen(unit, filename, readwrite, blocksize, status)
!   if (status > 0) then
!      !       call printerror(status)
!      print*,'FITS Error: ',status
!      call ftclos(unit, status)
!      return
!   endif
    
!   call ftthdu(unit, nhdu, status)
!   my_getnumext_fits = nhdu - 1

!   call ftclos(unit, status)
!   return
! end function my_getnumext_fits

!====================================================================
subroutine write_spice_header(headfits, nlh, nlmax, ncor, nsides, aboutmap)
  !====================================================================
  use healpix_types
#ifndef PIO
  use spice_parameters, only : decouple, version, thetamax, &
       & apodize, apodizesigma, apodizetype, map_present, map2_present, &
       masks_present, masks2_present, weights_present, weights2_present, &
       weightpower, weightpower2, &
       weightpowerp, weightpowerp2, &
       fwhm, fwhm2, correct_beam, correct_beam2, beam_present, beam2_present,&
       beam_file, beam_file2, &
       correct_transfer_function, subtract_average, &
       subtract_dipole, &
       pairsthresholding, npairsthreshold, tolerance
  use spice_common, only : fits_out, &
       mapfile, maskfile, weightfile, maskfilep, weightfilep, &
       mapfile2, maskfile2, weightfile2, maskfilep2, weightfilep2, &
       tf_file, listmapfile, listmapw8
  use head_fits, only: add_card
  implicit none
  integer(I4B),     intent(in)                  :: nlh, nlmax, ncor
  character(len=80), intent(inout), dimension(1:nlh) :: headfits
  integer(I4B),     intent(in), dimension(1:2)  :: nsides
  logical(LGT),     intent(in)                  :: aboutmap

  call add_card(headfits)
  call add_card(headfits,'COMMENT','---------------------------------')
  call add_card(headfits)
  call add_card(headfits,"CREATOR","Spice", "Software creating the FITS file")
  call add_card(headfits,"VERSION",version, "Version of the simulation software")
  call add_card(headfits)
  call add_card(headfits,"POLAR",  (ncor>1),"Polarisation included (True/False)")
  call add_card(headfits,"BCROSS", (ncor>4),"Magnetic cross terms included (True/False)")
  call add_card(headfits,"ASYMCL", (ncor>6),"Asymmetric pol cross terms (XY vs YX) included")
  call add_card(headfits,'APODIZE',apodize, "Apodization of Xi  (True/False)")
  call add_card(headfits,'DECOUPLE',decouple,"Decouple E and B polarization (True/False)")
  call add_card(headfits,'SUBAV',subtract_average .or. subtract_dipole &
       & ,"Remove average from T map (True/False)")
  call add_card(headfits,'SUBDIPOL',subtract_dipole,"Remove dipole from T map (True/False)")
  call add_card(headfits,'TF_CORRE',correct_transfer_function,'Tranfer funct. correction (True/False)')
  call add_card(headfits)
  call add_card(headfits,'MAX-LPOL',nlmax,  "Maximum L multipole order")
  call add_card(headfits,'NCOR',   ncor,    "Number of fields (1 or 4)")
  if (nsides(1) == nsides(2)) then
     call add_card(headfits,'NSIDE',  nsides(1),   "Resolution Parameter of Input Map")
  else
     call add_card(headfits,'NSIDE1',  nsides(1),   "Resolution Parameter of Input Map #1")
     call add_card(headfits,'NSIDE2',  nsides(2),   "Resolution Parameter of Input Map #2")
  endif
  call add_card(headfits,'THETAMAX',thetamax, "[Deg] Largest lag in angul. correl. fct Xi")
  call add_card(headfits,'TOLERANC',tolerance, "Relative tolerance on E/B decoupling integrals")
  if (apodize) then
     call add_card(headfits,'APOTYPE',apodizetype, "Apodization type: 0=Gaussian, 1=Cosine")
     call add_card(headfits,'APOSIGMA',apodizesigma, "[Deg] Apodization parameter")
  endif
  if (pairsthresholding) then
     call add_card(headfits,'NPTHRESH',npairsthreshold,"Minimal # of pixel pairs used in correlation")
  endif
  if (aboutmap) then
     if (map_present)  call add_card(headfits,'MAPFILE1',trim(listmapfile(1,1,1)),'Input map file 1')
     if (map2_present) call add_card(headfits,'MAPFILE2',trim(listmapfile(1,1,2)),'Input map file 2')
  endif
  if (masks_present) then
     call add_card(headfits,'MASKFIL1',trim(maskfile),'Input mask file 1')
     if (maskfilep /= maskfile) &
          & call add_card(headfits,'MASKF_P1',trim(maskfilep),'Input pol mask file 1')
  endif
  if (masks2_present) then
     call add_card(headfits,'MASKFIL2',trim(maskfile2),'Input mask file 2')
     if (maskfilep2 /= maskfile2) &
          & call add_card(headfits,'MASKF_P2',trim(maskfilep2),'Input pol mask file 2')
  endif
  if (weights_present) then
     call add_card(headfits,'W8FILE1',trim(weightfile),'Input weight file 1')
     call add_card(headfits,'W8POWER1',weightpower,'power index applied to weight file data 1')
     if (weightpower /= weightpowerp .or. weightfilep /= weightfile) then
        call add_card(headfits,'W8FIL_P1',trim(weightfilep),'Input pol weight file 1')
        call add_card(headfits,'W8POW_P1',weightpowerp,'power index applied to weight pol file data 1')
     endif
  endif
  if (weights2_present) then
     call add_card(headfits,'W8FILE2',trim(weightfile2),'Input weight file 2')
     call add_card(headfits,'W8POWER2',weightpower2,'power index applied to weight file data 2')
     if (weightpower2 /= weightpowerp2 .or. weightfilep2 /= weightfile2) then
        call add_card(headfits,'W8FIL_P2',trim(weightfilep),'Input pol weight file 2')
        call add_card(headfits,'W8POW_P2',weightpowerp2,'power index applied to weight pol file data 2')
     endif
  endif
  if (correct_beam) then
     if (beam_present) then 
        call add_card(headfits,'BEAMFIL1',trim(beam_file),'Beam file 1')
     else
        call add_card(headfits,'FWHM1',fwhm/60.0_dp,'[Deg] beam FWHM 1')
     endif
  endif
  if (correct_beam2) then
     if (beam2_present) then 
        call add_card(headfits,'BEAMFIL2',trim(beam_file2),'Beam file 2')
     else
        call add_card(headfits,'FWHM2',fwhm2/60.0_dp,'[Deg] beam FWHM 2')
     endif
  endif
  if (correct_transfer_function) then
     call add_card(headfits,'TF_FILE',trim(tf_file),'Transfer function file')
  endif
  ! missing:  maskfilep, weightfilep, ...
  call add_card(headfits)
  call add_card(headfits,'COMMENT','---------------------------------')
  call add_card(headfits)
#endif

  return
end subroutine write_spice_header


!=======================================================
subroutine get_fits_header(filename, header, hdu)
!=======================================================
  ! routine to read fits header only
  ! adapted from listhead.c
  use healpix_types
  use healpix_modules, only: fatal_error, printerror
  implicit none
  character(len=filenamelen),         intent(in) :: filename
  character(len=80), dimension(1:),   intent(out):: header
  integer(I4B),    optional,          intent(in) :: hdu


  integer(i4b) :: status, unit, readwrite, blocksize, nhdu
  integer(i4b) :: nkeys, keysadd, ik, myhdu, hdutype
    !-----------------------------------------------------------------------
  status         = 0
  unit           = 149
  readwrite      = 0 ! Read only
  call ftopen(unit, filename, readwrite, blocksize, status)
  if (status > 0) then
     call printerror(status)
     call ftclos(unit, status)
     return
  endif


  myhdu = 1
  if (present(hdu)) then
     myhdu = hdu
  endif

  ! total number of HDU
  call ftthdu(unit, nhdu, status)

  if (myhdu+1 > nhdu) then 
     print*,trim(filename)//' contains ',nhdu,' HDUs'
     print*,' = 1 primary + ',nhdu-1,' extensions'
     print*,'Asking for extension #',myhdu
     call fatal_error
  endif

  call ftmahd(unit, myhdu+1, hdutype, status) ! move to HDU (1=primary)

  ! number of keyword (expect 'END')
  call ftghsp(unit, nkeys, keysadd, status)

  if (size(header) < nkeys) then
     print*,'Header array is too short to read header of '//trim(filename)
     print*,'Size: ',size(header)
     print*,'Need: ',nkeys
     call fatal_error
  endif

  ! read header
  if (status > 0) call printerror(status)
  do ik=1,nkeys
     call ftgrec(unit, ik, header(ik), status)
  enddo
  if (status > 0) call printerror(status)

  ! close and leave
  call ftclos(unit, status)
  return

end subroutine get_fits_header

end module deal_with_files
