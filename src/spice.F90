! edited to run with Healpix 1.2
!=======================================================================
#ifdef LOCALCODE
function spice_lib(param_filename)
#else
program spice
#endif
!=======================================================================
  use spice_common
  use spice_parameters, only: version
  use spice_subs, only: compute
  use windows
  use misc_utils, only: brag_openmp, wall_clock_time
  use deal_with_files, only: check_headers, input_files, output_results, &
        & check_open_output_files
#ifdef LOCALCODE
  use piolibkind
#endif

  implicit none
#ifdef LOCALCODE
    character(len=DMCPIOSTRINGMAXLEN)              :: param_filename
    INTEGER(KIND=PIOINT)                           :: spice_lib
#endif

#ifdef MPI
  ! MPI is *NOT* used in the code.
  ! It is only invoked here in order to prompt Nersc's IPM monitoring
  include 'mpif.h'
  integer :: ierr, nproc, myid
#endif
  real (SP) :: ctime0, ctimef, ctime1
  real (SP) :: wtime0, wtimef, wtime1
 
#ifdef MPI
  call Mpi_Init(ierr)
  call MPI_Comm_size(MPI_COMM_WORLD,nproc,ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD,myid, ierr)
#endif
  call cpu_time(ctime0)
  call wall_clock_time(wtime0)

#ifdef PIO
!  print *,"init logger "
  ! create logger instance
  pio_log=-1
  pio_error = PIOInitLogger(pio_log)
  print*,pio_error, pio_log
  if (pio_error < 0) then
     write(0,*) "Spice: Cannot initialize logger, Aborting"
     call exit(pio_error)
  endif
  write(pio_logStr,*) '--------------------------'
  pio_error = PIOInfoLogger(pio_log,pio_logStr)
  write(pio_logStr,*) 'Starting Spice '//version
  pio_error = PIOInfoLogger(pio_log,pio_logStr)
  write(pio_logStr,*) '--------------------------'
  pio_error = PIOInfoLogger(pio_log,pio_logStr)
#endif

  listmapfile = ' '
  nlmaps      = 0
  listmapw8   = 1.0_dp
  extramapTQUfile_input = ' '

  ! read run parameters
  ! make sure that all internal flags are consistent 
  ! with each other and with run parameters
!  print *,"Deal With File"
#ifdef LOCALCODE
  call deal_with_options(param_filename)
#else
  call deal_with_options
#endif
  ! get relevant environment variables
  call init_system

  ! provide final filenames, check that files can be read or written
  call deal_with_filenames

  ! check header of existing files
  call check_headers

  ! reasses options that depend on actual filenames
  call reassess_options

#ifndef PIO
  call check_open_output_files
#endif
  ! summarize parameters and files before starting
  if (verbose) call summarize_all

  if (verbose) call brag_openmp
  ! compute correction due to apodization (if not reading precomputed one)
  !!!if (.not.windowinput .and. polarization) then
  !!!   fullkernels = kernelsoutput ! compute either all kernels (if required for output) or only TE one
  !!! if (.not.windowinput  .and. (kernelsoutput .or. windowoutput)) then
  if (.not.    (windowinput .or.tenorminput)  &
       & .and. (kernelsoutput .or. windowoutput .or. tenormoutput)) then
     fullkernels = kernelsoutput
     call compute_kernels
     if (.not. windowoutput) then
        if (allocated(kcross))   deallocate(kcross)
        if (.not. kernelsoutput .and. allocated(kernels)) deallocate(kernels)
     endif
  endif

  ! input data (heavy on memory, delay as much as possible)
  call input_files

  ! do actual calculation
  call compute

  ! output results
  call output_results

  ! deallocate global variables
  if (allocated(xi_final)) deallocate(xi_final)
  if (allocated(cl))       deallocate(cl)
  if (allocated(cl_mask))  deallocate(cl_mask)
  if (allocated(cl_final)) deallocate(cl_final)
  if (allocated(kcross))   deallocate(kcross)
  if (allocated(kernels))  deallocate(kernels)
  if (allocated(TEnorm))   deallocate(TEnorm)
  if (allocated(mu))       deallocate(mu)
  if (allocated(w))        deallocate(w)
  if (allocated(Pl))       deallocate(Pl)
  if (allocated(cov_matrix))  deallocate(cov_matrix)

  call cpu_time(ctimef)
  call wall_clock_time(wtimef)
  ctimef = ctimef - ctime0
  wtimef = wtimef - wtime0

#ifdef PIO
  write(pio_logStr,*) '--------------------------'
  pio_error = PIOInfoLogger(pio_log,pio_logStr)
  write(pio_logStr,*) 'Spice Finished Successfully'
  pio_error = PIOInfoLogger(pio_log,pio_logStr)
  write(pio_logStr,*) '--------------------------'
  pio_error = PIOInfoLogger(pio_log,pio_logStr)
  pio_error = PIOFreeLogger(pio_log)
#endif
  print*,"CPU  Time [s]: ",ctimef
  print*,"Wall Time [s]: ",wtimef

#ifdef MPI
  call Mpi_Finalize(ierr)
#endif

#ifdef LOCALCODE
  spice_lib=0
  end function spice_lib
#else
  stop
  end program spice
#endif
!=======================================================================
subroutine show_used_maps(jmap)
  !=======================================================================
  use spice_common, only: i4b, listmapfile, listmapw8, nlmaps
  use healpix_modules, only: string
  implicit none
  integer(i4b), intent(IN) :: jmap
  integer(i4b) :: kstokes, i
  character(len=1) :: sjmap
  character(len=1), dimension(1:3) :: sstokes
  character(len=12) :: form

  sjmap = string(jmap, format='(i1)')
  sstokes(1:3) = (/ 'T', 'Q', 'U' /)

  kstokes = 1
  form='(f7.3)'
  write(*,*)    'map file ('//sjmap//')......'
  do i=1,nlmaps(kstokes,jmap)
     write(*,*) '  '//trim(string(listmapw8(i,kstokes,jmap), format=form)) &
          &  //' *   '//trim(listmapfile(i,kstokes,jmap))
  enddo

#ifdef PIO
  do kstokes=2,3
     write(*,*)    'map '//sstokes(kstokes)//' file ('//sjmap//')....'
     do i=1,nlmaps(kstokes,jmap)
        write(*,*) '  '//trim(string(listmapw8(i,kstokes,jmap), format=form)) &
             &  //' *   '//trim(listmapfile(i,kstokes,jmap))
     enddo
  enddo
#endif
  return
end subroutine show_used_maps

!=======================================================================
subroutine summarize_all ! summary
!=======================================================================
  use spice_common
  implicit none
  character(len=FILENAMELEN) :: file_n,file_n2,file_n3,tmp
  integer(i4b):: k
  logical(lgt) :: flag

  write(*,*)
  write(*,*) '==========================================================='
  write(*,*) 'SpICE '//trim(version)//' has been called with the following attributes : '
  write(*,*) 
  if (nsides(1)==nsides(2)) then
     write(*,*) 'nside detected... ',nsides(1)
  else
     write(*,*) 'nsides detected... ',nsides(1),nsides(2)
  endif
  write(*,*) 'max multipole used',nlmax
  write(*,*) 'polarization..... ',trim(file_n3(polarization,default,none))
  !---------------------------
  call show_used_maps(1)

  flag = (masks_present .and. maskfilep_toread)
  if (polarization) tmp = 'mask file....(TP). '
  if (flag .or. .not. polarization) tmp = 'mask file....(T).. '
  write(*,*) trim(tmp),trim(file_n(masks_present,maskfile,none))
  if (flag) &
       & write(*,*) 'mask file (Pol).. ',trim(file_n(maskfilep_toread,maskfilep,none))

  flag = (weights_present .and. weightfilep_toread)
  if (polarization) tmp = 'weight file..(TP). '
  if (flag .or. .not. polarization) tmp = 'weight file..(T).. '
  write(*,*) trim(tmp),trim(file_n(weights_present,weightfile,none))
  if (flag) &
       & write(*,*) 'weight file.(Pol). ',trim(file_n(weightfilep_toread,weightfilep,none))

  write(*,*) 'weight power..... ',weightpower
  if (polarization) write(*,*) 'weight power (Pol)',weightpowerp
  write(*,*) 'dump alm(1)...... ',trim(file_n(alm1_dump, alm1_out_file, none))
  tmp=file_n2(correct_beam,fwhm,none) ! XLF does not accept a write in a write
  write(*,*) 'gauss beam arcmin ',trim(tmp)
  write(*,*) 'beam file........ ',trim(file_n(beam_present,beam_file,none))
  write(*,*) 'pix. window file. ',trim(file_n(correct_pix,pixwin_file_name,none))
  !---------------------------
  if (map2_present .or. masks2_present .or. weights2_present &
       & .or. correct_beam2 .or. weightpower/=weightpower2) then

     call show_used_maps(2)


     flag = (masks2_present .and. maskfilep2_toread)
     if (polarization) tmp = 'mask file...(TP2). '
     if (flag .or. .not. polarization) tmp = 'mask file...(T2).. '
     write(*,*) trim(tmp),trim(file_n(masks2_present,maskfile2,none))
     if (flag) &
       & write(*,*) 'mask file (Pol2). ',trim(file_n(maskfilep2_toread,maskfilep2,none))

     flag = (weights2_present .and. weightfilep2_toread)
     if (polarization) tmp = 'weight file.(TP2). '
     if (flag .or. .not. polarization) tmp = 'weight file..(T2). '
     write(*,*) trim(tmp),trim(file_n(weights2_present,weightfile2,none))
     if (flag) &
       & write(*,*) 'weight file (Pol2). ',trim(file_n(weightfilep2_toread,weightfilep2,none))

     write(*,*) 'weight power (2). ',weightpower2
     if (polarization) write(*,*) 'weight pow  (Pol2)',weightpowerp2
     write(*,*) 'dump alm(2)...... ',trim(file_n(alm2_dump, alm2_out_file, none))
     tmp=file_n2(correct_beam2,fwhm2,none) ! XLF does not accept a write in a write
     write(*,*) 'g.beam (2) arcmin ',trim(tmp)
     write(*,*) 'beam file (2).... ',trim(file_n(beam2_present,beam_file2,none))
     write(*,*) 'pix. wind. file(2)',trim(file_n(correct_pix2,pixwin_file_name2,none))
     write(*,*) 'symmetrized C(l)  ',trim(file_n3(symmetric_cl,default,none))
  endif
  !---------------------------
  write(*,*) 'subtract dipole.. ',trim(file_n(subtract_dipole,default,none))
  write(*,*) 'subtract average. ',trim(file_n(subtract_average,default,none))
  write(*,*) 'cor. file........ ',trim(file_n(coroutput,   cor_file_name,none))
  write(*,*) 'Cl file.......... ',trim(file_n(cloutput,    cl_file_name,none))
  write(*,*) 'Covariance file.. ',trim(file_n(do_cov,   cov_out_file,none))
  write(*,*) 'raw Cl map output ',trim(file_n(clmapoutput, cl_outmap_file,none))
  write(*,*) 'raw Cl map input  ',trim(file_n(clmapinput,  cl_inmap_file,none))
  write(*,*) 'raw Cl mask outp. ',trim(file_n(clmaskoutput,cl_outmask_file,none))
  write(*,*) 'raw Cl mask inpu. ',trim(file_n(clmaskinput, cl_inmask_file,none))
!  write(*,*) 'raw Cl file...... ',trim(file_n(cloutput,    cl_file_name,none))
  write(*,*) 'FITS output...... ',trim(file_n3(fits_out, default, none))
  tmp=file_n2(normalize,quad_uK,none)
  write(*,*) 'norm. factor..... ',trim(tmp)
!   tmp=file_n2(pairsthresholding,npairsthreshold,none)
!   write(*,*) 'pairs thresh val. ',trim(tmp)
  tmp=file_n2(apodize,apodizesigma,none)
  write(*,*) 'apodizing width.. ',trim(tmp)
  write(*,*) 'apodizing type... ',apodizetype
  write(*,*) 'decouple......... ',trim(file_n3(decouple,default,none))
  if (decouple) then
     write(*,*) 'tolerance........ ',tolerance
  endif
  write(*,*) 'thetamax......... ',thetamax
  write(*,*) 'kernels output... ',trim(file_n(kernelsoutput,kernels_out_file,none))
  write(*,*) 'window output.... ',trim(file_n(windowoutput,window_out_file,none))
  write(*,*) 'window input..... ',trim(file_n(windowinput, window_in_file, none))
  write(*,*) 'TE norm output... ',trim(file_n(tenormoutput,tenorm_out_file,none))
  write(*,*) 'TE norm input.... ',trim(file_n(tenorminput, tenorm_in_file, none))
  write(*,*) 'noise cor file... ',trim(file_n(subtract_noise_cor,noisecor_file_name,none))
  write(*,*) 'noise Cl file.... ',trim(file_n(subtract_noise_cl, noisecl_file_name, none))
  write(*,*) 'overwrite mode... ',trim(file_n(overwrite,default,none))
#ifndef PIO
  write(*,*) 'dry format....... ',trim(file_n3(dry,default,none))
#endif
#ifndef PIO
  write(*,*) 'output opt. file. ',trim(file_n(output_spice_rc,spice_rc_outfile,none))
  write(*,*) 'input opt. file.. ',trim(file_n(input_spice_rc, spice_rc_infile, none))
#endif
  write(*,*) '==========================================================='
  write(*,*)
end subroutine summarize_all

!=======================================================================
function file_n(present, file_name, none)
!=======================================================================
  use healpix_types
  implicit none
  logical,          intent(in) :: present
  character(len=*), intent(in) :: file_name, none
  character(len=FILENAMELEN) :: file_n

  if (present) then
     file_n=file_name
  else
     file_n=none
  endif
end function file_n

!=======================================================================
function file_n2(present, xnum, none)
!=======================================================================
  use healpix_types
  implicit none
  logical,          intent(in) :: present
  real(DP),         intent(in) :: xnum
  character(len=*), intent(in) :: none
  character(len=FILENAMELEN) :: file_n2
  if (present) then 
      write(file_n2,'(ES15.5)') xnum
  else 
     file_n2=none
  endif 
end function file_n2

!=======================================================================
function file_n3(present, default, none)
!=======================================================================
  use healpix_types
  implicit none
  logical,          intent(in) :: present
  character(len=*), intent(in) :: none,default
  character(len=FILENAMELEN)   :: file_n3
  if (present) then 
     file_n3=default
  else 
     file_n3=none
  endif 
end function file_n3
