! subroutine deal_with_options
! subroutine check_argument
! subroutine error_usage
! subroutine usage_reminder
! subroutine man
! subroutine about
! subroutine history
! subroutine make_option_active
! subroutine create_option
! subroutine check_options_choice
! subroutine check_compute
! subroutine create_spicerc
! subroutine read_spicerc
! subroutine reassess_options

!=======================================================================
#ifdef LOCALCODE
subroutine deal_with_options(param_file)
#else
subroutine deal_with_options
#endif
!=======================================================================
  use misc_utils, only: fatal_error
  use spice_common
#ifdef PIO
!   use pie_io, only: parContent, init_parContent
  use pie_io, only: init_parContent
  use piolib
#endif
  implicit none

#ifdef LOCALCODE
  character(len=DMCPIOSTRINGMAXLEN) :: param_file
#endif
  character(len=FILENAMELEN), allocatable, dimension(:) :: argument
  integer(I4B) :: i,narguments,nargu
  integer(I4B) :: ijunk=-1,jjunk=-1,iskipjunk=-1
!#ifndef LINUX
  integer(I4B) :: iargc
!#endif
  
  ! count the number of options available in code and put it in common noption
  allocate(argument(1))
  narguments=1
  call make_option_active(ijunk,jjunk,iskipjunk,argument,narguments,.true.,.true.)
  deallocate(argument)

#ifdef PIO
#ifdef LOCALCODE
  print *,'param_file ',param_file
  spice_rc_infile= param_file
#else
  nargu=iargc()
  if (nargu == 0) then
     spice_rc_infile=spice_rc_file_default
  elseif (nargu == 1) then
     call getarg(1,spice_rc_infile)
  elseif (nargu > 1) then
     write(*,*) 'USAGE : spice config_file.txt'
     CALL FATAL_ERROR
  endif
#endif
  narguments=noption-jstartnormal+1
  allocate(argument(narguments))
  call init_parContent(pie_param , spice_rc_infile, pio_log)
#else 
  narguments=iargc()
  allocate(argument(narguments))
  do i=1,narguments
     call getarg(i,argument(i))
  enddo
#endif

  call check_argument(argument,narguments)
  deallocate(argument)

  if (input_spice_rc) call read_spicerc
  call check_options_choice
  if (output_spice_rc) then
     call create_spicerc
     if (verbose) write(*,*) "Exiting."
     CALL FATAL_ERROR
  endif
  call check_compute

  deallocate(option_done)
  deallocate(option)
  deallocate(option_param)
  deallocate(iskip_option)

! #ifdef PIO
!   pio_error=PIOcloseparfile(mypar)
! #endif
end subroutine deal_with_options


!=======================================================================
subroutine check_argument(argument,narguments)
!=======================================================================
  use spice_common
#ifdef PIO
  use piolib
#endif
  implicit none

  integer(I4B) :: narguments
  character(len=*) :: argument(narguments)
  integer(I4B) :: j,i,iskip,ijunk=-1,iskipjunk=-1
  logical :: option_ok

  allocate(option(noption))
  allocate(option_param(noption))
  allocate(iskip_option(noption))  
  option = ''
  option_param = ''
  iskip_option = 0

  ! create lists containing options and their default value
  do j=1,noption
     call make_option_active(ijunk,j,iskipjunk,argument,narguments,.true.,.false.)
  enddo

#ifdef PIO
  argument(1:narguments)=option(jstartnormal:noption)
#endif

  allocate(option_done(noption))
  option_done=.false.

  ! parse argument list, to find actual value of options
  i=0
  iskip=1
  ! check each argument
  do while (i+iskip.le.narguments)
     i=i+iskip
     option_ok=.false.
     ! find matching code option
     do j=1,noption
        if (trim(argument(i)) == trim(option(j))) then
           if (.not.option_done(j)) then
              call make_option_active(i,j,iskip,argument,narguments,.false.,.false.)
              option_done(j)=.true.
              option_ok=.true.
           endif
        endif
     enddo
     if (.not.option_ok) then
        print*,'Unrecognised option: '//trim(argument(i))
        call error_usage('Mismatched option count')
     endif
#ifdef PIO
     if (i == narguments-1) iskip=2
#endif
  enddo
end subroutine check_argument

!=======================================================================
subroutine error_usage(message)
!=======================================================================
  ! issues error message, together with short usage reminder
  !---------------------------------------------------------------------
  use misc_utils, only: fatal_error
  use spice_common
#ifdef PIO
  use piolib
#endif
  implicit none
  character*(*) :: message

#ifdef PIO
  write(*,*) 'Something is wrong with the config file.'
  write(*,*) 'ERROR : '//trim(message)
  CALL FATAL_ERROR
#else
  write(*,*) 'ERROR : '//trim(message)
  write(*,*)
  call usage_reminder
  CALL FATAL_ERROR
#endif

end subroutine error_usage

!=======================================================================
subroutine usage_reminder
!=======================================================================
  ! issues short usage reminder
  !---------------------------------------------------------------------
  use spice_common
  implicit none

  write(*,9)
  write(*,9) 'USAGE : spice'
  write(*,9) '              -about'
  write(*,9) '              -alm1fileout [name] -alm2fileout [name] '
  write(*,9) '              -apodizesigma [dfloat] -apodizetype [0|1] '
  write(*,9) '              -beam [dfloat] -beam_file [name]'
  write(*,9) '              -beam2 [dfloat] -beam_file2 [name]'
  write(*,9) '              -clfile [name]'
  write(*,9) '              -cl_inmap_file [name]'
  write(*,9) '              -cl_inmask_file [name]'
  write(*,9) '              -cl_outmap_file [name]'
  write(*,9) '              -cl_outmask_file [name]'
  write(*,9) '              -corfile [name]'
  write(*,9) '              -covfileout [name]'
  write(*,9) '              -decouple ['//trim(default)//'|'//trim(none)//']'
  write(*,9) '              -dry ['//trim(default)//'|'//trim(none)//']'
  write(*,9) '              -extramapfile [name] -extramapfile2 [name]'
  write(*,9) '              -fits_out ['//trim(default)//'|'//trim(none)//']'
  write(*,9) '              -help'
  write(*,9) '              -history'
  write(*,9) '              -kernelsfileout [name]'
  write(*,9) '              -listmapfiles1_*    [name]  -listmapfiles2_*   [name]'
  write(*,9) '              -listmapweights1_* [dfloat] -listmapweights2_* [dfloat]'
  write(*,9) '              -mapfile [name] -mapfile2 [name]'
  write(*,9) '              -maskfile [name] -maskfile2 [name]'
  write(*,9) '              -maskfilep [name] -maskfilep2 [name]'
  write(*,9) '              -nlmax [integer]'
  write(*,9) '              -normfac [dfloat]'
!  write(*,9) '              -npairsthreshold [dfloat]'
  write(*,9) '              -noiseclfile [name]'
  write(*,9) '              -noisecorfile [name]'
  write(*,9) '              -optinfile [name] -optoutfile [name]'
  write(*,9) '              -overwrite ['//trim(default)//'|'//trim(none)//']'
  write(*,9) '              -pixelfile [name]'
  write(*,9) '              -pixelfile2 [name]'
  write(*,9) '              -polarization['//trim(default)//'|'//trim(none)//']' 
  write(*,9) '              -subav ['//trim(default)//'|'//trim(none)//']'
  write(*,9) '              -subdipole ['//trim(default)//'|'//trim(none)//']'
  write(*,9) '              -symmetric_cl ['//trim(default)//'|'//trim(none)//']'
  write(*,9) '              -tenormfilein [name] -tenormfileout [name]'
  write(*,9) '              -tf_file [name|'//trim(none)//']'
  write(*,9) '              -thetamax [dfloat|'//trim(none)//']'
  write(*,9) '              -tolerance [dfloat|'//trim(none)//']'
  write(*,9) '              -usage'
  write(*,9) '              -verbosity [0|1|2]'
  write(*,9) '              -version'
  write(*,9) '              -weightfile [name] -weightfile2 [name]' 
  write(*,9) '              -weightfilep [name] -weightfilep2 [name]' 
  write(*,9) '              -weightpower [dfloat] -weightpower2 [dfloat] '
  write(*,9) '              -weightpowerp [dfloat] -weightpowerp2 [dfloat] '
  write(*,9) '              -windowfilein [name] -windowfileout [name]'
9 format(a)

end subroutine usage_reminder

!=======================================================================
subroutine man
!=======================================================================
  use spice_common
  use misc_utils, only: string
  implicit none
  integer(I4B) :: j
  character(len=*), parameter :: hpx_fits_file = '(healpix SP fits file)'

  write(*,9)
  write(*,9) 'INSTALLATION : '
  write(*,9) '   See the file named INSTALL in the top directory '
  write(*,9) '  '
  write(*,9) 'STANDARD USAGE : '
  write(*,9) '   spice [-keyword(1) parameter(1) ... -keyword(N) parameter(N)] '
  write(*,9) '         [-optinfile parameter_file]'
  write(*,9) '  '
  write(*,9) '   All keywords are optional'
  write(*,9) '  '
#ifdef PIO
  write(*,9) 'PLANCK PIPELINE USAGE (when compiled with PIOLIB): '
  write(*,9) '   spice parameter_file'
  write(*,9) 
  write(*,9) 
#endif
  write(*,9) ' Format of parameter_file :'
  write(*,9) '              keyword = parameter'
  write(*,9) '       see example below'
  write(*,9)
  write(*,9) 'LIST OF KEYWORDS :     keyword [possible parameters](default)'
  write(*,9) 
  write(*,9) '-about'
  write(*,9) '   print copyright info and stop'
  write(*,9) 
  write(*,9) '-alm1fileout [name|'//trim(default)//'|'//trim(none)//']('//trim(none)//')'
  write(*,9) '   FITS file containing a_lm coefficients of (masked and weighted) map1'
  write(*,9) '   '//trim(none)//' desactivates output of the alm'
  write(*,9) '   '//trim(default)//' activates output of the alm into '//trim(alm1file_default)
  write(*,9) 
  write(*,9) '-alm2fileout [name|'//trim(default)//'|'//trim(none)//']('//trim(none)//')'
  write(*,9) '   FITS file containing a_lm coefficients of (masked and weighted) map2'
  write(*,9) '   '//trim(none)//' desactivates output of the alm'
  write(*,9) '   '//trim(default)//' activates output of the alm into '//trim(alm1file_default)
  write(*,9) 
  write(*,9) '-apodizesigma [dfloat|'//trim(none)//']('//trim(none)//')'
  write(*,9) '   scale factor in DEGREES of the correlation function tappering,'
  write(*,9) '   see option -apodizetype.'
  write(*,9) '   For better results, -apodizesigma should be close to -thetamax.'
  write(*,9) '   '//trim(none)//' desactivates apodization.'
  write(*,9) '   NOTE : if noise subtraction is active (options -noisecorfile or'
  write(*,9) '   -noiseclfile), the input noise data are assumed to be already apodized.'
  write(*,9)
  write(*,9) '-apodizetype [int](0)'
  write(*,9) '   type of apodization'
  write(*,9) '   0: the correlation function is multiplied by a gaussian window'
  write(*,9) '      (equal to 1    at theta=0,'
  write(*,9) '       equal to 0.5  at theta= -apodizesigma/2,'
  write(*,9) '       equal to 1/16 at theta= -apodizesigma).'
  write(*,9) '   1: the correlation function is multiplied by a cosine window '
  write(*,9) '      (equal to 1   at theta=0,'
  write(*,9) '       equal to 0.5 at theta= -apodizesigma/2,'
  write(*,9) '       equal to 0   at theta= -apodizesigma).'
  write(*,9)
  write(*,9) '-beam [dfloat|'//trim(none)//']('//trim(none)//')'
  write(*,9) '   gaussian beam FWHM (arcmin) for map 1.'
  write(*,9) '   If specified, output c_l and c_theta are corrected for it'
  write(*,9) '   NOTE : if noise subtraction is active (options -noisecorfile or'
  write(*,9) '   -noiseclfile), the input noise data are assumed to be already corrected'
  write(*,9) '   for the beam'
  write(*,9)
  write(*,9) '-beam_file [name|'//trim(none)//']('//trim(none)//')'
  write(*,9) '   beam window function for map 1.'
  write(*,9) '   If specified, output c_l and c_theta are corrected for it'
  write(*,9) '   The file should be plain ASCII with 2 columns : l and B(l) '
  write(*,9) '   or FITS (ASCII or BINARY) table containing B(l) starting at l=0'
  write(*,9) '   (B(l) is sqrt of beam power spectrum).'
  write(*,9) '   This over-rules the gaussian correction described by beam.'
  write(*,9) '   '//trim(none)//' will not read any beam_file (but a gaussian'
  write(*,9) '      correction may be performed if -beam is set).'
  write(*,9) '   NOTE : if noise subtraction is active (options -noisecorfile or'
  write(*,9) '   -noiseclfile), the input noise data are assumed to be already corrected'
  write(*,9) '   for the beam.'
  write(*,9) '   Starting with version 2.90 and if linked with HEALPix 3.10+'
  write(*,9) '   file.fits[extname] means extension ''extname'' of file.fits'
  write(*,9)
  write(*,9) '-beam2 [dfloat|'//trim(none)//']('//trim(none)//')'
  write(*,9) '   gaussian beam FWHM (arcmin) for map 2 (if applicable).'
  write(*,9) '   (see beam)'
  write(*,9) '  Note: if the second beam is the same as the first one,'
  write(*,9) '   then -beam2 should be explicitely given the same value as -beam'
  write(*,9)
  write(*,9) '-beam_file2 [name|'//trim(none)//']('//trim(none)//')'
  write(*,9) '   beam window function for map 2 (if applicable).'
  write(*,9) '   If specified, output c_l and c_theta are corrected for it'
  write(*,9) '   The file should be plain ASCII with 2 columns : l and B(l) '
  write(*,9) '   or FITS (ASCII or BINARY) table containing B(l) starting at l=0'
  write(*,9) '   (B(l) is sqrt of beam power spectrum).'
  write(*,9) '   This over-rules the gaussian correction described by -beam2.'
  write(*,9) '   '//trim(none)//' will not read any -beam_file2 (but a gaussian'
  write(*,9) '      correction may be performed if -beam2 is set.'
  write(*,9) '   NOTE : if noise subtraction is active (options -noisecorfile or'
  write(*,9) '   -noiseclfile), the input noise data are assumed to be already corrected'
  write(*,9) '   for the beam.'
  write(*,9) '  Note: if the second beam is the same as the first one,'
  write(*,9) '   then -beam_file2 should be explicitely given the same value as -beam_file'
  write(*,9)
  write(*,9) '-clfile [name|'//trim(default)//'|'//trim(none)//']('//trim(none)//')'
  write(*,9) '   name of output file for Cl'
  write(*,9) '   '//trim(default)//' activates Cl file output with name ' &
 &           //trim(cl_file_name_default)
  write(*,9) '   '//trim(none)//' desactivates output of Cl file'
  write(*,9) '   NOTE : the current format of this file is a two columns array'
  write(*,9) '   the first one corresponding to the value of l, the second one'
  write(*,9) '   to the temperature Cl and, if polarization is activated, five'
  write(*,9) '   additional spectra corresponding to EE, BB, TE, TB, EB '
  write(*,9) '   and optionally ET, BT, BE are also written.'
  write(*,9)
  write(*,9) '-cl_inmap_file [name|'//trim(default)//'|'//trim(none)//']('//trim(none)//')'
  write(*,9) '   name of input file containing the raw Cls of the map as obtained by'
  write(*,9) '   running previously SpICE with the option -cl_outmap_file.'
  write(*,9) '   If this option is active, it is not necessary anymore to read the input'
  write(*,9) '   map file, and option -cl_outmap_file is of course inactive.'
  write(*,9) '   '//trim(default)//' reads the map raw Cls from file ' &
 &           //trim(clmap_file_name_default)
  write(*,9) '   '//trim(none)//' computes the raw Cl directly from the map.'
  write(*,9) '   IMPORTANT NOTE : the cls computed that way can be affected by the nature'
  write(*,9) '      of the masks and/or the weights. If the masks/weights are changed,'
  write(*,9) '      the calculation will become inconsistent.'
  write(*,9)
  write(*,9) '-cl_inmask_file [name|'//trim(default)//'|'//trim(none)//']('//trim(none)//')'
  write(*,9) '   name of input file containing the raw Cl of the weight/masks as'
  write(*,9) '   obtained by running previously SpICE with the option -cl_outmask_file.'
  write(*,9) '   If option -cl_inmask_file is active, the Cls of the masks are not computed,'
  write(*,9) '   they are directly read in file specified by cl_inmask_file and option'
  write(*,9) '   -cl_outmask_file is of course inactive.'
  write(*,9) '   With option -cl_inmask_file active, it is in principle not necessary'
  write(*,9) '   to read masks nor weights, but if masks/weights are not read, it means'
  write(*,9) '   that the input map is supposed to have been processed prior to using spice,'
  write(*,9) '   i.e. it has been multiplied by the masks/weights already.'
  write(*,9) '   '//trim(default)//' reads the weights/mask raw Cls from file ' &
 &           //trim(clmask_file_name_default)
  write(*,9) '   '//trim(none)//' computes the raw Cl directly from the weights/masks.'
  write(*,9) '   IMPORTANT NOTE : the cls computed for the map are affected by the nature'
  write(*,9) '      of the masks and/or the weights. If the masks/weights are changed,'
  write(*,9) '      the calculation will become inconsistent.'
  write(*,9)
  write(*,9) '-cl_outmap_file [name|'//trim(default)//'|'//trim(none)//']('//trim(none)//')'
  write(*,9) '   name of output file containing the raw Cls of the map. This can be'
  write(*,9) '   convenient in order to run SpICE with different parameters, without'
  write(*,9) '   having to recompute these raw Cls.'
  write(*,9) '   '//trim(default)//' activates map raw Cls output with name ' &
 &           //trim(clmap_file_name_default)
  write(*,9) '   '//trim(none)//' desactivates output of map raw Cls file.'
  write(*,9)
  write(*,9) '-cl_outmask_file [name|'//trim(default)//'|'//trim(none)//']('//trim(none)//')'
  write(*,9) '   name of output file containing the raw Cls of the weights/masks.'
  write(*,9) '   This can be convenient in order to run SpICE with different parameters,'
  write(*,9) '   without having to recompute these raw Cls.'
  write(*,9) '   '//trim(default)//' activates mask raw Cls output with name ' &
 &           //trim(clmask_file_name_default)
  write(*,9) '   '//trim(none)//' desactivates output of mask raw Cls file.'
  write(*,9)
  write(*,9) '-corfile [name|'//trim(default)//'|'//trim(none)//']('//trim(none)//')'
  write(*,9) '   name of output correlation file'
  write(*,9) '   '//trim(default)//' activates correlation file output with name ' &
 &           //trim(cor_file_name_default)
  write(*,9) '   '//trim(none)//' desactivates output of correlation file'
  write(*,9) '   NOTE : the current format of this file is a three columns array'
  write(*,9) '   the first one corresponding to the separation theta in radians'
  write(*,9) '   the second one to cos(theta), the third one to the correlation'
  write(*,9) '   function of the temperature and, if polarization is activated,'
  write(*,9) '   five additional columns are written, corresponding respectively to'
  write(*,9) '   QQ, UU, TQ, TU and QU in the coupled case and to EE, BB, TQ, TU and QU'
  write(*,9) '   in the decoupled case (see option -decouple and -symmetric_cl).'
  write(*,9)
  write(*,9) '-covfileout [name|'//trim(default)//'|'//trim(none)//']('//trim(none)//')'
  write(*,9) '   Computes and ouputs the TT C(l) covariance matrices '
  write(*,9) '     These matrices are not needed by the code and are provided for '
  write(*,9) '   eg, cosmological interpretation of the result.'
  write(*,9) '   '//trim(default)//' activates output of the cov. matrices in default'
  write(*,9) '      FITS file '//trim(cov_file_name_default)
  write(*,9) '   '//trim(none)//' desactivates output of the matrices'
  write(*,9)
  write(*,9) '-decouple ['//trim(default)//'|'//trim(none)//']('//trim(none)//')'
  write(*,9) '   computes the decoupled correlation functions'
  write(*,9) '   by an integral over QQ and UU correlation functions in angle space.'
  write(*,9) '   NOTE : this option is active only if -polarization='//trim(default)//','
  write(*,9) '   otherwise it is ignored (see also -tolerance).'
  write(*,9) '   '//trim(default)//' activates decoupling :'
  write(*,9) '      Correlation functions computed are : TT, EE, BB, TQ, TU, QU, (QT, UT, UQ)'
  write(*,9) '      Cls computed are : TT, EE, BB, TE, TB, EB, (ET, BT, BE) '
  write(*,9) '      with no mixing between EE and BB'
  write(*,9) '   '//trim(none)//' desactivates decoupling :'
  write(*,9) '      Correlation functions computed are : TT, QQ, UU, TQ, TU, QU, (QT, UT, UQ)'
  write(*,9) '      Cls computed are : TT, EE, BB, TE, TB, EB, (ET, BT, BE) '
  write(*,9) '      with coupling between EE and BB'
  write(*,9) ' '
  write(*,9) '    In both cases, the TE (and TB) estimators are unbiased using window functions (see'
  write(*,9) '      options -tenormfilein, -tenormfileout, -windowfilein, -windowfileout).'
  write(*,9) 
  write(*,9) '-dry ['//trim(default)//'|'//trim(none)//']('//trim(none)//')'
  write(*,9) '   This unsafe option allows reading partly corrupted files, i.e. with'
  write(*,9) '   incomplete headers. However the size of the fits file must be correct,'
  write(*,9) '   i.e. it must correspond to a full sky map with 12*nside^2 pixels.'
  write(*,9)
  write(*,9) '-extramapfile [name|'//trim(none)//']('//trim(none)//')'
  write(*,9) '   name of the (optional) extra input map file '//hpx_fits_file
  write(*,9) '   to be added to mapfile before analysis. '
  write(*,9) '   Must be in the some Coordinates, Nside and Units as mapfile'
  write(*,9) '   (ordering can be different).'
  write(*,9)
  write(*,9) '-extramapfile2 [name|'//trim(none)//']('//trim(none)//')'
  write(*,9) '   name of the (optional) extra input map file '//hpx_fits_file
  write(*,9) '   to be added to mapfile2 before analysis. '
  write(*,9) '   Must be in the some Coordinates, Nside and Units as mapfile2'
  write(*,9) '   (ordering can be different).'
  write(*,9)
  write(*,9) '-fits_out ['//trim(default)//'|'//trim(none)//']('//trim(none)//')'
  write(*,9) '   If set, the Cl and correlation (-clfile, -corfile, -cl_outmap_file, -cl_outmask_file)'
  write(*,9) '   are written in FITS Ascii files, instead of plain text files.'
  write(*,9) '   This does not affect the input Cl files (-cl_inmap_file and -cl_inmask_file),'
  write(*,9) '   which will be read whatever is their format.'
  write(*,9)
  write(*,9) '-help'
  write(*,9) '   print this message and stop'
  write(*,9)
  write(*,9) '-history'
  write(*,9) '   print history of the program and stop'
  write(*,9) 
  write(*,9) '-kernelsfileout [name|'//trim(default)//'|'//trim(none)//']('//trim(none)//')'
  write(*,9) '   Compute and ouput the 4 kernels, each of size (lmax+1)*(2*lmax+1),'
  write(*,9) "   relating the Spice (polarized) estimator to the 'true' underlying CMB spectra."
  write(*,9) '     These kernels are not needed by the code and are provided for '
  write(*,9) '   eg, cosmological interpretation of the result.'
  write(*,9) '     Note that these kernels depend on the choice'
  write(*,9) '   of apodization (options -apodizetype and -apodizesigma) '
  write(*,9) '   and on -thetamax and -nlmax.'
  write(*,9) '     How to use them depends on the -decouple option '
  write(*,9) '   (see the FITS file header for more information).'
  write(*,9) '     This option is inactive if option -windowfilein is active.'
  write(*,9) '   '//trim(default)//' activates output of the kernels in default'
  write(*,9) '      FITS file '//trim(kernels_file_name_default)
  write(*,9) '   '//trim(none)//' desactivates output of the kernels'
  write(*,9)
  write(*,9) '-listmapfiles1_K [name|'//trim(none)//']('//trim(none)//')'
  write(*,9) '   with K integer in [1, '//trim(string(nefmax,format='(i2)'))//']'
  write(*,9) '   list of map files to be combined with weights -listmapweights1_K'
  write(*,9) '   to form the first map to analyze.'
  write(*,9) '   If this option is set, -mapfile (and -extramapfile) are ignored'
  write(*,9)
  write(*,9) '-listmapfiles2_K [name|'//trim(none)//']('//trim(none)//')'
  write(*,9) '   with K integer in [1, '//trim(string(nefmax,format='(i2)'))//']'
  write(*,9) '   list of map files to be combined with weights -listmapweights2_K'
  write(*,9) '   to form the optional second map to cross-correlate with the first one.'
  write(*,9) '   If this option is set, -mapfile2 (and -extramapfile2) are ignored'
  write(*,9)
  write(*,9) '-listmapweights1_K [dfloat|'//trim(none)//'](1.0)'
  write(*,9) '   with K integer in [1, '//trim(string(nefmax,format='(i2)'))//']'
  write(*,9) '   list of weights to be applied to -listmapfiles1_K above'
  write(*,9) '   to form the first map to analyze.'
  write(*,9) '   m1 = w1(1) * m1(1) + w1(2) * m1(2) + ... '
  write(*,9)
  write(*,9) '-listmapweights2_K [dfloat|'//trim(none)//'](1.0)'
  write(*,9) '   with K integer in [1, '//trim(string(nefmax,format='(i2)'))//']'
  write(*,9) '   list of weights to be applied to -listmapfiles2_K above'
  write(*,9) '   to form the second optional map to analyze.'
  write(*,9) '   m2 = w2(1) * m2(1) + w2(2) * m2(2) + ... '
  write(*,9)
  write(*,9) '-mapfile [name|'//trim(default)//']('//trim(default)//')'
  write(*,9) '   name of the input map file '//hpx_fits_file
  write(*,9) '   '//trim(default)//' is equivalent to take input file ' &
 &           //trim(mapfile_default)//'.fits'
  write(*,9) '   Starting with version 2.90 and if linked with HEALPix 3.10+'
  write(*,9) '   file.fits[extname] means extension ''extname'' of file.fits'
  write(*,9)
  write(*,9) '-mapfile2 [name|'//trim(default)//'|'//trim(none)//']('//trim(none)//')'
  write(*,9) '   name of the second (optional) input map file '//hpx_fits_file
  write(*,9) '   to cross-correlate with the first one (-mapfile)'
  write(*,9) '   If provided, ONLY cross-correlation is performed (no auto-correlation)'
  write(*,9) '   '//trim(default)//' is equivalent to take input file ' &
 &           //trim(mapfile2_default)//'.fits'
  write(*,9) '   '//trim(none)//' turns off map-map cross-correlation'
  write(*,9)
  write(*,9) '-maskfile [name|'//trim(default)//'|'//trim(none)//']('//trim(none)//')'
  write(*,9) '   name of the input mask file '//hpx_fits_file
  write(*,9) '   A mask file is composed of zeros (masked) and ones (in)'
  write(*,9) '   '//trim(default)//' activates maskfile with name ' &
 &           //trim(maskfile_default)//'.fits'
  write(*,9) '   '//trim(none)//' desactivates -maskfile'
  write(*,9) '   Starting with version 2.90 and if linked with HEALPix 3.10+'
  write(*,9) '   file.fits[extname] means extension ''extname'' of file.fits'
  write(*,9)
  write(*,9) '-maskfile2 [name|'//trim(default)//'|'//trim(none)//']('//trim(none)//')'
  write(*,9) '   name of the input mask file '//hpx_fits_file
  write(*,9) '   to be applied to the second map'
  write(*,9) '   '//trim(default)//' activates -maskfile2 with name ' &
 &           //trim(maskfile2_default)//'.fits'
  write(*,9) '   '//trim(none)//' desactivates -maskfile2'
  write(*,9) ' Notes:'
  write(*,9) '  1) if the you want the 2nd mask to be the same as the 1st one, '
  write(*,9) '    set -maskfile2 to the same name as -maskfile'
  write(*,9) '  2) the 2 masks should overlap at least partially'
  write(*,9) 
  write(*,9) '-maskfilep [name|'//trim(default)//']('//trim(default)//')'
  write(*,9) '   name of the input mask file '//hpx_fits_file
  write(*,9) '   to be applied to polarization.'
  write(*,9) '   By default, it is the same as -maskfile'
  write(*,9) '   Ignored if polarization is set to NO'
  write(*,9)
  write(*,9) '-maskfilep2 [name|'//trim(default)//']('//trim(default)//')'
  write(*,9) '   name of the input mask file '//hpx_fits_file
  write(*,9) '   to be applied to the second map polarization'
  write(*,9) '   By default, it is the same as -maskfile2'
  write(*,9) '   Ignored if polarization is set to NO'
  write(*,9) 
  write(*,9) '-nlmax [integer](3*Nside-1)'
  write(*,9) '   value of the largest multipole used in the analysis.'
  write(*,9) '   The default is 3*Nside-1, where Nside is the Healpix resolution parameter '
  write(*,9) '   of the map being analyzed.'
  write(*,9) '   Any value larger than 3*Nside-1 or <= 0 will be replaced by the default.'
  write(*,9) 
  write(*,9) '-noiseclfile [name|'//trim(default)//'|'//trim(none)//']('//trim(none)//')'
  write(*,9) '   name of input file for Cls of the noise'
  write(*,9) '   if -noiseclfile is specified, noise subtraction is applied'
  write(*,9) '   to both final Cls and correlation function (unless -noisecorfile is'
  write(*,9) '   specified for the latter)'
  write(*,9) '   '//trim(default)//' activates noise subtraction from file ' &
 &           //trim(noisecl_file_name_default)
  write(*,9) '   '//trim(none)//' desactivates noise subtraction from a Cl noise file'
  write(*,9)
  write(*,9) '-noisecorfile [name|'//trim(default)//'|'//trim(none)//']('//trim(none)//')'
  write(*,9) '   name of input file for correlation function of the noise'
  write(*,9) '   if -noisecorfile is specified, noise subtraction is applied'
  write(*,9) '   to both final correlation function and Cls (unless -noiseclfile is'
  write(*,9) '   specified for the latter)'
  write(*,9) '   '//trim(default)//' activates noise subtraction from file ' &
 &           //trim(noisecor_file_name_default)
  write(*,9) '   '//trim(none)//' desactivates noise subtraction from a correlation noise file'
  write(*,9)
  write(*,9) '-normfac [dfloat|'//trim(none)//']('//trim(none)//')'
  write(*,9) '   normalization factor. If normfac is specified, the correlation'
  write(*,9) '   function and the Cls are multiplied by normfac**2'
  write(*,9) '   '//trim(none)//' and -normfac=1.0 are equivalent : in that case'
  write(*,9) '   normalization is inactive to save calculations.'
  write(*,9) '   NOTE : if noise subtraction is active (options -noisecorfile or'
  write(*,9) '   -noiseclfile), the input noise data are assumed to be already'
  write(*,9) '   normalized appropriately'
  write(*,9) 
!   write(*,9) '-npairsthreshold [dfloat|'//trim(none)//']('//trim(none)//')'
!   write(*,9) '   thresholding value (>=0) for the available number npairs of pairs in a bin'
!   write(*,9) '   if npairs <= npairsthreshold in a bin then the correlation function'
!   write(*,9) '   is set to zero in this bin.'
!   write(*,9) '   '//trim(none)//' and npairsthreshold=0 are equivalent.'
!   write(*,9) '   In that case, the correlation function is set to zero only in the bins'
!   write(*,9) '   where xi_mask = 0 exactly to avoid numerical divergences'
!   write(*,9) '   and if that happens, a warning is written if verbosity=1 or 2.'
!   write(*,9) '   NOTE 1 : if noise subtraction is active (options -noisecorfile or'
!   write(*,9) '   -noiseclfile), the input noise data are assumed to be appropriately'
!   write(*,9) '   thresholded'
!   write(*,9) '   NOTE 2 : this option is active only if there is an input mask or weight'
!   write(*,9) '   file'
!   write(*,9)
  write(*,9) '-optinfile [name|'//trim(default)//'|'//trim(none)//']('//trim(none)//')'
  write(*,9) '   Read the options in an file (only those which are not specified'
  write(*,9) '   interactively are taken into account).'
  write(*,9) '   '//trim(default)//' reads default input option file ' &
 &               //trim(spice_rc_file_default)
  write(*,9) '   '//trim(none)//' desactivates the reading.'
  write(*,9) '   NOTE 1 : The file containing the options is a 2-column file,'
  write(*,9) '   the first column corresponds to the option (without the "-" at the'
  write(*,9) '   beginning) and the second column to the parameter (see below for'
  write(*,9) '   an example).'
  write(*,9) '   NOTE 2 : It is not necessary to specify all the options in this'
  write(*,9) '   file, and the order in which they appear does not matter. If an'
  write(*,9) '   option is specified neither in the option file nor interactively,'
  write(*,9) '   then the default value is assumed for this option.'
  write(*,9)
  write(*,9) '-optoutfile [name|'//trim(default)//'|'//trim(none)//']('//trim(none)//')'
  write(*,9) '   Output the options in an file (read with option -optinfile or specified'
  write(*,9) '   interactively. In that case, no calculation is performed, only the'
  write(*,9) '   file containing the options is created.'
  write(*,9) '   '//trim(default)//' output option file with default name ' &
 &               //trim(spice_rc_file_default)
  write(*,9) '   '//trim(none)//' desactivates the output'
  write(*,9) '   NOTE : see note of -optinfile.'
  write(*,9)
  write(*,9) '-overwrite ['//trim(default)//'|'//trim(none)//']('//trim(default)//')'
  write(*,9) '   Activates overwriting for all output files.'
  write(*,9) '   '//trim(default)//' overwrite all the output files generated previously.'
  write(*,9) '   '//trim(none)//' checks if any of the output files already exist. If it is the'
  write(*,9) '   case, SpICE stops without mercy.'
  write(*,9) '   To have more sophisticated handling of overwriting, the standard'
  write(*,9) '   Healpix syntax can be used : file names in the form !name means'
  write(*,9) '   that the name of the file to be used is name and the "!" explicitely allows'
  write(*,9) '   overwriting.'
  write(*,9)
  write(*,9) '-pixelfile [name|'//trim(default)//'|'//trim(none)//']('//trim(none)//')'
  write(*,9) '   name of input pixel correction file for 1st map'
  write(*,9) '   '//trim(default)//' activates pixel window correction with file '
  write(*,9) '   ${HEALPIXDATA}/'//trim(pixwin_file_name_default)//'xxxx.fits'
  write(*,9) '   or'
  write(*,9) '   ${HEALPIX}/data/'//trim(pixwin_file_name_default)//'xxxx.fits'
  write(*,9) '   where xxxx corresponds to the value of nside of input map file.'
  write(*,9) '   '//trim(none)//' desactivates pixel window correction.'
  write(*,9) '   NOTE : if noise subtraction is active (options -noisecorfile or'
  write(*,9) '   -noiseclfile), the input noise data are assumed to be already corrected'
  write(*,9) '   for pixel window'
  write(*,9)
  write(*,9) '-pixelfile2 [name|'//trim(default)//'|'//trim(none)//'](pixelfile)'
  write(*,9) '   name of input pixel correction file for 2nd map (for cross-analysis)'
  write(*,9) ' Default value is the one chosen for -pixelfile.'
  write(*,9) ' If set to '//trim(default)//' explicitly, or implicitly via -pixelfile, '
  write(*,9) ' the standard pixel window function adapted to the 2nd map Nside will be used.'
  write(*,9) ' It must be set to '//trim(none)//' to prevent pixel correction of 2nd map.'
  write(*,9)
  write(*,9) '-polarization ['//trim(default)//'|'//trim(none)//']('//trim(none)//')'
  write(*,9) '   to compute polarized power spectra/correlation functions'
  write(*,9) '   in addition to standard temperature Cls/correlation'
  write(*,9) '   function. The input map file is supposed to be in the standard'
  write(*,9) '   polarized Healpix fits format, so it contains 3 maps (T,Q,U)'
  write(*,9) '   instead of one.'
  write(*,9) '   The polarized Cls and correlation functions computed depend on the'
  write(*,9) '   choice of parameter decouple explained below (see option -decouple).'
  write(*,9) 
  write(*,9) '-subav ['//trim(default)//'|'//trim(none)//']('//trim(none)//')'
  write(*,9) '   in case the input map is expected to have an offset, its (weighted)'
  write(*,9) '   average can be subtracted from it prior to the calculations,'
  write(*,9) '   in order to insure zero average.'
  write(*,9) '   When there are masks/weights, W(i), this average is computed'
  write(*,9) "   simply as [sum_i W(i) map(i)]/[sum_i W(i)]."
  write(*,9) '   See options -maskfile, -weightfile, -weightpower for more information'
  write(*,9) '   on the masks/weights.'
  write(*,9) '   It is important to notice that the estimator of the correlation function'
  write(*,9) '   or the Cls is not expected anymore to be unbiased if option -subav'
  write(*,9) '   is active. It would be the case only if the offset would be exactly' 
  write(*,9) '   known and not estimated directly from the data.'
  write(*,9) '   '//trim(default)//' activates average subtraction.'
  write(*,9) '   '//trim(none)//' desactivates average subtraction.'
  write(*,9)
  write(*,9) '-subdipole ['//trim(default)//'|'//trim(none)//']('//trim(none)//')'
  write(*,9) '   the best fit monopole and dipole will regressed out from the'
  write(*,9) '   (masked and weighted) data before spectral analysis. '
  write(*,9) '   '//trim(default)//' activates dipole and monopole subtraction.'
  write(*,9) '   '//trim(none)//' desactivates dipole and monopole subtraction.'
  write(*,9) '   -subdipole=YES implies -subav=YES and has the same caveats.'
  write(*,9)
  write(*,9) '-symmetric_cl ['//trim(default)//'|'//trim(none)//']('//trim(none)//')'
  write(*,9) '   by default, when 2 maps are provided, and -polarization is set,'
  write(*,9) '   cross spectra such as TE are defined as T_map1 * E_map2 '
  write(*,9) '  (in which case it is recommended to use the least noisy map as map2)'
  write(*,9) '   and 9 spectra are produced.'
  write(*,9) '   If -symmetric_cl is set to '//trim(default)//', then'
  write(*,9) '   TE = (T_map1 * E_map2 + T_map2 * E_map1)/2'
  write(*,9) '   and 6 spectra are produced.'
  write(*,9) '  In the single map case, TE=ET and only 6 spectra are produced.'
  write(*,9)
  write(*,9) '-tenormfilein [name|'//trim(default)//'|'//trim(none)//']('//trim(none)//')'
  write(*,9) '   this option is ignored if -polarization is not active :'
  write(*,9) '   reads the normalisation function needed to compute correctly the TE (and TB) spectra'
  write(*,9) '   and correlation functions, instead of computing it. '
  write(*,9) '   The calculation of this term depends on the apodization chosen '
  write(*,9) '   (options -apodizetype and -apodizesigma).'
  write(*,9) '   '//trim(default)//' activates input of the normalisation function with default'
  write(*,9) '      input file name '//trim(tenorm_file_name_default)
  write(*,9) '   '//trim(none)//' desactivates input of the normalisation function : it is computed'
  write(*,9) '      on run time.'
  write(*,9)
  write(*,9) '-tenormfileout [name|'//trim(default)//'|'//trim(none)//']('//trim(none)//')'
  write(*,9) '   Compute and ouput the normalisation function needed to compute the  '
  write(*,9) '   TE (and TB) spectra. Note that this function depends on the choice'
  write(*,9) '    of apodization (options -apodizetype and -apodizesigma).'
  write(*,9) '   This option is inactive, of course, if option -tenormfilein is active, or'
  write(*,9) '   if option -polarization is inactive.'
  write(*,9) '   '//trim(default)//' activates output of the normalisation function in default'
  write(*,9) '      file '//trim(tenorm_file_name_default)
  write(*,9) '   '//trim(none)//' desactivates output of the normalisation function (but it is'
  write(*,9) '      computed if needed).'
  write(*,9)
  write(*,9) '-tf_file [name|'//trim(none)//']('//trim(none)//')'
  write(*,9) '   name of the input tranfer function file.'
  write(*,9) '   The output C(l) will be divided by the transfer function(s) read'
  write(*,9) '   from this file.'
  write(*,9) '   The transfer function is defined as the ratio of 2 power spectra.'
  write(*,9) '   The file should have the same format as the C(l) produced by this code.'
  write(*,9)
  write(*,9) '-thetamax [dfloat|'//trim(none)//']('//trim(none)//')'
  write(*,9) '   maximum value of angle theta used in the integrals to compute the'
  write(*,9) '   power spectra from the correlation functions.'
  write(*,9) '   '//trim(none)//' desactivates upper bound (sets it to 180 degrees).'
  write(*,9)
  write(*,9) '-tolerance [dfloat|'//trim(none)//']('//trim(none)//')'
  write(*,9) '   relative tolerance on convergence of QQ and UU integrals (see -decouple)'
  write(*,9) '   Active only when -decouple is set'
  write(*,9) '   Smaller values give more accurate EE and BB spectra, in longer times.'
  write(*,9) '   '//trim(none)//' leaves it to its default value (1.e-5).'
  write(*,9) 
  write(*,9) '-usage'
  write(*,9) '   print short usage reminder'
  write(*,9)
  write(*,9) '-verbosity [0|1|2|'//trim(default)//'|'//trim(none)//'](1)'
  write(*,9) '   verbose option :'
  write(*,9) '   0 or '//trim(none)//' : no verbose'
  write(*,9) '   1 or '//trim(default)//' : standard verbose'
  write(*,9) '   2 : full verbose'
  write(*,9) 
  write(*,9) '-version'
  write(*,9) '   prints version number and exits'
  write(*,9)
  write(*,9) '-weightfile [name|'//trim(default)//'|'//trim(none)//']('//trim(none)//')'
  write(*,9) '   name of the input weight file '//hpx_fits_file
  write(*,9) '   Weights should be always positive or equal to zero.'
  write(*,9) '   If masks are present, weights are multiplied by the masks.'
  write(*,9) '   If masks are not present, weights are supposed to be already masked,'
  write(*,9) '   i.e. are equal to zero for masked pixels'
  write(*,9) '   '//trim(default)//' activates -weightfile with name ' &
 &           //trim(weightfile_default)//'.fits'
  write(*,9) '   '//trim(none)//' desactivates -weightfile'
  write(*,9)
  write(*,9) '-weightfile2 [name|'//trim(default)//'|'//trim(none)//']('//trim(none)//')'
  write(*,9) '   name of the input weight file '//hpx_fits_file
  write(*,9) '   to be applied to the second map'
  write(*,9) '   '//trim(default)//' activates -weightfile2 with name ' &
 &           //trim(weightfile2_default)//'.fits'
  write(*,9) '   '//trim(none)//' desactivates -weightfile2'
  write(*,9)
  write(*,9) '-weightfilep [name|'//trim(default)//']('//trim(default)//')'
  write(*,9) '   name of the input weight file '//hpx_fits_file
  write(*,9) '   to be applied to polarization (Q and U Stokes parameters).'
  write(*,9) '   By default, it is the same as -weightfile'
  write(*,9) '   Ignored if polarization is set to NO'
  write(*,9)
  write(*,9) '-weightfilep2 [name|'//trim(default)//']('//trim(default)//')'
  write(*,9) '   name of the input weight file '//hpx_fits_file
  write(*,9) '   to be applied to the second map polarization'
  write(*,9) '   By default, it is the same as -weightfile2'
  write(*,9) '   Ignored if polarization is set to NO'
  write(*,9) 
  write(*,9) '-weightpower [dfloat|'//trim(none)//'](1.0)'
  write(*,9) '   the weights are taken to this power.'
  write(*,9) '   if -weightpower = 0.0, the weights are transformed into pure masks'
  write(*,9) '   if -weightpower < 0.0, only the non zero values of the weights are'
  write(*,9) '   modified.'
  write(*,9) '   '//trim(none)//' is equivalent to -weightpower=1.0.'
  write(*,9)
  write(*,9) '-weightpower2 [dfloat|'//trim(none)//'](1.0)'
  write(*,9) '   the weights of the second map are taken to this power.'
  write(*,9) '   '//trim(none)//' is equivalent to -weightpower2=1.0.'
  write(*,9)
  write(*,9) '-weightpowerp [dfloat|'//trim(default)//']('//trim(default)//')'
  write(*,9) '   the weights applied to 1st map (Q,U) are taken to this power.'
  write(*,9) '   Default: same as 1st Temperature weight power (see -weightpower)'
  write(*,9)
  write(*,9) '-weightpowerp2 [dfloat|'//trim(default)//']('//trim(default)//')'
  write(*,9) '   the weights applied to 2nd map (Q,U) are taken to this power.'
  write(*,9) '   Default: same as 2nd Temperature weight power (see -weightpower2)'
  write(*,9)
  write(*,9) '-windowfilein [name|'//trim(default)//'|'//trim(none)//']('//trim(none)//')'
  write(*,9) '   this option is ignored if polarization is not active :'
  write(*,9) '   reads the kcross term needed to compute correctly the TE (and TB) spectra'
  write(*,9) '   and correlation functions, instead of computing it. '
  write(*,9) '   The calculation of this term depends on the apodization chosen '
  write(*,9) '   (options -apodizetype and -apodizesigma).'
  write(*,9) '   '//trim(default)//' activates input of the window function with default'
  write(*,9) '      input file name '//trim(window_file_name_default)
  write(*,9) '   '//trim(none)//' desactivates input of the window function : it is computed'
  write(*,9) '      on run time.'
  write(*,9)
  write(*,9) '-windowfileout [name|'//trim(default)//'|'//trim(none)//']('//trim(none)//')'
  write(*,9) '   Compute and ouput the kcross term needed to compute the  '
  write(*,9) '   TE (and TB) spectra. Note that this kcross term depends on the choice'
  write(*,9) '    of apodization (options -apodizetype and -apodizesigma).'
  write(*,9) '   This option is inactive, of course, if option -windowfilein is active, or'
  !write(*,9) '   if options -polarization is inactive, or if -kernelsfileout is active.'
  write(*,9) '   if option -polarization is inactive.' ! 2017-08-09, allow windowfileout+kernelsfileout
  write(*,9) '   '//trim(default)//' activates output of the window function in default'
  write(*,9) '      file '//trim(window_file_name_default)
  write(*,9) '   '//trim(none)//' desactivates output of the window function (but it is'
  write(*,9) '      computed if needed).'
  write(*,9)
  write(*,9)
  write(*,9) 'IMPORTANT NOTE'
  write(*,9)
  write(*,9) 'The order of the operations is not the same whether mask/weights are present'
  write(*,9) 'or not.'
  write(*,9)
  write(*,9) 'If there is no mask/weight, spice is nearly equivalent to anafast, and'
  write(*,9) 'the order of operations is the following :'
  write(*,9) '   (1a) map average subtraction, if required'
  write(*,9) '   (2a) calculation of Cl'
  write(*,9) '   (3a) pixel window correction, beam correction, transfer function correction,'
  write(*,9) '       renormalization, if required'
  write(*,9) '   (4a) calculation of the 2pt correlation function, xi, if needed'
  write(*,9) '   (5a) apodization of xi and calculation of the final Cl, if  needed'
  write(*,9) '   (6a) subtraction of noise bias on Cl and/or xi, if applicable'
  write(*,9) 
  write(*,9) 'If there is any mask/weight, the order of operations is the following'
  write(*,9) '   (1b) multiplication of the map by the masks/weights'
  write(*,9) '   (2b) subtraction of the (weighted) average'
  write(*,9) '   (3b) calculation of the Cl of the map and of the masks/weights'
  write(*,9) '   (4b) calculation of the correlation function of the map and of the'
  write(*,9) '       masks/weights, and the corresponding ratio, xi'
!  write(*,9) '   (5b) pair count thresholding, if required'
  write(*,9) '   (5b) decoupling if needed (polarization case)'
  write(*,9) '   (6b) apodization of xi, if required'
  write(*,9) '   (7b) renormalization, if required'
  write(*,9) '   (8b) calculation of Cl from xi, if needed'
  write(*,9) '   (9b) pixel window correction, beam correction, transfer function correction,'
  write(*,9) '       if required, to obtain the final Cl'
  write(*,9) '  (10b) calculation of final xi, if needed'
  write(*,9) '  (11b) subtraction of noise bias on Cl and/or xi, if applicable'
  write(*,9)
  write(*,9) 'Note, therefore, that if there are masks/weights, beam correction and'
  write(*,9) 'pixel window correction are performed AFTER apodization, at variance with'
  write(*,9) 'the case without mask/weight. Such ordering in the procedure is expected'
  write(*,9) 'to introduce a small bias in the final results. Ideally, beam and pixel'
  write(*,9) 'window correction should be performed before apodization, but this is'
  write(*,9) 'possible only prior to masking/weighting on a FULL SKY map, i.e. before'
  write(*,9) 'operation (1b).'
  write(*,9)
  write(*,9) 'EXAMPLES'
  write(*,9)
  write(*,9) 'spice'
  write(*,9) 'spice -mapfile '//trim(default)//' -corfile '//trim(default)// &
 &            ' -clfile '//trim(default)
  write(*,9) '   reads default input map, computes correlation function'
  write(*,9) '   and Cl and output them on default output files'
  write(*,9) '   no masks/weights are taken into account'
  write(*,9)
  write(*,9) 'spice -mapfile mymap.fits -weightfile weights.fits -clfile '//trim(none)
  write(*,9) '   reads input map mymap.fits, input file weights.fits for weights,'
  write(*,9) '   output correlation functions only'
  write(*,9) 
  write(*,9) 'spice -mapfile mymap.fits -maskfile mymasks.fits '
  write(*,9) '      -weightfile myweights.fits -beam 18.d0'
  write(*,9) '   reads input map mymap.fits, input mask file mymasks.fits,'
  write(*,9) '   input weight file myweights.fits, computes Cl and' 
  write(*,9) '   correlation function with correction for a beam of'
  write(*,9) '   18 arcmin, output them on default output files. If both'
  write(*,9) '   mask and weight files are specified, the weights are'
  write(*,9) '   multiplied by the masks, allowing one to read a full unmasked'
  write(*,9) '   weight file'
  write(*,9)
  write(*,9) 'spice -mapfile mymap.fits -weightfile myweights.fits  '
  write(*,9) '      -noisecorfile mynoise.cor -noiseclfile mynoise.cl'
  write(*,9) '   reads input map mymap.fits, input weight file myweights.fits,'
  write(*,9) '   computes Cls and correlation function, subtracts from Cls'
  write(*,9) '   the Cls of the noise from file mynoise.cl, subtracts from'
  write(*,9) '   the correlation function the correlation function of the noise'
  write(*,9) '   from mynoise.cor, output the results for the Cls and the'
  write(*,9) '   correlation function on default output files.'
  write(*,9) '   Since only the weights are specified, they are supposed'
  write(*,9) '   to be appropriately masked.'
  write(*,9) '   Even if only one noise file is specified, either from the'
  write(*,9) '   option -noisecorfile or from the option -noiseclfile,'
  write(*,9) '   noise correction is applied to both Cls and the the'
  write(*,9) '   correlation function. But this latter way is less accurate'
  write(*,9) '   because of the supplementary Legendre transform it involves.'
  write(*,9) 
  write(*,9) 'spice -mapfile mymap.fits -polarization '//trim(default)//' \ '
  write(*,9) '      -maskfile mymasks.fits'
  write(*,9) '   reads multiple input mapfile mymap.fits (T,Q,U), computes'
  write(*,9) '   the temperature and polarization Cls (TT, GG, CC, TG, TC and GC) and'
  write(*,9) '   correlation functions (TT, QQ, UU, TQ, TU, QU), and writes the'
  write(*,9) '   results in default output files. All the calculations use'
  write(*,9) '   appropriate mask corrections.'
  write(*,9)
  write(*,9) 'spice -optinfile YES'
  write(*,9) 'spice -optinfile '//trim(spice_rc_file_default)
  write(*,9) '   Read options and parameters in default file '//trim(spice_rc_file_default)
  write(*,9) '   and run spice accordingly.'
  write(*,9) 
  write(*,9) 'EXAMPLE OF OPTIONS INPUT FILE ('//trim(spice_rc_file_default)//') :'
  do j=jstartnormal,noption
     write(*,9) trim(option(j)(2:))//' = '//trim(option_param(j))
  enddo
  write(*,9) 
9 format(a)
  
  STOP
end subroutine man

!=======================================================================
subroutine print_version
!=======================================================================
  use spice_parameters, only: version
  print*,trim(version)
  stop
end subroutine print_version
!=======================================================================
subroutine about
!=======================================================================
  use spice_parameters
  implicit none
  write(*,9)
  write(*,9) 'SpICE : Spatially Inhomogeneous Correlation Estimator'
  write(*,9)
  write(*,9) '                       Version '//trim(version)
  write(*,9)
  write(*,9) '(c) 2001 I. Szapudi [3], S. Prunet [1,*], S. Colombi [1,*]'
  write(*,9) '    No polarization'
  write(*,9)
  write(*,9) '(c) 2002 S. Prunet [1,*], E. Hivon [1,2], ' 
  write(*,9) '         S. Colombi [1,*], I. Szapudi [3]'
  write(*,9) '         A. Challinor [4]. G. Chon [4]'
  write(*,9) '    Polarization included'
  write(*,9) 
  write(*,9) '    email : prunet@iap.fr'
  write(*,9) '            hivon@iap.fr'
  write(*,9) '            colombi@iap.fr'
  write(*,9) '            szapudi@IfA.Hawaii.Edu'
  write(*,9) '            a.d.challinor@mrao.cam.ac.uk'
  write(*,9) '            g.chong@mrao.cam.ac.uk'
  write(*,9) 
  write(*,9) '[1] Institut d''Astrophysique de Paris'
  write(*,9) '    98 bis boulevard Arago, F-75014, Paris, France'
  write(*,9)
  write(*,9) '[2] IPAC/Caltech'
  write(*,9) '    Mail Code 100-22, 770 S. Wilson Av., Pasadena, CA 91125, USA'
  write(*,9)
  write(*,9) '[3] Institude For Astronomy, University of Hawaii'
  write(*,9) '    2680 Woodlawn Dr. Honolulu, HI 96822'
  write(*,9)
  write(*,9) '[4] Astrophysics Group, Cavendish Laboratory'
  write(*,9) '    Madingley Road, CAmbridge CB3 0HE, United Kingdom'
  write(*,9)
  write(*,9) '[*] NIC (Numerical Investigations in Cosmology), CNRS'
  write(*,9)
  write(*,9)
  write(*,9) 'The software can be downloaded from http://www2.iap.fr/users/hivon/software/PolSpice'
  write(*,9)
  write(*,9)
  write(*,9) 'LICENSE INFORMATION'
  write(*,9) 
  write(*,9) 'This sofware is an extension of the HEALPix package'
  write(*,9) '(c) 2000 by K.M. Gorski, B.D. Wandelt, E. Hivon, F. Hansen,'
  write(*,9) '    A.J. Banday, M. Bartelmann'
  write(*,9) '(c) 2003 by K.M. Gorski, E. Hivon, M. Reinecke, A.J. Banday,' 
  write(*,9) '    B.D. Wandelt, and M. Bartelmann'  
  write(*,9) '(c) 2005 by K.M. Gorski, E. Hivon, A.J. Banday, B.D. Wandelt,'
  write(*,9) '    F. Hansen, M. Reinecke, M. Bartelmann  (ApJ 622 759)'
  write(*,9) 'and is thus subject to the same restrictions.'
  write(*,9) 'For more information on the HEALPix package see'
  write(*,9) 'http://healpix.sourceforge.net'
  write(*,9)
  write(*,9)
  write(*,9) 'DISCLAIMER'  
  write(*,9) 
  write(*,9) 'The SpICE software is provided "as is", and the authors'
  write(*,9) '(SC, EH, SP, IS, GC, AC) make no representations or warranties,'
  write(*,9) 'expressed or implied, of merchantability or fitness for any particular'
  write(*,9) 'purpose, and disclaim all liability for direct or consequential damages'
  write(*,9) 'resulting from anyone''s use of the program.'
9 format(a)
   
  STOP
end subroutine about

!=======================================================================
subroutine history
!=======================================================================
  write(*,9) 'HISTORY :'
  write(*,9) '+++++++++'
  write(*,9)
  write(*,9) 'UNPOLARIZED SPICE :'
  write(*,9) 
  write(*,9) '07/01/01/SP&IS       : version 0.0'
  write(*,9) '                       - Basic version of the code, xi and Cl'
  write(*,9) '07/15/01/SC&SP&IS    : version 0.1'
  write(*,9) '                       - Improve input/output'
  write(*,9) '                       - Build friendly interface'
  write(*,9) '                       - pixel window and beam window correction'
  write(*,9) '11/16/01/SC&SP       : version 0.2'
  write(*,9) '                       - add pair count thresholding, option'
  write(*,9) '                         npairsthreshold'
  write(*,9) '                       - add apodizationg of the 2pt correlation'
  write(*,9) '                         function, option apodizesigma'
  write(*,9) '                       - add unsafe reading option, option dry'
  write(*,9) '22/05/02/SC&SP       : version 1.0'
  write(*,9) '                       - correct for a bug : the average from'
  write(*,9) '                         the whole map was always subtracted'
  write(*,9) '                         from it. Instead, an option subav'
  write(*,9) '                         is added, where there the (weighted)'
  write(*,9) '                         average is subtracted from the map'
  write(*,9) '                         if required'
  write(*,9) '07/06/02/SC&SP       : version 1.01'
  write(*,9) '                       - correct for a bug in the handling'
  write(*,9) '                         of the flags coroutput, apodize,'
  write(*,9) '                         verbose and megaverbose in subroutine'
  write(*,9) '                         compute.'
  write(*,9) '                       - correct for a bug in the wrapper'
  write(*,9) '                         normfac option was here twice.'
  write(*,9) 
  write(*,9) 'POLARIZED SPICE :'
  write(*,9) 
  write(*,9) '01/04/02/SC&EH&SP&IS : version 2.0'
  write(*,9) '                       - add option polarization to compute'
  write(*,9) '                         polarization Cls and correlation functions'
  write(*,9) '07/06/02/EH&SP       : version 2.01'
  write(*,9) '                       - add decouple option to compute pure'
  write(*,9) '                         E and B modes correlation functions'
  write(*,9) '                         and power spectra'
  write(*,9) '07/06/02/SC&SP       : version 2.02'
  write(*,9) '                       - clean up wrapper.'
  write(*,9) '07/03/03/EH          : version 2.1'
  write(*,9) '                       - increase compatibility with Healpix 1.20'
  write(*,9) '                       - introduce nlmax, beam_file, apodizetype'
  write(*,9) '10/15/04/SC&SP       : version 2.2'
  write(*,9) '                       - clean up the code'
  write(*,9) '                       - introduce cl_outmap_file, cl_inmap_file,'
  write(*,9) '                         cl_outmask_file, cl_inmask_file,'
  write(*,9) '                         windowfileout, windowfilein.'
  write(*,9) '                       - PIOLIB compatibility              [Planck]'
  write(*,9) '12/04/06/EH&SP       : version 2.21'
  write(*,9)'                        - bug correction in pointer boundaries'
  write(*,9)'                        - bug correction in PIOLIB output   [Planck]'
  write(*,9)'                        - slight edition of documentation (about noise)'
  write(*,9) '05/10/06/EH          : version 2.22'
  write(*,9)'                        - maps cross-correlation'
  write(*,9) '??      /EH          : version 2.30'
  write(*,9)'                        - different beam for each map'
  write(*,9) '02/02/07/EH          : version 2.40'
  write(*,9)'                        - fixed default behavior of second beam option'
  write(*,9)'                        - reorganized help message'
  write(*,9) '2007-03-08/EH        : version 2.45'
  write(*,9)'                        - first stab at PIE-ization        [Planck]'
  write(*,9)'                        - allow different mask for I and (Q,U)'
  write(*,9) '2007-04-04/EH        : version 2.46'
  write(*,9)'                        - simpler ASCII output format for g95'
  write(*,9)'                        - accepts polarized cut sky maps'
  write(*,9)'                        - crashing bug correction in compute windows'
  write(*,9) '2007-04-18/EH        : version 2.47'
  write(*,9)'                        - some more bug correction (nmask=1 by default)'
  write(*,9) '2007-04-27/EH        : version 2.48'
  write(*,9) '                       - PIE-ization done                 [Planck]'
  write(*,9)'                        +++ Delouis test with thinC        [Planck]'
  write(*,9) '2007-04-27/EH        : version 2.49'
  write(*,9) '                       - bug correction in Legendre for x=-1'
!   write(*,9) '2007-07-06/EH        : version ?? '
!   write(*,9) '                       - output in CL group               [Planck]'
!   write(*,9) '                       - added fits_out for FITS C(l) output      '
!   write(*,9) '2007-09-10/EH        : version 2.49_b '
!   write(*,9) '                       - bug correction in argument parsing'
!   write(*,9) '                       - output in VECT group             [Planck]'
!   write(*,9) '                           added -version option'
  write(*,9) '2007-09-10/EH        : version 2.50'
  write(*,9) '                       - set bad pixels to 0                [Planck]'
  write(*,9) '                       - added fits_out for FITS C(l) output      '
  write(*,9) '                       - added version option                     '
  write(*,9) '2007-12-07/EH        : version 2.52'
  write(*,9) '                       - added transfer function                  '
  write(*,9) '                       - Healpix version 2.* or more is required  '
  write(*,9) '2007-12-27/EH        : version 2.53'
  write(*,9) '                       - correctly account for DMC beam file [Planck] '
  write(*,9) '2008-02-29/EH        : version 2.54'
  write(*,9) '                       - copy run parameters in FITS header'
  write(*,9) '2008-04-23/EH        : version 2.55'
  write(*,9) '                       - allow for different weight files and powers '
  write(*,9) '                        for temp. and pol. maps'
  write(*,9) '2008-05-29/EH        : version 2.560'
  write(*,9) '                       - compute TB and EB spectra'
  write(*,9) '2008-06-02/EH        : version 2.57'
  write(*,9) '                       - correct TB and EB from window functions, update documentation'
  write(*,9) '2008-07-01/EH        : version 2.57a'
  write(*,9) '                       - corrected QU correlation: QU_new = -2 * QU_old'
  write(*,9) '2008-08-25/EH        : version 2.57b'
  write(*,9) '                       - add subdipole option (-> corrects bug on monopole subtraction introduced in v2.55)'
  write(*,9) '2008-09-02/EH        : version 2.59'
  write(*,9) '                       - use DMCPIOSTRINGMAXLEN instead of hardcoded value [Planck]'
  write(*,9) '                       - automatic generation of IO routines with PIE [Planck]'
  write(*,9) '2008-09-11/EH        : version 2.60'
  write(*,9) '                       - compiles on CC-in2p3 [Planck]'
  write(*,9) '2008-10-01/EH        : version 2.61'
  write(*,9) '                       - correction of crashing bug when dealing with cl_inmask_file AND polarized mask'
  write(*,9) '2008-12-05/EH        : version 2.62'
  write(*,9) '                       - edited 3-J symbol routine (rec3jj.f90) to fix problem affecting large |m| values'
  write(*,9) '                         (not relevant for CMB polarization)'
  write(*,9) '2009-01-22/EH        : version 2.63'
  write(*,9) '                       - output C(l) in CL tuples [Planck]'
  write(*,9) '2009-01-23/EH        : version 2.63a'
  write(*,9) '                       - can turn code into library [Planck]'
  write(*,9) '                       - renormalize TE (and TB) in all cases, not only decoupling'
  write(*,9) '2009-01-29/EH        : version 2.64'
  write(*,9) '                       - compute and output kernels'
  write(*,9) '2009-02-02/EH        : version 2.65'
  write(*,9) '                       - deallocate global variables before exiting'
  write(*,9) '                       - output kernels in TAB3D objects [Planck]'
  write(*,9) '2009-02-17/EH        : version 2.66'
  write(*,9) '                       - some cleanup'
  write(*,9) '2009-03-23/EH        : version 2.66a'
  write(*,9) '                       - corrected typos in history'
  write(*,9) '2009-05-29/EH        : version 2.68'
  write(*,9) '                       - compiles on magique3 [Planck]'
  write(*,9) '2009-06-10/EH        : version 2.69'
  write(*,9) '                       - npairsthreshold removed (was actually not implemented)'
  write(*,9) '2010-03-12/EH        : version 2.70'
  write(*,9) '                       - started importing cov_temp (TT covariance matrix)'
  write(*,9) '2010-06-04/EH        : version 2.73'
  write(*,9) '                       - debugged TT covariance matrix'
  write(*,9) '2010-06-18/EH        : version 2.75'
  write(*,9) '                       - added symmetric_cl option'
  write(*,9) '2010-08-25/EH        : version 2.76'
  write(*,9) '                       - issues a warning if input maps contain NaN-valued pixels'
  write(*,9) '2011-01-31/EH        : version 2.80'
  write(*,9) '                       - public release'
  write(*,9) '2011-02-15/EH        : version 2.81'
  write(*,9) '                       - fixed bug with string keywords in FITS writing'
  write(*,9) '                       - split a very long line.'
  write(*,9) '2012-03-22/EH        : version 2.83'
  write(*,9) '                       - cosmetic changes [Planck]'
  write(*,9) '2012-05-25/EH        : version 2.84'
  write(*,9) '                       - minor bug and typo corrections'
  write(*,9) '                       - public release'
  write(*,9) '2012-07-13/EH        : version 2.85'
  write(*,9) '                       - parallel calculation of the C(l)xC(l) TT covariance matrix'
  write(*,9) '2012-11-26/EH        : version 2.86'
  write(*,9) '                       - correction of a bug preventing covariance calculation for single apodized map'
  write(*,9) '2013-03-29/EH        : version 2.90'
  write(*,9) '                       - makes use of CFTISIO Extended File Name features available'
  write(*,9) '                          in HEALPIX 3.10, allowing the reading of '
  write(*,9) '                          arbitrary FITS extensions in beam_file*, mapfile* and maskfile*'
  write(*,9) '2013-10-07/EH        : version 2.91'
  write(*,9) '                       - better handling of discrepant masks for cross-maps TE within option symmetric_cl'
  write(*,9) '2014-04-03/EH        : version 3.00'
  write(*,9) '                       - added listmapfiles*_* and listmapweights*_* options to allow linear combination '
  write(*,9) '                           of input maps'
  write(*,9) '                       - computes correctly cross-power spectrum of 2 maps when '
  write(*,9) '                           the first one is masked or weighted (maskfile and/or weightfile are set),'
  write(*,9) '                           and the second one is not (maskfile2 and weightfile2 not set)'
  write(*,9) '                       - turned-off monopole and dipole removal on (Q,U) maps in mask (and weight) free cases'
  write(*,9) '                       - correct default value of extramap* in example parameter file (--help option)'
  write(*,9) '2014-04-08/EH        : version 3.01'
  write(*,9) '                       - correction of a bug introduced in 3.00 (for autospectra)'
  write(*,9) '2015-03-02/EH        : version 3.03'
  write(*,9) '                       - fix un-initialized map variables (creating problems at small Nsides)'
  write(*,9) '                       - all lines shorter than 132 characters (gfortran default limitation)'
  write(*,9) '2015-08-26/EH        : version 3.10'
  write(*,9) '                       - when 2 polarized maps are cross-correlated, 9 or 6 spectra are computed, '
  write(*,9) '                            depending on SYMMETRIC_CL setting'
  write(*,9) '                       - introduced TOLERANCE keyword for E/B decoupling'
  write(*,9) '                       - minor bugs correction in keyword parsing'
  write(*,9) '                       - new python wrapper (ispice.py)'
  write(*,9) '                       - bug correction in IDL wrapper (ispice.pro)'
  write(*,9) '2015-10-01/EH        : version 3.13'
  write(*,9) '                       - add SUBDIPOL in output FITS header'
  write(*,9) '                       - replace rewind (ignored by some compilers) by close+open'
  write(*,9) '2015-10-08/EH        : version 3.14'
  write(*,9) '                       - fixed DMC specific bug in C(l) writing [Planck]'
  write(*,9) '2015-10-12/EH        : version 3.16   (v03-01-06)'
  write(*,9) '                       - fixed DMC specific bug in map name parsing [Planck]'
  write(*,9) '                       - fixed bug in listmapweights2'
  write(*,9) '                       - improved  ispice.pro  and ispice.py'
  write(*,9) '2016-08-09/EH        : version 3.20   (v03-02-00)'
  write(*,9) '                       - fixed error in reading cut sky polarized FITS files'
  write(*,9) '                       - improved  read_spice.pro  and ispice.py'
  write(*,9) '2016-12-01/EH        : version 3.30'
  write(*,9) '                       - introduced bin_llcl.py'
  write(*,9) '                       - first attempt at supporting correlation of maps '
  write(*,9) '                             with discrepant pixel sizes (Nside)'
  write(*,9) '2017-02-21/EH        : version 3.31   (v03-03-01)'
  write(*,9) '                       - slight edits to documentation'
  write(*,9) '2017-04-05/EH        : version 3.32   (v03-03-02)'
  write(*,9) '                       - kernel calculation, relating the output C(l)s to the true ones,'
  write(*,9) '                         is now possible even when polarization is not set.'
  write(*,9) '2017-08-04/EH        : version 3.33   (v03-03-03)'
  write(*,9) '                       - more sensible trigger of kernel calculation'
  write(*,9) '                       - faster map2alm step by ignoring empty Healpix rings'
  write(*,9) '2017-08-10/EH        : version 3.34   (v03-03-04)'
  write(*,9) '                       - added tenormfilein and tenormfileout (smaller file size that windowfile*)'
  write(*,9) '2017-08-28/EH        : version 3.40   (v03-04-00)'
  write(*,9) '                       - faster calculation of kernels (and tenormfileout)'
  write(*,9) '2018-01-16/EH        : version 3.51   (v03-05-01)'
  write(*,9) '                       - support for double precision maps and alms'
  write(*,9) '                       - support for cmake '
  write(*,9) '2019-06-06/EH        : version 3.52   (v03-05-02)'
  write(*,9) '                       - ispice.py runs in python 2 and 3 '
  write(*,9) '2019-12-20/EH        : version 3.61   (v03-06-01)'
  write(*,9) '                       - smaller memory footprint, in particular for cut-sky format input maps'
  write(*,9) '                         with fsky < 25%'
  write(*,9) '                       - requires HEALPix 3.60 or more'
  write(*,9) '2020-01-30/EH        : version 3.62   (v03-06-02)'
  write(*,9) '                       - fixed bug in reading of mask(s) and weight(s) for polarisation;'
  write(*,9) '                       - can read partial sky IQU FITS files as generated by healpy'
  write(*,9) '                       - improvements to ispice.py'
  write(*,9) '2020-02-04/EH        : version 3.64   (v03-06-04)'
  write(*,9) '                       - compliant with gfortran'
  write(*,9) '                       - debugging in mask and weight reading'
  write(*,9) '2020-04-08/EH        : version 3.65   (v03-06-05)'
  write(*,9) '                       - can mix full sky masks and partial sky maps'
  write(*,9) '2020-04-10/EH        : version 3.66   (v03-06-06)'
  write(*,9) '                       - debugged parsing of C(l) text files'
  write(*,9) '2020-05-22/EH        : version 3.67   (v03-06-07)'
  write(*,9) '                       - fixed segfault in ClCl covariance matrix,'
  write(*,9) '                       - added omp_num_threads, kernelsfileout, noiseclfile, noisecorfile'
  write(*,9) '                         and normfac to ispice.py'
  write(*,9) '2020-05-26/EH        : version 3.68   (v03-06-08)'
  write(*,9) '                       - make sure partial maps and masks are correctly sorted on pixel index'
  write(*,9) '2020-09-22/EH        : version 3.69   (v03-06-09)'
  write(*,9) '                       - allow dumping of alm of (masked and weighted) maps'
  write(*,9) '2020-10-16/EH        : version 3.70   (v03-07-00)'
  write(*,9) '                       - fixed parsing bugs'
9 format(a)
  STOP
end subroutine history

!======================================
subroutine process_listmap(inc, i, j, narguments, count, create, none, default, &
     & argument, nefmax, &
     & nlmaps,  extramapTQUfile_input, listmapw8)
!======================================
  use healpix_types
  use misc_utils, only:string
#ifdef PIO
  use piolib
  use spice_common, only: parContent, pie_param
#endif
  use spice_common, only: option, option_param

  implicit none
  integer(i4b),     intent(INOUT):: inc
  integer(i4b),     intent(IN)   :: i, j, narguments, nefmax
  logical(lgt),     intent(IN)   :: count, create
  character(len=*), intent(IN)   :: none, default
  character(len=*), intent(IN),    dimension(1:narguments) :: argument
  integer(i4b),     intent(INOUT), dimension(1:3,1:2)      :: nlmaps
  character(len=*), intent(INOUT), dimension(1:nefmax,1:3,1:2) :: extramapTQUfile_input
  real(DP),         intent(INOUT), dimension(1:nefmax,1:3,1:2) :: listmapw8
  integer(I4B)      :: imap, jmap, kstokes
  character(len=10) :: simap, sjmap, skstokes
  character(len=30) :: myopt


  ! First (and optional second) list(s) of maps
  do jmap=1,2 
     sjmap = adjustl(string(jmap,format='(i1)'))
     do imap=1,nefmax
        simap = adjustl(string(imap,format='(i3)'))
        inc=inc+1
        if ((.not.count) .and. j == inc) then
           if (create) then
              myopt = trim('-listmapfiles'//trim(sjmap)//'_'//trim(simap))
              call create_option(j,2,myopt,none)
           else
              if (i < narguments) then
#ifdef PIO
                 if (jmap==1 .and. imap <= pie_param%n_listmapfiles1_) &
                      & option_param(j)=pie_param%listmapfiles1_(imap)
                 if (jmap==2 .and. imap <= pie_param%n_listmapfiles2_) &
                      & option_param(j)=pie_param%listmapfiles2_(imap)
#else
                 option_param(j)=argument(i+1)
#endif
                 if (trim(option_param(j)) /= none) then
                    extramapTQUfile_input(imap,1,jmap) = option_param(j)
                    nlmaps                    (1,jmap) = nlmaps(1,jmap)+1
                 endif
              else
                 call error_usage(option(j))
              endif
           endif
        endif
     enddo
  enddo

  ! First (and optional second) list(s) of weights
  do jmap=1,2
     sjmap = adjustl(string(jmap,format='(i1)'))
     do imap=1,nefmax
        simap = adjustl(string(imap,format='(i3)'))
        inc=inc+1
        if ((.not.count) .and. j == inc) then
           if (create) then
              myopt = trim('-listmapweights'//trim(sjmap)//'_'//trim(simap))
              call create_option(j,2,myopt,string(listmapw8(imap,1,jmap)))
           else
              if (i < narguments) then
#ifdef PIO
                 if (jmap==1 .and. imap <= pie_param%n_listmapweights1_) &
                      & option_param(j) = string(pie_param%listmapweights1_(imap))
                 ! bug correction 2015-10-12
                 if (jmap==2 .and. imap <= pie_param%n_listmapweights2_) &
                      & option_param(j) = string(pie_param%listmapweights2_(imap))
#else
                 option_param(j) = argument(i+1)
#endif
                 if (trim(option_param(j)) /= none .and. trim(option_param(j)) /= default) then
                    read(option_param(j),*,err=1) listmapw8(imap,1,jmap)
                    listmapw8(imap,2:3,jmap) = listmapw8(imap,1,jmap) ! same weights for Q,U as for T
                 endif
              else
                 call error_usage(option(j))
              endif
           endif
        endif
     enddo
  enddo

  ! First (and optional second) list(s) of Q and U maps (only for DMC/PIOLIB)
#ifdef PIO
  do kstokes=2,3
     if (kstokes==2) skstokes='Q'
     if (kstokes==3) skstokes='U'
     do jmap=1,2
        sjmap = adjustl(string(jmap,format='(i1)'))
        do imap=1,nefmax
           simap = adjustl(string(imap,format='(i3)'))

           inc=inc+1
           if ((.not.count) .and. j == inc) then
              if (create) then
                 myopt = trim('-listmapfiles'//trim(skstokes)//trim(sjmap)//'_'//trim(simap))
                 call create_option(j,2,myopt,none)
              else
                 if (i < narguments) then
                    if (kstokes==2) then
                       if (jmap==1 .and.    imap <= pie_param%n_listmapfilesQ1_) &
                            & option_param(j)=pie_param%listmapfilesQ1_(imap)
                       if (jmap==2 .and.    imap <= pie_param%n_listmapfilesQ2_) &
                            & option_param(j)=pie_param%listmapfilesQ2_(imap)
                    else                                                                                    
                       if (jmap==1 .and.    imap <= pie_param%n_listmapfilesU1_) &
                            & option_param(j)=pie_param%listmapfilesU1_(imap)
                       if (jmap==2 .and.    imap <= pie_param%n_listmapfilesU2_) &
                            & option_param(j)=pie_param%listmapfilesU2_(imap)
                    endif
                    if (trim(option_param(j)) /= none) then
                       extramapTQUfile_input(imap,kstokes,jmap) = option_param(j)
                       nlmaps                    (kstokes,jmap) = nlmaps(kstokes,jmap)+1
                    endif
                 else
                    call error_usage(option(j))
                 endif
              endif
           endif
        enddo ! loop on imap
     enddo ! loop on jmap
  enddo ! loop on kstokes

#endif

  return
1 call error_usage(option(j))
end subroutine process_listmap

!=======================================================================
subroutine make_option_active(i,j,iskip,argument,narguments,create,count)
!=======================================================================
  ! if count==T : count number of options offered by code
  ! if count==F and create==T : create list of options with (useless) default value
  ! if count==F and create==F : read actual value chosen by user with
  !   either argument(i) or a pie variable
  ! j = location in option list
  ! i = index of argument
  use spice_common
  use misc_utils, only:string, fatal_error
#ifdef PIO
!   use pie_io, only: parContent
  use piolib
#endif
  implicit none

  integer(I4B),     intent(IN)                           :: narguments, i, j
  integer(I4B),     intent(OUT)                          :: iskip
  character(len=*), intent(IN), dimension(1:narguments)  :: argument
  logical,          intent(IN)                           :: create,count
  integer(I4B)                                           :: inc, nverbose

#ifdef PIO
  character(len=FILENAMELEN) :: set_optional_dmc_keyword
#endif


  if ((.not.count) .and. (.not.create)) then
#ifdef PIO
     iskip=1
#else
     iskip=iskip_option(j)
#endif
  endif
  
  if (count) noption=0
  inc=0

!=======================================================================
! Special options
!=======================================================================

  inc=inc+1
  if ((.not.count) .and. j == inc) then
     if (create) then
        call create_option(j,1,'-about',default)
     else
        call about
     endif
  endif

  inc=inc+1
  if ((.not.count) .and. j == inc) then
     if (create) then
        call create_option(j,1,'-help',default)
     else
        call man
     endif
  endif

  inc=inc+1
  if ((.not.count) .and. j == inc) then
     if (create) then
        call create_option(j,1,'-history',default)
     else
        call history
     endif
  endif

  inc=inc+1
  if ((.not.count) .and. j == inc) then
     if (create) then
        call create_option(j,2,'-optoutfile',default)
     else
        if (i < narguments) then
           spice_rc_outfile_input=argument(i+1)
           if (spice_rc_outfile_input /= default) then
              if (spice_rc_outfile_input /= none) then
                 default_spice_rc_out=.false.
                 output_spice_rc=.true.
              endif
           else
              output_spice_rc=.true.
           endif
        else
           call error_usage(option(j))
        endif
     endif
  endif

  inc=inc+1
  if ((.not.count) .and. j == inc) then
     if (create) then
        call create_option(j,2,'-optinfile',default)
     else
        if (i < narguments) then
           spice_rc_infile_input=argument(i+1)
           if (spice_rc_infile_input /= default) then
              if (spice_rc_infile_input /= none) then
                 default_spice_rc_in=.false.
                 input_spice_rc=.true.
              endif
           else
              input_spice_rc=.true.
           endif
        else
           call error_usage(option(j))
        endif
     endif
  endif

  inc=inc+1
  if ((.not.count) .and. j == inc) then
     if (create) then
        call create_option(j,1,'-usage',default)
     else
        call usage_reminder
        !call fatal_error
        stop
     endif
  endif

  inc=inc+1
  if ((.not.count) .and. j == inc) then
     if (create) then
        call create_option(j,1,'-version',default)
     else
        call print_version
     endif
  endif

!=======================================================================
! Normal options
!=======================================================================

  inc=inc+1
  jstartnormal=inc ! Special number 

  if ((.not.count) .and. j == inc) then
     if (create) then
        call create_option(j,2,'-alm1fileout',none)
     else
        if (i < narguments) then
#ifdef PIO
           if (pie_param%flag_alm1fileout) option_param(j)     = pie_param%alm1fileout
#else
           option_param(j)=trim(argument(i+1))
#endif
           if (trim(option_param(j)) == none) then
              alm1_out_file = ''
              alm1_dump     = .false.
           elseif (trim(option_param(j)) == default) then
              alm1_out_file = trim(alm1file_default)
              alm1_dump     = .true.
           else
              alm1_out_file = trim(option_param(j))
              alm1_dump     = .true.
           endif
        else
           call error_usage(option(j))
        endif
     endif
  endif

  inc=inc+1
  if ((.not.count) .and. j == inc) then
     if (create) then
        call create_option(j,2,'-alm2fileout',none)
     else
        if (i < narguments) then
#ifdef PIO
           if (pie_param%flag_alm2fileout) option_param(j)     = pie_param%alm2fileout
#else
           option_param(j)=trim(argument(i+1))
#endif
           if (trim(option_param(j)) == none) then
              alm2_out_file = ''
              alm2_dump     = .false.
           elseif (trim(option_param(j)) == default) then
              alm2_out_file = trim(alm2file_default)
              alm2_dump     = .true.
           else
              alm2_out_file = trim(option_param(j))
              alm2_dump     = .true.
           endif
        else
           call error_usage(option(j))
        endif
     endif
  endif

  inc=inc+1
  if ((.not.count) .and. j == inc) then
     if (create) then
        call create_option(j,2,'-apodizesigma',none)
     else
        if (i < narguments) then
#ifdef PIO
           option_param(j) = pie_param%apodizesigma
#else
           option_param(j)=argument(i+1)
#endif
           if (option_param(j) /= none) then
              read(option_param(j),*,err=1) apodizesigma
              apodize=.true.
              if (apodizesigma <= 0.d0) then
                 write(*,*) 'Improper apodizesigma value.'
                 CALL FATAL_ERROR
              endif
           endif
        else
           call error_usage(option(j))
        endif
     endif
  endif

  inc=inc+1
  if ((.not.count) .and. j == inc) then
     if (create) then
        call create_option(j,2,'-apodizetype',none)
     else
        if (i < narguments) then
#ifdef PIO
           option_param(j) = string(int(pie_param%apodizetype))
#else
           option_param(j)=argument(i+1)
#endif
           if (option_param(j) /= none) then
              read(option_param(j),*,err=1) apodizetype
              if (apodizetype /= 0 .and. apodizetype /= 1 ) then
                 write(*,*) 'Improper apodizetype value.'
                 CALL FATAL_ERROR
              endif
           endif
        else
           call error_usage(option(j))
        endif
     endif
  endif

  inc=inc+1
  if ((.not.count) .and. j == inc) then
     if (create) then
        call create_option(j,2,'-beam',none)
     else
        if (i < narguments) then
#ifdef PIO
           option_param(j) = pie_param%beam
#else
           option_param(j)=argument(i+1)
#endif
           if (option_param(j) /= none) then
              read(option_param(j),*,err=1) fwhm
              correct_beam=.true.
           endif
        else
           call error_usage(option(j))
        endif
     endif
  endif

  inc=inc+1
  if ((.not.count) .and. j == inc) then
     if (create) then
        call create_option(j,2,'-beam_file',none) 
     else
        if (i < narguments) then
#ifdef PIO
           if (pie_param%flag_beam_file) option_param(j)     = pie_param%beam_file
           beam_file=trim(option_param(j))
#else
           option_param(j)=argument(i+1)
           beam_file=trim(option_param(j))
#endif
           if (trim(beam_file) /= none) then
              correct_beam=.true.
              beam_present=.true.
           endif
        else
           call error_usage(option(j))
        endif
     endif
  endif

  inc=inc+1
  if ((.not.count) .and. j == inc) then
     if (create) then
        call create_option(j,2,'-beam2',none)
     else
        if (i < narguments) then
#ifdef PIO
           option_param(j) = pie_param%beam2
#else
           option_param(j)=argument(i+1)
#endif
           if (option_param(j) /= none) then
              read(option_param(j),*,err=1) fwhm2
              correct_beam2=.true.
           endif
        else
           call error_usage(option(j))
        endif
     endif
  endif

  inc=inc+1
  if ((.not.count) .and. j == inc) then
     if (create) then
        call create_option(j,2,'-beam_file2',none) 
     else
        if (i < narguments) then
#ifdef PIO
           if (pie_param%flag_beam_file2) option_param(j)     = pie_param%beam_file2
           beam_file2=trim(option_param(j))
#else
           option_param(j)=argument(i+1)
           beam_file2=trim(option_param(j))
#endif
           if (trim(beam_file2) /= none) then
              correct_beam2=.true.
              beam2_present=.true.
           endif
        else
           call error_usage(option(j))
        endif
     endif
  endif

  inc=inc+1
  if ((.not.count) .and. j == inc) then
     if (create) then
        call create_option(j,2,'-clfile',none) 
     else
        if (i < narguments) then
#ifdef PIO
            if (pie_param%flag_clfile) option_param(j)     = pie_param%clfile
#else
           option_param(j)=argument(i+1)
#endif
           cl_file_name_input=trim(option_param(j))
           if (cl_file_name_input /= none) then
              cloutput=.true.
              if (cl_file_name_input == default) then
                 option_param(j)=cl_file_name_default
              else
                 default_cl_file=.false.
              endif
           endif
        else
           call error_usage(option(j))
        endif
     endif
  endif

  inc=inc+1
  if ((.not.count) .and. j == inc) then
     if (create) then
        call create_option(j,2,'-cl_outmap_file',none)
     else
        if (i < narguments) then
#ifdef PIO
           if (pie_param%flag_cl_outmap_file) option_param(j)     = pie_param%cl_outmap_file
#else
           option_param(j)=argument(i+1)
#endif
           cl_outmap_file_input=trim(option_param(j))
           if (cl_outmap_file_input == default) then
              clmapoutput=.true.
           elseif (cl_outmap_file_input /= none) then
              default_clmapout_file=.false.
              clmapoutput=.true.
           endif
        else
           call error_usage(option(j))
        endif
     endif
  endif

  inc=inc+1
  if ((.not.count) .and. j == inc) then
     if (create) then
        call create_option(j,2,'-cl_inmap_file',none)
     else
        if (i < narguments) then
#ifdef PIO
           if (pie_param%flag_cl_inmap_file) option_param(j)     = pie_param%cl_inmap_file
#else
           option_param(j)=argument(i+1)
#endif
           cl_inmap_file_input=trim(option_param(j))
           if (cl_inmap_file_input == default) then
              clmapinput=.true.
           elseif (cl_inmap_file_input /= none) then
              default_clmapin_file=.false.
              clmapinput=.true.
           endif
        else
           call error_usage(option(j))
        endif
     endif
  endif

  inc=inc+1
  if ((.not.count) .and. j == inc) then
     if (create) then
        call create_option(j,2,'-cl_outmask_file',none)
     else
        if (i < narguments) then
#ifdef PIO
           if (pie_param%flag_cl_outmask_file) option_param(j)     = pie_param%cl_outmask_file
#else
           option_param(j)=argument(i+1)
#endif
           cl_outmask_file_input=trim(option_param(j))
           if (cl_outmask_file_input == default) then
              clmaskoutput=.true.
           elseif (cl_outmask_file_input /= none) then
              default_clmaskout_file=.false.
              clmaskoutput=.true.
           endif
        else
           call error_usage(option(j))
        endif
     endif
  endif

  inc=inc+1
  if ((.not.count) .and. j == inc) then
     if (create) then
        call create_option(j,2,'-cl_inmask_file',none)
     else
        if (i < narguments) then
#ifdef PIO
           if (pie_param%flag_cl_inmask_file) option_param(j)     = pie_param%cl_inmask_file
#else
           option_param(j)=argument(i+1)
#endif
           cl_inmask_file_input=trim(option_param(j))
           if (cl_inmask_file_input == default) then
              clmaskinput=.true.
           elseif (cl_inmask_file_input /= none) then
              default_clmaskin_file=.false.
              clmaskinput=.true.
           endif
        else
           call error_usage(option(j))
        endif
     endif
  endif

  inc=inc+1
  if ((.not.count) .and. j == inc) then
     if (create) then
        call create_option(j,2,'-corfile',none)
     else
        if (i < narguments) then
#ifdef PIO
           if (pie_param%flag_corfile) option_param(j)     = pie_param%corfile
#else
           option_param(j)=argument(i+1)
#endif
           cor_file_name_input=trim(option_param(j))
           if (cor_file_name_input /= none) then
              coroutput=.true.
              if (cor_file_name_input == default) then                 
                 option_param(j)=cor_file_name_default
              else
                 default_cor_file=.false.
              endif
           endif
        else
           call error_usage(option(j))
        endif
     endif
  endif

  inc=inc+1
  if ((.not.count) .and. j == inc) then
     if (create) then
        call create_option(j,2,'-covfileout',none)
     else
        if (i < narguments) then
#ifdef PIO
           if (pie_param%flag_covfileout) option_param(j)     = pie_param%covfileout
#else
           option_param(j)=argument(i+1)
#endif
           cov_out_file_input=trim(option_param(j))
           if (cov_out_file_input == default) then
              do_cov=.true.
           elseif (cov_out_file_input /= none) then
              default_covout_file=.false.
              do_cov=.true.
           endif
        else
           call error_usage(option(j))
        endif
     endif
  endif

  inc=inc+1
  if ((.not.count) .and. j == inc) then
     if (create) then
        call create_option(j,2,'-decouple',none)
     else
        if (i < narguments) then
#ifdef PIO
           option_param(j) = pie_param%decouple
#else
           option_param(j)=argument(i+1)
#endif
           if (option_param(j) == default) then
              decouple=.true.
           elseif (option_param(j) /= none) then
              write(*,*) 'Improper decouple option.'
              CALL FATAL_ERROR
           endif
        else
           call error_usage(option(j))
        endif
     endif
  endif

#ifndef PIO
  inc=inc+1
  if ((.not.count) .and. j == inc) then
     if (create) then
        call create_option(j,2,'-dry',none)
     else
        if (i < narguments) then
           option_param(j)=argument(i+1)
           if (option_param(j) == default) then
              dry=.true.
           elseif (option_param(j) /= none) then
              write(*,*) 'Improper dry option.'
              CALL FATAL_ERROR
           endif
        else
           call error_usage(option(j))
        endif
     endif
  endif
#endif

  inc=inc+1
  if ((.not.count) .and. j == inc) then
     if (create) then
        call create_option(j,2,'-extramapfile',none)
     else
        if (i < narguments) then
#ifdef PIO
           if (pie_param%flag_extramapfile) option_param(j) = pie_param%extramapfile
           extramap1_present            = pie_param%flag_extramapfile
#else
           option_param(j)              = argument(i+1)
           extramap1_present            = (trim(option_param(j)) /= default &
                &                    .and. trim(option_param(j)) /= none)
#endif
           extramapTQUfile_input(2,1,1) = trim(option_param(j))
        else
           call error_usage(option(j))
        endif
     endif
  endif

  inc=inc+1
  if ((.not.count) .and. j == inc) then
     if (create) then
        call create_option(j,2,'-extramapfile2',none)
     else
        if (i < narguments) then
#ifdef PIO
           if (pie_param%flag_extramapfile2) option_param(j) = pie_param%extramapfile2
           extramap2_present            = pie_param%flag_extramapfile2
#else
           option_param(j)              = argument(i+1)
           extramap2_present            = (trim(option_param(j)) /= default &
                &                    .and. trim(option_param(j)) /= none)
#endif
           extramapTQUfile_input(2,1,2) = trim(option_param(j))
        else
           call error_usage(option(j))
        endif
     endif
  endif

#ifdef PIO
  inc=inc+1
  if ((.not.count) .and. j == inc) then
     if (create) then
        call create_option(j,2,'-extramapQfile',none)
     else
        if (i < narguments) then
           if (pie_param%flag_extramapQfile) option_param(j) = pie_param%extramapQfile
           readextraQfile1               = pie_param%flag_extramapQfile
           extramapTQUfile_input(2,2,1)  = trim(option_param(j))
        else
           call error_usage(option(j))
        endif
     endif
  endif

  inc=inc+1
  if ((.not.count) .and. j == inc) then
     if (create) then
        call create_option(j,2,'-extramapUfile',none)
     else
        if (i < narguments) then
           if (pie_param%flag_extramapUfile) option_param(j) = pie_param%extramapUfile
           readextraUfile1               = pie_param%flag_extramapUfile
           extramapTQUfile_input(2,3,1)  = trim(option_param(j))
        else
           call error_usage(option(j))
        endif
     endif
  endif

  inc=inc+1
  if ((.not.count) .and. j == inc) then
     if (create) then
        call create_option(j,2,'-extramapQfile2',none)
     else
        if (i < narguments) then
           if (pie_param%flag_extramapQfile2) option_param(j)  = pie_param%extramapQfile2
           readextraQfile2               = pie_param%flag_extramapQfile2
           extramapTQUfile_input(2,2,2)  = trim(option_param(j))
        else
           call error_usage(option(j))
        endif
     endif
  endif

  inc=inc+1
  if ((.not.count) .and. j == inc) then
     if (create) then
        call create_option(j,2,'-extramapUfile2',none)
     else
        if (i < narguments) then
           if (pie_param%flag_extramapUfile2) option_param(j)  = pie_param%extramapUfile2
           readextraUfile2               = pie_param%flag_extramapUfile2
           extramapTQUfile_input(2,3,2)  = trim(option_param(j))
        else
           call error_usage(option(j))
        endif
     endif
  endif
#endif


#ifndef PIO
  inc=inc+1
  if ((.not.count) .and. j == inc) then
     if (create) then
        call create_option(j,2,'-fits_out',none)
     else
        if (i < narguments) then
           option_param(j)=argument(i+1)
           if (option_param(j) == default) then
              fits_out=.true.
           elseif (option_param(j) /= none) then
              write(*,*) 'Improper fits_out option.'
              CALL FATAL_ERROR
           endif
        else
           call error_usage(option(j))
        endif
     endif
  endif
#endif

  inc=inc+1
  if ((.not.count) .and. j == inc) then
     if (create) then
        call create_option(j,2,'-kernelsfileout',none)
     else
        if (i < narguments) then
#ifdef PIO
           if (pie_param%flag_kernelsfileout) option_param(j)     = pie_param%kernelsfileout
#else
           option_param(j)=argument(i+1)
#endif
           kernels_out_file_input=trim(option_param(j))
           if (kernels_out_file_input == default) then
              kernelsoutput=.true.
           elseif (kernels_out_file_input /= none) then
              default_kernelsout_file=.false.
              kernelsoutput=.true.
           endif
        else
           call error_usage(option(j))
        endif
     endif
  endif

!-------------------------- list map
  call process_listmap(inc, i, j, narguments, count, create, none, default, &
       &  argument, nefmax, &
       &  nlmaps,  extramapTQUfile_input, listmapw8)


  inc=inc+1
  if ((.not.count) .and. j == inc) then
     if (create) then
        call create_option(j,2,'-mapfile',default)
        option_param(j)=trim(mapfile_default)//'.fits'
     else
        if (i < narguments) then
#ifdef PIO
           if (pie_param%flag_mapfile) option_param(j)=pie_param%mapfile
#else
           option_param(j)=argument(i+1)
#endif
           mapfile_input=trim(option_param(j))
           if (mapfile_input /= default) then
              default_map_file=.false.
           else
              option_param(j)=trim(mapfile_default)//'.fits'
           endif
        else
           call error_usage(option(j))
        endif
     endif
  endif

  inc=inc+1
  if ((.not.count) .and. j == inc) then
     if (create) then
        call create_option(j,2,'-mapfile2',none)
     else
        if (i < narguments) then
#ifdef PIO
           if (pie_param%flag_mapfile2) option_param(j)=pie_param%mapfile2
#else
           option_param(j)=argument(i+1)
#endif
           mapfile2_input=trim(option_param(j))
           if (mapfile2_input /= default) then
              if (mapfile2_input /= none) then
                 default_map_file2=.false.
                 map2_present=.true.
              endif
           else
              map2_present=.true.
              option_param(j)=trim(mapfile2_default)//'.fits'
           endif
        else
           call error_usage(option(j))
        endif
     endif
  endif

#ifdef PIO
  inc=inc+1
  if ((.not.count) .and. j == inc) then
     if (create) then
        call create_option(j,2,'-mapQfile',none)
     else
        if (i < narguments) then
           if (pie_param%flag_mapQfile) option_param(j) = pie_param%mapQfile
           mapQfile_input  = trim(option_param(j))
           readQfile       = pie_param%flag_mapQfile
        else
           call error_usage(option(j))
        endif
     endif
  endif

  inc=inc+1
  if ((.not.count) .and. j == inc) then
     if (create) then
        call create_option(j,2,'-mapUfile',none)
     else
        if (i < narguments) then
           if (pie_param%flag_mapUfile) option_param(j) = pie_param%mapUfile
           mapUfile_input  = trim(option_param(j))
           readUfile       = pie_param%flag_mapUfile
        else
           call error_usage(option(j))
        endif
     endif
  endif
#endif

#ifdef PIO
  inc=inc+1
  if ((.not.count) .and. j == inc) then
     if (create) then
        call create_option(j,2,'-mapQfile2',none)
     else
        if (i < narguments) then
           if (pie_param%flag_mapQfile2) option_param(j)  = pie_param%mapQfile2
           mapQfile2_input  = trim(option_param(j))
           readQfile2       = pie_param%flag_mapQfile2
        else
           call error_usage(option(j))
        endif
     endif
  endif

  inc=inc+1
  if ((.not.count) .and. j == inc) then
     if (create) then
        call create_option(j,2,'-mapUfile2',none)
     else
        if (i < narguments) then
           if (pie_param%flag_mapUfile2) option_param(j)  = pie_param%mapUfile2
           mapUfile2_input  = trim(option_param(j))
           readUfile2       = pie_param%flag_mapUfile2
        else
           call error_usage(option(j))
        endif
     endif
  endif
#endif


  inc=inc+1
  if ((.not.count) .and. j == inc) then
     if (create) then
        call create_option(j,2,'-maskfile',none)
     else
        if (i < narguments) then
#ifdef PIO
           if (pie_param%flag_maskfile) option_param(j)   = pie_param%maskfile
           maskfile_input    = trim(option_param(j))
           masks_present     = pie_param%flag_maskfile
           default_mask_file = .false.
#else
           option_param(j)=argument(i+1)
           maskfile_input=trim(option_param(j))
           if (maskfile_input /= default) then
              if (maskfile_input /= none) then
                 default_mask_file=.false.
                 masks_present=.true.
              endif
           else
              masks_present=.true.
              option_param(j)=trim(maskfile_default)//'.fits'
              maskfile_input = trim(option_param(j)) ! added by EH
           endif
#endif
        else
           call error_usage(option(j))
        endif
     endif
  endif

  !---
  inc=inc+1
  if ((.not.count) .and. j == inc) then
     if (create) then
        call create_option(j,2,'-maskfile2',none)
     else
        if (i < narguments) then
#ifdef PIO
           if (pie_param%flag_maskfile2) option_param(j)   = pie_param%maskfile2
#else
           option_param(j)=argument(i+1)
#endif
           maskfile2_input=trim(option_param(j))
           if (maskfile2_input /= default) then
              if (maskfile2_input /= none) then
                 default_mask_file2=.false.
                 masks2_present=.true.
              endif
           else
              masks2_present=.true.
              option_param(j)=trim(maskfile2_default)//'.fits'
           endif
        else
           call error_usage(option(j))
        endif
     endif
  endif

  inc=inc+1
  if ((.not.count) .and. j == inc) then
     if (create) then
        call create_option(j,2,'-maskfilep',default)
     else
        if (i < narguments) then
#ifdef PIO
           if (pie_param%flag_maskfilep) option_param(j)   = pie_param%maskfilep
           maskfilep_input   = trim(option_param(j))
           maskfilep_toread  = pie_param%flag_maskfilep
#else
           option_param(j)   = argument(i+1)
           maskfilep_input   = trim(option_param(j))
           maskfilep_toread  = (maskfilep_input /= default)
#endif
!            if (trim(option_param(j)) == none) then
!               print*,trim(none)//' is not a valid value for '//trim(option(j))
!               call fatal_error
!            endif
        else
           call error_usage(option(j))
        endif
     endif
  endif

  inc=inc+1
  if ((.not.count) .and. j == inc) then
     if (create) then
        call create_option(j,2,'-maskfilep2',default)
     else
        if (i < narguments) then
#ifdef PIO
           if (pie_param%flag_maskfilep2) option_param(j)    = pie_param%maskfilep2
           maskfilep2_input   = trim(option_param(j))
           maskfilep2_toread  = pie_param%flag_maskfilep2
#else
           option_param(j)    = argument(i+1)
           maskfilep2_input   = trim(option_param(j))
           maskfilep2_toread  = (maskfilep2_input /= default)
#endif
!            if (trim(option_param(j)) == none) then
!               print*,trim(none)//' is not a valid value for '//trim(option(j))
!               call fatal_error
!            endif
        else
           call error_usage(option(j))
        endif
     endif
  endif

  inc=inc+1
  if ((.not.count) .and. j == inc) then
     if (create) then
        call create_option(j,2,'-nlmax',none)
     else
        if (i < narguments) then
#ifdef PIO
           option_param(j)=string(int(pie_param%nlmax))
#else
           option_param(j)=argument(i+1)
#endif
           if (option_param(j) /= none) then
              read(option_param(j),*,err=1) nlmax
           endif
        else
           call error_usage(option(j))
        endif
     endif
  endif

  inc=inc+1
  if ((.not.count) .and. j == inc) then
     if (create) then
        call create_option(j,2,'-normfac',none)
     else
        if (i < narguments) then
#ifdef PIO
           option_param(j) = string(pie_param%normfac)
#else
           option_param(j)=argument(i+1)
#endif
           if (option_param(j) /= none) then
              read(option_param(j),*,err=1) quad_uK
              if (quad_uK /= 1.0_DP) normalize=.true.
           endif
        else
           call error_usage(option(j))
        endif
     endif
  endif

  inc=inc+1
  if ((.not.count) .and. j == inc) then
     if (create) then
        call create_option(j,2,'-npairsthreshold',none)
     else
        if (i < narguments) then
#ifdef PIO
           option_param(j) = string(pie_param%npairsthreshold)
#else
           option_param(j)=argument(i+1)
#endif
           if (option_param(j) /= none) then
!               write(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
!               write(*,*) 'WARNING: npairsthreshold is deprecated'
!               write(*,*) 'keyword value will be ignored.'
!               write(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
              read(option_param(j),*,err=1) npairsthreshold
              if (npairsthreshold /= 0.d0) then
                 pairsthresholding=.true.
!               elseif (npairsthreshold < 0.d0) then
!                  write(*,*) 'Improper npairsthreshold value.'
!                  CALL FATAL_ERROR
              endif
           endif
        else
           call error_usage(option(j))
        endif
     endif
  endif

  inc=inc+1
  if ((.not.count) .and. j == inc) then
     if (create) then
        call create_option(j,2,'-noisecorfile',none)
     else
        if (i < narguments) then
#ifdef PIO
           if (pie_param%flag_noisecorfile) option_param(j)     = pie_param%noisecorfile
#else
           option_param(j)=argument(i+1)
#endif
           noisecor_file_name_input=trim(option_param(j))
           if (noisecor_file_name_input /= default) then
              if (noisecor_file_name_input /= none) then
                 subtract_noise_cor=.true.
                 default_noisecor_file=.false.
              endif
           else
              subtract_noise_cor=.true.
              option_param(j)=noisecor_file_name_default
           endif
        else
           call error_usage(option(j))
        endif
     endif
  endif

  inc=inc+1
  if ((.not.count) .and. j == inc) then
     if (create) then
        call create_option(j,2,'-noiseclfile',none)
     else
        if (i < narguments) then
#ifdef PIO
           if (pie_param%flag_noiseclfile) option_param(j)     = pie_param%noiseclfile
#else
           option_param(j)=argument(i+1)
#endif
           noisecl_file_name_input=trim(option_param(j))
           if (noisecl_file_name_input /= default) then
              if (noisecl_file_name_input /= none) then
                 subtract_noise_cl=.true.
                 default_noisecl_file=.false.
              endif
           else
              subtract_noise_cl=.true.
              option_param(j)=noisecl_file_name_default
           endif
        else
           call error_usage(option(j))
        endif
     endif
  endif

#ifdef PIO
  inc=inc+1
  if ((.not.count) .and. j == inc) then
     if (create) then
        call create_option(j,2,'-out_cl_grp',none)
     else
        if (i < narguments) then
           option_param(j) = pie_param%out_cl_grp
           out_cl_grp      = trim(option_param(j))
           know_out_cl_grp = (out_cl_grp /= none)
        else
           call error_usage(option(j))
        endif
     endif
  endif
#endif

  inc=inc+1
  if ((.not.count) .and. j == inc) then
     if (create) then
        call create_option(j,2,'-overwrite',default) 
     else
        if (i < narguments) then
#ifdef PIO
           option_param(j) = pie_param%overwrite
#else
           option_param(j)=argument(i+1)
#endif
           if (option_param(j) == none) then
              overwrite=.false.
           elseif (option_param(j) /= default) then
              write(*,*) 'Improper overwrite option'
              CALL FATAL_ERROR
           endif
        else
           call error_usage(option(j))
        endif
     endif
  endif

  inc=inc+1
  if ((.not.count) .and. j == inc) then
     if (create) then
        call create_option(j,2,'-pixelfile',none) 
     else
        if (i < narguments) then
#ifdef PIO
           option_param(j)  = pie_param%pixelfile
#else
           option_param(j)  = argument(i+1)
#endif
           !print*,'pixelfile =  ', trim(argument(i+1))//', '//trim(option_param(j))
           pixwin_file_name_input=trim(option_param(j))
           if (pixwin_file_name_input == default) then ! YES -> correct, with default file
              correct_pix=.true.
           elseif (pixwin_file_name_input /= none) then ! filename -> correct, with filename
              default_pix_file=.false.
              correct_pix=.true.
           endif ! otherwise (NO) -> do not correct
        else
           call error_usage(option(j))
        endif
     endif
  endif

  inc=inc+1
  if ((.not.count) .and. j == inc) then
     if (create) then
        call create_option(j,2,'-pixelfile2',none) 
     else
        if (i < narguments) then
#ifdef PIO
           option_param(j)  = pie_param%pixelfile2
#else
           option_param(j)  = argument(i+1)
#endif
           pwf2_set = .true.
           !print*,'pixelfile2=  ', trim(argument(i+1))//', '//trim(option_param(j))
           pixwin_file_name_input2=trim(option_param(j))
           if (pixwin_file_name_input2 == default) then
              correct_pix2=.true.
           elseif (pixwin_file_name_input2 /= none) then
              default_pix_file2=.false.
              correct_pix2=.true.
           endif
        else
           call error_usage(option(j))
        endif
     endif
  endif

  inc=inc+1
  if ((.not.count) .and. j == inc) then
     if (create) then
        call create_option(j,2,'-polarization',none)
     else
        if (i < narguments) then
#ifdef PIO
           option_param(j) = trim(pie_param%polarization)
#else
           option_param(j)=argument(i+1)
#endif
           if (option_param(j) == default) then
              polarization=.true.
           elseif (option_param(j) /= none) then
              write(*,*) 'Improper polarization option: '//trim(option_param(j))//trim(none)
              print*,option_param(j) == none
              print*,trim(option_param(j)) == trim(none)
              CALL FATAL_ERROR
           endif
        else
           call error_usage(option(j))
        endif
     endif
  endif

  inc=inc+1
  if ((.not.count) .and. j == inc) then
     if (create) then
        call create_option(j,2,'-subav',none)
     else
        if (i < narguments) then
#ifdef PIO
           if (pie_param%flag_subav) then
              option_param(j) = pie_param%subav
           else
              option_param(j) = none
           endif
           option_param(j) = set_optional_dmc_keyword(pie_param%flag_subav, pie_param%subav, none)
#else
           option_param(j)=argument(i+1)
#endif
           if (option_param(j) == default) then
              subtract_average=.true.
           elseif (option_param(j) /= none) then
              write(*,*) 'Improper subav option.'
              CALL FATAL_ERROR
           endif
        else
           call error_usage(option(j))
        endif
     endif
  endif

  inc=inc+1
  if ((.not.count) .and. j == inc) then
     if (create) then
        call create_option(j,2,'-subdipole',none)
     else
        if (i < narguments) then
#ifdef PIO
           option_param(j) = set_optional_dmc_keyword(pie_param%flag_subdipole, pie_param%subdipole, none)
#else
           option_param(j)=argument(i+1)
#endif
           if (option_param(j) == default) then
              subtract_dipole=.true.
           elseif (option_param(j) /= none) then
              write(*,*) 'Improper subdipole option.'
              CALL FATAL_ERROR
           endif
        else
           call error_usage(option(j))
        endif
     endif
  endif

  inc=inc+1
  if ((.not.count) .and. j == inc) then
     if (create) then
        call create_option(j,2,'-symmetric_cl',none)
     else
        if (i < narguments) then
#ifdef PIO
           option_param(j) = set_optional_dmc_keyword(pie_param%flag_symmetric_cl, pie_param%symmetric_cl, none)
#else
           option_param(j)=argument(i+1)
#endif
           if (option_param(j) == default) then
              symmetric_cl=.true.
           elseif (option_param(j) /= none) then
              write(*,*) 'Improper symmetric_cl option.'
              CALL FATAL_ERROR
           endif
        else
           call error_usage(option(j))
        endif
     endif
  endif

  inc=inc+1
  if ((.not.count) .and. j == inc) then
     if (create) then
        call create_option(j,2,'-tenormfilein',none)
     else
        if (i < narguments) then
#ifdef PIO
           if (pie_param%flag_tenormfilein) option_param(j)     = pie_param%tenormfilein
#else
           option_param(j)=argument(i+1)
#endif
           tenorm_in_file_input=trim(option_param(j))
           if (tenorm_in_file_input == default) then
              tenorminput=.true.
           elseif (tenorm_in_file_input /= none) then
              default_tenormin_file=.false.
              tenorminput=.true.
           endif
        else
           call error_usage(option(j))
        endif
     endif
  endif

  inc=inc+1
  if ((.not.count) .and. j == inc) then
     if (create) then
        call create_option(j,2,'-tenormfileout',none)
     else
        if (i < narguments) then
#ifdef PIO
           if (pie_param%flag_tenormfileout) option_param(j)     = pie_param%tenormfileout
#else
           option_param(j)=argument(i+1)
#endif
           tenorm_out_file_input=trim(option_param(j))
           if (tenorm_out_file_input == default) then
              tenormoutput=.true.
           elseif (tenorm_out_file_input /= none) then
              default_tenormout_file=.false.
              tenormoutput=.true.
           endif
        else
           call error_usage(option(j))
        endif
     endif
  endif
  !---

  inc=inc+1
  if ((.not.count) .and. j == inc) then
     if (create) then
        call create_option(j,2,'-tf_file',none) 
     else
        if (i < narguments) then
#ifdef PIO
           if (pie_param%flag_tf_file) option_param(j)     = pie_param%tf_file
           tf_file=trim(option_param(j))
           correct_transfer_function=pie_param%flag_tf_file
#else
           option_param(j)=argument(i+1)
           tf_file=trim(option_param(j))
           if (trim(tf_file) /= none) then
              correct_transfer_function=.true.
           endif
#endif
        else
           call error_usage(option(j))
        endif
     endif
  endif

  inc=inc+1
  if ((.not.count) .and. j == inc) then
     if (create) then
        call create_option(j,2,'-thetamax',none)
     else
        if (i < narguments) then
#ifdef PIO
           option_param(j) = pie_param%thetamax
#else
           option_param(j)=argument(i+1)
#endif
           if (option_param(j) /= none) then
              read(option_param(j),*,err=1) thetamax
              if ((thetamax <= 0.0_DP).or.(thetamax > 180.0_DP)) then
                 write(*,*) 'Improper thetamax value.'
                 CALL FATAL_ERROR
              endif
           endif
        else
           call error_usage(option(j))
        endif
     endif
  endif

  inc=inc+1
  if ((.not.count) .and. j == inc) then
     if (create) then
        call create_option(j,2,'-tolerance',none)
     else
        if (i < narguments) then
#ifdef PIO
           option_param(j) = pie_param%tolerance
#else
           option_param(j)=argument(i+1)
#endif
           if (option_param(j) /= none) then
              read(option_param(j),*,err=1) tolerance
              if ((tolerance <= 0.0_DP).or.(tolerance > 1.0_DP)) then
                 write(*,*) 'Improper tolerance value.'
                 CALL FATAL_ERROR
              endif
           endif
        else
           call error_usage(option(j))
        endif
     endif
  endif

  inc=inc+1
  if ((.not.count) .and. j == inc) then
     if (create) then
        call create_option(j,2,'-verbosity',default)
        option_param(j)='1'
     else
        if (i < narguments) then
#ifdef PIO
           ! mapping PIE verbosity (0 or 1) -> Spice verbosity (1 or 2)
           pie_param%verbosity = min(pie_param%verbosity + 1, 2)
           option_param(j)=string(int(pie_param%verbosity))
#else
           option_param(j)=argument(i+1)
#endif
           if (option_param(j) == default) then
              nverbose=1
           elseif (option_param(j) == none) then
              nverbose=0
           else
              read(option_param(j),*,err=1) nverbose
           endif
           if (nverbose.eq.0) then
              verbose=.false.
           elseif (nverbose.eq.2) then
              verbose=.true.
              megaverbose=.true.
           elseif (nverbose.ne.1) then
              write(*,*) 'Improper verbose option.'
              CALL FATAL_ERROR
           endif
        else
           call error_usage(option(j))
        endif
     endif
  endif

  inc=inc+1
  if ((.not.count) .and. j == inc) then
     if (create) then
        call create_option(j,2,'-weightfile',none)
     else
        if (i < narguments) then
#ifdef PIO
           if (pie_param%flag_weightfile) option_param(j)     = pie_param%weightfile
           weightfile_input    = trim(option_param(j))
           weights_present     = pie_param%flag_weightfile
           default_weight_file = .false.
#else
           option_param(j)=argument(i+1)
           weightfile_input=trim(option_param(j))
           if (weightfile_input /= default) then
              if (weightfile_input /= none) then
                 default_weight_file=.false.
                 weights_present=.true.
              endif
           else 
              weights_present=.true.
              option_param(j)=trim(weightfile_default)//'.fits'
           endif
#endif
        else
           call error_usage(option(j))
        endif
     endif
  endif

  inc=inc+1
  if ((.not.count) .and. j == inc) then
     if (create) then
        call create_option(j,2,'-weightfilep',default)
     else
        if (i < narguments) then
#ifdef PIO
           if (pie_param%flag_weightfilep) option_param(j)   = pie_param%weightfilep
           weightfilep_input   = trim(option_param(j))
           weightfilep_toread  = pie_param%flag_weightfilep
#else
           option_param(j)   = argument(i+1)
           weightfilep_input   = trim(option_param(j))
           weightfilep_toread  = (weightfilep_input /= default)
#endif
           !print*,'weightfilep=',trim(argument(i+1)),trim(weightfilep_input),weightfilep_toread
!            if (trim(option_param(j)) == none) then
!               print*,trim(none)//' is not a valid value for '//trim(option(j))
!               call fatal_error
!            endif
        else
           call error_usage(option(j))
        endif
     endif
  endif

  inc=inc+1
  if ((.not.count) .and. j == inc) then
     if (create) then
        call create_option(j,2,'-weightfile2',none)
     else
        if (i < narguments) then
#ifdef PIO
           if (pie_param%flag_weightfile2) option_param(j)     = pie_param%weightfile2
#else
           option_param(j)=argument(i+1)
#endif
           weightfile2_input=trim(option_param(j))
           if (weightfile2_input /= default) then
              if (weightfile2_input /= none) then
                 default_weight_file2=.false.
                 weights2_present=.true.
              endif
           else 
              weights2_present=.true.
              option_param(j)=trim(weightfile2_default)//'.fits'
           endif
        else
           call error_usage(option(j))
        endif
     endif
  endif

  inc=inc+1
  if ((.not.count) .and. j == inc) then
     if (create) then
        call create_option(j,2,'-weightfilep2',default)
     else
        if (i < narguments) then
#ifdef PIO
           if (pie_param%flag_weightfilep2) option_param(j)    = pie_param%weightfilep2
           weightfilep2_input   = trim(option_param(j))
           weightfilep2_toread  = pie_param%flag_weightfilep2
#else
           option_param(j)    = argument(i+1)
           weightfilep2_input   = trim(option_param(j))
           weightfilep2_toread  = (weightfilep2_input /= default)
#endif
!            if (trim(option_param(j)) == none) then
!               print*,trim(none)//' is not a valid value for '//trim(option(j))
!               call fatal_error
!            endif
        else
           call error_usage(option(j))
        endif
     endif
  endif

  inc=inc+1
  if ((.not.count) .and. j == inc) then
     if (create) then
        call create_option(j,2,'-weightpower',string(weightpower))
     else
        if (i < narguments) then
#ifdef PIO
           option_param(j) = string(pie_param%weightpower)
#else
           option_param(j) = argument(i+1)
#endif
           if (option_param(j) /= none) then
              read(option_param(j),*,err=1) weightpower
           else
              option_param(j) = string(weightpower)
           endif
        else
           call error_usage(option(j))
        endif
     endif
  endif

  inc=inc+1
  if ((.not.count) .and. j == inc) then
     if (create) then
        call create_option(j,2,'-weightpower2',string(weightpower2))
     else
        if (i < narguments) then
#ifdef PIO
           option_param(j) = string(pie_param%weightpower2)
#else
           option_param(j) = argument(i+1)
#endif
           if (option_param(j) /= none) then
              read(option_param(j),*,err=1) weightpower2
           else
              option_param(j) = string(weightpower2)
           endif
        else
           call error_usage(option(j))
        endif
     endif
  endif

  inc=inc+1
  if ((.not.count) .and. j == inc) then
     if (create) then
        call create_option(j,2,'-weightpowerp',string(weightpowerp))
     else
        if (i < narguments) then
#ifdef PIO
           option_param(j) = string(pie_param%weightpowerp)
           set_wp_p1 = (pie_param%weightpowerp /= 1.0_dp)
#else
           option_param(j) = argument(i+1)
           set_wp_p1 = (option_param(j) /= default)
#endif
           if (set_wp_p1) read(option_param(j),*,err=1) weightpowerp
        else
           call error_usage(option(j))
        endif
     endif
  endif

  inc=inc+1
  if ((.not.count) .and. j == inc) then
     if (create) then
        call create_option(j,2,'-weightpowerp2',string(weightpowerp2))
     else
        if (i < narguments) then
#ifdef PIO
           option_param(j) = string(pie_param%weightpowerp2)
           set_wp_p2 = (pie_param%weightpowerp2 /= 1.0_dp)
#else
           option_param(j) = argument(i+1)
           set_wp_p2 = (option_param(j) /= default)
#endif
           if (set_wp_p2) read(option_param(j),*,err=1) weightpowerp2
        else
           call error_usage(option(j))
        endif
     endif
  endif

  inc=inc+1
  if ((.not.count) .and. j == inc) then
     if (create) then
        call create_option(j,2,'-windowfilein',none)
     else
        if (i < narguments) then
#ifdef PIO
           if (pie_param%flag_windowfilein) option_param(j)     = pie_param%windowfilein
#else
           option_param(j)=argument(i+1)
#endif
           window_in_file_input=trim(option_param(j))
           if (window_in_file_input == default) then
              windowinput=.true.
           elseif (window_in_file_input /= none) then
              default_windowin_file=.false.
              windowinput=.true.
           endif
        else
           call error_usage(option(j))
        endif
     endif
  endif

  inc=inc+1
  if ((.not.count) .and. j == inc) then
     if (create) then
        call create_option(j,2,'-windowfileout',none)
     else
        if (i < narguments) then
#ifdef PIO
           if (pie_param%flag_windowfileout) option_param(j)     = pie_param%windowfileout
#else
           option_param(j)=argument(i+1)
#endif
           window_out_file_input=trim(option_param(j))
           if (window_out_file_input == default) then
              windowoutput=.true.
           elseif (window_out_file_input /= none) then
              default_windowout_file=.false.
              windowoutput=.true.
           endif
        else
           call error_usage(option(j))
        endif
     endif
  endif
  !---


  if (count) noption=inc

  return
1 call error_usage(option(j))
end subroutine make_option_active

!=======================================================================
subroutine create_option(j,iskip,command,param)
!=======================================================================
  use spice_common, only: I4B, option, option_param, iskip_option
  implicit none

  character(len=*), intent(IN) :: command,param
  integer(I4B),     intent(IN) :: iskip,j
  
  option(j)=command
  option_param(j)=param
  iskip_option(j)=iskip
  
end subroutine create_option

!=======================================================================
subroutine check_options_choice
!=======================================================================
  use misc_utils, only: fatal_error
  use spice_common
  implicit none
#ifdef PIO
  logical test_t, test_q, test_u
#endif

  if (.not. alm1_dump)   alm1_out_file=''
  if (.not. alm2_dump)   alm2_out_file=''
  if (.not. alm1_weight) alm1_weight_file=''
  if (.not. alm2_weight) alm2_weight_file=''

  if (decouple .and. (.not.polarization) .and. megaverbose) then
     write(*,*) 'WARNING in check_options_choice :'
     write(*,*) 'Polarization is off and decouple parameter is activated.'
     write(*,*) 'Decoupling will be ignored by default.'
  endif
  decouple=decouple.and.polarization

!   if (.not.decouple.and.windowinput.and.megaverbose) then
!      write(*,*) 'WARNING in check_options_choice :'
!      write(*,*) 'Decoupling is off and input of window function is activated.'
!      write(*,*) 'Input of window function will be ignored by default.'
!   endif
!   windowinput=windowinput.and.decouple
!   windowoutput=windowoutput.and.(.not.windowinput).and.decouple
  if (windowinput .and. (.not.polarization) .and. megaverbose) then
     write(*,*) 'WARNING in check_options_choice :'
     write(*,*) 'Polarization is off and input of window function is activated.'
     write(*,*) 'Input of TE window function will be ignored.'
  endif
  windowinput = windowinput .and. polarization

  if (windowoutput .and. (.not.polarization) .and. megaverbose) then
     write(*,*) 'WARNING in check_options_choice :'
     write(*,*) 'Polarization is off and output of window function is activated.'
     write(*,*) 'Output of TE window function will be ignored.'
  endif
  windowoutput= windowoutput .and. (.not.windowinput) .and. polarization

  if (tenorminput .and. (.not.polarization) .and. megaverbose) then
     write(*,*) 'WARNING in check_options_choice :'
     write(*,*) 'Polarization is off and input of tenorm function is activated.'
     write(*,*) 'Input of TE tenorm function will be ignored.'
  endif
  tenorminput = tenorminput .and. polarization

  if (tenormoutput .and. (.not.polarization) .and. megaverbose) then
     write(*,*) 'WARNING in check_options_choice :'
     write(*,*) 'Polarization is off and output of tenorm function is activated.'
     write(*,*) 'Output of TE tenorm function will be ignored.'
  endif
  tenormoutput= tenormoutput .and. (.not.tenorminput) .and. polarization

! commented out on 2017-04-05
!   if (kernelsoutput .and. (.not.polarization) .and. megaverbose) then
!      write(*,*) 'WARNING in check_options_choice :'
!      write(*,*) 'Polarization is off and output of kernels is activated.'
!      write(*,*) 'Output of kernels will be ignored.'
!   endif
!   kernelsoutput= kernelsoutput .and. polarization

  if (kernelsoutput .and. windowinput  .and. megaverbose) then
     write(*,*) 'WARNING in check_options_choice :'
     write(*,*) 'output of kernels is ON and input of window is ON.'
     write(*,*) 'Calculation and Output of kernels will be ignored.'
  endif
  kernelsoutput= kernelsoutput .and. (.not.windowinput)

  if (kernelsoutput .and. tenorminput .and. megaverbose) then
     write(*,*) 'WARNING in check_options_choice :'
     write(*,*) 'output of kernels is ON and input of tenorm is ON.'
     write(*,*) 'Calculation and Output of kernels will be ignored.'
  endif
  kernelsoutput= kernelsoutput .and. (.not.tenorminput)

!  commented out 2017-08-09, allow windowfileout+kernelsfileout
!   if (kernelsoutput .and. windowoutput .and. megaverbose) then
!      write(*,*) 'WARNING in check_options_choice :'
!      write(*,*) 'output of kernels is ON and output of window is ON.'
!      write(*,*) 'Window output will be ignored.'
!   endif
!   windowoutput= windowoutput .and. (.not.kernelsoutput)

  if (pairsthresholding .and. megaverbose) then
     write(*,*) 'WARNING in check_options_choice:'
     write(*,*) 'npairsthreshold is now deprecated,'
     write(*,*) 'keyword value will be ignored.'
  endif
  pairsthresholding = .false.
  npairsthreshold = 0.0_dp

  if (clmapinput) then
     map_present=.false.
     readQfile=.false. !SP
     readUfile=.false. !SP
  endif


  clmapoutput  = clmapoutput  .and. .not. clmapinput
  clmaskoutput = clmaskoutput .and. .not. clmaskinput

  if (.not.masks_present.and..not.weights_present) clmaskoutput=.false.

  if (polarization) then
     nmap=3 ! I,Q,U
     ncor=6 ! TT, EE, BB, TE, TB, EB
     if (.not.symmetric_cl .and. (map2_present .or. extramap2_present .or. nlmaps(1,2)> 0)) then
        ncor=9 ! TT, EE, BB, TE, TB, EB, ET, BT, BE
     endif
  else 
     nmap=1 ! I
     ncor=1 ! TT
  endif

  ! set weightpower for polarization to same value as T ones, 
  ! unless they have been set by the user
  if (polarization) then
     if (.not. set_wp_p1) weightpowerp  = weightpower
     if (.not. set_wp_p2) weightpowerp2 = weightpower2
!   else
!      if (set_wp_p1 .or. set_wp_p2) then
!      endif
  endif


#ifdef PIO
  test_t = (map_present .or. (nlmaps(1,1)>0)) ! either mapfile  or listmapfiles1_*
  test_q = (readQfile   .or. (nlmaps(2,1)>0)) ! either mapQfile or listmapfilesQ1_*
  test_u = (readUfile   .or. (nlmaps(3,1)>0)) ! either mapUfile or listmapfilesU1_*
  !!print*,map_present,readQfile,readUfile,nlmaps(1:3,1)
  !!print*,test_t, test_q, test_u
  if (.not.test_t) then
     if ((test_q.or.test_u).and.polarization) then
        write(*,*) 'WARNING in check_options_choice :'
        write(*,*) 'there is no T map and one of the Q/U file'
        write(*,*) 'is specified. It will be ignored'        
        readQfile=.false.
        readUfile=.false.
        nlmaps(2:3,1) = 0
     endif
  endif

  if (polarization.and.(.not.clmapinput)) then !SP
     if (.not.test_q.or..not.test_u) then
        write(*,*) 'ERROR in check_options_choice :'
        write(*,*) 'Polarization is active and one of the Q/U file'
        write(*,*) 'is not specified'
        CALL FATAL_ERROR
     endif
  else
     if ((test_q.or.test_u).and.megaverbose) then
        write(*,*) 'WARNING in check_options_choice :'
        write(*,*) 'Polarization is inactive and one of the Q/U file'
        write(*,*) 'is specified. It will be ignored'
        readQfile=.false.
        readUfile=.false.
        nlmaps(2:3,1) = 0
     endif
  endif
#endif
end subroutine check_options_choice

!=======================================================================
subroutine check_compute
!=======================================================================
  use spice_common
  implicit none

  if ( (.not.cloutput) .and. (.not.coroutput) .and. (.not.clmapoutput) &
 &     .and.(.not.clmaskoutput).and.(.not.output_spice_rc) &
 &     .and.(.not.clmapoutput).and.(.not.clmaskoutput) &
 &     .and.(.not.windowoutput) ) then
     write(*,*) '    Spice '//version
     write(*,*) 'No output required ? No calculation required. ;0)'
#ifndef PIO
     write(*,*) 'Type '
     write(*,*) '  spice  -usage'
     write(*,*) 'for a short help, or'
     write(*,*) '  spice  -help'
     write(*,*) 'for extended help.'
#endif
     STOP
  endif
  
end subroutine check_compute

!=======================================================================
subroutine create_spicerc
!=======================================================================
  use spice_common
  implicit none
  integer(I4B) :: lunit=10
  integer(I4B) :: j

  call get_generic_filename(spice_rc_file_default, &
 &                          spice_rc_outfile_input, &
 &                          spice_rc_outfile,default_spice_rc_out,.false.,.false.,overwrite,stringlength,0)

  if (verbose) write(*,*) 'Create options file ',trim(spice_rc_outfile)

  open(unit=lunit,file=spice_rc_outfile,form='formatted',status='unknown')

  do j=jstartnormal,noption
     write(lunit,'(a)') trim(option(j)(2:))//'  '//trim(option_param(j))
  enddo
  
  close(lunit)

end subroutine create_spicerc

!=======================================================================
subroutine read_spicerc
!=======================================================================
! routine to parse parameter file
! modified (and commented) by EH, July 2003, 
! to emulate more closely Healpix parameter files
!  - ignore lines starting with #
!  - accept syntax : key [=] value
! Original wrapper works for multiple parameters options. This modified 
! version works only for single parameter options.
!=======================================================================
  use misc_utils, only: fatal_error
  use spice_common
  implicit none
  integer(I4B) :: lunit=10
  logical :: eof,option_ok
  logical, allocatable, dimension(:) :: option_done_in
  character(len=FILENAMELEN), allocatable, dimension(:) :: argument
  integer(I4B) :: narguments,i,iskip,j,k,idx
  character(len=FILENAMELEN) :: inargument

  call get_generic_filename(spice_rc_file_default, &
 &                          spice_rc_infile_input, &
 &                          spice_rc_infile,default_spice_rc_in,.false.,.true.,.false.,stringlength,0)

  open(unit=lunit,file=spice_rc_infile,form='formatted',status='old')
  
  ! count the number of words in file
  narguments=0
  do 
     read(lunit,'(a)',end=1) inargument
     inargument = adjustl(inargument) ! remove leading blanks
     if (inargument(1:1) /= '#') then ! ignore lines starting with '#'
        idx = index(inargument,'=')
        if (idx > 0) then
           j = 2
        else
           j=0
           do i=2,len_trim(inargument)+1
              if (inargument(i:i) == ' '.and.inargument(i-1:i-1) /= ' ') &
                   &      j=j+1
           enddo
        endif
        narguments=narguments+j
     endif
  enddo
  
!!!! 1 rewind(lunit) ! ignored by some compilers (ifort 15.0.4 ?)
1 continue
  close(lunit)
  open(unit=lunit,file=spice_rc_infile,form='formatted',status='old')

  ! read argument key and value into 'argument' array
  allocate(argument(narguments))
  narguments=0
  do 
     read(lunit,'(a)',end=2) inargument
     inargument = adjustl(inargument) ! remove leading blanks
     if (inargument(1:1) /= '#') then ! ignore lines starting with '#'
        idx = index(inargument,'=')
        if (idx > 0) then
           argument(narguments+1) = '-'//trim(adjustl(inargument(1:idx-1)))
           argument(narguments+2) =      &
                & trim(adjustl(inargument(idx+1:len_trim(inargument))))
           narguments = narguments + 2
        else
           k=1
           j=0
           do i=2,len_trim(inargument)+1
              if (inargument(i:i) == ' '.and.inargument(i-1:i-1) /= ' ') then
                 narguments=narguments+1
                 j=j+1
                 argument(narguments)=inargument(k:i-1)
              endif
              if (inargument(i:i) /= ' '.and.inargument(i-1:i-1) == ' ') k=i
           enddo
           ! replace key by -key
           if (j > 0) argument(narguments-j+1)= &
                & '-'//trim(argument(narguments-j+1))
        endif
     endif
  enddo


2 close(lunit)

  allocate(option_done_in(noption))
  option_done_in=.false.
  i=0
  iskip=1
  do while (i+iskip.le.narguments) ! loop on keywords in user-generated parameter file
     i=i+iskip
     option_ok=.false.
     do j=1,noption ! loop on all available keywords
        if (trim(argument(i)) == trim(option(j))) then
           if (.not.option_done(j)) then
              if (.not.option_done_in(j)) then
                 call make_option_active(i,j,iskip,argument,narguments, &
 &                                       .false.,.false.)
                 option_done_in(j)=.true.
                 option_ok=.true.
              endif
           else 
              if (.not.option_done_in(j)) then
                 iskip=iskip_option(j)
                 option_done_in(j)=.true.
                 option_ok=.true.
              endif
           endif
        endif
     enddo
     if (.not.option_ok) then
        write(*,*) 'ERROR : unknown or duplicate option read : '//trim(argument(i)(2:))
        CALL FATAL_ERROR
     endif
  enddo

  do j=1,noption
     option_done(j)=option_done(j).or.option_done_in(j)
  enddo

  deallocate(option_done_in)
  deallocate(argument)
  if (megaverbose) write(*,*) 'Input option file was read : ',trim(spice_rc_infile)

end subroutine read_spicerc
!------------------------------------------------------
subroutine reassess_options
  use spice_common
  use misc_utils, only: fatal_error

  nmask1 = 1 ! default
  nmask2 = 1 ! default
  nweight1 = 1
  nweight2 = 1
  nmask = 1 ! default
  nweight = 1
  if (polarization) then
     ! check dimension of mask map and C(l)
     if (clmaskinput) then
        ! ncmask has been determined from file containing mask C(l)
        if (ncmask /= 1 .and. ncmask /= 3 .and. ncmask /=4) then
           print*,' invalid values of ncmask, expected 1, 3 or 4, got:', ncmask
           call fatal_error
        endif
        if (ncmask >= 3) then
           nmask = 2 ! should not be necessary
           nweight = 2 ! should not be necessary
           nmask1 = 2 ! added 2008-10-01
           nmask2 = 2 ! added 2008-10-01
        endif
     else
        ncmask = 1
        ! different mask for temperature and polarization
        if (masks_present .and. maskfilep_toread)   nmask1 = 2
        if (masks2_present .and. maskfilep2_toread) nmask2 = 2
        nmask = max(nmask1, nmask2) ! 1 or 2
        ! different weight for temperature and polarization
        if (weights_present .and. weightfilep_toread)   nweight1 = 2
        if (weights2_present .and. weightfilep2_toread) nweight2 = 2
        nweight = max(nweight1, nweight2) ! 1 or 2
        ! mask*weight C(l) must account for polarization
        if (nweight > 1 .or. nmask > 1) then
           ncmask = 3
           if (ncor > 6) ncmask = 4 ! treat both (and separately TE and ET)
        endif
     endif
  else
     ! if no polarization, ignore polarization specific mask
     ncmask = 1
     maskfilep_toread = .false.
     maskfilep2_toread = .false.
     weightfilep_toread = .false.
     weightfilep2_toread = .false.
  endif

!  print*,'nmask, ncmask',nmask, ncmask
  !-------------------------------

  return
end subroutine reassess_options

#ifdef PIO
function set_optional_dmc_keyword(flag, value, default)
  use healpix_types
  ! returns value if flag=True, and default otherwise
  ! cf file_n3
  character(len=FILENAMELEN)   :: set_optional_dmc_keyword
  logical(LGT), intent(in)     :: flag
  character(len=*), intent(in) :: value
  character(len=*), intent(in) :: default

  if (flag) then
     set_optional_dmc_keyword = trim(value)
  else
     set_optional_dmc_keyword = trim(default)
  endif

  return
end function set_optional_dmc_keyword
#endif
