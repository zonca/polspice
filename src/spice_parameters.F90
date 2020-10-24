module spice_parameters
  use healpix_types
  implicit none

#if DOUBLE
  integer(I4B),     parameter          :: KMAP  =   DP
  character(len=*), parameter, private :: vprec = '_DP'
#else
  integer(I4B),     parameter          :: KMAP  =   SP
  character(len=*), parameter, private :: vprec = '_SP'
#endif
  integer(I4B), parameter :: KMAPC = KMAP

  character(len=FILENAMELEN) :: healpix_data
  !character(len=*), parameter :: version='2.9.1'
  !character(len=*), parameter :: version='2.9.2' ! multi extra files
  !character(len=*), parameter :: version='3.0.1' ! many bug corrections
  !character(len=*), parameter :: version='3.0.3' ! fix un-initialized variables
  !character(len=*), parameter :: version='3.1.0' ! computes TE and ET, let user pick tolerance
  !character(len=*), parameter :: version='3.1.3' ! added SUBDIPOL in output FITS file,replace rewind by close+open
  !character(len=*), parameter :: version='3.1.4' ! fixed DMC specific bug
  !character(len=*), parameter :: version='3.1.5' ! fixed another DMC specific bug (mapQfile2)
  !character(len=*), parameter :: version='3.1.6' ! fixed listmapweights2 bug
  !character(len=*), parameter :: version='3.2.0' ! read correctly polarized cut sky map
  !character(len=*), parameter :: version='3.3.0' ! different Nside(s)
  !character(len=*), parameter :: version='3.3.1' ! slight edits to documentation on -tolerance and -decouple
  !character(len=*), parameter :: version='3.3.2' ! computes kernel in T only case as well
  !character(len=*), parameter :: version='3.3.3' ! skip empty Healpix rings
  !character(len=*), parameter :: version='3.3.4' ! added tenormfilein and tenormfileout
  !character(len=*), parameter :: version='3.4.0' ! faster computation of kernels (and TEnorm)
  !character(len=*), parameter :: version='3.4.1' ! more portable Makefile, cosmetic editions
  !character(len=*), parameter :: version='3.5.0'//vprec ! allow use of double precision Ylm transforms
  !character(len=*), parameter :: version='3.5.1'//vprec ! cleaned Makefile and CMakeLists.txt
  !character(len=*), parameter :: version='3.5.2'//vprec ! ispice.py runs in python 2 and 3
  !character(len=*), parameter :: version='3.6.1'//vprec ! reduced memory footprint on maps, requires Healpix 3.60
  !character(len=*), parameter :: version='3.6.2'//vprec ! bug correction for polarization mask and weight
  !character(len=*), parameter :: version='3.6.3'//vprec ! compliant with gfortran
  !character(len=*), parameter :: version='3.6.4'//vprec ! debugged mask and weight reading
  !character(len=*), parameter :: version='3.6.5'//vprec ! can mix full sky and partial sky maps
  !character(len=*), parameter :: version='3.6.6'//vprec ! debugged text files parsing
  !character(len=*), parameter :: version='3.6.7'//vprec ! fixed segfault in ClCl covariance matrix, improved ispice.py
  !character(len=*), parameter :: version='3.6.8'//vprec ! sort pixels of partial maps (and test for duplicates)
  !character(len=*), parameter :: version='3.6.9'//vprec ! dump alm's of (masked and weighted) maps info FITS files
  character(len=*), parameter :: version='3.7.0'//vprec ! fixed bugs with alm*fileout


  character(len=*),parameter :: mapfile_default='map'
  character(len=*),parameter :: mapfile2_default='map2'
  character(len=*),parameter :: alm1file_default='alm1.fits'
  character(len=*),parameter :: alm2file_default='alm2.fits'
  character(len=*),parameter :: alm1wf_default  ='almw1.fits'
  character(len=*),parameter :: alm2wf_default  ='almw2.fits'
  character(len=*),parameter :: maskfile_default='mask'
  character(len=*),parameter :: maskfile2_default='mask2'
  character(len=*),parameter :: weightfile_default='weight'
  character(len=*),parameter :: weightfile2_default='weight2'
  character(len=*),parameter :: cor_file_name_default='spice.cor'
  character(len=*),parameter :: cl_file_name_default='spice.cl'
  character(len=*),parameter :: clmap_file_name_default='spice.clrawmap'
  character(len=*),parameter :: clmask_file_name_default='spice.clrawmasks' 
  character(len=*),parameter :: window_file_name_default='spice.window'
  character(len=*),parameter :: tenorm_file_name_default='spice.tenorm'
  character(len=*),parameter :: kernels_file_name_default='spice.kernels'
  character(len=*),parameter :: cov_file_name_default='spice.covariance'
  character(len=*),parameter :: pixwin_file_name_default='pixel_window_n'
  character(len=*),parameter :: pixwin_file_name_default2='pixel_window_n'
  character(len=*),parameter :: noisecor_file_name_default='noise.cor'
  character(len=*),parameter :: noisecl_file_name_default='noise.cl'
  character(len=*),parameter :: spice_rc_file_default='.spicerc'
  character(len=*),parameter :: default='YES'
  character(len=*),parameter :: none='NO'
  integer(I4B), parameter :: stringlength=256
  integer(I4B), parameter :: mytype_junk=0, mytype_table=1, mytype_table3=10
  integer(I4B), parameter :: mytype_vect=2, mytype_cl=3, mytype_corr=4
  integer(I4B), parameter :: mytype_vect_tuple=12, mytype_cl_tuple=13, mytype_corr_tuple=14
  integer(I4B), parameter :: nefmax = 10 ! max number of extra mapfiles to add 
  logical :: alm1_dump=.false.
  logical :: alm2_dump=.false.
  logical :: alm1_weight=.false.
  logical :: alm2_weight=.false.
  logical :: output_spice_rc=.false.
  logical :: default_spice_rc_out=.true.
  logical :: input_spice_rc=.false.
  logical :: default_spice_rc_in=.true.
  logical :: verbose=.true.
  logical :: megaverbose=.false.
  logical :: map_present=.true.
  logical :: default_map_file=.true.
  logical :: map2_present=.false.
  logical :: default_map_file2=.true.

  logical :: masks_present=.false.
  logical :: maskfilep_toread=.false.
  logical :: maskfilep2_toread=.false.
  logical :: default_mask_file=.true.
  logical :: masks2_present=.false.
  logical :: default_mask_file2=.true.
  logical :: pwf2_set=.false.

  logical ::      extramap1_present =.false., extramap2_present =.false.
  integer(i4b) :: extramap1_number  =0,       extramap2_number  =0
  logical ::      readextraQfile1=.false.,    readextraUfile1=.false.
  logical ::      readextraQfile2=.false.,    readextraUfile2=.false.
  integer(i4b) :: extramapQ1_number  =0,       extramapU1_number  =0
  integer(i4b) :: extramapQ2_number  =0,       extramapU2_number  =0

  logical :: weights_present=.false.
  logical :: weightfilep_toread=.false.
  logical :: weightfilep2_toread=.false.
  logical :: default_weight_file=.true.
  logical :: weights2_present=.false.
  logical :: default_weight_file2=.true.

  logical :: coroutput=.false.
  logical :: default_cor_file=.true.
  logical :: cloutput=.false.
  logical :: default_cl_file=.true.
  logical :: clmapoutput=.false.
  logical :: default_clmapout_file=.true.
  logical :: clmapinput=.false.
  logical :: default_clmapin_file=.true.
  logical :: clmaskoutput=.false.
  logical :: default_clmaskout_file=.true.
  logical :: clmaskinput=.false.
  logical :: default_clmaskin_file=.true.

  logical :: kernelsoutput=.false.
  logical :: default_kernelsout_file=.true.
  logical :: fullkernels = .false.
  logical :: windowoutput=.false.
  logical :: default_windowout_file=.true.
  logical :: windowinput=.false.
  logical :: default_windowin_file=.true.
  logical :: tenormoutput=.false.
  logical :: default_tenormout_file=.true.
  logical :: tenorminput=.false.
  logical :: default_tenormin_file=.true.

  logical :: correct_pix=.false.
  logical :: correct_pix2=.false.
  logical :: default_pix_file=.true.
  logical :: default_pix_file2=.true.
  logical :: correct_beam=.false.
  logical :: correct_beam2=.false.
  logical :: normalize=.false.
  logical :: subtract_average=.false.
  logical :: subtract_dipole=.false.
  logical :: subtract_noise_cor=.false.
  logical :: default_noisecor_file=.true.
  logical :: subtract_noise_cl=.false.
  logical :: default_noisecl_file=.true.

  real(DP) :: weightpower=1.d0
  real(DP) :: weightpower2=1.d0
  real(DP) :: weightpowerp=1.d0
  real(DP) :: weightpowerp2=1.d0
  logical  :: set_wp_p1=.false.
  logical  :: set_wp_p2=.false.

  logical ::   pairsthresholding=.false.
  real(DP) :: npairsthreshold=0.d0
  logical :: apodize=.false.
  real(DP) :: apodizesigma=-1.d0
  logical :: decouple=.false.
  real(DP) :: thetamax = 180.d0
  logical :: dry=.false.
  logical :: polarization=.false.
  integer(I4B) :: nlmax = -1
  character(len=FILENAMELEN) :: beam_file=none, beam_file2=none
  logical :: beam_present=.false.
  logical :: beam2_present=.false.
  integer(I4B) :: apodizetype = 0
  logical :: overwrite=.true.
  logical :: readQfile=.false.
  logical :: readUfile=.false.
  logical :: readQfile2=.false.
  logical :: readUfile2=.false.
  logical :: know_out_cl_grp=.false.
  real(DP) :: fwhm = 0.0d0
  real(DP) :: fwhm2 = 0.0d0
  logical :: fits_out = .false.
  logical :: correct_transfer_function = .false.
  real(DP) :: tolerance = 1.d-5

  logical :: do_cov = .false.
  logical :: default_covout_file=.true.
  logical :: symmetric_cl=.false.

  type maptype
     real(KMAP),   allocatable, dimension(:,:) :: map
     integer(I4B), allocatable, dimension(:)   :: pixel
     integer(I8B), allocatable, dimension(:)   :: pixel8
     integer(I8B)                              :: npix
     integer(I4B)                              :: nside, ordering
     integer(I4B)                              :: type, nmaps
     character(len=FILENAMELEN)                :: coordsys
     real(KMAP)                                :: fmissval
     real(DP), dimension(2) :: zbounds
     !real(DP)     :: THRESHOLD = 0.01_DP
     real(DP)     :: THRESHOLD = 0.25_DP
     integer(I4B) :: RING=1, NESTED=2
     integer(I4B) :: VOID=-1, FULL=0, CUT4B=10, CUT8B=11
  end type maptype

  type almtype
     complex(KMAPC), allocatable, dimension(:,:,:) :: almtab
     complex(KMAPC), allocatable, dimension(:,:)   :: almvec
     integer(I4B) :: lmax, mmax
     integer(I4B) :: type, npol
     integer(I4b) :: TAB=2, VECT=1, VOID=-1
  end type almtype

  contains
    subroutine init_system
      use extension, only: getEnvironment
      logical(LGT) :: found 
      integer :: i_try
      character(len=*), parameter :: prompt = '===> '
#ifdef PIO
    ! for DPC: use preprocessing variable HEALPIXDATA
!!!!!!     character(len=*), parameter :: hpxdir = 'HEALPIXDIR'
#ifdef HEALPIXDATA
     healpix_data = HEALPIXDATA
     healpix_data = trim(healpix_data) // '/'
#else
     print*,prompt//'Preprocessing variable HEALPIXDATA not defined'
     print*,prompt//'Aborting.'
     call exit(-1)
#endif

#else
     ! outside DPC: use 
     !       1)pre-processing variable HEALPIXDATA  
     !       2)environment variables HEALPIXDATA
     !       3)environment variables HEALPIX
     character(len=*), parameter :: hpxdir = 'HEALPIX'

#ifdef HEALPIXDATA
     healpix_data = HEALPIXDATA
#else
      call getEnvironment('HEALPIXDATA', healpix_data)
#endif
      if (trim(healpix_data) == '') then
         print*,prompt//'HEALPIXDATA not defined, trying '//hpxdir//'...'
         call getEnvironment(hpxdir, healpix_data)
         if (trim(healpix_data) /= '') then
            healpix_data = trim(healpix_data)//'/data/'
         else
            print*,prompt//hpxdir//' not defined, trying default values ...'
            found = .false.

            i_try = 1
            do 
               if (.not. found) then
                  if (i_try == 1) healpix_data = '/usr/local/src/Healpix/data'
                  if (i_try == 2) healpix_data = '/usr/local/src/Healpix_2.00/data'
                  if (i_try == 3) healpix_data = '/usr/local/lib/Healpix/data'
                  if (i_try > 3) exit
                  print*,trim(healpix_data)
                  inquire(file=trim(healpix_data)//'/pixel_window_n0002.fits',exist=found)
                  if (found) exit
                  i_try = i_try + 1
               endif
            enddo

            if (.not. found)  then
               print*,prompt//'Failed to find Healpix directory. Abort'
               stop
            else
               healpix_data = trim(healpix_data) // '/'
            endif
         endif
      endif
#endif
      print*,prompt//'Using HEALPIXDATA='//trim(healpix_data)

    end subroutine init_system

end module spice_parameters
