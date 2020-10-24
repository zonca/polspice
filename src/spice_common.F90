module spice_common
! edited to run with Healpix 1.2 2003-06
  use healpix_types
  use spice_parameters
#ifdef PIO
  use pie_io, only : parContent
  use piolib
#endif
  implicit none

#ifdef PIO
  ! define PIE structure containing all parameters read from parameter file
  type(parContent)                              :: pie_param
  INTEGER(KIND=PIOLONG)                         :: pio_error ! error code
  INTEGER(KIND=PIOLONG)                         :: pio_log   ! pointer toward logger
  character(LEN=1024)                           :: pio_logStr ! string to pass on to logger
#endif

  integer(I4B), dimension(1:2) :: nsides
  integer(I4B), dimension(1:2) :: npixtots
  integer(I4B), dimension(1:2) :: obsnpixs, obsnpixs_m, obsnpixs_w
  integer(I4B) :: nmap,ncor
  real(KMAP), allocatable, dimension(:,:)   :: mask_map,  weight_map ! edited
  real(KMAP), allocatable, dimension(:,:)   :: mask2_map, weight2_map ! edited
  real(KMAP), allocatable, dimension(:,:)   :: map_in, map2_in
  type(maptype) :: tmap1_in,     tmap2_in
  type(maptype) :: tmask1_map,   tmask2_map
  type(maptype) :: tweight1_map, tweight2_map
  real(KMAP), allocatable, dimension(:,:)   :: extramap1_in, extramap2_in
  real(DP), allocatable, dimension(:,:)   :: xi, xi_final, xi_noise, xi_mask
  real(DP), allocatable, dimension(:,:)   :: cl, cl_final, cl_noise, cl_mask
  real(DP), allocatable, dimension(:,:)   :: TEnorm
  real(DP), allocatable, dimension(:,:)   :: kcross
  real(DP), allocatable, dimension(:,:)   :: cl_map_precomp, cl_mask_precomp
  real(DP), allocatable, dimension(:)     :: mu,w
  real(DP), allocatable, dimension(:,:)   :: wlpix1,wlpix2,gb,gb2
  real(DP), allocatable, dimension(:,:,:) :: Pl
  real(DP), allocatable, dimension(:,:)   :: transfer_function
  real(DP), allocatable, dimension(:,:,:)   :: kernels
  real(DP), allocatable, dimension(:,:,:)   :: cov_matrix

  character(len=FILENAMELEN) :: spice_rc_outfile,spice_rc_infile
  character(len=FILENAMELEN) :: spice_rc_outfile_input,spice_rc_infile_input

  character(len=FILENAMELEN) :: maskfile, weightfile
  character(len=FILENAMELEN) :: maskfile2,weightfile2
  character(len=FILENAMELEN) :: pixwin_file_name, pixwin_file_name2
  character(len=FILENAMELEN) :: noisecor_file_name,noisecl_file_name
  character(len=FILENAMELEN) :: cor_file_name,cl_file_name

  character(len=FILENAMELEN) :: maskfile_input, weightfile_input
  character(len=FILENAMELEN) :: maskfile2_input,weightfile2_input
  character(len=FILENAMELEN) :: pixwin_file_name_input, pixwin_file_name_input2
  character(len=FILENAMELEN) :: noisecor_file_name_input,noisecl_file_name_input
  character(len=FILENAMELEN) :: cor_file_name_input,cl_file_name_input

  character(len=FILENAMELEN) :: cl_outmap_file, cl_outmask_file
  character(len=FILENAMELEN) :: cl_inmap_file, cl_inmask_file
  character(len=FILENAMELEN) :: cl_outmap_file_input, cl_outmask_file_input
  character(len=FILENAMELEN) :: cl_inmap_file_input, cl_inmask_file_input
  character(len=FILENAMELEN) :: tf_file

  character(len=FILENAMELEN) :: window_out_file, window_out_file_input
  character(len=FILENAMELEN) :: window_in_file,  window_in_file_input
  character(len=FILENAMELEN) :: tenorm_out_file, tenorm_out_file_input
  character(len=FILENAMELEN) :: tenorm_in_file,  tenorm_in_file_input

  character(len=FILENAMELEN) :: kernels_out_file, kernels_out_file_input
  character(len=FILENAMELEN) :: cov_out_file,     cov_out_file_input

  character(len=FILENAMELEN) :: alm1_out_file,    alm2_out_file
  character(len=FILENAMELEN) :: alm1_weight_file, alm2_weight_file

  !
  character(len=FILENAMELEN), dimension(1:nefmax,1:3, 1:2) :: &
       &                   listmapfile
  real(DP),     dimension(1:nefmax,1:3, 1:2) :: listmapw8
  integer(i4b), dimension(         1:3, 1:2) :: nlmaps


  character(len=FILENAMELEN) :: mapfile,  mapfile2
  character(len=FILENAMELEN) :: mapQfile, mapQfile2
  character(len=FILENAMELEN) :: mapUfile, mapUfile2

  character(len=FILENAMELEN), dimension(1:nefmax, 1:3, 1:2) :: &
       &        extramapTQUfile_input

  character(len=FILENAMELEN) :: mapfile_input,  mapfile2_input
  character(len=FILENAMELEN) :: mapQfile_input, mapQfile2_input
  character(len=FILENAMELEN) :: mapUfile_input, mapUfile2_input

  character(len=FILENAMELEN) :: maskfilep, maskfilep_input
  character(len=FILENAMELEN) :: maskfilep2, maskfilep2_input
  character(len=FILENAMELEN) :: weightfilep, weightfilep_input
  character(len=FILENAMELEN) :: weightfilep2, weightfilep2_input

  integer(i4b), dimension(1:nefmax, 3, 2) :: ordering_tqu

  integer(I4B) :: ordering_t1, ordering_t2
  integer(I4B) :: ordering_et1, ordering_et2
!  integer(I4B) :: ordering_q1, ordering_q2
!  integer(I4B) :: ordering_u1, ordering_u2
  integer(I4B) :: ordering_m1, ordering_m2
  integer(I4B) :: ordering_w1, ordering_w2
  integer(I4B) :: windowfile_nx, windowfile_ny
  integer(I4B) :: nmask1, nmask2, nmask  ! 1 or 2
  integer(I4B) :: ncmask ! 1 or 3
  integer(I4B) :: nweight1, nweight2, nweight ! 1 or 2
  integer(I4B) :: nkernels, ncov

  character(len=FILENAMELEN) :: out_cl_grp

  integer :: noption,jstartnormal
  character(len=FILENAMELEN), allocatable, dimension(:) :: option,option_param
  integer(I4B), allocatable, dimension(:) :: iskip_option
  logical, allocatable, dimension(:) :: option_done(:)

!  real(DP) :: fwhm,fwhm2,quad_uK,decouplethreshold
  real(DP) :: quad_uK,decouplethreshold

  integer(I8B) :: mypar
end module spice_common
