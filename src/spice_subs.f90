module spice_subs


  use misc_utils, only: assert_alloc, wall_clock_time
!  implicit none

  private
  public :: compute

contains

subroutine apply_w8_power(map, npix, nmap, weightpower, verbose, comment)
  use healpix_types
  use spice_parameters, only: KMAP
  implicit none
  integer(i4b),                          intent(in) :: npix, nmap
  real(KMAP), dimension(0:npix-1, 1:nmap), intent(inout) :: map
  real(DP),                              intent(in) ::  weightpower
  logical(LGT),                          intent(in) :: verbose
!   character(len=*), optional,            intent(in) :: comment
  character(len=*),                      intent(in) :: comment
  character(len=40) :: mycomment

  mycomment = '         '

! can not deal correctly with optional comment (ifort 10.0.016 on MacOSx)
!  if (present(comment)) mycomment = adjustl(trim(comment))
!  print*,'mycomment:',mycomment
  mycomment = adjustl(trim(comment))

  if (weightpower == 0.d0) then
     if (verbose) write(*,*) 'Put all the '//trim(mycomment)//' weights to unity (pure mask)'
     where (map > 0.) map = 1.
  elseif (weightpower > 0.d0) then
     if (weightpower /= 1.d0) then
        if (verbose) write(*,*) 'Change the '//trim(mycomment)//' weights w to w**',weightpower
        ! dont do anything if weightpower=1. to avoid unnecessary calculation
        map = map**weightpower
     endif
  elseif (weightpower < 0.d0) then
     if (verbose) write(*,*) 'Change the non zero '//trim(mycomment)//' weights w to w**',weightpower
     where (map > 0.) map = map**weightpower
  endif

  return
end subroutine apply_w8_power

!============================================================================
!subroutine compute_V_matrix(ClMap1, ClMap2, ClMask_sq, V_matrix)
subroutine compute_V_matrix(V_matrix, cl_local, cl_mask_sq, cl_auto1, cl_auto2)
  !============================================================================
  !
  ! compute V matrix
  ! V(l1,l2) = 2 C(l1) C(l2) sum_l3 (2*l3+1)/(4Pi) W_l3 J(l1, l2, l3; 0,0,0)^2 
  !
  ! 2012-07-13: parallel implementation

  use healpix_types, only: i4b, lgt, sp, dp, FOURPI
  use spice_common,  only: nlmax
  use rec3jj_module, only: rec3jj
!  use misc_utils

  implicit none

  real(dp), intent(out), dimension(0:nlmax,0:nlmax,1:1) :: V_matrix
  real(dp), intent(in),  dimension(0:nlmax, 1:1) :: Cl_local
  real(dp), intent(in),  dimension(0:nlmax, 1:1), optional :: Cl_mask_sq
  real(dp), intent(in),  dimension(0:nlmax, 1:1), optional :: Cl_auto1, Cl_auto2

  real(dp), allocatable, dimension(:) :: jmsymbol
  real(dp) :: l3min, l3max, sum, factor, dm1, dm2, dl1, dl2, tmp
  character(len=*), parameter :: code = 'compute_V_matrix (Parallel)'
  integer(i4b) :: l1, l2, l3, m1, m2, status, info
  logical(lgt) :: do_cross
  real(sp) :: wt1, wt2
!!  integer(i4b) :: i, j, ios
  


  call wall_clock_time(wt1)
  do_cross = (present(Cl_auto1) .and. present(Cl_auto2))
  V_matrix = 0.0_dp

  if (present(Cl_mask_sq)) then
     ! cut sky -> compute full matrix
     dm1 = 0.0_dp
     dm2 = 0.0_dp
!$OMP PARALLEL DEFAULT(NONE) &
!$OMP SHARED(V_matrix, cl_mask_sq, Cl_auto1, Cl_auto2, Cl_local, &
!$OMP        nlmax, dm1, dm2, do_cross) &
!$OMP PRIVATE(jmsymbol, status, factor, sum, tmp, l1, l2, l3, l3min, l3max, dl1, info)
     allocate(jmsymbol(0:2*nlmax),stat=status)
     call assert_alloc(status, code, 'can not allocate JMsymbol')

!$OMP DO schedule(dynamic, 16)
     do l1=0, nlmax
        dl1 = real(l1, kind=DP)
        do l2=l1, nlmax ! only half of matrix
           call rec3jj(jmsymbol, dl1, real(l2,kind=dp), dm1, dm2, l3min, l3max, 2*nlmax+1,info)
           sum = 0.0_dp
           do l3=nint(l3min), min(nint(l3max), nlmax)
              !sum =  sum + (2.0_dp*real(l3,kind=dp)+1.0_dp) * cl_mask_sq(l3,1) * jmsymbol(l3-nint(l3min))**2
              sum =  sum + (2.0_dp*l3 +1.0_dp) * cl_mask_sq(l3,1) * jmsymbol(l3-nint(l3min))**2
           end do
           if (do_cross) then
              tmp = Cl_auto1(l1,1)*Cl_auto1(l2,1)*Cl_auto2(l1,1)*Cl_auto2(l2,1)
              if (tmp < 0.) tmp = 0.d0
              factor = Cl_local(l1,1)*Cl_local(l2,1) + sqrt(abs(tmp))
           else
              factor = 2.0_dp * Cl_local(l1,1) * Cl_local(l2,1)
           endif
           V_matrix(l1, l2, 1) =  factor * sum / FOURPI
           if (l2 /= l1) V_matrix(l2, l1, 1) = V_matrix(l1, l2, 1) ! other half by symmetry
        end do
     end do

!$OMP END DO
     deallocate(jmsymbol)
!$OMP END PARALLEL

  else
     ! full sky -> compute diagonal of matrix only
     do l1=0, nlmax
        factor = 2.0_dp * l1 + 1.0_dp
        if (do_cross) then
           V_matrix(l1, l1, 1) = (Cl_local(l1,1)**2 + Cl_auto1(l1,1)*Cl_auto2(l1,1)) / factor
        else
           V_matrix(l1, l1, 1) = 2.0_dp * Cl_local(l1,1)**2 / factor
        endif
     enddo
  endif
  call wall_clock_time(wt2)
  print*,code,wt2-wt1,' [s]'

!   print*,'V matrix minmax ', minval(V_matrix),      maxval(V_matrix)

  return

end subroutine compute_V_matrix


!============================================================================
 subroutine Compute_Cov_matrix(V_matrix, only_non_trivial, xi_mask)
   !============================================================================
  ! 2012-07-13: parallel implementation

  use healpix_types
  use spice_common, only: mu, w, Pl
  use spice_parameters, only: apodizesigma, thetamax, apodizetype, nlmax
!  use misc_utils
  use apodize_mod, only: apodizefunction
  use deal_with_xi_and_cl, only: do_xi_from_cl, do_cl_from_xi, correct_final_cl

  implicit none

  real(dp),     intent(inout), dimension(0:nlmax,0:nlmax,1:1) :: V_matrix
  logical,      intent(in)                                    :: only_non_trivial
  real(dp),     intent(in),    dimension(0:nlmax,1:1), optional :: Xi_mask


  real(dp), allocatable, dimension(:,:) :: XiV
  real(dp), allocatable, dimension(:,:) :: Vll
  real(dp), allocatable, dimension(:)   :: apfunction
  integer(i4b) :: ll, i, status, ios, idim, ncor_local, status_private
  character(len=40) :: message=''
  character(len=*), parameter :: code = 'compute_Cov_matrix (Parallel)'
  logical :: correct_mask, verbose_local, polarization_local, decouple_local
  real(sp) :: wt1, wt2

!  print*,'Entering '//code
  call wall_clock_time(wt1)

  verbose_local=.false.
  polarization_local=.false.
  decouple_local=.false.
  ncor_local = 1
  correct_mask = present(Xi_mask)

  allocate(apfunction(0:nlmax),stat=status)
  call assert_alloc(status, code, 'can not allocate apfunction')

  apfunction = 0.0_dp
  do i=0,nlmax
     call apodizefunction(mu(i+1), apodizesigma, thetamax, apfunction(i), type=apodizetype)
  enddo
  
  do idim = 1, 2
!$OMP PARALLEL DEFAULT(NONE) &
!$OMP SHARED(V_matrix, Pl, Xi_mask, apfunction, w, &
!$OMP        nlmax, idim, message, correct_mask, only_non_trivial,  &
!$OMP        ncor_local, polarization_local, decouple_local, verbose_local) &
!$OMP PRIVATE(Vll, XiV, ll, i, status_private)
     allocate(XiV(0:nlmax,1:1),stat=status_private)
     call assert_alloc(status_private, code, 'can not allocate XiV')
     allocate(Vll(0:nlmax,1:1),stat=status_private)
     call assert_alloc(status_private, code, 'can not allocate Vll')
!$OMP DO schedule(dynamic, 16)

     do ll=0, nlmax
        ! select row (or column)
        if (idim == 1)     Vll(:,1)=V_matrix(:,ll,1)
        if (idim == 2)     Vll(:,1)=V_matrix(ll,:,1)

        ! C(l) -> Xi (theta)
        call do_xi_from_cl(XiV, Vll, Pl(:,:,0:ncor_local-1), nlmax, verbose_local, message, ncor_local, polarization_local)

        ! correct Xi(theta)
        if (correct_mask) then
           do i=0, nlmax
              if (Xi_mask(i,1) .ne. 0.0_dp) then
                 XiV(i,1)=XiV(i,1)/Xi_mask(i,1)*apfunction(i)
              else
                 XiV(i,1)=0.0_dp
              end if
           end do
        else
           XiV(:,1) = XiV(:,1) * apfunction(:)
        endif

        ! Xi (theta) -> C(l)
        Vll(:,1)=0.0_dp
        call do_cl_from_xi(XiV, Vll, w, Pl(:,:,0:ncor_local-1), nlmax, &
             & verbose_local, message, ncor_local, polarization_local, decouple_local)

        ! correct C(l)
        call correct_final_cl(Vll(:,1:1), only_non_trivial, verbose_local)

        ! put back row (or column) in matrix
        if (idim == 1) V_matrix(:,ll,1)=Vll(:,1)
        if (idim == 2) V_matrix(ll,:,1)=Vll(:,1)

     end do
!$OMP END DO
     deallocate(XiV, Vll)
!$OMP END PARALLEL
  enddo

  deallocate(apfunction)

  call wall_clock_time(wt2)
  print*,code,wt2-wt1,' [s]'

  return

end subroutine Compute_Cov_matrix

!=======================================================================
subroutine compute
!=======================================================================
  use spice_parameters !, only: I4B, SP, DP, PI, maptype, &
!        & masks_present, weights_present, masks2_present, weights2_present, &
!        & clmaskinput, cloutput, correct_beam, correct_pix, normalize, &
!        & subtract_noise_cl, subtract_noise_cor, subtract_average, subtract_dipole, &
!        & thetamax, nlmax, do_cov, map_present, map2_present, &
!        & verbose, megaverbose, symmetric_cl, clmapoutput, polarization, &
!        & coroutput, apodize, decouple, apodizetype, apodizesigma, &
!        & weightpower, weightpower2, weightpowerp, weightpowerp2
  use spice_common !, only: nsides, npixtots, nmap, ncor, ncov, ncmask, &
!        & tmap1_in, tmap2_in, tmask1_map, tmask2_map, tweight1_map, tweight2_map, &
!        & cl, cl_mask, cl_noise, cl_final, cl_map_precomp, cl_mask_precomp, &
!        & xi, xi_mask, xi_noise, xi_final, cov_matrix, &
!        & mu, w, Pl, transfer_function, wlpix1, wlpix2, gb, gb2, &
!        & nmask1, nmask2, nweight1, nweight2, obsnpixs
  use misc_utils, only: fatal_error
  use deal_with_xi_and_cl, ONLY : do_cl_with_alm, &
       & do_xi_from_cl, do_cl_from_xi, &
       & correct_final_cl, correct_final_xi, correct_xi_from_mask
  use map_type, only : alloc_map_type, dealloc_map_type, &
       mult_map_type, copy_map_type, remove_dipole_map_type
  use apodize_mod, only : apodizefunction
  use my_cumul, only : cumul
  use do_legendre_mod, only: do_legendre
  use pix_tools, only : convert_nest2ring
  implicit none

  integer(I4B) :: i,j,l, nm_eff1, nm_eff2, nm_eff, k1, k2, np, npl, on1, on2
  logical :: mask_correction,nontrivial_corrections,trivial_corrections
  logical :: do_mask, do_mask2
  logical :: subtract_noise,warning, only_non_trivial
!   real(SP),pointer,dimension(:) :: m_map, m_map2
!   real(KMAP), allocatable, dimension(:,:) :: m_map, m_map2
!   real(KMAP), allocatable, dimension(:,:) :: m_map_sq, m_map2_sq
  type(maptype) :: tm_map1,    tm_map2
  type(maptype) :: tm_map1_sq, tm_map2_sq
  real(DP) :: average,sumweights, product, tmp_prod
!  real(DP), dimension(1:3) :: top_w8
  real(DP) :: maskthreshold, mumin, fact1
  !character(len=80), dimension(1:100) :: header
  integer(I4B) :: lb, lb2, kmask, kmask2
  character(len=40) :: message
!!  real(DP), allocatable, dimension(:) :: sintheta
!  logical(LGT) :: symmetric = .true.
  character(len=8) :: coordsys = '' ! place holder

  real(DP) :: tempo
  real(DP), dimension(4) :: lowmultipoles
  real(DP), dimension(2) :: zbounds = (/ -1.d0, 1.d0 /)
  real(DP), dimension(:,:), allocatable :: cl_auto1, cl_auto2, cl_mask_sq
  real(DP), dimension(:,:), allocatable :: xi_a1, xi_a2, xi_f1, xi_f2, cl_for_cov
  integer(i4b) :: degree_remove = 0
  character(len=80), dimension(1:60) :: header
  real(SP) :: ct0,ct1,wt0, wt1
  real(KMAP), parameter :: ONE = 1.0_KMAP

  do_mask  = masks_present .or.weights_present
  do_mask2 = masks2_present.or.weights2_present
  mask_correction=do_mask.or.clmaskinput.or.do_mask2

  nontrivial_corrections=correct_beam.or.correct_pix
  trivial_corrections=normalize
  subtract_noise=subtract_noise_cl.or.subtract_noise_cor

  if (subtract_average) degree_remove = 1
  if (subtract_dipole)  degree_remove = 2
  
  ! Allocate files
  allocate(xi(0:nlmax,1:ncor),xi_final(0:nlmax,1:ncor))
  allocate(cl(0:nlmax,1:ncor),cl_final(0:nlmax,1:ncor))
  if (mask_correction) then
     allocate(xi_mask(0:nlmax,1:ncmask))
     allocate(cl_mask(0:nlmax,1:ncmask))
     if (do_cov) then
        allocate(cl_mask_sq(0:nlmax,1:ncmask))
        cl_mask_sq = 0.0_dp
     endif
  endif
  if (map2_present) then
     allocate(cl_auto1(0:nlmax,1:ncor))
     allocate(cl_auto2(0:nlmax,1:ncor))
  endif

!!!!EH  if (correct_beam.and..not.beam_present) allocate(gb(0:nlmax,1))


  ! Compute lower limit of integration in cos(theta)
  ! Decouple threshold is now the maximum value in degrees of the integration 
  mumin = cos(thetamax *PI/180.d0) 

  !=====================================================================
  ! CASE WITH NO MASK
  !=====================================================================
  if (.not.mask_correction) then
     only_non_trivial = .false. ! flag for final C(l) correction (does gain correction)

     ! start analysis from map
     if (map_present) then
        ! reordering and combination done in input_files
        ! Subtract monopole or monopole & dipole
        if (degree_remove == 1 .or. degree_remove == 2) then
           if (megaverbose .and. degree_remove == 1) write(*,*) 'Subtract monopole from map'
           if (megaverbose .and. degree_remove == 2) write(*,*) 'Subtract monopole and dipole from map'
           do j=1, 1 ! nmap  2014-03-27: turn-off monopole and dipole removal for Q,U maps
!               call remove_dipole(nsides(1), map_in(:,j), 1, degree_remove, lowmultipoles, zbounds, silent=.not.megaverbose)
              call remove_dipole_map_type(tmap1_in, j, degree_remove, lowmultipoles, silent=.not.megaverbose)
           enddo
        endif

        ! cross correlate 2 maps
        if (map2_present) then
           ! Subtract monopole and/or dipole
           if (degree_remove == 1 .or. degree_remove == 2) then
              if (megaverbose .and. degree_remove == 1) write(*,*) 'Subtract monopole from map 2'
              if (megaverbose .and. degree_remove == 2)  write(*,*) 'Subtract monopole and dipole from map 2'
              do j=1, 1 ! nmap  2014-03-27: turn-off monopole and dipole removal for Q,U maps
!                  call remove_dipole(nsides(2), map2_in(:,j), 1, degree_remove, lowmultipoles, zbounds, silent=.not.megaverbose)
                 call remove_dipole_map_type(tmap2_in, j, degree_remove, lowmultipoles, silent=.not.megaverbose)
              enddo
           endif
           ! Compute the cross-Cls of the maps
           call do_cl_with_alm(nsides,npixtots,nlmax,tmap1_in,cl_final, &
                &                        megaverbose,'(mapXmap)',nmap,ncor,polarization, &
                &                        tmap2=tmap2_in, cl_auto1=cl_auto1, cl_auto2=cl_auto2, &
                &                        symmetric=symmetric_cl, &
                &                        dump_alm1=alm1_out_file, dump_alm2=alm2_out_file, &
                &                        edit_alm1='',            edit_alm2='')
        else
           ! Compute the auto Cls of the map
           call do_cl_with_alm(nsides,npixtots,nlmax,tmap1_in,cl_final, &
                &                        megaverbose,'(map)',nmap,ncor,polarization, &
                &                        dump_alm1=alm1_out_file, &
                &                        edit_alm1='')
        endif

        if (clmapoutput) cl(0:nlmax,1:ncor)=cl_final(0:nlmax,1:ncor)
     else ! start directly from precomputed power spectrum
        cl_final(0:nlmax,1:ncor) = cl_map_precomp(0:nlmax,1:ncor)
     endif

     ! deallocate global variables
     ! if (allocated(map_in))  deallocate(map_in)
     ! if (allocated(map2_in)) deallocate(map2_in)

     if (do_cov) then
        ncov = 1
        allocate(mu(1:nlmax+1),w(1:nlmax+1),Pl(0:nlmax,0:nlmax,0:ncov-1)) ! 2020-05-19
        call do_legendre(nlmax, mumin, mu, w, Pl, megaverbose, ncov) ! 2020-05-19, compute mu, w and Pl needed by cov
        allocate(cov_matrix(0:nlmax,0:nlmax,1:ncov))
        if (map2_present) then
           ! 2 maps: provide cross-spectrum and 2 auto-spectra
           call compute_V_matrix(cov_matrix, cl_final(:,1:1), cl_auto1=cl_auto1(:,1:1), cl_auto2=cl_auto2(:,1:1))
        else                                                                                     
           ! single map: provide (auto-)spectrum
           call compute_V_matrix(cov_matrix, cl_final(:,1:1))
        endif
        call compute_Cov_matrix(cov_matrix, only_non_trivial)
        deallocate(mu, w, Pl) ! 2020-05-19
     endif

     ! Apply various corrections
     call correct_final_cl(cl_final, only_non_trivial, megaverbose)

     ! Compute legendre polynomials/roots (sometimes unnecessary but simpler
     ! that way and anyway negligible in CPU)
     ! The roots are rescaled from mumin to 1.d0
     npl = min(ncor, 5)
     allocate(mu(1:nlmax+1),w(1:nlmax+1),Pl(0:nlmax,0:nlmax,0:npl-1))
     call do_legendre(nlmax, mumin, mu, w, Pl, megaverbose, npl) 

     ! Compute xi if required
     if (coroutput.or.apodize.or.decouple) then    
        ! Compute xi from the Cls of the map
        call do_xi_from_cl(xi_final,cl_final,Pl,nlmax,megaverbose, &
&                          '(map)',ncor,polarization)
     endif

     if (decouple) then
        call cumul(nlmax,mu,xi_final,cl_final,ncor,megaverbose, &
&                  '(map)',thetamax,mask_correction)
     endif

     ! Apodize the correlation function !
     if (apodize) then
        if (megaverbose) write(*,*) 'Apodize xi'
        do l=0,nlmax
           call apodizefunction(mu(l+1),apodizesigma,thetamax,tempo,type=apodizetype)
           xi_final(l,:) = xi_final(l,:) * tempo
        enddo
        if (cloutput) then
           call do_cl_from_xi(xi_final,cl_final,w,Pl,nlmax,megaverbose, &
 &                            '(final)',ncor,polarization,decouple,mu=mu, &
 &                             apodizesigma=apodizesigma,thetamax=thetamax)
        endif
     endif
     
  !=====================================================================
  ! CASE WITH MASKS
  !=====================================================================
  elseif (mask_correction) then
     call cpu_time(ct0)
     call wall_clock_time(wt0)
     only_non_trivial = .true. ! flag for final C(l) correction (does not do gain correction)

     nm_eff1 = max(nmask1,nweight1) ! either 1 or 2
     nm_eff2 = max(nmask2,nweight2) ! either 1 or 2
     nm_eff  = max(nm_eff1, nm_eff2) ! make sure the 2 mask*weight have same size for cross-correlation
     on1     = max(obsnpixs_m(1), obsnpixs_w(1)) ! 2020-04-08
     ! ----- merge mask and weight ---------
     if (megaverbose) write(*,*) 'Create the masked/weighted map'
     call alloc_map_type(tm_map1, nsides(1), nm_eff, npixtots(1), on1, 1, coordsys, value=ONE) ! 2019-10
     !tm_map1%map = 1.0 ! very important
     if (masks_present) then
        ! conversion to RING of tmask1_map (if required) already done in input_files
        ! tm_map1 = tmask1_map
        call mult_map_type(tm_map1, tmask1_map)
        call dealloc_map_type(tmask1_map)
     endif
     if (weights_present) then
        ! conversion to RING of tweight1_map (if required) already done in input_files
        ! tm_map1 *= tweight1_map
        call mult_map_type(tm_map1, tweight1_map)
        call dealloc_map_type(tweight1_map)
     endif

     if (map2_present) then ! second map -> second mask required
        on2     = max(obsnpixs_m(2), obsnpixs_w(2)) ! 2020-04-08
        call alloc_map_type(tm_map2, nsides(2), nm_eff, npixtots(2), on2, 1, coordsys, value=ONE) ! 2019-10
        !tm_map2%map = 1.0 ! very important
        if (do_mask2) then ! either mask2 or weight2 provided by user
           if (megaverbose) write(*,*) 'Create the masked/weighted map (2)'
           if (masks2_present) then
              ! conversion to RING of tmask2_map (if required) already done in input_files
              ! tm_map2 = tmask2_map
              call mult_map_type(tm_map2, tmask2_map)
              call dealloc_map_type(tmask2_map)
           endif
           if (weights2_present) then
              ! conversion to RING of tweight2_map (if required) already done in input_files
              ! tm_map2 *= tweight2_map
              call mult_map_type(tm_map2, tweight2_map)
              call dealloc_map_type(tweight2_map)
           endif
        else
           if (megaverbose) write(*,*) 'No mask nor weights for second map'
        endif
     endif


     ! ----- deal with weightpower ---------
     if (weights_present) then
        np = int(tm_map1%npix, I4B)
        if (nm_eff > 1 .and. weightpower /= weightpowerp) then
           call apply_w8_power(tm_map1%map(:,1:1), np, 1, weightpower,  megaverbose, '(temp. 1)')
           call apply_w8_power(tm_map1%map(:,2:2), np, 1, weightpowerp, megaverbose, '(pol. 1)')
        else
           call apply_w8_power(tm_map1%map, np, nm_eff, weightpower, megaverbose, '(map 1)')
        endif
     endif
     if (map2_present .and. weights2_present) then
        np = int(tm_map2%npix, I4B)
        if (nm_eff > 1 .and. weightpower2 /= weightpowerp2) then
           call apply_w8_power(tm_map2%map(:,1:1), np, 1, weightpower2,  megaverbose, '(temp. 2)')
           call apply_w8_power(tm_map2%map(:,2:2), np, 1, weightpowerp2, megaverbose, '(pol. 2)')
        else
           call apply_w8_power(tm_map2%map, np, nm_eff, weightpower2, megaverbose, '(map 2)')
        endif
     endif
        
     ! --- check that the 2 mask*weights are not totally disjoint ---
     ! --- and find out largest weight (required for pair thresholding)
     if (map2_present .and. do_mask2 .and. nsides(1)==nsides(2)  .and. tm_map1%type == tm_map1%FULL) then
        if (megaverbose) write(*,*) 'Checking consistency of the 2 masks and/or weights'
        lb  = lbound(tm_map1%map,1)
        lb2 = lbound(tm_map2%map,1)
        do kmask = 1, ncmask ! ncmask = 1 or 3 or 4
           k1 = kmask ; k2 = kmask
           if (kmask == 3) then
              k1 = 1
              k2 = 2
           endif
           if (kmask == 4) then
              k1 = 2
              k2 = 1
           endif
           product = 0.0_dp
           do i=0,npixtots(1)-1
              tmp_prod = tm_map1%map(i+lb,k1) * tm_map2%map(i+lb2,k2)
              product   =  product + tmp_prod
!              top_w8(kmask) =  max(top_w8(kmask),  tmp_prod)
           enddo
           if (product < 1.d-50) then
              if (kmask == 1) write(*,*) 'ERROR in spice : the 2 masks (TxT) appear to be disjoint'
              if (kmask == 2) write(*,*) 'ERROR in spice : the 2 masks (PxP) appear to be disjoint'
              if (kmask == 3) write(*,*) 'ERROR in spice : the 2 masks (TxP) appear to be disjoint'
              call fatal_error
           endif
        enddo
     endif
     ! ---  Subtract monopole and/or dipole -----
     if (map_present) then
        ! subtract monopole and/or dipole
        if (degree_remove == 1 .or. degree_remove == 2) then
           if (megaverbose .and. degree_remove == 1) write(*,*) 'Subtract monopole from map'
           if (megaverbose .and. degree_remove == 2) write(*,*) 'Subtract monopole and dipole from map'
!            call remove_dipole(nsides(1), map_in(:,1), 1, degree_remove, lowmultipoles, zbounds, &
!                 & mask = m_map(:,1), weights = m_map(:,1), silent=.not.megaverbose)
           call remove_dipole_map_type(tmap1_in, 1, degree_remove, lowmultipoles, &
                & mask = tm_map1, weights = tm_map1, silent=.not.megaverbose)
        endif
        ! apply mask*weight
        if (do_mask) then
           if (megaverbose) write(*,*) 'Multiply the map1 by the mask/weight map'
           call mult_map_type(tmap1_in, tm_map1)
        endif
     endif
     if (map2_present) then
        ! subtract monopole and/or dipole
        if (degree_remove == 1 .or. degree_remove == 2) then
           if (megaverbose .and. degree_remove == 1) write(*,*) 'Subtract monopole from map2'
           if (megaverbose .and. degree_remove == 2) write(*,*) 'Subtract monopole and dipole from map2'
           if (do_mask2) then 
              call remove_dipole_map_type(tmap2_in, 1, degree_remove, lowmultipoles, &
                   & mask = tm_map2, weights = tm_map2, silent=.not.megaverbose)
           else
              call remove_dipole_map_type(tmap2_in, 1, degree_remove, lowmultipoles, &
                   & silent=.not.megaverbose)
           endif
        endif
        ! apply mask*weight
        if (do_mask2) then
           if (megaverbose) write(*,*) 'Multiply the map2 by the mask/weight map'
           call mult_map_type(tmap2_in, tm_map2)
        endif
     endif
     call cpu_time(ct1)
     call wall_clock_time(wt1)
     if (verbose) then
        print*,'Mask & map (CPU,wall) time [s]: ',ct1-ct0,wt1-wt0
     endif

     ! --- C(l) calculations ----
     call cpu_time(ct0)
     call wall_clock_time(wt0)
     if (map_present) then
        ! Compute the Cls of the map
        if (map2_present) then ! cross-correlate 2 maps
           call do_cl_with_alm(nsides,npixtots,nlmax,tmap1_in, cl, &
                &              megaverbose,'(mapXmap)',nmap,ncor,polarization, &
                &              tmap2=tmap2_in, cl_auto1=cl_auto1, cl_auto2=cl_auto2, &
                &              symmetric=symmetric_cl, &
                &              dump_alm1=alm1_out_file, dump_alm2=alm2_out_file, &
                &              edit_alm1='',            edit_alm2='')
        else ! auto-correlate 1 map
           call do_cl_with_alm(nsides,npixtots,nlmax,tmap1_in, cl, &
                &              megaverbose,'(map)',nmap,ncor,polarization, &
                &              dump_alm1=alm1_out_file, &
                &              edit_alm1='')
        endif
     else
        ! Start directly from precomputed power-spectra
        cl(0:nlmax,1:ncor) = cl_map_precomp(0:nlmax,1:ncor)
     endif

     ! deallocate global variables
     call dealloc_map_type(tmap1_in)
     call dealloc_map_type(tmap2_in)


     ! Compute the Cls of the mask/weight map
     if (clmaskinput) then
        ! Start directly from precomputed power-spectra
        cl_mask(0:nlmax,1:ncmask) = cl_mask_precomp(0:nlmax,1:ncmask)           
     else
        ! Compute the Cls of the weight map or pure mask map
        if (map2_present) then ! 2 maps -> assume the 2 masks to be different
           message = '(masksXmasks)'
           if (weights_present .or. weights2_present) message = '(weightsXweights)'
           call do_cl_with_alm(nsides, npixtots, nlmax, tm_map1, cl_mask, &
                &              megaverbose,message,nm_eff,ncmask,.false.,&
                &              tmap2=tm_map2,symmetric=symmetric_cl)
        else ! auto-correlate 1 mask
           message = '(masks)'
           if (weights_present) message = '(weights)'
           call do_cl_with_alm(nsides, npixtots, nlmax, tm_map1, cl_mask, &
                &              megaverbose,message,nm_eff,ncmask,.false.)
        endif
        if (do_cov) then
           if (weights_present) then
              ! compute the C(l) of weight^2 for covariance (only for non-binary weight)
              ! allocate(m_map_sq (0:npixtots(1)-1,1:nm_eff))
              ! m_map_sq  = m_map**2
              !call alloc_map_type(tm_map1_sq, nsides(1), nm_eff, npixtots(1), tm_map1%npix, 1, coordsys) ! 2019-10
              !tm_map1_sq%map = (tm_map1%map)**2
              call copy_map_type(tm_map1_sq, tm_map1)
              call mult_map_type(tm_map1_sq, tm_map1_sq)
              if (map2_present) then ! 2 maps -> assume the 2 masks to be different
                 ! allocate(m_map2_sq(0:npixtots(2)-1,1:nm_eff))
                 ! m_map2_sq = m_map2**2
                 !call alloc_map_type(tm_map2_sq, nsides(2), nm_eff, npixtots(2), tm_map2%npix, 1, coordsys) ! 2019-10
                 !tm_map2_sq%map = (tm_map2%map)**2
                 call copy_map_type(tm_map2_sq, tm_map2)
                 call mult_map_type(tm_map2_sq, tm_map2_sq)
                 message = '(weights^2Xweights^2)'
                 call do_cl_with_alm(nsides, npixtots, nlmax, tm_map1_sq, cl_mask_sq, &
                      &              megaverbose,message,nm_eff,ncmask,.false.,&
                      &              tmap2=tm_map2_sq,symmetric=symmetric_cl)
                 !deallocate(m_map2_sq)
                 call dealloc_map_type(tm_map2_sq)
              else
                 message = '(weights^2)'
                 call do_cl_with_alm(nsides, npixtots, nlmax, tm_map1_sq, cl_mask_sq, &
                      &              megaverbose,message,nm_eff,ncmask,.false.)
              endif
              !deallocate(m_map_sq)
              call dealloc_map_type(tm_map1_sq)
           else
              cl_mask_sq(:,1:ncmask) = cl_mask(:,1:ncmask)
           endif
        endif
     endif ! clmaskinput
!      print*,'cl_mask  ',cl_mask(0:9,1)
!      print*,'cl_mask2 ',cl_mask_sq(0:9,1)
     call cpu_time(ct1)
     call wall_clock_time(wt1)
     if (verbose .and. (map_present .or. .not. clmaskinput)) then
        print*,'Map to C(l) (CPU,wall) time [s]: ',ct1-ct0,wt1-wt0
     endif

!      deallocate(m_map)
     call dealloc_map_type(tm_map1)
     call dealloc_map_type(tm_map2)

     ! Compute legendre polynomials/roots 
     ! The roots are rescaled from mumin to 1.d0
     npl = min(ncor, 5)
     allocate(mu(1:nlmax+1),w(1:nlmax+1),Pl(0:nlmax,0:nlmax,0:npl-1))
     call do_legendre(nlmax, mumin, mu, w, Pl, megaverbose, npl) 

     ! Compute xi of the masks from the Cls of the masks
     do kmask = 1, ncmask
        message = '(masks)'
        if (ncmask > 1) message = '(mask TxT)'
        if (kmask == 2) message = '(mask PxP)'
        if (kmask == 3) message = '(mask TxP)'
        if (kmask == 4) message = '(mask PxT)'
        call do_xi_from_cl(xi_mask(:,kmask),cl_mask(:,kmask), &
             &             Pl(:,:,0:0),nlmax, megaverbose,message,1,.false.)
     enddo

     ! Compute xi from the Cls of the map (cl -> xi)
     call do_xi_from_cl(xi, cl, Pl, nlmax, megaverbose, '(map)', ncor, polarization)
     ! Correct for the masks/weights effects on xi (xi -> xi_final)
     call correct_xi_from_mask(nlmax, ncor, ncmask, xi, xi_mask, xi_final, verbose)

     if (map2_present .and. do_cov) then
        allocate(xi_a1(0:nlmax, 1:ncor))
        allocate(xi_a2(0:nlmax, 1:ncor))
        allocate(xi_f1(0:nlmax, 1:ncor))
        allocate(xi_f2(0:nlmax, 1:ncor))
        call do_xi_from_cl(xi_a1, cl_auto1, Pl, nlmax, megaverbose, '(map)', ncor, polarization)
        call do_xi_from_cl(xi_a2, cl_auto2, Pl, nlmax, megaverbose, '(map)', ncor, polarization)

        call correct_xi_from_mask(nlmax, ncor, ncmask, xi_a1, xi_mask, xi_f1, verbose)
        call correct_xi_from_mask(nlmax, ncor, ncmask, xi_a2, xi_mask, xi_f2, verbose)
        deallocate(xi_a1, xi_a2)
     endif


     if (decouple) then
        ! C(l) -> xi_final (for E and B only) (cl -> xi_final)
        ! this requires cl_mask
        call cumul(nlmax, mu, xi_final, cl, ncor, megaverbose,'(final)', thetamax, mask_correction)
!
        !   TODO : DEAL WITH COVARIANCE OF CROSS-SPECTRA AND POLARIZED MAPS
!
!        print*,'cumul done'
     endif

     ! Apodize the correlation function !
     if (apodize) then
        ! xi_final -> xi_final
        if (megaverbose) write(*,*) 'Apodize final xi'
        do l=0,nlmax
           call apodizefunction(mu(l+1), apodizesigma, thetamax, tempo, type=apodizetype)
           xi_final(l,:) = xi_final(l,:) * tempo
           if (map2_present .and. do_cov) then
              xi_f1(l,:) = xi_f1(l,:) * tempo
              xi_f2(l,:) = xi_f2(l,:) * tempo
           endif
        enddo
     endif

     ! Apply trivial corrections
     if (trivial_corrections) then
        call correct_final_xi(xi_final)
        if (map2_present .and. do_cov) then
           call correct_final_xi(xi_f1)
           call correct_final_xi(xi_f2)
        endif
     endif
     
     ! Compute Cls if required (xi_final -> cl)
     if (cloutput .or. nontrivial_corrections .or. do_cov) then
        call do_cl_from_xi(xi_final, cl_final, w, Pl, nlmax, megaverbose, &
 &                         '(final)', ncor, polarization, decouple, mu=mu, &
 &                         apodizesigma=apodizesigma, thetamax=thetamax)
        if (map2_present .and. do_cov) then
           call do_cl_from_xi(xi_f1, cl_auto1, w, Pl, nlmax, megaverbose, &
                &                         '(final)', ncor, polarization, decouple, mu=mu, &
                &                         apodizesigma=apodizesigma, thetamax=thetamax)
           call do_cl_from_xi(xi_f2, cl_auto2, w, Pl, nlmax, megaverbose, &
                &                         '(final)', ncor, polarization, decouple, mu=mu, &
                &                         apodizesigma=apodizesigma, thetamax=thetamax)  
           deallocate(xi_f1, xi_f2)
        endif

        ! the C(l) used for the covariance are not corrected from the beam, pixel, ...
        if (do_cov) then
           allocate(cl_for_cov(0:nlmax, 1:ncor))
           cl_for_cov = cl_final
        endif

        ! Apply nontrivial corrections (beam, pixel, transfer function, ...)
        if (nontrivial_corrections) then
           call correct_final_cl(cl_final, only_non_trivial, megaverbose)

           ! Recompute xi_final if required
           if (coroutput) then 
              call do_xi_from_cl(xi_final,cl_final,Pl,nlmax, &
 &                               megaverbose,'(final result)',ncor,polarization) 
           endif
        endif
     endif


     ! ---------------------------
     ! deal with covariance matrix
     ! ---------------------------
     if (do_cov) then
        ncov = 1
        allocate(cov_matrix(0:nlmax,0:nlmax,1:ncov))
        ! compute pseudo-C(l) covariance
        if (map2_present) then
           ! 2 maps: provide (final) cross-spectrum and 2 auto-spectra
           call compute_V_matrix(cov_matrix, cl_for_cov(:,1:1), cl_mask_sq=cl_mask_sq(:,1:1), &
                cl_auto1=cl_auto1(:,1:1), cl_auto2=cl_auto2(:,1:1))
        else                                                                                     
           ! single map: provide (auto-)spectrum
           call compute_V_matrix(cov_matrix, cl_for_cov(:,1:1), cl_mask_sq=cl_mask_sq(:,1:1))
        endif
        ! correct pseudo-C(l) covariance from mask, beam and pixel to get final C(l) covariance
        call compute_Cov_matrix(cov_matrix, only_non_trivial, xi_mask=xi_mask(:,1:1))
     endif
!     print*,'finished with Cov matrix'


  endif ! end of test on mask_correction
  !=====================================================================
  ! SUBTRACT NOISE
  !=====================================================================
  if (subtract_noise) then
     if (cloutput) then
        if (megaverbose) write(*,*) 'Subtract the noise from final Cls'
        if (.not.subtract_noise_cl) then
           if (megaverbose) write(*,*) 'No Cl noise file provided : use xi noise file'
           call do_cl_from_xi(xi_noise,cl_noise,w,Pl,nlmax,megaverbose, &
 &                            '(noise)',ncor,polarization,decouple, mu=mu, &
 &                            apodizesigma=apodizesigma,thetamax=thetamax)
        endif
        cl_final=cl_final-cl_noise
     endif
     if (coroutput) then
        if (megaverbose) write(*,*) 'Subtract the noise from final xi'
        if (.not.subtract_noise_cor) then
           if (megaverbose) write(*,*) 'No xi noise file provided : use Cl noise file'
           call do_xi_from_cl(xi_noise,cl_noise,Pl,nlmax,megaverbose, &
 &                            '(noise)',ncor,polarization)
        endif
        xi_final=xi_final-xi_noise
     endif
  endif
  ! deallocate global variables
  if (allocated(xi_noise)) deallocate(xi_noise)
  if (allocated(cl_noise)) deallocate(cl_noise)
  if (allocated(xi_mask))  deallocate(xi_mask)
  if (allocated(xi))       deallocate(xi)
  if (allocated(cl_map_precomp))  deallocate(cl_map_precomp)
  if (allocated(cl_mask_precomp)) deallocate(cl_mask_precomp)
  if (allocated(transfer_function)) deallocate(transfer_function)
  if (allocated(wlpix1)) deallocate(wlpix1)
  if (allocated(wlpix2)) deallocate(wlpix2)
  if (allocated(gb))  deallocate(gb)
  if (allocated(gb2)) deallocate(gb2)

  ! deallocate local variables
  if (allocated(cl_auto1)) deallocate(cl_auto1)
  if (allocated(cl_auto2)) deallocate(cl_auto2)
  if (allocated(cl_for_cov)) deallocate(cl_for_cov)
  if (allocated(cl_mask_sq)) deallocate(cl_mask_sq)

end subroutine compute

end module spice_subs

