!
! gcc -I. -I${HEALPIX}/include sharp_interface.c -c -o tsharp.o ; ifort -I${HEALPIX}/include_ifort test_sharp.f90 -L${HEALPIX}/lib -lsharp -L${HEALPIX}/lib_ifort -lhealpix -fopenmp -o test_sharp tsharp.o
!


module healpix_sharp_f90_cut
  use, intrinsic :: iso_c_binding
  implicit none
  public !:: sharp_hpcut_map2alm, sharp_hpcut_map2alm_pol

!   interface shap_hpcut_map2alm
!      subroutine sharp_hpcut_map2alm(nside, lmax, mmax, nrings, rings, np, map, alm, wgt)
!        integer(i4b), intent(in) :: nside, lmax, mmax, nrings
!        integer(i8b), intent(in) :: np
!        integer(i4b), intent(in), dimension(:)  :: rings
!        real(sp),     intent(in), dimension(:)  :: map
!        complex(spc), intent(out), dimension(:) :: alm
  
  interface !sharp_hpcut_map2alm
     subroutine sharp_hpcut_map2alm_x_s(nside, lmax, mmax, nrings, rings, np, map, alm, wgt) &
          bind(C, name="sharps_map2alm_cut")
       use iso_c_binding
       integer(c_int32_t),  intent(in), value :: nside, lmax, mmax, nrings
       integer(c_int64_t),  intent(in), value :: np
       integer(c_int32_t),  intent(in)        :: rings(*)
       real   (c_float),    intent(in)        :: map(*)
       complex(c_float),    intent(out)       :: alm(*)
       real   (c_double),   intent(in)        :: wgt(*)
     end subroutine sharp_hpcut_map2alm_x_s
     
     subroutine sharp_hpcut_map2alm_x_d(nside, lmax, mmax, nrings, rings, np, map, alm, wgt) &
          bind(C, name="sharpd_map2alm_cut")
       use iso_c_binding
       integer(c_int32_t),  intent(in), value :: nside, lmax, mmax, nrings
       integer(c_int64_t),  intent(in), value :: np
       integer(c_int32_t),  intent(in)        :: rings(*)
       real   (c_double),   intent(in)        :: map(*)
       complex(c_double),   intent(out)       :: alm(*)
       real   (c_double),   intent(in)        :: wgt(*)
     end subroutine sharp_hpcut_map2alm_x_d
!   end interface

!   interface sharp_hpcut_map2alm_pol
     subroutine sharp_hpcut_map2alm_pol_x_s(nside, lmax, mmax, nrings, rings, np, map, alm, wgt) &
          bind(C, name="sharps_map2alm_pol_cut")
       use iso_c_binding
       integer(c_int32_t),  intent(in), value :: nside, lmax, mmax, nrings
       integer(c_int64_t),  intent(in), value :: np
       integer(c_int32_t),  intent(in)        :: rings(*)
       real   (c_float),    intent(in)        :: map(*)
       complex(c_float),    intent(out)       :: alm(*)
       real   (c_double),   intent(in)        :: wgt(*)
     end subroutine sharp_hpcut_map2alm_pol_x_s
     
     subroutine sharp_hpcut_map2alm_pol_x_d(nside, lmax, mmax, nrings, rings, np, map, alm, wgt) &
          bind(C, name="sharpd_map2alm_pol_cut")
       use iso_c_binding
       integer(c_int32_t),  intent(in), value :: nside, lmax, mmax, nrings
       integer(c_int64_t),  intent(in), value :: np
       integer(c_int32_t),  intent(in)        :: rings(*)
       real   (c_double),   intent(in)        :: map(*)
       complex(c_double),   intent(out)       :: alm(*)
       real   (c_double),   intent(in)        :: wgt(*)
     end subroutine sharp_hpcut_map2alm_pol_x_d
  end interface

end module healpix_sharp_f90_cut




program test
  use healpix_types
  use healpix_modules
  use healpix_sharp_f90_cut
  implicit none

  integer(i4b) :: nside, lmax, mmax, nrings
  integer(i4b) :: i, l, m
  integer(i8b) :: np
  integer(i4b), allocatable, dimension(:)     :: rings
  real(sp),     allocatable, dimension(:,:)   :: map
  complex(spc), allocatable, dimension(:,:,:) :: alm
  real(dp),     allocatable, dimension(:)     :: w

  

  nside = 4
  lmax = 2*nside
  mmax = lmax
  np   = 6 * nside * nside
  nrings = 2*nside

  allocate(w(1:4*nside-1)) ! weights for full sky map
  w(:) = 1.0_dp

  allocate(map(0:np-1,1:1)) ! partial map
  map(:,:) = 1.0

  allocate(rings(1:nrings)) ! list of rings covered by partial map
  do i=1,nrings
     rings(i) = i
  enddo

  allocate(alm(1:1, 0:lmax, 0:mmax))

  call sharp_hpcut_map2alm_x_s(nside, lmax, mmax, nrings, rings, np, map, alm, w)
  !call sharp_hpcut_map2alm_pol_x_s(nside, lmax, mmax, nrings, rings, np, map, alm, w)
  !call sharp_hpcut_map2alm(nside, lmax, mmax, nrings, rings, np, map, alm, w)
  do l=0,4
     do m=0,l
        print*,l,m,alm(1,l,m)
     enddo
     print*
  enddo


end program test
