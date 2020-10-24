
module healpix_sharp_f90_cut
  use, intrinsic :: iso_c_binding
  implicit none
  public !:: sharp_hpcut_map2alm, sharp_hpcut_map2alm_pol
  
  interface !sharp_hpcut_map2alm
     subroutine sharp_hpcut_map2alm_x_s(nside, lmax, mmax, nda, nrings, rings, np, map, alm, wgt) &
          bind(C, name="sharps_map2alm_cut")
       use iso_c_binding
       integer(c_int32_t),  intent(in), value :: nside, lmax, nda, mmax, nrings
       integer(c_int64_t),  intent(in), value :: np
       integer(c_int32_t),  intent(in)        :: rings(*)
       real   (c_float),    intent(in)        :: map(*)
       complex(c_float),    intent(out)       :: alm(*)
       real   (c_double),   intent(in)        :: wgt(*)
     end subroutine sharp_hpcut_map2alm_x_s
     
     subroutine sharp_hpcut_map2alm_x_d(nside, lmax, mmax, nda, nrings, rings, np, map, alm, wgt) &
          bind(C, name="sharpd_map2alm_cut")
       use iso_c_binding
       integer(c_int32_t),  intent(in), value :: nside, lmax, nda, mmax, nrings
       integer(c_int64_t),  intent(in), value :: np
       integer(c_int32_t),  intent(in)        :: rings(*)
       real   (c_double),   intent(in)        :: map(*)
       complex(c_double),   intent(out)       :: alm(*)
       real   (c_double),   intent(in)        :: wgt(*)
     end subroutine sharp_hpcut_map2alm_x_d
!   end interface

!   interface sharp_hpcut_map2alm_pol
     subroutine sharp_hpcut_map2alm_pol_x_s(nside, lmax, mmax, nda, nrings, rings, np, map, alm, wgt) &
          bind(C, name="sharps_map2alm_pol_cut")
       use iso_c_binding
       integer(c_int32_t),  intent(in), value :: nside, lmax, nda, mmax, nrings
       integer(c_int64_t),  intent(in), value :: np
       integer(c_int32_t),  intent(in)        :: rings(*)
       real   (c_float),    intent(in)        :: map(*)
       complex(c_float),    intent(out)       :: alm(*)
       real   (c_double),   intent(in)        :: wgt(*)
     end subroutine sharp_hpcut_map2alm_pol_x_s
     
     subroutine sharp_hpcut_map2alm_pol_x_d(nside, lmax, mmax, nda, nrings, rings, np, map, alm, wgt) &
          bind(C, name="sharpd_map2alm_pol_cut")
       use iso_c_binding
       integer(c_int32_t),  intent(in), value :: nside, lmax, nda, mmax, nrings
       integer(c_int64_t),  intent(in), value :: np
       integer(c_int32_t),  intent(in)        :: rings(*)
       real   (c_double),   intent(in)        :: map(*)
       complex(c_double),   intent(out)       :: alm(*)
       real   (c_double),   intent(in)        :: wgt(*)
     end subroutine sharp_hpcut_map2alm_pol_x_d
  end interface
end module healpix_sharp_f90_cut




