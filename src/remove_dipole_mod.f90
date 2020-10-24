module remove_dipole_mod

  USE healpix_types
  USE misc_utils
  use pix_tools, only: nside2npix, npix2nside, pix2vec_ring, pix2vec_nest, query_strip
  IMPLICIT none

! emulate remove_dipole as contained in Healpix 2.10 (ie, with 'weights' optional argument)
!
  interface remove_dipole
     module procedure remove_dipole_real, remove_dipole_double
  end interface

contains
  !==============================================================
  subroutine remove_dipole_real( nside, map, ordering, degree, multipoles, zbounds, &
       & fmissval, mask, weights, silent)
    !============================================================
    use num_rec, only : dsvdcmp, dsvbksb
    ! parameters
    character(len=*),      parameter :: code = "REMOVE_DIPOLE_REAL"
    integer(kind=I4B),     parameter :: KMAP = SP
    real   (kind=KMAP),    parameter :: fbad_value = -1.6375e30_KMAP
    !
    include 'remove_dipole_inc_test.f90'
  end subroutine remove_dipole_real
  !==============================================================
  subroutine remove_dipole_double( nside, map, ordering, degree, multipoles, zbounds, &
       & fmissval, mask, weights, silent)
    !============================================================
    use num_rec, only : dsvdcmp, dsvbksb
    ! parameters
    character(len=*),      parameter :: code = "REMOVE_DIPOLE_DOUBLE"
    integer(kind=I4B),     parameter :: KMAP = DP
    real   (kind=KMAP),    parameter :: fbad_value = -1.6375e30_KMAP
    !
    include 'remove_dipole_inc_test.f90'
  end subroutine remove_dipole_double


end module remove_dipole_mod
