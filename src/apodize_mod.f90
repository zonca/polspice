module apodize_mod
  public :: apodizefunction

contains
  !===========================================================
  subroutine apodizefunction(mu,fwhm_degrees,thetamax,resultat,type)
    !===========================================================
    ! This routine computes the value of the apodization function
    ! as a fontion of mu=cos(theta) and apodizesigma

    use healpix_types
    implicit none

    real(DP),     intent(IN)           :: mu, fwhm_degrees, thetamax
    real(DP),     intent(OUT)          :: resultat
    integer(i4b), intent(IN), optional :: type
    ! Local variables
    real(DP) :: argument, apodizesigma
    integer(i4b) :: typef
    typef = 0

    if (present(type)) typef = type

    ! special case for no apodization and decouple=T
    if (fwhm_degrees <= 0.0_dp) then
       resultat = 1.0_dp
       return
    endif

    if (typef == 0) then
       ! First get sigma in radians from fwhm value in degrees
       apodizesigma = fwhm_degrees*PI/180.d0/sqrt(8.d0*log(2.d0))

       if (acos(mu)*180.d0/PI < thetamax) then
          resultat = exp(-(acos(mu)/apodizesigma)**2/2.d0)
       else
          resultat = 0.d0
       endif
    endif
    if (typef == 1) then
       ! This is the half-cosine case
       apodizesigma = fwhm_degrees*PI/180.d0
       argument = max(acos(mu)*180.0/PI/thetamax , acos(mu)/apodizesigma)
       if (argument < 1.0d0) then
          resultat = ( 1.d0 + cos(argument*PI) )/2.d0
       else
          resultat = 0.d0
       endif
    endif

  end subroutine apodizefunction


end module apodize_mod
