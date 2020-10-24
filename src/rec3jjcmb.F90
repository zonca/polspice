module rec3jjcmb_module
  !-------------------------------
  ! module for computation of Wigner 3J symbols J(l1,l2,l2;m1,m2,m3)
  ! needed for CMB Temperature + Polarization power spectra analysis
  ! ie (m1,m2,m3) = (0,0,0) and (-2,2,0) exclusively.
  !
  ! Adapted from a code by F. Elsner used for power spectra covariance matrix calculation
  !
  ! It is faster (by 30 to 50%) than the more general computation
  ! based on SLATEC routines available in rec3jj.f90
  !
  ! 2017-08-24 / IAP / Eric Hivon
  !-------------------------------
  

  use healpix_types
  use healpix_modules
  implicit none

  private
  public :: w3j_220, w3j_000, init_w3j

  integer(I4B), parameter         :: ncof = 6, g0 = 5
  real(DP),     dimension(0:ncof) :: cof
  real(DP)                        :: stp 
  logical(LGT)                    :: init_done = .false.
  real(DP),     dimension(:), allocatable :: ramp
  integer(I4B) :: l_top
  
  interface gammln
     module procedure gammln_s, gammln_v
  end interface

contains

  subroutine init_w3j(lmax)
    integer(I4B), optional, intent(in) :: lmax
    
    if (.not. init_done) then

       ! Lanczos coefficient of LogGamma 
       ! (see NumRec and https://mrob.com/pub/ries/lanczos-gamma.html#code)
       cof(0:ncof) =  (/ & 
            1.000000000190015d0, &
            76.18009172947146d0, -86.50532032941677d0,  &
            24.01409824083091d0, -1.231739572450155d0,  &
            .1208650973866179d-2, -.5395239384953d-5 /)
       stp = 2.5066282746310005d0 ! sqrt(2 Pi)
       
       l_top = 2**20 ! 1.e6
       if (present(lmax)) then
          if (lmax > 0) l_top = lmax + 2 ! add 2 to avoid high-ell boundary effect in EXTEND_W3J_mm0
       endif
       !print*,'l_top = ',l_top

       init_done = .true.
    endif
    return

  end subroutine init_w3j

!-----------------------------------------------------------------------
! Compute Wigner 3j (l1,l2,l3; m1,m2,m3) 
!  where (m1,m2,m3) = (0,0,0), (-1,1,0) and (-2,2,0)
!-----------------------------------------------------------------------

  SUBROUTINE W3J_220(w3jarr, l1, l2, l3_min, l3_max)

    INTEGER(I4B),                           INTENT(IN)    :: l1
    INTEGER(I4B),                           INTENT(IN)    :: l2
    real(DP),             dimension(0:,0:), intent(out)   :: w3jarr
    integer(I4B), intent(out) :: l3_min, l3_max

    REAL(DP)                       :: l1_dp, l2_dp, l3_dp
    REAL(DP)                       :: p12, s12, f12, t12, l1_sq, l2_sq
    integer(I4B) :: l3


    ! J(l1,l2,l3;0,0,0)
    CALL W3J_000(w3jarr(0:l1+l2,0), l1, l2, l3_min, l3_max)

    if ((l1 < 2) .or. (l2 < 2)) then
       w3jarr(:,1:2)= 0.0_dp
       return
    endif

    ! J(l1,l2,l3;-1,1,0), even l1+l2+l3
    l1_dp = REAL(l1, KIND=DP)
    l2_dp = REAL(l2, KIND=DP)
    l1_sq = l1_dp*(l1_dp+1.0_DP)
    l2_sq = l2_dp*(l2_dp+1.0_DP)
    s12   = l1_sq + l2_sq
    p12   = DSQRT( l1_sq * l2_sq )
    f12   = 0.5_dp / p12
    t12   = DSQRT((l1_dp-1.0_DP)*(l1_dp+2.0_DP)  &
         &      * (l2_dp-1.0_DP)*(l2_dp+2.0_DP))

    do l3 = l3_min, l3_max, 2
       w3jarr(l3,1) =  f12 * ( l3*(l3+1) - s12 ) * w3jarr(l3,0)
    enddo

    ! J(l1,l2,l3;-1,1,0), odd l1+l2+l3
    CALL EXTEND_W3J_mm0(w3jarr(0:l3_max,1), -1.0_DP, l3_min, l3_max, l1, l2)

    ! J(l1,l2,l3;-2,2,0)
    do l3 = l3_min, l3_max
       w3jarr(l3,2) = ( (l3*(l3+1)+2 - s12) * w3jarr(l3,1)   &
            &          - p12                * w3jarr(l3,0) ) &
            &         / t12
    enddo


  END SUBROUTINE W3J_220

!-----------------------------------------------------------------------


  SUBROUTINE EXTEND_W3J_mm0(w3_110, m1_dp, l3_min, l3_max, l1, l2)
    !-----------------------------------------------------------------------
    ! Extend Wigner 3j (m1,m2,m3)=(m,-m,0) with m/=0 to odd l1 + l2 + l3
    !-----------------------------------------------------------------------

    REAL(DP),          DIMENSION(0:),       INTENT(INOUT) :: w3_110
    REAL(DP),                               INTENT(IN)    :: m1_dp
    INTEGER(I4B),                           INTENT(IN)    :: l3_min
    INTEGER(I4B),                           INTENT(IN)    :: l3_max
    INTEGER(I4B),                           INTENT(IN)    :: l1
    INTEGER(I4B),                           INTENT(IN)    :: l2

    INTEGER(I4B)                                          :: l3
    REAL(DP)                                              :: l1_dp
    REAL(DP)                                              :: l2_dp, m2_dp
    REAL(DP)                                              :: l3_dp, l3_dp1
    REAL(DP)                                              :: l3sq, l3sq1
    REAL(DP)                                              :: xs, xd, dm
    REAL(DP)                                              :: factor1
    REAL(DP)                                              :: factor2
    REAL(DP)                                              :: denominator

    m2_dp = -m1_dp
    l1_dp = REAL(l1, KIND=DP)
    l2_dp = REAL(l2, KIND=DP)
    l3_dp = REAL(l3_min+1, KIND=DP)
    xd     = (l1_dp - l2_dp         )**2
    xs     = (l1_dp + l2_dp + 1.0_DP)**2
    dm     =  m1_dp - m2_dp 

    DO l3 = l3_min+1,l3_max-1,2

       l3_dp  = real(l3, kind = dp)
       l3_dp1 = l3_dp  + 1.0_DP
       l3sq   = l3_dp  * l3_dp
       l3sq1  = l3_dp1 * l3_dp1

       factor1 = l3_dp1 * DSQRT(  (l3sq  - xd) * (xs - l3sq ) * l3sq  )
       factor2 = l3_dp  * DSQRT(  (l3sq1 - xd) * (xs - l3sq1) * l3sq1 )
       denominator =  (l3_dp+l3_dp1) * l3_dp*l3_dp1 * dm

       w3_110(l3) = (factor1 * w3_110(l3-1) + factor2 *w3_110(l3+1))&
            &       / denominator

    ENDDO


  END SUBROUTINE EXTEND_W3J_MM0

!-----------------------------------------------------------------------

  SUBROUTINE W3J_000(w3jarr, l1, l2, l3_min, l3_max)

    INTEGER(I4B),                           INTENT(IN)    :: l1
    INTEGER(I4B),                           INTENT(IN)    :: l2
    real(DP), dimension(0:), intent(out) :: w3jarr
    integer(I4B),                           intent(out)   :: l3_min, l3_max

    INTEGER(I4B)                   :: l
    INTEGER(I4B)                   :: l3
    REAL(DP)                       :: l_dp, hl_dp
    REAL(DP)                       :: l1_dp
    REAL(DP)                       :: l2_dp
    REAL(DP)                       :: l3_dp
    REAL(DP)                       :: sign
    REAL(DP)                       :: x(8)
    REAL(DP)                       :: g_x(8)
    REAL(DP)                       :: two_l1, two_l2, two_l3
    real(DP)                       :: x1, x2, x3


    if (.not. init_done) then
       call fatal_error('Initialisation not done')
    endif

    w3jarr = 0.0_DP

    l3_min = ABS(l1 - l2)
    l3_max = min(l1 + l2, l_top)

    ! start recursion:  (l1, l2, |l1-l2| ; 0,0,0)
    l3     = l3_min
    l      = l1 + l2 + l3

    l_dp  = REAL(l,  KIND=DP)
    l1_dp = REAL(l1, KIND=DP)
    l2_dp = REAL(l2, KIND=DP)
    l3_dp = REAL(l3_min, KIND=DP)

    x(1) = 1.0_DP + l_dp - 2.0_DP*l1_dp
    x(2) = 1.0_DP + l_dp - 2.0_DP*l2_dp
    x(3) = 1.0_DP + l_dp - 2.0_DP*l3_dp
    x(4) = 1.0_DP + l_dp + 1.0_DP
    x(5) = 1.0_DP + 0.5_DP*l_dp
    x(6) = 1.0_DP + 0.5_DP*l_dp - l1_dp
    x(7) = 1.0_DP + 0.5_DP*l_dp - l2_dp
    x(8) = 1.0_DP + 0.5_DP*l_dp - l3_dp

#ifdef F2008
    g_x = log_gamma(x)
#else
    g_x = gammln(x)
#endif

    sign = (-1.0_DP)**(l/2)

    w3jarr(l3) = sign&
         & * DEXP(0.5_DP * (g_x(1) + g_x(2) + g_x(3) - g_x(4))&
         &                 +g_x(5) - g_x(6) - g_x(7) - g_x(8) )


    ! other terms:   (l1, l2, l3 ; 0,0,0)
    two_l1 = 2.0_DP*l1_dp
    two_l2 = 2.0_DP*l2_dp

    DO l3 = l3_min+2, l3_max, 2

       l_dp  = l_dp  + 2.0_DP
       l3_dp = l3_dp + 2.0_DP
       hl_dp = 0.5_DP * l_dp
       two_l3 = 2.0_DP*l3_dp
       x1    = l_dp - two_l1
       x2    = l_dp - two_l2
       x3    = l_dp - two_l3

       w3jarr(l3) = - w3jarr(l3-2)&
            & * DSQRT(  &
            &   ( x1 - 1.0_DP ) * x1 &
            &  *( x2 - 1.0_DP ) * x2 &
            &  /((x3 + 1.0_DP ) *(x3 + 2.0_DP)&
            &    *l_dp * ( l_dp + 1.0_DP)  )  )  &
            & *  hl_dp * (hl_dp + 1.0_DP - l3_dp)&
            & / ((hl_dp - l1_dp) * (hl_dp - l2_dp))

    ENDDO


  END SUBROUTINE W3J_000


  !=======================================================================
  !=======================================================================
  ! routines for the computation of Log(Gamma)
  ! useless if log_gamma is available on compiler 
  ! (as is the case in F2008 compilers)
  !=======================================================================
  !=======================================================================

  subroutine vdlgamma_old(n, x, y)
    integer(I4B), intent(in) :: n
    real(DP),     intent(in),  dimension(1:n) :: x
    real(DP),     intent(out), dimension(1:n) :: y
    integer(I4B) :: i

    do i=1, n
       y(i) = gammln(x(i))
    enddo

  end subroutine vdlgamma_old
    
  !=======================================================================
  FUNCTION gammln_s(xx) result (gammln)
    !=======================================================================
    !     computes log of Gamma(xx)
    !     from NumRec changed to double precision
    !=======================================================================
    real(DP):: gammln,xx
    INTEGER(I4B):: j
    real(DP):: ser,tmp,x,y
    !=======================================================================
    
    x=xx
    y=x
    tmp=x+0.5d0 + g0
    tmp=(x+0.5d0)*log(tmp)-tmp
    ser=cof(0)
    do j=1,ncof
       y=y+1.d0
       ser=ser+cof(j)/y
    enddo
    gammln=tmp+log(stp*ser/x)
    return
  END FUNCTION gammln_s
  !=======================================================================
  function gammln_v(xin) result (yout)
    !=======================================================================
    !     computes log of Gamma(xx)
    !     from NumRec changed to double precision, works on array
    !=======================================================================
    real(DP),     intent(in),  dimension(1:) :: xin
    real(DP),                  dimension(1:size(xin)) :: yout
    INTEGER(I4B):: j
    real(DP),                  dimension(1:size(xin)) :: ser,tmp,y
    !=======================================================================
    
    y   = xin
    tmp =  xin+0.5d0 + g0
    tmp = (xin+0.5d0)*log(tmp)-tmp
    ser = cof(0)
    do j=1,ncof
       y   = y+1.d0
       ser = ser+cof(j)/y
    enddo
    yout = tmp+log(stp*ser/xin)
    return
  end function gammln_v
  


end module rec3jjcmb_module
