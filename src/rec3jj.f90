module rec3jj_module

  public :: rec3jj

  contains

! ======================================================================
subroutine rec3jj(thrcof, l2, l3, m2, m3, l1min, l1max, ndim, ier)
  use healpix_types
  implicit none

  real(DP),     intent(in) :: l2,l3, m2,m3
  integer(i4b), intent(in) :: ndim
  real(DP),     intent(out), dimension(ndim) :: thrcof
  integer(i4b), intent(out) :: ier
  real(DP),     intent(out) :: l1min,l1max

  real(DP) :: newfac, l1, m1

  real(DP) :: a1, a2, a1s, a2s, c1, c1old, c2, cnorm, denom, dv
  real(DP) :: oldfac, ratio, sign1, sign2, sumbac, sumfor, sumuni, sum1, sum2, thresh
  real(DP) :: x, x1, x2, x3, y, y1, y2, y3
  integer(i4b) :: i, index, lstep, n, nfin, nlim, nstep2, nfinp1, nfinp2, nfinp3

  real(DP),parameter :: zero=0.00000000000000_DP
  real(DP),parameter :: eps =0.01000000000000_DP
  real(DP),parameter :: one =1.00000000000000_DP
  real(DP),parameter :: two =2.00000000000000_DP
!
! some compilers (eg, ifort 9) do not accept SQRT() in parameter definition
!
!   real(DP),parameter :: megahuge=sqrt(MAX_DP/20.0_dp)
!   real(DP),parameter :: srhuge  =sqrt(megahuge)
  real(DP),parameter :: megahuge= 1.e100_dp
  real(DP),parameter :: srhuge  = 1.e50_dp
!                   ----------
  real(DP),parameter :: tiny    =one/megahuge
  real(DP),parameter :: srtiny  =one/srhuge

! ---------------------------------------------------------------------------------
! adapted from Slatec's DRC3JJ.F (http://www.netlib.org/slatec/src/)
!   (Gordon & Schulten, 1975)
!
! modified by Gayoung CHON, 2002?
! re-edited by Eric Hivon, 2008-12-04: 
!  - reset to original values the parameters MEGAHUGE, TINY...  to achieve better precision 
!   on small 3J's (more relevant for large abs(m))
!  - implicit none
!  -  1.0 -> ONE
!  - removed useless variables (HALF, l1cmin, l1cmax, lmatch)
! ---------------------------------------------------------------------------------
!
! routine to generate set of 3j-coeffs (l1,l2,l3\\ m1,m2,m3)

! by recursion from l1min = max(abs(l2-l3),abs(m1)) 
!                to l1max = l2+l3
! the resulting 3j-coeffs are stored as thrcof(l1-l1min+1)

! to achieve the numerical stability, the recursion will proceed
! simultaneously forwards and backwards, starting from l1min and l1max
! respectively.
!
! lmatch is the l1-value at which forward and backward recursion are matched.
!
! ndim is the length of the array thrcof
!
! ier = -1 for all 3j vanish(l2-abs(m2)<0, l3-abs(m3)<0 or not integer)
! ier = -2 if possible 3j's exceed ndim
! ier >= 0 otherwise
!

!  lmatch = zero
  m1 = -(m2+m3)

! check relative magnitude of l and m values
  
  if (l2-abs(m2)) 5,1,1
1 if (l3-abs(m3)) 5,2,2
2 if (mod(l2+abs(m2)+eps,one).ge.eps+eps) go to 5
  if (mod(l3+abs(m3)+eps,one).ge.eps+eps) go to 5

! limits for l1
  l1min = max(abs(l2-l3),abs(m1))
  l1max = l2+l3
  if (l1min.lt.l1max-eps) go to 20
  if (l1min.eq.l1max) go to 10
 
! reached if l2-/m2/ and l3-/m3/ less than zero or not integer
5 ier = -1
  return

! reached if l1 can take only one value, i.e.l1min=l1max
10 ier = 0
  thrcof(1) = (-1)**nint(abs(l2+m2-l3+m3))/sqrt(l1min+l2+l3+one)
  return

20 ier = 0
  nfin = nint(l1max-l1min+one)
  if (ndim-nfin) 21,23,23

! the dimension of thrcof is not large enough to hold all the allowed values
! of l1

21 ier = -2
  write(*,*) 'ndim < nfin'
  return

! starting forward recursion from l1min taking nstep1 steps
23 l1 = l1min
  newfac = zero ! added by EH
  c1 = zero     ! added by EH
  thrcof(1) = srtiny
  sum1 = (two * l1 + one)*tiny

  lstep = 1
30 lstep = lstep+1
  l1 = l1 + one

  oldfac = newfac
  a1 = (l1+l2+l3 + one)*(l1-l2+l3)*(l1+l2-l3)*(-l1+l2+l3 + one)
  a2 = (l1+m1)*(l1-m1)
  newfac = sqrt(a1*a2)
  if (l1.lt.one+eps) go to 40

  dv = -l2*(l2 + one)*m1 + l3*(l3 + one)*m1 + l1*(l1 - one)*(m3-m2)
  denom = (l1 - one)*newfac

  if (lstep-2) 32,32,31
  
31 c1old = abs(c1)
32 c1 = -(two * l1 - one)*dv/denom
  go to 50

! if l1=1, (l1-1) has to be factored out of dv, so

40 c1 = -(two * l1 - one)*l1*(m3-m2)/newfac

50 if (lstep.gt.2) go to 60

! if l1=l1min+1 the third term in the recursion eqn vanishes, hence
  x = srtiny*c1
  thrcof(2) = x
  sum1 = sum1+tiny*(two * l1 + one)*c1*c1
  if(lstep.eq.nfin) go to 220
  go to 30

60 c2 = -l1*oldfac/denom

! recursion to the next 3j-coeff x  
  x = c1*thrcof(lstep-1) + c2*thrcof(lstep-2)
  thrcof(lstep) = x
  sumfor = sum1
  sum1 = sum1 + (two * l1 + one)*x*x
  if (lstep.eq.nfin) go to 100

! see if last unnormalised 3j-coeff exceeds srhuge
  if (abs(x).lt.srhuge) go to 80

! REACHED IF LAST 3J-COEFFICIENT LARGER THAN SRHUGE
! SO THAT THE RECURSION SERIES THRCOF(1), ... , THRCOF(LSTEP)
! HAS TO BE RESCALED TO PREVENT OVERFLOW

  ier = ier+1
  do i = 1, lstep
     if (abs(thrcof(i)).lt.srtiny) thrcof(i)= zero
     thrcof(i) = thrcof(i)/srhuge
  end do

  sum1 = sum1/megahuge
  sumfor = sumfor/megahuge
  x = x/srhuge

! as long as abs(c1) is decreasing, the recursion proceeds towards increasing
! 3j-valuse and so is numerically stable. Once an increase of abs(c1) is 
! detected, the recursion direction is reversed.

80 if (c1old-abs(c1)) 100,100,30

! keep three 3j-coeffs around lmatch for comparison with backward recursion

100 continue
!  lmatch = l1 - one
  x1 = x
  x2 = thrcof(lstep-1)
  x3 = thrcof(lstep-2)
  nstep2 = nfin-lstep+3

! --------------------------------------------------------------------------
!
! starting backward recursion from l1max taking nstep2 stpes, so that
! forward and backward recursion overlap at 3 points 
! l1 = lmatch-1, lmatch, lmatch+1

  nfinp1 = nfin+1
  nfinp2 = nfin+2
  nfinp3 = nfin+3
  l1 = l1max
  thrcof(nfin) = srtiny
  sum2 = tiny*(two * l1 + one)
 
  l1 = l1 + two
  lstep=1
110 lstep = lstep + 1
  l1= l1 - one

  oldfac = newfac
  a1s = (l1+l2+l3)*(l1-l2+l3 - one)*(l1+l2-l3 - one)*(-l1+l2+l3 + two)
  a2s = (l1+m1 - one)*(l1-m1 - one)
  newfac = sqrt(a1s*a2s)

  dv = -l2*(l2 + one)*m1 + l3*(l3 + one)*m1 +l1*(l1 - one)*(m3-m2)

  denom = l1*newfac
  c1 = -(two * l1 - one)*dv/denom
  if (lstep.gt.2) go to 120

! if l2=l2max+1, the third term in the recursion vanishes

  y = srtiny*c1
  thrcof(nfin-1) = y
  sumbac = sum2
  sum2 = sum2 + tiny*(two * l1-3)*c1*c1

  go to 110

120 c2 = -(l1 - one)*oldfac/denom

! recursion to the next 3j-coeff y
  y = c1*thrcof(nfinp2-lstep)+c2*thrcof(nfinp3-lstep)

  if (lstep.eq.nstep2) go to 200
  
  thrcof(nfinp1-lstep) = y
  sumbac = sum2
  sum2 = sum2+(two * l1-3)*y*y

! see if last unnormalised 3j-coeff exceeds srhuge
  if (abs(y).lt.srhuge) go to 110

! reached if 3j-coeff larger than srhuge so that the recursion series
! thrcof(nfin),..., thrcof(nfin-lstep+1) has to be rescaled to prevent overflow

  ier=ier+1
  do i = 1, lstep
     index=nfin-i+1
     if (abs(thrcof(index)).lt.srtiny) thrcof(index)=zero
     thrcof(index) = thrcof(index)/srhuge
  end do

  sum2=sum2/megahuge
  sumbac=sumbac/megahuge

  go to 110

! the forward recursion 3j-coeffs x1, x2, x3 are to be matched with the 
! corresponding backward recursion vals y1, y2, y3



200 y3 = y
  y2 = thrcof(nfinp2-lstep)
  y1 = thrcof(nfinp3-lstep)

! determine now ratio such that yi=ratio*xi (i=1,2,3) holds with minimal error

  ratio = (x1*y1+x2*y2+x3*y3)/(x1*x1+x2*x2+x3*x3)
  nlim = nfin-nstep2+1

  if (abs(ratio).lt. one) go to 211

  do n = 1, nlim
     thrcof(n) = ratio*thrcof(n)
  end do
  
  sumuni = ratio*ratio*sumfor + sumbac
  go to 230

211 nlim = nlim+1
  ratio = one/ratio
  do n = nlim, nfin
     thrcof(n) = ratio*thrcof(n)
  end do
  sumuni = sumfor + ratio*ratio*sumbac
  go to 230

220 sumuni=sum1

! normalise 3j-coeffs


230 cnorm = one/sqrt(sumuni)

! sign convention for last 3j-coeff determines overall phase

  sign1 = sign(one,thrcof(nfin))
  sign2 = (-1)**nint(abs(l2+m2-l3+m3))
  if (sign1*sign2) 235,235,236
235 cnorm = -cnorm

236 if (abs(cnorm).lt.one) go to 250

  do n = 1, nfin
     thrcof(n) = cnorm*thrcof(n)
  end do
  return

250 thresh = tiny/abs(cnorm)

  do n = 1, nfin
     if (abs(thrcof(n)).lt.thresh) thrcof(n) = zero
     thrcof(n) = cnorm*thrcof(n)
  end do
  return 

end subroutine rec3jj

end module rec3jj_module
