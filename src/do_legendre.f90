module do_legendre_mod

private
public :: do_legendre, Legendre

contains

!=======================================================================
!subroutine do_legendre(nside, nlmax, mumin, mu, w, Pl, verbose, ncor) ! 2010-03-12
!subroutine do_legendre(nlmax, mumin, mu, w, Pl, verbose, ncor) ! 2019-12-06
subroutine do_legendre(nlmax, mumin, mu, w, Pl, verbose, npl)
!=======================================================================
! This subroutine precalculate quantities related to Legendre polynomials
! for SpICE.
! 2010-03-12: removed unused nside
!=======================================================================
  use healpix_types
  implicit none

  integer(I4B), intent(in) :: nlmax, npl
  real(DP), intent(in)  :: mumin
  real(DP), intent(out) :: mu(nlmax+1), w(nlmax+1)
  real(DP), intent(out) :: Pl(0:nlmax,0:nlmax,0:npl-1)
  logical, intent(in)   :: verbose

  integer(I4B) itheta
  real(DP) :: secPl(0:nlmax,0:npl-1)

  if (verbose) write(*,*) 'Compute Legendre polynomials data'

  call gauleg_double(mumin, 1.0d0, mu, w, nlmax+1)

  do itheta = 0, nlmax
     call Legendre(nlmax, mu(itheta+1), secPl, npl)
     Pl(itheta,:,:) = secPl
  end do

end subroutine do_legendre

!=======================================================================
subroutine Legendre(nlmax, x, Pl, npl)
!=======================================================================
! Compute the Legendre polynomials and the polarized Legendre polynomials
! G_{l2}^{+} et G_{l2}^{-} (voir Kamionkowski et al. 97) at the roots of 
! the legendre polynomial corresponding to order nlmax+1.
!=======================================================================
  use healpix_types
  implicit none

  integer(I4B), intent(in)  :: nlmax, npl
  real(DP),     intent(in)  :: x
  real(DP),     intent(out) :: Pl(0:nlmax,0:npl-1)
  !
  integer(I4B) :: l
  real(DP) :: dl, s, m, sign
  real(DP) :: t1, t2, t3, t4
  real(DP), allocatable, dimension(:) :: Pl2, Nl !m=2
  real(DP) :: nfact = 2.0_dp !SZ(cmbfast), 1 for KKS
  logical :: dopol, bcoupling

  if (npl > 1) then
     dopol=.true.
     ! bcoupling = (npl >= 6) ! 2019-12
     bcoupling = .true.
  else
     dopol=.false.
     bcoupling = .false.
  endif

  if (dopol) allocate(Pl2(0:nlmax), Nl(0:nlmax))

  ! Initialize to zero  
  Pl = 0.0_dp

  ! First compute the regular Legendre polynomials P_{l0}(x);
  ! Border cases  
  if (x==1.0_dp) then
     if (dopol) then
        Pl(:,0) = 1.0_dp
        Pl(2:nlmax,1) = 0.250_dp
        Pl(2:nlmax,2) = 0.250_dp
        Pl(:,3) = 0.0_dp
     else
        Pl(:,0)=1.0_dp
     endif
     return
  endif
  if (x==-1.0_dp) then
     if (dopol) then
        sign=1.0_dp
        Pl(0,0) = 1.0_dp
        Pl(1,0) = -1.0_dp
        do l=2,nlmax
           Pl(l,0) = sign
           Pl(l,1) = sign*0.250_dp
           Pl(l,2) = -sign*0.250_dp
           sign = -sign
        enddo
        Pl(:,3) = 0.0_dp
     else
        Pl(:,0)=1.0_dp
        ! bug correction, 2007-04-27
        do l=1,nlmax,2
           Pl(l,0) = -1.0_dp
        enddo
     endif
     return
  endif
  
  Pl(0,0) = 1.0_dp
  Pl(1,0) = x
  do l=2,nlmax
     dl = real(l,kind=DP)
     Pl(l,0) = ( (2.0_dp*dl-1.0_dp)*x*Pl(l-1,0)-(dl-1.0_dp)*Pl(l-2,0) )/dl
  enddo

  ! Now do the polarized part
  if (dopol) then
     Nl(0:1)=0.0_dp
     do l=2,nlmax
        dl=real(l,kind=DP)
        Nl(l) = sqrt(2.0_dp/nfact)/sqrt((dl+2.0_dp)*(dl+1.0_dp)*dl*(dl-1.0_dp)) !TB used later
     enddo
  
     ! Now compute P_{l2}(x) and N_{l}P_{l2}(x) for l>=2
     Pl2(0:1) = 0.0_dp
     Pl2(2) = 3.0_dp*(1-x*x) 
     Pl(2,3) = Nl(2)*Pl2(2)
     do l=3,nlmax
        dl = real(l,kind=DP)
        Pl2(l) = ( (2.0_dp*dl-1.0_dp)*x*Pl2(l-1) - (dl+1.0_dp)*Pl2(l-2) ) / (dl-2.0_dp)
        Pl(l,3) = Nl(l)*Pl2(l)
     enddo

     ! Finaly compute the G^{\pm}_{l2}(x)    
     s = sin(acos(x))
     m = 2.0_dp

     do l=2,nlmax
        dl = real(l,kind=DP)
        t1 = (dl+m) * x/(s*s)                        * Pl2(l-1)
        t2 = ( (dl-m*m)/(s*s) + (dl-1.0_dp)*dl*0.50_dp ) * Pl2(l)
        Pl(l,1) = Nl(l)**2 * (t1-t2) ! N_l^2 G^{+}_{l2}
        t3 = (dl-1.0_dp) * x * Pl2(l)
        t4 = (dl+m)        * Pl2(l-1)
        Pl(l,2) = Nl(l)**2 * (t3-t4) * m/(s*s) ! N_l^2 G^{-}_{l2}
     enddo
  
     if (bcoupling) then
        Pl(:,4) = Pl(:,3)
        ! Pl(:,5) = Pl(:,2) ! 2019-12
     endif

     deallocate(Pl2, Nl)

  endif

end subroutine Legendre

!=======================================================================
subroutine gauleg_double(x1,x2,x,w,n)
!=======================================================================
! From Numerical Recipes in FORTRAN, Second Edition,
! Press, Teukolsky, Vetterling & Flannery, 1992
! Cambridge University Press
!
! This subroutine has been modified compared to the original
! one (everything is set in double precision and pi=3.14 is computed
! in double precision).
! 2009-06-26: reduced eps from 3.d-15 to 3*EPSILON(1.d0)=6.66e-16
!=======================================================================
  use healpix_types
  implicit none
  integer(I4B), intent(in)  :: n
  real(DP),     intent(in)  :: x1,x2
  real(DP),     intent(out) :: x(n),w(n)
  !
  real(DP), parameter :: EPS = 3*EPSILON(1.0_dp)
  integer(I4B) :: i,j,m
  real(DP) :: p1,p2,p3,pp,xl,xm,z,z1
  m=(n+1)/2
  xm=0.50_dp*(x2+x1)
  xl=0.50_dp*(x2-x1)
  do i=1,m
     z=cos(pi*(i-.250_dp)/(n+.50_dp))
1    continue
     p1=1.0_dp
     p2=0.0_dp
     do j=1,n
        p3=p2
        p2=p1
        p1=((2.0_dp*j-1.0_dp)*z*p2-(j-1.0_dp)*p3)/j
     enddo
     pp=n*(z*p1-p2)/(z*z-1.0_dp)
     z1=z
     z=z1-p1/pp
     if(abs(z-z1).gt.EPS) goto 1
     x(i)=xm-xl*z
     x(n+1-i)=xm+xl*z
     w(i)=2.0_dp*xl/((1.0_dp-z*z)*pp*pp)
     w(n+1-i)=w(i)
  enddo
end subroutine gauleg_double

end module do_legendre_mod
