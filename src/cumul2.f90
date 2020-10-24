module my_cumul

! requires spice_parameters and spice_common
private
public :: cumul

contains


!=======================================================================
subroutine trapzd(integrand,thetamin,thetamax,j,cl,nlmax,ncor,s,mask_correction)
!=======================================================================
! This computes the jth refinement of the trapezoidal
! extended rule
! From Numerical Recipes in FORTRAN, Second Edition,
! Press, Teukolsky, Vetterling & Flannery, 1992
! Cambridge University Press
!======================================================================= 
  use healpix_types
  implicit none

  real(DP), external :: integrand
  real(DP),                            intent(IN)           :: thetamin,thetamax
  real(DP),                            intent(INOUT)        :: s
  integer(I4B),                        intent(IN)           :: j, nlmax, ncor
  real(DP), dimension(0:nlmax,1:ncor), intent(IN)           :: cl
  logical(LGT),                        intent(IN)           :: mask_correction
  ! local
  real(DP) :: del,x,somme
  integer(I4B) :: it,jj

  if (j==1) then
     s = 0.50_DP*(thetamax-thetamin)*(integrand(thetamin,cl,nlmax &
          &,ncor,mask_correction) + integrand(thetamax,cl,nlmax &
          &,ncor,mask_correction))
     return
  else
     it=2**(j-2)
     del = (thetamax-thetamin)/dble(it)
     x = thetamin+0.5*del
     somme=0.0_DP

     do jj=1,it
        somme = somme + integrand(x,cl,nlmax,ncor,mask_correction)
        x = x + del
     enddo
     s = 0.50_DP * (s + (thetamax-thetamin)*somme/dble(it))
  endif

end subroutine trapzd

!=======================================================================
function simpson2(integrand,thetamin,thetamax,cl,nlmax,ncor,nmax&
     &,verbose,mask_correction)
!=======================================================================
! Extended Simpson integrator of the integrand function
! From Numerical Recipes in FORTRAN, Second Edition,
! Press, Teukolsky, Vetterling & Flannery, 1992
! Cambridge University Press
!=======================================================================
  use healpix_types
  use spice_parameters, only: tolerance
  implicit none

  integer(I4B),                        intent(IN)           :: nlmax,ncor
  integer(I4B),                        intent(OUT)          :: nmax
  real(DP), dimension(0:nlmax,1:ncor), intent(IN)           :: cl
  logical(LGT),                        intent(IN)           :: verbose, mask_correction
  real(DP),                            intent(IN)           :: thetamin, thetamax
  real(DP), external                                        :: integrand
  real(DP)                                                  :: simpson2, residual
  ! local
  integer(I4B) :: j, jmax
  real(DP) :: st, oldsimpson, oldst !, tolerance

  tolerance = abs(tolerance)
  jmax = 20
  if (tolerance <= 1.d-8) jmax = 25
  if (tolerance <= 1.d-9) jmax = 30
  
  oldst = -1d30
  oldsimpson = -1d30
  residual = 1d30

  st = 0.0_DP
  do j=1,JMAX
     if (verbose) write(*,'(a, i3, i3, 1pe10.2, 1pe10.2)') &
          & 'Simpson2: iter, max, residual, tolerance: ', &
          & j-1, JMAX, residual, tolerance
     nmax=j
     call trapzd(integrand,thetamin,thetamax,j,cl,nlmax,ncor,st,mask_correction)
     simpson2 = (4.0_DP*st-oldst)/3.0_DP
     residual = simpson2-oldsimpson
     if ((abs(residual) < tolerance*abs(oldsimpson)).and. j>5) then
        if (verbose) write(*,'(a, i3, 3x, 1pe10.2)') &
            &  'Simpson2: Convergence reached:            ',&
            & j, abs(residual / oldsimpson)
        return
     endif
     residual = abs(residual / oldsimpson)
     oldsimpson = simpson2
     oldst = st
  enddo
  if (verbose) then
     write(*,*) 'Simpson integrator did not converge to required relative accuracy: ',tolerance
     write(*,*) 'After ',JMAX,' iterations, one gets: ',residual
  endif

end function simpson2


!=======================================================================
function fsub_first(beta,cl,nlmax,ncor,mask_correction)
!=======================================================================
! This will be the first integrand
!======================================================================= 
  use healpix_types
  implicit none

  integer(I4B),                        intent(IN)          :: nlmax, ncor
  real(DP),                            intent(IN)          :: beta
  real(DP), dimension(0:nlmax,1:ncor), intent(IN)          :: cl
  real(DP)                                                 :: fsub_first
  logical(LGT),                        intent(IN)          :: mask_correction
  ! local
  real(DP) :: c_plus

  c_plus = cplus(cos(beta),cl,nlmax,ncor,mask_correction)
  fsub_first = sin(beta)/cos(beta/2.0_DP)**4*c_plus
end function fsub_first

!=======================================================================
function fsub_second(beta,cl,nlmax,ncor,mask_correction)
!=======================================================================
! This will be the second integrand
!=======================================================================
  use healpix_types
  implicit none

  integer(I4B),                        intent(IN)          :: nlmax, ncor
  real(DP),                            intent(IN)          :: beta
  real(DP), dimension(0:nlmax,1:ncor), intent(IN)          :: cl
  logical(LGT),                        intent(IN)          :: mask_correction
  real(DP)                                                 :: fsub_second
  ! local
  real(DP) :: c_plus

  c_plus = cplus(cos(beta),cl,nlmax,ncor,mask_correction)
  fsub_second = tan(beta/2.0_DP)**3*c_plus
end function fsub_second

!========================================================================
function cplus(x,cl,nlmax,ncor,mask_correction)
!========================================================================
! This function computes C+(x)=<QQ>+<UU>(x) for given x and power spectra
! Corrects for edge effects
!========================================================================
  use healpix_types
  use spice_common, only : cl_mask, ncmask
  use do_legendre_mod, only: Legendre
  implicit none

  integer(I4B),                        intent(IN)           :: nlmax, ncor
  real(DP),                            intent(IN)           :: x
  real(DP), dimension(0:nlmax,1:ncor), intent(IN)           :: cl
  logical(LGT),                        intent(IN)           :: mask_correction
  real(DP)                                                  :: cplus
  ! local
  integer(I4B) :: l
  real(DP), allocatable, dimension(:,:) :: Pl
  real(DP) :: cmask
  integer :: kmask

  if (ncor < 4) then
     print*,'No polarisation involved, you should not be in here.'
     return
  endif
  allocate(Pl(0:nlmax,0:ncor-1))
  
  call Legendre(nlmax,x,Pl,ncor)
  
  kmask = 1
  if (ncmask >= 3) kmask = 2 ! PP part of mask C(l)

  cplus=0.0_DP
  cmask=0.0_DP
  if (mask_correction) then
     do l=0,nlmax
        cplus = cplus + 2.0_DP*(cl(l,2)+cl(l,3))*(Pl(l,1)+Pl(l,2))*dble(2*l+1)
        cmask = cmask + cl_mask(l,kmask)*Pl(l,0)*dble(2*l+1)
     enddo
     if (cmask > 0.0_DP) then
        cplus = cplus/cmask
     else
        cplus = 0.0_DP
     endif
  else
     do l=0,nlmax
        cplus = cplus + 2.0_DP*(cl(l,2)+cl(l,3))*(Pl(l,1)+Pl(l,2))*dble(2*l+1)
     enddo
     cplus = cplus/4.0_DP/PI
  endif
     
  deallocate(Pl)
  return
end function cplus

!=======================================================================
subroutine cumul_simpson(integrand,theta,sum,cl,nlmax,ncor,nmax,thetamax,mask_correction)
!=======================================================================
! This computes the cumulative integral from 0 to theta(j) of 
! the integrand function, with a fixed set of function 
! evaluations, with an extended Simpson rule.
! Care is taken of the end-point positions (theta(j)) and a 
! specific simple Simpson rule is dedicated to them.
! The integral is calculated up to thetamax = theta(jmin)
!=======================================================================
  use healpix_types
  implicit none

  real(DP), external :: integrand
  integer(I4B), intent(IN)                             :: nlmax, ncor, nmax
  real(DP), dimension(nlmax+1),   intent(IN)           :: theta
  real(DP), dimension(nlmax+1),   intent(OUT)          :: sum
  real(DP), dimension(0:nlmax,1:ncor), intent(IN)      :: cl
  real(DP),                       intent(IN)           :: thetamax
  logical(LGT),                   intent(IN)           :: mask_correction
  ! local
  real(DP) ,allocatable,  dimension(:) :: ftemp
  integer(I4B)                         :: i, j, npoints
  integer(I4B), dimension(nlmax+1)     :: lastn
  real(DP)                             :: step, step2, endpoint, temp, xsti

  npoints = 2**(nmax-1) ! Fine-grid sampling of [0,thetamax]

  allocate(ftemp(0:npoints))
  step = thetamax/dble(npoints)

  ! Precompute the integrand on the fine grid
!$OMP PARALLEL &
!$OMP   DEFAULT(shared) &
!$OMP   PRIVATE(i,xsti)
!!$OMP   SHARED(npoints, ftemp, step, cl, nlmax, ncor, mask_correction) &
!$OMP DO SCHEDULE(DYNAMIC, 5)
  do i=0,npoints
     xsti = dble(step*i)
     ftemp(i) = integrand(xsti,cl,nlmax,ncor,mask_correction)
  enddo
!$OMP END DO
!$OMP END PARALLEL


  sum=0
  lastn=-1

!$OMP PARALLEL &
!$OMP   DEFAULT(SHARED) &
!$OMP   PRIVATE(j,temp,i,endpoint,step2)
!$OMP DO SCHEDULE(DYNAMIC, 5)
  do j=1,nlmax+1
     ! Last even index of the fine grid smaller than theta(j)
     ! This is the endpoint of the extended Simpson integration
     lastn(j) = 2*nint(theta(j)/(2.0_DP*step))

     if (lastn(j)==0) then
        sum(j)=0
     else
        sum(j) = (ftemp(0)+ftemp(lastn(j)))/3.0_DP ! First and last term => special coeff
        temp=0.0_DP
        do i=1,lastn(j),2
           temp = temp + ftemp(i)
        enddo
        sum(j) = sum(j) + temp * 4.0_DP/3.0_DP !odd terms
        temp=0.0_DP
        do i=2,lastn(j)-2,2
           temp = temp + ftemp(i)
        enddo
        sum(j) = sum(j) + 2.0_DP/3.0_DP*temp ! Adding even terms
        sum(j) = sum(j)*step
     endif
     ! Now correct for end point by adding a last Simpson rule
     step2 = (theta(j)-dble(step*lastn(j)))/2.0_DP
     endpoint = step2*(ftemp(lastn(j))+4.0_DP * &
          integrand(theta(j)-step2,cl,nlmax,ncor,mask_correction) &
          &   + integrand(theta(j),cl,nlmax,ncor,mask_correction))/3.0_DP
     sum(j) = sum(j) + endpoint
  enddo
!$OMP END DO
!$OMP END PARALLEL
     
  deallocate(ftemp)

end subroutine cumul_simpson

!=======================================================================
subroutine cumul(nlmax,mu,xi,cl,ncor,verbose,message,thetamax_degrees,mask_correction)
!=======================================================================
! Computes XiE and XiB by cumulative integration
! xi contains the edge-corrected values at the (rescaled) Pl roots
! cl contains the non-edge-corrected values of the cls
!=======================================================================
  use healpix_types
  use healpix_modules, only: wall_clock_time
  implicit none
   
  integer(I4B),                        intent(IN)           :: nlmax,ncor
  real(DP), dimension(nlmax+1),        intent(IN)           :: mu
  real(DP), dimension(0:nlmax,1:ncor), intent(INOUT)        :: xi
  real(DP), dimension(0:nlmax,1:ncor), intent(IN)           :: cl
  real(DP),                            intent(IN)           :: thetamax_degrees
  character(len=*),                    intent(IN)           :: message
  logical(LGT),                        intent(IN)           :: verbose, mask_correction
  ! local
  integer(I4B) :: j,nmax,n
  real(DP), dimension(nlmax+1) :: sum1, sum2, c_beta, c_plus, theta
  real(DP) :: t1, t2, x2, x3, thetamax
  real(SP) :: ct0,ct1,wt0, wt1


  call cpu_time(ct0)
  call wall_clock_time(wt0)
  if (ncor < 4) return

  if (verbose) write(*,*) 'Compute pure E and B modes correlation function ' &
 &                            //trim(message)

  ! Convert thetamax to radians
  thetamax = thetamax_degrees *PI/180.d0

  t1 = simpson2(fsub_first,  0.d0,thetamax,cl,nlmax,ncor,nmax,verbose,mask_correction)
  t2 = simpson2(fsub_second, 0.d0,thetamax,cl,nlmax,ncor,   n,verbose,mask_correction)

  nmax=max(n,nmax)

  ! Now compute the cumulative distribution by an extended Simpson rule
  ! plus a regular Simpson rule for the (unevenly sampled) end point.

  do j=1,nlmax+1
     theta(j) = dacos(mu(j))
  enddo

  sum1 = 0
  sum2 = 0
  c_beta = 0

  call cumul_simpson(fsub_first,  theta,sum1,cl,nlmax,ncor,nmax,thetamax,mask_correction)
  call cumul_simpson(fsub_second ,theta,sum2,cl,nlmax,ncor,nmax,thetamax,mask_correction)

!$OMP PARALLEL DO &
!$OMP DEFAULT(SHARED) &
!$OMP PRIVATE(j,x2,x3)
  do j=1,nlmax+1
     c_plus(j) = cplus(mu(j),cl,nlmax,ncor,mask_correction)
     c_beta(j) = c_plus(j)
     c_beta(j) = c_beta(j) + sum1(j)/sin(theta(j)/2.0_DP)**2 &
          & -sum2(j)*2.d0*(2.d0+cos(theta(j)))/sin(theta(j)/2.0_DP)**4
     x2 = xi(j-1,2)
     x3 = xi(j-1,3)
     xi(j-1,2) = 0.5_DP * (c_beta(j)+x2-x3) ! 1/2[C(\beta)+ReC_{-}(\beta)]
     xi(j-1,3) = 0.5_DP * (c_beta(j)-x2+x3) ! 1/2[C(\beta)-ReC_{-}(\beta)]
  enddo
!$OMP END PARALLEL DO


  call cpu_time(ct1)
  call wall_clock_time(wt1)
  if (verbose) print*,'E/B decoupling (CPU,wall) time [s]: ',ct1-ct0,wt1-wt0

end subroutine cumul

end module my_cumul
