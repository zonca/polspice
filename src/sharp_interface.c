#include "libsharp/sharp.h"
#include "libsharp/sharp_geomhelpers.h"
#include "libsharp/sharp_almhelpers.h"
#include "c_utils.h"

#include <stdio.h>

/* definition of the routines

map2alm_cut(int nside, int lmax, int mmax, int nrings, int *rings, ptrdiff_t npix, 
            FLT *map, FLT *alm, double *wgt)
  where map is T map, and has size npix, corresponding to subset of nrings (full) rings

map2alm_pol_cut(int nside, int lmax, int mmax, int nrings, int *rings, ptrdiff_t npix, 
                FLT *map, FLT *alm, double *wgt)
  where map is TQU map, and has size 3*npix, corresponding to subset of nrings (full) rings

*/


void sharp_ringsubset_healpix_geom_info (int nside, int stride, int nrings,
  const int *rings, const double *weight, sharp_geom_info **geom_info)
{
  const double pi=3.141592653589793238462643383279502884197;
  ptrdiff_t npix=(ptrdiff_t)nside*nside*12;
  ptrdiff_t ncap=2*(ptrdiff_t)nside*(nside-1);

  double *theta=RALLOC(double,nrings);
  double *weight_=RALLOC(double,nrings);
  int *nph=RALLOC(int,nrings);
  double *phi0=RALLOC(double,nrings);
  ptrdiff_t *ofs=RALLOC(ptrdiff_t,nrings);
  int *stride_=RALLOC(int,nrings);
  ptrdiff_t curofs=0, checkofs; /* checkofs used for assertion introduced when adding rings arg */
  for (int m=0; m<nrings; ++m)
    {
    int ring = (rings==NULL)? (m+1) : rings[m];
    ptrdiff_t northring = (ring>2*nside) ? 4*nside-ring : ring;
    stride_[m] = stride;
    if (northring < nside)
      {
      theta[m] = 2*asin(northring/(sqrt(6.)*nside));
      nph[m] = 4*northring;
      phi0[m] = pi/nph[m];
      checkofs = 2*northring*(northring-1)*stride;
      }
    else
      {
      double fact1 = (8.*nside)/npix;
      double costheta = (2*nside-northring)*fact1;
      theta[m] = acos(costheta);
      nph[m] = 4*nside;
      if ((northring-nside) & 1)
        phi0[m] = 0;
      else
        phi0[m] = pi/nph[m];
      checkofs = (ncap + (northring-nside)*nph[m])*stride;
      // ofs[m] = curofs;
      }
    if (northring != ring) /* southern hemisphere */
      {
      theta[m] = pi-theta[m];
      checkofs = (npix - nph[m])*stride - checkofs;
      // ofs[m] = curofs;
      }
    weight_[m]=4.*pi/npix*((weight==NULL) ? 1. : weight[northring-1]);
    if (rings==NULL) {
        UTIL_ASSERT(curofs==checkofs, "Bug in computing ofs[m]");
    }
    ofs[m] = curofs;
    curofs+=nph[m];
    //printf("%d %d %d\n",m,rings[m],ofs[m]);
  }

  sharp_make_geom_info (nrings, nph, ofs, stride_, phi0, theta, weight_,
    geom_info);

  DEALLOC(theta);
  DEALLOC(weight_);
  DEALLOC(nph);
  DEALLOC(phi0);
  DEALLOC(ofs);
  DEALLOC(stride_);
  }




#define CONCAT(a,b) a ## b 

#define FLT double 
#define FLAG SHARP_DP 
#define X(arg) CONCAT(sharpd_,arg) 
#include "sharp_interface_inc.c" 
#undef FLT 
#undef FLAG 
#undef X 

#define FLT float 
#define FLAG 0 
#define X(arg) CONCAT(sharps_,arg) 
#include "sharp_interface_inc.c" 
#undef FLT 
#undef X 

#undef CONCAT 


// void main() {
//   int nside=8;
//   int stride=1;
//   //  int nrings=4*nside-1;
//   int nrings = 2;
//   int *rings=RALLOC(int,nrings);
//   rings[0] = 4;
//   rings[1] = 6;
//   sharp_geom_info *ginfo;
//   sharp_ringsubset_healpix_geom_info(8, 1, nrings, rings, NULL, &ginfo);
// }
