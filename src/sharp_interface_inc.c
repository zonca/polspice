void X(map2alm_cut) (int nside, int lmax, int mmax, int nda, int nrings, int *rings, ptrdiff_t npix, 
		     FLT *map, FLT *alm, double *wgt)
  {
  sharp_geom_info *ginfo;
  sharp_alm_info  *ainfo;

  CHECK_STACK_ALIGN(8);
  sharp_ringsubset_healpix_geom_info(nside, 1, nrings, rings, wgt, &ginfo);
  if (nda == 2) {
    sharp_make_rectangular_alm_info (lmax,mmax,1,&ainfo);
  } else {
    sharp_make_triangular_alm_info (lmax,mmax,1,&ainfo);
  }

  sharp_execute(SHARP_MAP2ALM,0,&alm,&map,ginfo,ainfo,FLAG,NULL,NULL);

  sharp_destroy_alm_info(ainfo);
  sharp_destroy_geom_info(ginfo);
  }

void X(map2alm_pol_cut) (int nside, int lmax, int mmax, int nda, int nrings, int *rings, ptrdiff_t npix, 
			 FLT *map, FLT *alm, double *wgt)
  {
  sharp_geom_info *ginfo;
  sharp_alm_info  *ainfo;
  void *mapptr[2], *almptr[2];

  CHECK_STACK_ALIGN(8);
  sharp_ringsubset_healpix_geom_info(nside, 1, nrings, rings, wgt, &ginfo);
  if (nda == 2) {
    sharp_make_rectangular_alm_info (lmax,mmax,3,&ainfo);
  } else {
    sharp_make_triangular_alm_info (lmax,mmax,3,&ainfo);
  }

  sharp_execute(SHARP_MAP2ALM,0,&alm,&map,ginfo,ainfo,FLAG,0,0);
  mapptr[0]=map+npix; mapptr[1]=map+2*npix;
  almptr[0]=alm+2; almptr[1]=alm+4;
  sharp_execute(SHARP_MAP2ALM,2,&almptr[0],&mapptr[0],ginfo,ainfo,FLAG,0,0);

  sharp_destroy_alm_info(ainfo);
  sharp_destroy_geom_info(ginfo);
  }


