;pro test_cov


; tests:
; -----
; (single map)
; no mask       OK
; mask = 1      OK
; binary mask   OK
; arbitrary mask    not sure
init_healpix
loadct, 0
!P.color = 255

nside = 128
lmax  = 2*nside
npix = nside2npix(nside)
if undefined(type) then type=2
mask = replicate(0.1,npix)
case type of
    0: ; uniform mask
    1: mask[0:3] = 0
    2: mask [0:npix/3] = 0
    3: mask = findgen(npix)/npix
    else: message,'unvalid choice'
endcase
fwhm = 60.0d0
bin = 40

fsky = mean(mask gt 0,/double)
f2 = mean(mask^2,/double)
f4 = mean(mask^4,/double)

cl_in = !healpix.path.test+'cl.fits'
snside = string(nside, form='(i4.4)')
slmax  = string(lmax,  form='(i4.4)')

grp_map = !DMCDATA+'/tmp_efh_MAP_'+snside
grp_cl  = !DMCDATA+'/tmp_efh_CL_'+slmax
grp_cov = !DMCDATA+'/tmp_efh_T3_'+slmax
obj_map = grp_map+'/in_map'
obj_mask = grp_map+'/in_mask'
obj_cl  = grp_cl + '/out_cl'
obj_cov = grp_cov + '/out_cov'

;;;;;;;;;fits2cl, cl0, cl_in, mult=l, llfact=fl
; isynfast, cl_in, map, nlmax = lmax, fwhm = fwhm, simu=1, /silent

; mypiowrite, map, obj_map, /map
mypiowrite, mask, obj_mask, /map
mollview,mask, title='Mask'
loadct,0

ispice, /dmc, obj_map, nlmax = lmax, obj_cl, covfileout = obj_cov, fwhm1 = fwhm, /pixel, binpath='/home/hivon/wrk/tmp/InstallArea/Linux-x86_64/bin/spice',/keep, /subav, weightfile1=obj_mask, mapfile2=obj_map, weightfile2=obj_mask ; maskfile1=obj_mask


cov = mypioread(obj_cov)
cl  = mypioread(obj_cl, mult=l, llfact=fl)
d = diag_matrix(cov[*,*,0])

bin_cl_cov_matrix, cov[*,*,0], bin, l_binned, cov_binned, /flatten, /uniform
bin_llcl,          cl[*,0],    bin, l_b2,     c_binned,   /flatten, /uniform
db = diag_matrix(cov_binned)

title = 'Fsky = '+string(fsky)
window,/free
l0 = 105
tek_color
plot, l, fl*fl[l0]*cov[*,l0],title='Cov(l, '+string(l0,form='(i4)')+'), '+title
oplot, l_binned, cov_binned[*,l0/bin], lines=2, col=2

window,/free
y = sqrt(abs(d))/cl*sqrt((l+0.5)*f2*f2/f4)
plot, title=title, l,y  ;, xr=[0,lmax] ;,/ylog

window,/free
yb = sqrt(abs(db))/( c_binned/sqrt((l_binned+0.5)*f2*f2/f4*bin) )
plot, title=title, l_binned,yb  ;, xr=[0,lmax] ;,/ylog

plottv,/free,alog10(abs(cov)),/iso,/scale,col=39,title=title
plottv,/free,alog10(abs(cov_binned)),/iso,/scale,col=39,title=title
print,fsky,f2,f4
print,mean(y[20:200]),mean(yb[20/bin:200/bin])


;return
end
