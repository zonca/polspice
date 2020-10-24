

fmap0 = '/tmp/map_iqu0.fits'
fmap1 = '/tmp/map_iqu1.fits'
fcl0 = '/tmp/cl_spice0.fits'
fcl1 = '/tmp/cl_spice1.fits'
fcor0 = '/tmp/cor_spice0.fits'
fcor1 = '/tmp/cor_spice1.fits'
fcl2 = '/tmp/cl_anafast1.fits'
fmask = '/tmp/mask.fits'
host=getenv('HOST') ; only on IDL 6.1
spawn,'pwd',/sh,cwd
junk = where(host eq ['n21','n22','n23','n24'], nmagique)
if (nmagique eq 1) then begin
; on magique 2
    binpath = cwd+'/spice'
    healpixdir = '~/src/Healpix'
endif else begin
;
    binpath = '~/power_spectrum/spice/HL2_Spice/src/spice'
    binpath = '~/src/HL2_Spice/src/spice'
    healpixdir = '~/Healpix'
endelse
;
nside = 256
lmax = 1.5*nside
fwhm = 30

; ; input map
; isynfast, healpixdir+'/test/cl.fits',fmap0,nside=nside,lmax=lmax,sim=2,/silent,fwhm=fwhm
; ; mollview,fmap0
; ; mollview,fmap0,2

; rotated map
read_tqu,fmap0,tqu0
; a= 22.5 * !dtor
; a= -22.5 * !dtor
; a= -45 * !dtor
if undefined(adeg) then begin
    read,prompt='Enter rotation in Deg',adeg
endif
a=adeg * !dtor
c = cos (2*a) & s = sin(2*a)
rot = [[1, 0, 0],[0, c, s],[ 0, -s, c]] 
; ; ; matrix: (1 0 0)
; ; ;         (0 c s)
; ; ;         (0-s c)
tqu = rot ## tqu0
write_tqu,fmap1,tqu,/ring

; ; mask
; npix = nside2npix(nside)
; ip = lindgen(npix)
; pix2ang_ring, nside, ip, theta, phi
; dec = 90. - theta * !radeg
; mask = (abs(dec) gt 15)
; ;mollview,mask
; print,total(mask)/npix
; write_fits_map, fmask, mask, /ring

;ispice, fmap1, fcl1, /pol, /subav, bin=binpath, show=1,nlmax=lmax,decou=1,maskfile1=fmask,apodizesigma=90.,corfile=fcor1,fwhm1=fwhm

;ispice, fmap0, fcl0, /pol, /subav, bin=binpath, /show
; ispice, fmap1, fcl1, /pol, /subav, bin=binpath, show=0,nlmax=lmax,decou=0,maskfile1=fmask,apodizesigma=90.,corfile=fcor1

; ;ispice, fmap0, fcl1, /pol, /subav, bin=binpath, show=0,
; nlmax=lmax,decou=1,corfile=fcor1

;ispice, fmap0, corfile=fcor0, /pol, /subav, bin=binpath,nlmax=nlmax,fwhm1=fwhm, /silent
;read_spice, fcor0, cor0, theta=theta

ispice, fmap1, corfile=fcor1, /pol, /subav, bin=binpath,nlmax=nlmax,fwhm1=fwhm,/silent ;/keep ;, /silent
read_spice, fcor1, cor1, theta=theta

r00=rot[0,0]
; rotation of (TT, QQ, UU, TQ, TU, QU) from rotation of (T,Q,U)
;          TT     QQ                 UU                 TQ            TU           QU                                     
bigrot = [[r00^2, 0,                 0,                 0,            0,            0                                     ], $  ; TT
          [0,     rot[1,1]^2,        rot[2,1]^2,        0,            0,            rot[1,1]*rot[2,1] + rot[2,1]*rot[1,1] ], $  ; QQ
          [0,     rot[1,2]^2,        rot[2,2]^2,        0,            0,            rot[1,2]*rot[2,2] + rot[2,2]*rot[1,2] ], $  ; UU
          [0,     0,                 0,                 r00*rot[1,1], r00*rot[2,1], 0                                     ], $  ; TQ
          [0,     0,                 0,                 r00*rot[1,2], r00*rot[2,2], 0                                     ], $  ; TU
          [0,     rot[1,1]*rot[1,2], rot[2,1]*rot[2,2], 0,            0,            rot[1,1]*rot[2,2] + rot[1,2]*rot[2,1] ] $  ; QU
         ]

print,bigrot
cora = bigrot ## cor0
;;;help,cor0,bigrot,cora
xp = theta ;*0 + 1.
tek_color

junk=min(theta, minloc)
print,transpose(cor0[minloc[0],*])
print,transpose(cor1[minloc[0],*])

var0 = dblarr(6) & var1 = var0
for i=0,2 do var0[i] = mean(tqu0[*,i]^2,/double)
var0[3] = mean(tqu0[*,0]*tqu0[*,1],/double)
var0[4] = mean(tqu0[*,0]*tqu0[*,2],/double)
var0[5] = mean(tqu0[*,1]*tqu0[*,2],/double)
print,transpose(cor0[minloc[0],*])/var0

tqu1 = tqu
for i=0,2 do var1[i] = mean(tqu1[*,i]^2,/double)
var1[3] = mean(tqu1[*,0]*tqu1[*,1],/double)
var1[4] = mean(tqu1[*,0]*tqu1[*,2],/double)
var1[5] = mean(tqu1[*,1]*tqu1[*,2],/double)
print,transpose(cor1[minloc[0],*])/var1

window,/free,xs=800,ys=900
!p.multi=[0,2,3]
titles = ['!6TT','QQ','UU','TQ','TU','QU']
for i=0,5 do begin
    fact = 1
    ; if (i eq 5) then fact = -0.5
    plot,/nodata,[0,!pi],minmax([cor0[*,i]*xp,cor1[*,i]*xp]),title=titles[i]+string(adeg),charsize=2
    oplot,theta,cor0[*,i]*xp,thick=1 ;charsize=2
    oplot,theta,cor1[*,i]*xp,thick=1,col=2 ;charsize=2
    oplot,theta,cora[*,i]*xp*fact,thick=1,col=3,lines=3
    oplot,[0,!pi],[0,0],lines=2,col=4
endfor

; window,/free,xs=800,ys=800
; !p.multi=[0,2,2]
; titles = ['!6QQ+UU','TQ+TU','TQ-TU','QU']
; for i=0,3 do begin
;     case i of
;         0: begin ; QQ + UU
;             y0 = cor0[*,1]+cor0[*,2]
;             y1 = cor1[*,1]+cor1[*,2]
;         end
;         1: begin ; TQ + TU
;             y0 = cor0[*,3]+cor0[*,4]
;             y1 = cor1[*,3]+cor1[*,4]
;         end
;         2: begin ; TQ - TU
;             y0 = cor0[*,3]-cor0[*,4]
;             y1 = cor1[*,3]-cor1[*,4]
;         end
;         3: begin ; QU
;             y0 = cor0[*,5]
;             y1 = cor1[*,5]
;         end
;     endcase
;     plot,/nodata,[0,!pi],minmax([y0*xp,y1*xp]),title=titles[i],charsize=2
;     oplot,theta,y0*xp,thick=1 ;charsize=2
;     oplot,theta,y1*xp,thick=1,col=2 ;charsize=2
;     oplot,[0,!pi],[0,0],lines=2,col=4
; endfor
;ianafast, fmap1, fcl2,/pol, /show,/silent,nlmax=lmax

!p.multi=0

end
