pro paramfile2thincscript, paramfile, module, thinc_script, omp=omp

jobname = 'myJob'
nickname = module

openw, lunit2, thinc_script, /get_lun
printf, lunit2, 'from thinc import *'
printf, lunit2, ' '
printf, lunit2, 'BeginPipe()'
printf, lunit2, jobname+' = NewJob("'+module+'", label = "'+nickname+'", cast = True)'

; nl = numlines(paramfile) for IDL before 5.6
nl = file_lines(paramfile)
openr, lunit1, paramfile,    /get_lun
line_in = ''
line_out = ''
for i=0, nl-1 do begin
    readf, lunit1, line_in
    line_in = strtrim(line_in,2)
    if (line_in ne '' && strmid(line_in,0,1) ne '#' ) then begin
        parse = str_sep(line_in,'=')
        kwd = strtrim(parse[0],2)
        val = strtrim(parse[1],2)

        line_out = 'Set('+jobname+', "'+kwd+'", "'+val+'")'
        ;;;;;;print, thinc_script, line_out
        printf, lunit2, line_out
    endif
endfor

extra = ''
if keyword_set(omp) then begin
    extra += ', exportEnviron={"OMP_NUM_THREADS":"'+strtrim(omp,2)+'"}'
endif
printf,lunit2, 'Submit('+jobname+extra+')'
printf,lunit2, 'EndPipe()'
free_lun, lunit2

return
end
;===========================================================
function spice_yes_no, kwd
; kwd not set -> NO
; kwd="NO"    -> NO
; kwd="YES"   -> YES
if keyword_set(kwd) then begin ; kwd set, test for NO
    yes_no = (strupcase(kwd) eq 'NO')? 'NO' : 'YES'
endif else begin ; kwd not set -> NO
    yes_no = 'NO'
endelse
return, yes_no
end
;===========================================================

pro ispice, mapfile1, clfile $
            , alm1fileout        =  alm1fileout_usr  $
            , alm2fileout        =  alm2fileout_usr  $
            , apodizesigma	 =  apodizesigma $
            , apodizetype	 =  apodizetype $
            , binpath            = binpath $
            , fwhm1		 =  fwhm1 $
            , beam_file1	 =  beam_file1 $
            , fwhm2		 =  fwhm2 $
            , beam_file2	 =  beam_file2 $
            , cl_outmap_file	 =  cl_outmap_file $
            , cl_inmap_file	 =  cl_inmap_file $
            , cl_outmask_file	 =  cl_outmask_file $
            , cl_inmask_file	 =  cl_inmask_file $
            , corfile		 =  corfile $
            , covfileout	 =  covfileout $
            , decouple		 =  decou_usr $
;            , dry		 =  dry $
            , fits_out		 =  fits_out_usr $
            , dmc                = dmc $
            , help               = help $
            , keep_tmp_files     = keep_tmp_files $
            , kernelsfileout	 =  kernelsfileout $
            , lmw1               = lmw1 $  
            , lmw2               = lmw2 $  
            , mapfile2		 =  mapfile2 $
            , maskfile1		 =  maskfile1 $
            , maskfile2		 =  maskfile2 $
            , maskfilep1	 =  maskfilep1 $
            , maskfilep2	 =  maskfilep2 $
            , nlmax	         =  nlmax $
            , normfac		 =  normfac $
            , npairsthreshold	 =  npairsthreshold $
            , noisecorfile	 =  noisecorfile $
            , noiseclfile	 =  noiseclfile $
;            , overwrite		 =  overwrite $
            , pixelfile		 =  pixelfile_usr $
;             , pixelfile1         =  pixelfile1_usr $
;             , pixelfile2	 =  pixelfile2_usr $
            , pf2                = pixelfile2_usr $
            , polarization	 =  polar_usr $
            , show_cl            = show_cl $
            , silent             = silent $
            , subav		 =  subav_usr $
            , subdipole		 =  subdipole_usr $
            , symmetric_cl       =  symm_usr $
            , tenormfilein	 =  tenormfilein $
            , tenormfileout	 =  tenormfileout $
            , tf_file		 =  tf_file $
            , thetamax		 =  thetamax $
            , tmpdir             =  tmpdir $
            , tolerance		 =  tolerance $
;;            , verbosity		 =  verbosity $
            , weightfile1	 =  weightfile1 $
            , weightfilep1	 =  weightfilep1 $
            , weightfile2	 =  weightfile2 $
            , weightfilep2	 =  weightfilep2 $
            , weightpower1	 =  weightpower1 $
            , weightpowerp1	 =  weightpowerp1 $
            , weightpower2	 =  weightpower2 $
            , weightpowerp2	 =  weightpowerp2 $
            , windowfilein	 =  windowfilein $
            , windowfileout	 =  windowfileout $
            , xtramapfile1       =  xtramapfile1 $
            , xtramapfile2       =  xtramapfile2 
;+
; NAME:
;        ISPICE
;
; PURPOSE:
;        Interface to Spice C(l) (angular power spectrum) estimator
;
; CATEGORY:
;
;
; CALLING SEQUENCE:
;        ispice, map [ clfile , 
;    alm1fileout=, alm2fileout=, apodizesigma=, apodizetype=, binpath=, dmc=, fwhm1=,
;    beam_file1=, fwhm2=, beam_file2=, 
;    cl_outmap_file=, cl_inmap_file=, cl_outmask_file=, cl_inmask_file=, corfile=, covfileout=
;    decouple=, dmc=, fits_out =, help=, keep_tmp_files =, kernelsfileout =,    lmw1=, lmw2=, 
;    mapfile2=, maskfile1=, maskfile2=, maskfilep1=, maskfilep2=, 
;    nlmax=, normfac=, npairsthreshold=, noisecorfile=, noiseclfile=, 
;    pixelfile=, pf2=, polarization=, show_cl=, silent=, subav=, subdipole=, symmetric_cl=,
;    tenormfilein=, tenormfileout=,
;    tf_file=, thetamax=, tmpdir=, tolerance=, 
;    weightfile1=, weightfilep1=, weightfile2=, weightfilep2=, 
;    weightpower1=, weightpowerp1=, weightpower2=, weightpowerp2=, windowfilein=, windowfileout=,
;    xtramapfile1=, xtramapfile2= ]
;
;
;     The inputs and keywords are very similar to those of the F90 Spice code and described in 
;   http://www2.iap.fr/users/hivon/software/PolSpice/README.html
;   or via     spice -help
;   Only major discrepancies are detailed below.
;
; INPUTS:
;
;    map:    single map file, 
;            or array of multiple map files, to be added together with the optional weights LMW1
;            Spice will compute the spectra of this map, or its cross-spectra with mapfile2
;
;
; OPTIONAL INPUTS:
;
;    clfile: Memory array or FITS or ASCII file that will contain the computed auto or cross C(l)
;
;
; KEYWORD PARAMETERS:
;
;    binary keywords can be set by "YES" or 1, and unset by "NO" or 0.
;    These include: DECOUPLE, FITS_OUT, POLARIZATION, SUBAV, SUBDIPOLE, SYMMETRIC_CL
;
;    ALM1FILEOUT: FITS file in which to dump alm of 1st map
;
;    ALM2FILEOUT: FITS file in which to dump alm of 2nd map (if applicable)
;
;    BINPATH: full path to spice executable
;
;    DMC: to work in Planck Data Management environment. Soon to be deprecated.
;
;    FITS_OUT: output files are in FITS format instead of plain ASCII. Turned on by default
;
;    FWHM1, FWHM2: beam FHWM for maps #1 and #2 (interface to Spice's "beam" and "beam2" keywords)
;
;    HELP: prints out this help header and exits
;
;    KEEP_TMP_FILES: keeps temporary files instead of removing them when leaving.
;
;    LMW1, LMW2: map weights (interface to Spice's "listmapweights1_*" and "listmapweights2_*")
;
;    MAPFILE2:  single map file, 
;               or array of multiple map files, to be added together with the optional weights LMW2
;
;    PIXELFILE, PF2: pixel files for maps #1 and #2 (interface to Spice's "pixelfile" and "pixelfile2").
;      Must be set explicitely to "NO" or 0 to prevent correction of final C(l) 
;        from relevant pixel window function(s).
;
;    SHOW_CL: if set, the resulting C(l) is plotted.
;
;    SILENT: if set, the code runs silently. Otherwise use Spice's verbosity=2 setting.
;
;    TMPDIR:  directory in which are written temporary files 
;         [default: IDL_TMPDIR (see IDL documentation about IDL_TMPDIR)]
;
;    XTRAMAPFILE, XTRAMAPFILE2: addition map files, made redundant by map and MAPFILE2 now accepting 
;          arrays of maps
;
;    see   spice -help 
;    for more information on Spice native keywords
;
;
; SIDE EFFECTS:
;     Runs the F90 code 'spice' found locally or in BINPATH.
;     Creates temporay files in IDL_TMPDIR and removes them when leaving.
;
; RESTRICTIONS:
;     Requires HEALPIX (F90 and IDL) and Spice
;
; PROCEDURE:
;     See side effects
;
; EXAMPLE:
;
;
;
; MODIFICATION HISTORY:
;
;  creation: 2008-03-02, EH, IAP
;  2009-01-27: outputs kernels
;  2009-07-23: pixelfile='NO' will not correct for pixel window 
;  2010-01-08: bug correction of the above
;  2010-01-11: added subdipole, 
;              can process Planck-HFI data in DMC format (DMC= keyword)
;  2010-05: added covfileout
;  2012-03: added xtramapfile1 and xtramapfile2
;  2014-04-01: support linear combination of input maps on the fly:
;     allow map1_in and mapfile2 to be vector of FITS files or objects
;     added lmw1 and lmw2 (wrappers for listmapweights1_* and listmapweights2_)
; 2014-04-16: keywords DECOUPLE, POLARIZATION, SUBAV, SUBDIPOLE, SYMMETRIC_CL  now accept
;            "YES" and "NO" values (as well as 1 and 0)
; 2015-03-02: bug correction on lmw1 and lmw2
; 2015-08-26: added TOLERANCE; 
;             set explicitely FITS_OUT=0 to prevent FITS output
; 2016-12-02: added PF2 (in order to support cross-analysis of maps with different Nsides)
;             and TMPDIR
; 2020-09-18: added -alm1fileout and -alm2fileout do dump alm
;
;-

do_dmc = keyword_set(dmc)
options = do_dmc ? '' : '-optinfile'
local = {routine: 'ispice', exe: 'spice', options: options}
syntax = [local.routine+', map1_in, [clfile ,' $
,'    alm1fileout=, alm2fileout=, apodizesigma=, apodizetype=, binpath=, dmc=, fwhm1=,' $
,'    beam_file1=, fwhm2=, beam_file2=, ' $
,'    cl_outmap_file=, cl_inmap_file=, cl_outmask_file=, cl_inmask_file=, corfile=, covfileout=,' $
,'    decouple=, fits_out=, help=, keep_tmp_files =, kernelsfileout =,       lmw1 =, lmw2 =, ' $
,'    mapfile2=, maskfile1=, maskfile2=, maskfilep1=, maskfilep2=, ' $
,'    nlmax=, normfac=, npairsthreshold=, noisecorfile=, noiseclfile=, ' $
,'    pixelfile=, polarization=, show_cl=, silent=, subav=, subdipole=, symmetric_cl=, ' $
,'    tenormfilein=, tenormfileout=, ' $
,'    tf_file=, thetamax=, tolerance=, weightfile1=, weightfilep1=, weightfile2=, weightfilep2=, ' $
,'    weightpower1=, weightpowerp1=, weightpower2=, weightpowerp2=, windowfilein=, windowfileout=' $
,'    xtramapfile1=, xtramapfile2= ]' ]

if keyword_set(help) then begin
    doc_library,local.routine
    return
endif

if (n_params() eq 0) then begin
    print,syntax,form='(a)'
    print
    print,local.routine+',/help     for extended help'
    return
endif
if (n_params() gt 2) then begin
    print,syntax,form='(a)'
    message,'wrong number of arguments'
endif

if (size(mapfile1,/tname) ne 'STRING') then begin
    message,'First argument (input map) must be a FITS file',/info
    print
    print,syntax
    return
endif
multimap1 = (n_elements(mapfile1) gt 1)
multimap2 = (n_elements(mapfile2) gt 1)

; polarization = keyword_set(polar_usr) ? 'YES' : 'NO'
; decouple     = keyword_set(decou_usr) ? 'YES' : 'NO'
; symmetric_cl = keyword_set(symm_usr)  ? 'YES' : 'NO'
; subav        = keyword_set(subav_usr) ? 'YES' : 'NO'
; subdipole    = keyword_set(subdipole_usr) ? 'YES' : 'NO'
polarization = spice_yes_no(polar_usr) 
decouple     = spice_yes_no(decou_usr) 
symmetric_cl = spice_yes_no(symm_usr)  
subav        = spice_yes_no(subav_usr) 
subdipole    = spice_yes_no(subdipole_usr) 
pixelfile    = (keyword_set(pixelfile_usr) || ~(keyword_set(noisecorfile) || keyword_set(noiseclfile)) ) ? 'YES' : 'NO'
if (defined(pixelfile_usr) && spice_yes_no(pixelfile_usr) eq 'NO') then pixelfile = 'NO'
; pixelfile1    = (keyword_set(pixelfile1_usr) || ~(keyword_set(noisecorfile) || keyword_set(noiseclfile)) ) ? 'YES' : 'NO'
; if (defined(pixelfile1_usr) && strupcase(pixelfile1_usr) eq 'NO') then pixelfile1 = 'NO'
pixelfile2    = (keyword_set(pixelfile2_usr) || ~(keyword_set(noisecorfile) || keyword_set(noiseclfile)) ) ? 'YES' : 'NO'
if (defined(pixelfile2_usr) && spice_yes_no(pixelfile2_usr) eq 'NO') then pixelfile2 = 'NO'
fits_out     = defined(fits_out_usr) ? spice_yes_no(fits_out_usr) : "YES"

alm1fileout = size(alm1fileout_usr,/tname) eq 'STRING' ? alm1fileout_usr : spice_yes_no(alm1fileout_usr)
alm2fileout = size(alm2fileout_usr,/tname) eq 'STRING' ? alm2fileout_usr : spice_yes_no(alm2fileout_usr)
;-------------------
hpx_xface_generic, fullpath, tmp_par_file, binpath, init=local, tmpdir=tmpdir
NoFile = keyword_set(cxx) ? " " : " '' "

; deal with online data
; tmp_clfile   = hpx_mem2file((arg_present(clfile) || defined(clfile)) ? clfile : NoFile, /out)
tmp_clfile   = hpx_mem2file((arg_present(clfile) || defined(clfile)) ? clfile : 'NO', /out)

; writes parameter file
openw,lunit,tmp_par_file, /get_lun
printf,lunit,'# parameter file for IDL interface to '+fullpath
printf,lunit,'# written: '+systime()+' by '+local.routine
printf,lunit,' '

printf,lunit,hpx_add_parameter('alm1fileout',	alm1fileout,    /skip_if_not_set, default="NO")
printf,lunit,hpx_add_parameter('alm2fileout',	alm2fileout,    /skip_if_not_set, default="NO")
printf,lunit,hpx_add_parameter('apodizesigma',	apodizesigma,	skip_if_not_set=~do_dmc, default="NO") ; let code choose default
printf,lunit,hpx_add_parameter('apodizetype' ,	apodizetype,	skip_if_not_set=~do_dmc, default=0)
printf,lunit,hpx_add_parameter('beam',          fwhm1   ,	skip_if_not_set=~do_dmc, default="NO")
printf,lunit,hpx_add_parameter('beam_file',	beam_file1,	/skip_if_not_set) 
printf,lunit,hpx_add_parameter('beam2',         fwhm2,   	skip_if_not_set=~do_dmc, default="NO")
printf,lunit,hpx_add_parameter('beam_file2',    beam_file2,	/skip_if_not_set)
printf,lunit,hpx_add_parameter('clfile',        tmp_clfile,  	/skip_if_not_set)
printf,lunit,hpx_add_parameter('cl_outmap_file',  cl_outmap_file,	      /skip_if_not_set)
printf,lunit,hpx_add_parameter('cl_inmap_file',   cl_inmap_file,	      /skip_if_not_set)
printf,lunit,hpx_add_parameter('cl_outmask_file', cl_outmask_file,	      /skip_if_not_set)
printf,lunit,hpx_add_parameter('cl_inmask_file',  cl_inmask_file,	      /skip_if_not_set)
printf,lunit,hpx_add_parameter('corfile',	  corfile,  	              /skip_if_not_set)
printf,lunit,hpx_add_parameter('covfileout',      covfileout,	              /skip_if_not_set)
printf,lunit,hpx_add_parameter('decouple',	  decouple,     skip_if_not_set=~do_dmc, default="NO")
if (~do_dmc) then begin
    printf,lunit,hpx_add_parameter('dry',   	  'NO')
    ;printf,lunit,hpx_add_parameter('fits_out',	  'YES')
    printf,lunit,hpx_add_parameter('fits_out',	  fits_out)
endif
printf,lunit,hpx_add_parameter('kernelsfileout',  kernelsfileout,	     /skip_if_not_set)
;
if (multimap1) then begin
    if (do_dmc && keyword_set(polar_usr)) then begin
        if (n_elements(mapfile1) mod 3 ne 0) then message,'Expected (I,Q,U) triplet in map'
        for i=0,n_elements(mapfile1)-1,3 do begin
            k = i/3+1
            printf,lunit,hpx_add_parameter('listmapfiles1_' +strtrim(k,2), mapfile1[i],   /expand)
            printf,lunit,hpx_add_parameter('listmapfilesQ1_'+strtrim(k,2), mapfile1[i+1], /expand)
            printf,lunit,hpx_add_parameter('listmapfilesU1_'+strtrim(k,2), mapfile1[i+2], /expand)
        endfor        
    endif else begin
        for i=0,n_elements(mapfile1)-1 do begin
            ;print,'listmapfiles1_'+strtrim(i+1,2)+' = '+mapfile1[i]
            printf,lunit,hpx_add_parameter('listmapfiles1_'+strtrim(i+1,2), mapfile1[i], /expand)
        endfor
    endelse
endif else begin
    printf,lunit,hpx_add_parameter('mapfile',         mapfile1,      /expand)
    printf,lunit,hpx_add_parameter('extramapfile',    xtramapfile1,          /skip_if_not_set)
endelse
if (multimap2) then begin
    if (do_dmc && keyword_set(polar_usr)) then begin
        if (n_elements(mapfile2) mod 3 ne 0) then message,'Expected (I,Q,U) triplet in mapfile2'
        for i=0,n_elements(mapfile2)-1,3 do begin
            k = i/3+1
            printf,lunit,hpx_add_parameter('listmapfiles2_' +strtrim(k,2), mapfile2[i],   /expand)
            printf,lunit,hpx_add_parameter('listmapfilesQ2_'+strtrim(k,2), mapfile2[i+1], /expand)
            printf,lunit,hpx_add_parameter('listmapfilesU2_'+strtrim(k,2), mapfile2[i+2], /expand)
        endfor        
    endif else begin
        for i=0,n_elements(mapfile2)-1 do begin
            ;print,'listmapfiles2_'+strtrim(i+1,2)+' = '+mapfile2[i]
            printf,lunit,hpx_add_parameter('listmapfiles2_'+strtrim(i+1,2), mapfile2[i], /expand)
        endfor
    endelse
endif else begin
    printf,lunit,hpx_add_parameter('mapfile2',	      mapfile2,		     /skip_if_not_set)
    printf,lunit,hpx_add_parameter('extramapfile2',   xtramapfile2,	     /skip_if_not_set)
endelse
if defined(lmw1) then begin
;    for i=0, n_elements(lmw1)-1 do printf,lunit,hpx_add_parameter('listmapweights1_'+strtrim(i+1,2), lmw1[i+1],   /skip_if_not_set, default=1.0)
    for i=0, n_elements(lmw1)-1 do printf,lunit,hpx_add_parameter('listmapweights1_'+strtrim(i+1,2), lmw1[i],   /skip_if_not_set, default=1.0)
endif
if defined(lmw2) then begin
;    for i=0, n_elements(lmw2)-1 do printf,lunit,hpx_add_parameter('listmapweights2_'+strtrim(i+1,2), lmw2[i+1],   /skip_if_not_set, default=1.0)
    for i=0, n_elements(lmw2)-1 do printf,lunit,hpx_add_parameter('listmapweights2_'+strtrim(i+1,2), lmw2[i],   /skip_if_not_set, default=1.0)
endif
;
printf,lunit,hpx_add_parameter('maskfile',	  maskfile1,		     /skip_if_not_set)
printf,lunit,hpx_add_parameter('maskfile2',	  maskfile2,		     /skip_if_not_set)
printf,lunit,hpx_add_parameter('maskfilep',	  maskfilep1,	             /skip_if_not_set)
printf,lunit,hpx_add_parameter('maskfilep2',	  maskfilep2,	             /skip_if_not_set)
printf,lunit,hpx_add_parameter('nlmax',		  nlmax,	   skip_if_not_set=~do_dmc, default=-1)
printf,lunit,hpx_add_parameter('normfac',	  normfac,         skip_if_not_set=~do_dmc, default=1.0)
;printf,lunit,hpx_add_parameter('npairsthreshold', npairsthreshold, skip_if_not_set=~do_dmc, default="NO")
printf,lunit,hpx_add_parameter('npairsthreshold', npairsthreshold, skip_if_not_set=~do_dmc, default=0.0)
printf,lunit,hpx_add_parameter('noisecorfile',	  noisecorfile,	     /skip_if_not_set)
printf,lunit,hpx_add_parameter('noiseclfile',	  noiseclfile,	     /skip_if_not_set)
printf,lunit,hpx_add_parameter('overwrite',	  'YES')
printf,lunit,hpx_add_parameter('pixelfile',	  pixelfile,     skip_if_not_set=~do_dmc, default="YES")
; printf,lunit,hpx_add_parameter('pixelfile1',	  pixelfile1,     skip_if_not_set=~do_dmc, default="YES")
printf,lunit,hpx_add_parameter('pixelfile2',	  pixelfile2,     skip_if_not_set=~do_dmc, default="YES")
printf,lunit,hpx_add_parameter('polarization',	  polarization,	     /skip_if_not_set)
printf,lunit,hpx_add_parameter('subav',		  subav,         skip_if_not_set=~do_dmc, default="NO")
printf,lunit,hpx_add_parameter('subdipole',	  subdipole,	 skip_if_not_set=~do_dmc, default="NO")
printf,lunit,hpx_add_parameter('symmetric_cl',	  symmetric_cl,  skip_if_not_set=~do_dmc, default="NO")
printf,lunit,hpx_add_parameter('tenormfilein',	  tenormfilein,	     /skip_if_not_set)
printf,lunit,hpx_add_parameter('tenormfileout',   tenormfileout,     /skip_if_not_set)
printf,lunit,hpx_add_parameter('tf_file',	  tf_file,		     /skip_if_not_set)
printf,lunit,hpx_add_parameter('thetamax',	  thetamax,      skip_if_not_set=~do_dmc, default="NO")
printf,lunit,hpx_add_parameter('tolerance',	  tolerance,     skip_if_not_set=~do_dmc, default="NO")
printf,lunit,hpx_add_parameter('verbosity',	  keyword_set(silent)?0:2)
printf,lunit,hpx_add_parameter('weightfile',	  weightfile1,	     /skip_if_not_set)
printf,lunit,hpx_add_parameter('weightfile2',	  weightfile2,	     /skip_if_not_set)
printf,lunit,hpx_add_parameter('weightfilep',	  weightfilep1,	     /skip_if_not_set)
printf,lunit,hpx_add_parameter('weightfilep2',	  weightfilep2,	     /skip_if_not_set)
printf,lunit,hpx_add_parameter('weightpower',	  weightpower1,   skip_if_not_set=~do_dmc, default=1.0)
printf,lunit,hpx_add_parameter('weightpower2',	  weightpower2,   skip_if_not_set=~do_dmc, default=1.0)
printf,lunit,hpx_add_parameter('weightpowerp',	  weightpowerp1,  skip_if_not_set=~do_dmc, default=1.0)
printf,lunit,hpx_add_parameter('weightpowerp2',   weightpowerp2,  skip_if_not_set=~do_dmc, default=1.0)
printf,lunit,hpx_add_parameter('windowfilein',	  windowfilein,	     /skip_if_not_set)
printf,lunit,hpx_add_parameter('windowfileout',   windowfileout,     /skip_if_not_set)
free_lun, lunit

; execute command
if (do_dmc) then begin
    ;thinc_file = '/tmp/hivon_thinc.py'
    thinc_file = getenv('IDL_TMPDIR')+'spice_from_idl_'+string(long(systime(1)*100 mod 1.e8), form='(i8.8)')+'.py'
    paramfile2thincscript, tmp_par_file, 'spice', thinc_file, omp=8
    spawn, 'python '+thinc_file
endif else begin
    hpx_xface_generic, /run, fullpath, tmp_par_file, silent=silent
endelse

; deal with online data
if (arg_present(clfile) || defined(clfile)) then hpx_file2mem, tmp_clfile, clfile, /cl, show_cl = show_cl

; to_remove
hpx_xface_generic, clean = ~keyword_set(keep_tmp_files)


return
end
