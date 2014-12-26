;+
; Decompose spetra given a list of objects
; objs includes z, plate, fiber, mjd
; Please check jhusdss_get_tags() to see what tags should be used
;
; to-do: check if file exists before running the decomposition...
;-

pro jhusdss_lrg_decompose_spec, objs, nmfver=nmfver, overwrite=overwrite, $
       qaonly=qaonly, silent=silent, boss=boss

;; Let's not deal with n(objs) > 10^4 at this moment
nobjs = n_elements(objs)
if (nobjs gt 10000L) then $
    message, "Please don't be too ambitious. Give me at most 10^4 objects."

if (n_elements(nmfver) eq 0) then message, "nmfver required"

; path
if (n_elements(path) eq 0) then begin
    if (keyword_set(boss)) then begin
        path = jhusdss_get_path(/boss_spectra)
    endif else begin
        path = jhusdss_get_path(/spectra)
    endelse
endif
finalpath = path

if (~keyword_set(outpath)) then begin
   if (keyword_set(boss)) then begin
       outpath=jhusdss_get_path(/fitlrg)+'/'+$
          string(nmfver, format='(I3.3)')+'/Decompose_BOSS'
   endif else begin
       outpath=jhusdss_get_path(/fitlrg)+'/'+$
          string(nmfver, format='(I3.3)')+'/Decompose'
   endelse
endif

;; tags
if (n_elements(ztag) eq 0) then ztag = jhusdss_get_tags(/ztag)
zindex = tag_indx(objs[0], ztag)
if (zindex eq -1) then message, "Z (redshift) tag doesn't exist!"
if (n_elements(platetag) eq 0) then platetag = jhusdss_get_tags(/platetag)
pindex = tag_indx(objs[0], platetag)
if (pindex eq -1) then message, "PLATE ID tag doesn't exist!"
if (n_elements(fibertag) eq 0) then fibertag = jhusdss_get_tags(/fibertag)
findex = tag_indx(objs[0], fibertag)
if (findex eq -1) then message, "FIBER ID tag doesn't exist!"
if (n_elements(mjdtag) eq 0) then mjdtag = jhusdss_get_tags(/mjdtag)
mindex = tag_indx(objs[0], mjdtag)
if (mindex eq -1) then message, "MJD tag doesn't exist!"

splog, "Perfomring BC03 decomposition for ", nobjs, " LRGs,"
    
for i=0L, nobjs-1L do begin

    ;; check if file exists, see jhusdss_lrg_decompose_writeout.pro
    outfile = jhusdss_lrg_decompose_filename(objs[i].(pindex), objs[i].(findex), objs[i].(mindex), boss=boss)
    subpath = outpath+'/'+string(objs[i].(pindex), format='(i4.4)')
    if (jhusdss_direxist(subpath) eq 0) then spawn, 'mkdir -p '+subpath
    filename = subpath+'/'+outfile
    if (file_test(filename) eq 1) then begin
       if (~keyword_set(overwrite)) then begin
          splog, filename, ' already exists. Not overwriting ...'
          continue
       endif
    endif

    counter, i+1L, nobjs
    ;; read out the spectrum
    if (keyword_set(boss)) then finalpath = path+'/'+string(objs[i].(pindex), format='(i4.4)')
;   stop
    readspec, objs[i].(pindex), objs[i].(findex), mjd=objs[i].(mindex), path=finalpath, $
       wave=wave, flux=flux, invvar=ivar, andmask=andmask, ormask=ormask

    tmpwave = wave/(1.+objs[i].(zindex))
    jhusdss_lrg_fit, tmpwave, flux, ivar, continuum=nmf_continuum, lflux=lflux, $
          err_lflux=err_lflux, ew=ew, err_ew=err_ew, names=names, lines=lines, $
          linemodel=linemodel, sigma=sigma, err_sigma=err_sigma, $
          vdisp=objs[i].vdisp, clevel=clevel

    nmf_continuum = nmf_continuum+linemodel
    nmf_residual = flux/nmf_continuum
    tmpivar = ivar*nmf_continuum^2
    ;; lrg median filtering iteration

    mask = (ivar eq 0.)
    npix = n_elements(mask)
    filter_sizes=[91, 163]
    jhusdss_median_filter, reform(nmf_residual, 1, npix), reform(tmpivar, 1, npix), $
       mask=reform(mask, 1, npix), continuum=tmp_continuum, residual=tmp_residual, filter_sizes=filter_sizes

    filter_sizes=[143, 71]
    mask = (ivar eq 0.) or (abs(reform(tmp_continuum)-nmf_residual)*sqrt(tmpivar) gt 1.5)
    jhusdss_median_filter, reform(nmf_residual, 1, npix), reform(tmpivar, 1, npix), $
       mask=reform(mask, 1, npix),  continuum=tmp_continuum, residual=tmp_residual, filter_sizes=filter_sizes

    filter_sizes=[143, 71]
    mask = (ivar eq 0.) or (abs(reform(tmp_continuum)-nmf_residual)*sqrt(tmpivar) gt 1.5)
    jhusdss_median_filter, reform(nmf_residual, 1, npix), reform(tmpivar, 1, npix), $
       mask=reform(mask, 1, npix),  continuum=tmp_continuum, residual=tmp_residual, filter_sizes=filter_sizes
    med_continuum = reform(tmp_continuum)
    residual = reform(tmp_residual)

    loglam = alog10(tmpwave)
    jhusdss_lrg_decompose_writeout, objs[i], nmfver=nmfver, $
          loglam=loglam, nmf_continuum=reform(nmf_continuum, 1, npix), $
          med_continuum=reform(med_continuum, 1, npix), residual=reform(residual, 1, npix), $
          flux=reform(flux, 1, npix), ivar=reform(ivar, 1, npix), overwrite=overwrite, boss=boss

endfor

end
