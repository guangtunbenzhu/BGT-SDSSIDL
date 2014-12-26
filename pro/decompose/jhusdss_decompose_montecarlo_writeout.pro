;+
; Documentation needed!
;-

pro jhusdss_decompose_montecarlo_writeout, objs, nmfver=nmfver, basisfile=basisfile, $
       loglam=loglam, nmf_continuum=nmf_continuum, med_continuum=med_continuum, $
       residual=residual, flux=flux, ivar=ivar, eigen_values=eigen_values, $
       path=path, overwrite=overwrite, boss=boss

;; path
if (~keyword_set(path)) then begin
   if (keyword_set(boss)) then begin
       path=jhusdss_get_path(/nmfqso)+'/MC_'+$
          string(nmfver, format='(I3.3)')+'/Decompose_BOSS'
   endif else begin
       path=jhusdss_get_path(/nmfqso)+'/MC_'+$
          string(nmfver, format='(I3.3)')+'/Decompose'
   endelse
endif

if (jhusdss_direxist(path) eq 0) then spawn, 'mkdir -p '+path
if (n_elements(basisfile) eq 0) then $
    message, "What basis did you use?"

;if (n_elements(ztag) eq 0) then ztag = jhusdss_get_tags(/ztag)
zindex = tag_indx(objs[0], 'ZQSO')
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

strtmp = {ra:0.D, dec:0.D, plate:objs[0].(pindex), fiber:objs[0].(findex), $
          mjd:objs[0].(mindex), z:objs[0].(zindex), $
          wave:10.^loglam, nmf_continuum:fltarr((size(nmf_continuum))[2]), $
          med_continuum:fltarr((size(med_continuum))[2]), $
          residual:fltarr((size(residual))[2]), $
          flux:fltarr((size(ivar))[2]), ivar:fltarr((size(ivar))[2]), $
          basisfile:basisfile, nmfver:nmfver, eigen_values:reform(eigen_values[0,*])}

for i=0L, n_elements(objs)-1L do begin
;   zpath = path+'/'+jhusdss_get_zpath(objs[i].(zindex))
;   if (jhusdss_direxist(zpath) eq 0) then spawn, 'mkdir -p '+zpath
    outfile = jhusdss_decompose_name(objs[i].(pindex), objs[i].(findex))
    subpath = path+'/'+string(objs[i].(pindex), format='(i4.4)')
    if (jhusdss_direxist(subpath) eq 0) then spawn, 'mkdir -p '+subpath
    filename = subpath+'/'+outfile
    if (file_test(filename) eq 1) then begin
       if (~keyword_set(overwrite)) then begin
          splog, filename, ' already exists. Not overwriting ...'
          continue
       endif
    endif

    strtmp.ra = objs[i].ra
    strtmp.dec = objs[i].dec
    strtmp.plate = objs[i].(pindex)
    strtmp.fiber = objs[i].(findex)
    strtmp.mjd = objs[i].(mindex)
    strtmp.z = objs[i].(zindex)

    strtmp.nmf_continuum = nmf_continuum[i,*]
    strtmp.med_continuum = med_continuum[i,*]
    strtmp.residual = residual[i,*]
    strtmp.flux = flux[i,*]
    strtmp.ivar = ivar[i,*]
    strtmp.eigen_values = reform(eigen_values[i,*])

    mwrfits, strtmp, filename, /create
endfor

end
