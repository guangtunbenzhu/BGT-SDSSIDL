;+
; Obsolete?
;-
pro jhusdss_decompose_allinone, nmfver, overwrite=overwrite, $
       boss=boss

;if (~keyword_set(nmfver)) then nmfver=jhusdss_get_nmf_version()
if (n_elements(nmfver) eq 0) then $
    message, "nmf version required so that you won't overwrite files unintentionally."

if (~keyword_set(path)) then begin
   if (keyword_set(boss)) then begin
       path=jhusdss_get_path(/nmfqso)+'/'+$
          string(nmfver, format='(I3.3)')+'/Decompose_BOSS'
   endif else begin
       path=jhusdss_get_path(/nmfqso)+'/'+$
          string(nmfver, format='(I3.3)')+'/Decompose'
   endelse
endif

flux_file = path+'/'+jhusdss_decompose_allinone_filename(nmfver, boss=boss, /flux)
cont_file = path+'/'+jhusdss_decompose_allinone_filename(nmfver, boss=boss, /continuum)
residual_file = path+'/'+jhusdss_decompose_allinone_filename(nmfver, boss=boss, /residual)

if ((file_test(flux_file) or file_test(cont_file) or file_test(residual_file)) $
    and (~keyword_set(overwrite))) then begin
    splog, 'Files already exist. Not overwriting.'
    return
endif else begin
    splog, 'Will write into these files'
    splog, flux_file
    splog, cont_file
    splog, residual_file
endelse

path = jhusdss_get_path(/qso)
if (keyword_set(boss)) then begin
   filename = jhusdss_boss_qsofile()
endif else begin
   filename =  jhusdss_dr7_qsofile()
endelse

infile = path+'/'+filename
splog, 'reading '+infile
objs0 = mrdfits(infile, 1)
nspec = n_elements(objs0)

sdsswave = jhusdss_sdsswave_minmax(boss=boss)
loglam = jhusdss_get_loglam(minwave=sdsswave[0], maxwave=sdsswave[1])
nwave = n_elements(loglam)

;; should have used [nspec, nwave]!!!!!!!!!!!!
str_residual = {isitdecomposed:bytarr(nspec), wave:10^loglam, residual:fltarr(nwave, nspec), ivar:fltarr(nwave, nspec)}
str_flux = {isitdecomposed:bytarr(nspec), wave:10^loglam, flux:fltarr(nwave, nspec), ivar:fltarr(nwave, nspec)}
str_continuum = {isitdecomposed:bytarr(nspec), wave:10^loglam, nmf_continuum:fltarr(nwave, nspec), med_continuum:fltarr(nwave, nspec)}

;for i=0L, nspec-1L do begin
 for i=0L, nspec-1L do begin
    counter, i+1, nspec
    spec = jhusdss_decompose_loadspec(objs0[i].plate, objs0[i].fiber, nmfver, boss=boss, error=error)
    if (error) then continue
    
    wave = spec.wave*(1.+spec.z)
    flux = spec.flux/(1.+spec.z)
    ivar = spec.ivar*(1.+spec.z)^2
    nmf_continuum = spec.nmf_continuum/(1.+spec.z)
    med_continuum = spec.med_continuum
    residual = spec.residual
    ivar_residual = ivar*spec.med_continuum^2*spec.nmf_continuum^2

    ;; Don't know if mask works. Give it a try
    mask = (ivar eq 0.)
    curr_loglam = alog10(wave)

    combine1fiber, curr_loglam, flux, ivar, newloglam=loglam, $
                   newflux=finterp, newivar=iinterp, maxiter=0, $
                   finalmask=mask, andmask=maskinterp
    str_flux.isitdecomposed[i] = 1b
    str_flux.flux[*,i] = finterp
    str_flux.ivar[*,i] = iinterp

    combine1fiber, curr_loglam, residual, ivar_residual, newloglam=loglam, $
                   newflux=finterp, newivar=iinterp, maxiter=0, $
                   finalmask=mask, andmask=maskinterp
    str_residual.isitdecomposed[i] = 1b
    str_residual.residual[*,i] = finterp
    str_residual.ivar[*,i] = iinterp

    combine1fiber, curr_loglam, nmf_continuum, ivar, newloglam=loglam, $
                   newflux=finterp, newivar=iinterp, maxiter=0, $
                   finalmask=mask, andmask=maskinterp
    str_continuum.isitdecomposed[i] = 1b
    str_continuum.nmf_continuum[*,i] = finterp

    combine1fiber, curr_loglam, med_continuum, ivar, newloglam=loglam, $
                   newflux=finterp, newivar=iinterp, maxiter=0, $
                   finalmask=mask, andmask=maskinterp
    str_continuum.med_continuum[*,i] = finterp
endfor

mwrfits, str_flux, flux_file, /create
mwrfits, str_continuum, cont_file, /create
mwrfits, str_residual, residual_file, /create

end
