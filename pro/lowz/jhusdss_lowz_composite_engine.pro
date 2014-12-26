;; mask out detected lines, 02/29/2012 (Leap Year!) -- Guangtun
;; jhusdss_absorbers_composite_engine
function jhusdss_lowz_composite_engine, absorbers, loglam=loglam, $
         nmfver=nmfver, boss=boss, minwave=minwave, maxwave=maxwave

common com_spec, all_residual

if (n_elements(nmfver) eq 0) then nmfver = jhusdss_get_nmf_version()
if (n_elements(minwave) eq 0)then minwave = 1100.d0
if (n_elements(maxwave) eq 0) then maxwave = 6500.d0
if (n_elements(loglam) eq 0) then loglam = jhusdss_get_loglam(minwave=minwave, maxwave=maxwave)

nabs = n_elements(absorbers)

nwave = n_elements(loglam)
allflux = fltarr(nabs, nwave)
allivar = fltarr(nabs, nwave)

;; read in all lowz
if (n_elements(all_spec) eq 0) then begin
   if (keyword_set(boss)) then begin
       path=jhusdss_get_path(/nmfqso)+'/'+$
          string(nmfver, format='(I3.3)')+'/Decompose_BOSS'
   endif else begin
       path=jhusdss_get_path(/nmfqso)+'/'+$
          string(nmfver, format='(I3.3)')+'/Decompose'
   endelse

   residual_file = path+'/'+jhusdss_decompose_allinone_filename(nmfver, boss=boss, /residual)
   all_residual = mrdfits(residual_file, 1)
   inwave = all_residual.wave
endif

for i=0L, nabs-1L do begin
    counter, i+1, nabs
    indx = absorbers[i].index
;   spec = jhusdss_decompose_loadspec(absorbers[i].plate, absorbers[i].fiber, nmfver, boss=boss)
;;  use  NMF residuals
;   ;; mask out absorption lines and recalculate residuals
;   residual = jhusdss_redecompose_spec(spec, absorbers[i].zabs, new_med_continuum=new_med_continuum)
;   spec.med_continuum = new_med_continuum
;   spec.residual = spec.flux/spec.nmf_continuum
;   outspec = jhusdss_decompose_interp_spec(spec, absorbers[i].zabs, loglam)

    ;; interpolate 
    flux = reform(all_residual.residual[*, indx])
    ivar = reform(all_residual.ivar[*, indx])

    ;; Don't know if mask works. Give it a try
    mask = (ivar eq 0.)
    curr_loglam = alog10(inwave/(1.+absorbers[i].zabs))
    combine1fiber, curr_loglam, flux, ivar, newloglam=loglam, $
                   newflux=finterp, newivar=iinterp, maxiter=0, $
                   finalmask=mask, andmask=maskinterp

    allflux[i,*]=finterp
    allivar[i,*]=iinterp
endfor

jhusdss_composite_engine, allflux, allivar, /ivarweight, fmean=fmean, fmedian=fmedian, $
   fgeomean=fgeomean, nobjuse=nobjuse

;; reject weird ones?
stack = {nabs:nabs, zabs:median(absorbers.zabs), wave:10.^loglam, fluxmean:fmean, fluxmedian:fmedian, $
         fluxgeomean:fgeomean, nobjuse:nobjuse}

return, stack
       
end
