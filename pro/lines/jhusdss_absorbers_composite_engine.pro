;; mask out detected lines, 02/29/2012 (Leap Year!) -- Guangtun
function jhusdss_absorbers_composite_engine, absorbers, loglam=loglam, $
         nmfver=nmfver, boss=boss, minwave=minwave, maxwave=maxwave

if (n_elements(nmfver) eq 0) then nmfver = jhusdss_get_nmf_version()
if (n_elements(minwave) eq 0)then minwave = 1100.d0
if (n_elements(maxwave) eq 0) then maxwave = 6500.d0
if (n_elements(loglam) eq 0) then loglam = jhusdss_get_loglam(minwave=minwave, maxwave=maxwave)

nabs = n_elements(absorbers)

nwave = n_elements(loglam)
allflux = fltarr(nabs, nwave)
allivar = fltarr(nabs, nwave)

for i=0L, nabs-1L do begin
    counter, i+1, nabs
    spec = jhusdss_decompose_loadspec(absorbers[i].plate, absorbers[i].fiber, nmfver, boss=boss)
    ;; mask out absorption lines and recalculate residuals
    residual = jhusdss_redecompose_spec(spec, absorbers[i].zabs, new_med_continuum=new_med_continuum)
    spec.med_continuum = new_med_continuum
    spec.residual = residual
    outspec = jhusdss_decompose_interp_spec(spec, absorbers[i].zabs, loglam)
    allflux[i,*]=outspec.flux
    allivar[i,*]=outspec.ivar
endfor

jhusdss_composite_engine, allflux, allivar, /ivarweight, fmean=fmean, fmedian=fmedian, $
   fgeomean=fgeomean, nobjuse=nobjuse

;; reject weird ones?
stack = {nabs:nabs, zabs:median(absorbers.zabs), wave:10.^loglam, fluxmean:fmean, fluxmedian:fmedian, $
         fluxgeomean:fgeomean, nobjuse:nobjuse}

return, stack
       
end
