;; if you don't have to interpolate,
;; use jhusdss_composite_engine alone would be enough
function jhusdss_qso_composite_engine, inloglam, influx, inivar, z, outloglam, ivarweight=ivarweight

if (n_elements(inloglam) eq 0) then message, 'INLOGLAM required'
if (n_elements(influx) eq 0) then message, 'INFLUX required'
if (n_elements(inivar) eq 0) then message, 'INIVAR required'
if (n_elements(z) eq 0) then message, 'Z required'
if (n_elements(outloglam) eq 0) then message, 'OUTLOGLAM required'

nqso = (size(inloglam))[1]

nwave = n_elements(outloglam)
allflux = fltarr(nqso, nwave)
allivar = fltarr(nqso, nwave)

for i=0L, nqso-1L do begin
    counter, i+1L, nqso
    tmploglam = reform(inloglam[i,*])
    tmpflux = reform(influx[i,*])
    tmpivar = reform(inivar[i,*])
    tmpmask = (tmpivar le 0.)

    combine1fiber, tmploglam, tmpflux, tmpivar, newloglam=outloglam, $
       newflux=finterp, newivar=iinterp, maxiter=0, $
       finalmask=tmpmask, andmask=maskinterp
    allflux[i,*] = finterp
    allivar[i,*] = iinterp
endfor

jhusdss_composite_engine, allflux, allivar, fmean=fmean, fmedian=fmedian, $
   fgeomean=fgeomean, nobjuse=nobjuse

;; reject weird ones?
stack = {nqso:nqso, zqso:median(z), wave:10.^outloglam, fluxmean:fmean, fluxmedian:fmedian, $
         fluxgeomean:fgeomean, nobjuse:nobjuse}

return, stack

end
