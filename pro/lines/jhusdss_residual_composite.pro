pro jhusdss_residual_composite, nmfver, boss=boss, overwrite=overwrite

if (n_elements(nmfver) eq 0) then message, 'nmfver required'

;; output 
if (~keyword_set(boss)) then begin
   path=jhusdss_get_path(/nmfqso)+'/'+string(nmfver, format='(I3.3)')+'/Composite'
endif else begin
   path=jhusdss_get_path(/nmfqso)+'/'+string(nmfver, format='(I3.3)')+'/Composite_BOSS'
endelse
if (jhusdss_direxist(path) eq 0) then message, "Can't find the directory."

if (~keyword_set(boss)) then begin
   outfile = path+'/Residual_Composite_'+string(nmfver, format='(i3.3)') + '.fits'
endif else begin
   outfile = path+'/Residual_Composite_'+string(nmfver, format='(i3.3)') + '_BOSS.fits'
endelse

if (file_test(outfile) and ~keyword_set(overwrite)) then begin
   splog, 'File already exists. Use /overwrite if you want to overwrite it.'
   return
endif else begin
   splog, 'Will write the absorber catalog into this file: '
   print, outfile
endelse

qsopath = jhusdss_get_path(/qso)
if (keyword_set(boss)) then begin
   infile = jhusdss_boss_qsofile()
endif else begin
   infile =  jhusdss_dr7_qsofile()
endelse

mgiifile = qsopath+'/'+infile
splog, 'reading '+mgiifile
mgii = mrdfits(mgiifile, 1)
nspec = n_elements(mgii)

;; randomly select ~10^4 objects (there may be repetition)
;nran = 32000L
;iran = floor(randomu(seed, nran)*nspec)
;iran = iran[uniq(iran[bsort(iran)])]
;nran = n_elements(iran)
nran = nspec
iran = lindgen(nran)

;; wavelength grid
minwave=3800d0
maxwave=9200d0
loglam = jhusdss_get_loglam(minwave=minwave, maxwave=maxwave)
nwave = n_elements(loglam)

allflux = fltarr(nran, nwave)
allivar = fltarr(nran, nwave)

for i=0L, nran-1L do begin
    counter, i+1, nran
    is = iran[i]
    spec = jhusdss_decompose_loadspec(mgii[is].plate, mgii[is].fiber, nmfver, boss=boss, error=error)
    if (error) then begin
       splog, "Can't find the decomposed continuum"
       continue
    endif
    wave = spec.wave*(1.+spec.z)
    flux = spec.residual
    ivar = spec.ivar*spec.med_continuum^2*spec.nmf_continuum^2
    mask = (ivar eq 0.)
    curr_loglam = alog10(wave)
    combine1fiber, curr_loglam, flux, ivar, newloglam=loglam, $
       newflux=finterp, newivar=iinterp, maxiter=0, $
       finalmask=mask, andmask=maskinterp

    allflux[i,*] = finterp
    allivar[i,*] = iinterp
endfor

jhusdss_composite_engine, allflux, allivar, /ivarweight, fmean=fmean, fmedian=fmedian, $
   fgeomean=fgeomean, nobjuse=nobjuse

stack = {wave:10.^loglam, fluxmean:fmean, fluxmedian:fmedian, $
         fluxgeomean:fgeomean, nobjuse:nobjuse}

mwrfits, stack, outfile, /create

end
