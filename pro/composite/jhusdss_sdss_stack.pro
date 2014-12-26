;+
; NAME:
;   jhusdss_sdss_stack()
; PURPOSE:
;   Given a list of SDSS objects with redshift, calculate the composite spectra in the rest frame
; CALLING SEQUENCE:
;   stack = jhusdss_sdss_stack(objs, weight=weight, normoption=1, /ivarweight)
; INPUTS: 
;   objs - a list of Nobj SDSS objects, including as least the following tags:
;           .PLATE
;           .FIBERID
;           .MJD
;           .Z
; OPTIONAL INPUTS:
;   weight - weight of each spectra, by default 1 for all spectra: fltarr(nobj)+1
;   normoption - normalization option. 
;                0 - no normalization, default
;                1 - median flux between 5200 AA and 5800 AA
;                2 - median flux of the whole spectrum
;                3 - mean flux of the whole spectrum in the unmasked region
;                4 - use the best fit template between 5200 AA and 5800 AA, 
;                    does not exist yet.
;                otherwise - not defined yet, no normalization for now
;   minwave - minimum wavelength in the rest frame for stacking, by default 3000 AA
;   maxwave - maximum wavelength in the rest frame for stacking, by default 8000 AA
;   normminwave - minimum wavelength in the rest frame for normalization, by default 5200 AA
;   normmaxwave - maximum wavelength in the rest frame for normalization, by default 5800 AA
;   ztag - tag name for redshift, by default 'z'
;   platetag - tag name for plate id, by default 'plate'
;   fibertag - tag name for fiber id, by default 'fiberid'
;   mjdtag - tag name for mjd, by default 'mjd'
;   path - path to the individual spectra, by default '$RAW_DATA/SDSS/Spectra'
;   loglam - output wavelength vector, in log scale, by default dloglam = 1E-4, see the code
; KEYWORDS:
;   ivarweight - weight with ivar when stacking?
;   nomoremask - unless given, by default we chuck more pixels than andmask given 
;              - by the original mask based on in ivar
;              - see the code, likely need to be revisited
; OUTPUT:
;   stack - ['nobj':nobj, 'zmean': zmean, 'wave':wave, 
;            'fmedian':fmedian, 'fmean':fmean, 'fgeomean':fgeomean]
; ROUTINES CALLED:
;   readspec -- read in individual spectra
;   combine1fiber -- Cubic B-spline interpolation
; COMMENTS:
;   -- Not generalized enough
;   -- No errors being output -- Do a jackknife or bootstrap to test for ensemble errors
;   -- A simplified version of gt_sdss_lrgstack.pro
; REVISION HISTORY:
;   10-Nov-2011  Guangtun Zhu, JHU
;-
function jhusdss_sdss_stack, objs, weight=weight, normoption=normoption, $
         ivarweight=ivarweight, minwave=minwave, maxwave=maxwave, $
         normminwave=normminwave, normmaxwave=normmaxwave, $
         ztag=ztag, platetag=platetag, fibertag=fibertag, mjdtag=mjdtag, $
         path=path, loglam=loglam, nomoremask=nomoremask

if (n_elements(objs) eq 0) then begin
   doc_libray, 'jhusdss_sdss_stack'
   return, 0
endif

nobj = n_elements(objs)

;; set default parameters
if (n_elements(weight) eq 0) then weight = fltarr(nobj)+1.
if (n_elements(normoption) eq 0) then normoption = 0
if (n_elements(minwave) eq 0) then minwave = 3000.D
if (n_elements(maxwave) eq 0) then maxwave = 8000.D
if (n_elements(normminwave) eq 0) then normminwave = 5200.D
if (n_elements(normmaxwave) eq 0) then normmaxwave = 5500.D

if (n_elements(ztag) eq 0) then ztag = 'z'
zindex = tag_indx(objs[0], ztag)
if (zindex eq -1) then message, "Z (redshift) tag doesn't exist!"
if (n_elements(platetag) eq 0) then platetag = 'plate'
pindex = tag_indx(objs[0], platetag)
if (pindex eq -1) then message, "PLATE ID tag doesn't exist!"
if (n_elements(fibertag) eq 0) then fibertag = 'fiberid'
findex = tag_indx(objs[0], fibertag)
if (findex eq -1) then message, "FIBER ID tag doesn't exist!"
if (n_elements(mjdtag) eq 0) then mjdtag = 'mjd'
mindex = tag_indx(objs[0], mjdtag)
if (mindex eq -1) then message, "MJD tag doesn't exist!"

if (n_elements(path) eq 0) then path = getenv('RAW_DATA')+'/SDSS/Spectra'

;; initialization wavelength grid, in log space
if (n_elements(loglam) eq 0) then begin
   dloglam = 1.D-4
   rloglam = (alog10(maxwave)-alog10(minwave))
   nwave = long(rloglam/dloglam)
   loglam = alog10(minwave)+(dindgen(nwave)+0.5)*dloglam
endif 
nwave = n_elements(loglam)

;; save all the projected flux into these vectors
allflux = fltarr(nobj, nwave)
allivar = fltarr(nobj, nwave)
;; good one are 0, bad ones are 0! -- as opposed to mask/andmask/ormask
allmask = bytarr(nobj, nwave)
nuse = lonarr(nwave)

;; save mean, median, and geometric mean.
fmean = fltarr(nwave)
fmedian = fltarr(nwave)
fgeomean = fltarr(nwave)

tweight = total(weight)

for i=0L, nobj-1L do begin

    print, strtrim(string(i+1),2) + '/' +strtrim(string(nobj),2)

    ;; read out the spectrum
    readspec, objs[i].(pindex), objs[i].(findex), mjd=objs[i].(mindex), path=path, $
       wave=wave, flux=flux, invvar=ivar, andmask=andmask, ormask=ormask

    ;; blueshift the spectrum back to rest frame
    wave = wave/(1.+objs[i].(zindex))
    flux = flux*(1.+objs[i].(zindex))
    ivar = ivar/(1.+objs[i].(zindex))^2

    ;; more pixels masked than andmask?
    mask = lonarr(n_elements(wave))
    if (not keyword_set(nomoremask)) then begin
        masktmp = double((ivar lt 1.e-3) or (ivar ge 1.d+30))
        masktmp = smooth(masktmp, 5)
        mask = (masktmp ge 1.e-4)
    endif
    
    usemask = (andmask or mask)

    badmask = where(usemask ne 0, comp=goodmask, nbad)
    if (nbad gt 0) then ivar[badmask] = 0.
    if (n_elements(goodmask) lt 5) then begin
       print, "I found less than 5 good pixels. I am not using this spectrum!"
       continue
    endif
    ;; normalization
    case normoption of
       0: begin
          flux_norm = 1.
          end
       ;; won't work if there are too many masked pixels
       1: begin
          inorm = where(wave ge normminwave and wave le normmaxwave, nnorm)
          flux_norm = median(smooth(flux[inorm],5))
          end
       ;; won't work if there are too many masked pixels
       2: begin
          flux_norm = median(smooth(flux, 5))
          end
       3: begin
          flux_norm = mean(flux[goodmask])
          end
       else: begin
          flux_norm = 1.
          end
    endcase
    flux = flux/flux_norm

    ;; interpolate the spectrum, combine1fiber
    ;; use cubic b-spline, because kernel width close to the psf width of sdss spectrograph ???
    ;; need to understand the output ivar better
    curr_loglam = alog10(wave)
    combine1fiber, curr_loglam, flux, ivar, newloglam=loglam, $
                   newflux=finterp, newivar=iinterp, maxiter=0, $
                   finalmask=usemask, andmask=maskinterp

    ;; Is the following necessary?
    
    allflux[i,*] = finterp*(~maskinterp)
    ;; use ~(Logical negation) instead of NOT (bitwise)
    allivar[i,*] = iinterp*(~maskinterp)
    allmask[i,*] = ~(allivar[i,*] eq 0.)
endfor

;; median
;; mask?
for i=0L, nwave-1L do fmedian[i] = median(allflux[*,i])
;; number of used object at each wavelength, allmask is in byte type but can be totalled...
for i=0L, nwave-1L do nuse[i] = total(allmask[*,i])
;stop

;; mean and geometric mean
if (not keyword_set(ivarweight)) then begin
   for i=0L, nwave-1L do begin
       fmean[i] = total(allflux[*,i]*allmask[*,i]*weight)/(total(weight*allmask[*,i]) + (total(weight*allmask[*,i]) eq 0.))
       fgeomean[i] = jhusdss_geo_mean(allflux[*,i], weight=allmask[*,i]*weight)
   endfor
endif else begin
   for i=0L, nwave-1L do begin
       fmean[i] = total(allflux[*,i]*allivar[*,i]*weight)/(total(weight*allivar[*,i]) + (total(weight*allivar[*,i]) eq 0.))
       fgeomean[i] = jhusdss_geo_mean(allflux[*,i], weight=allivar[*,i]*weight)
   endfor
endelse

stack = {nobj:nobj, zmean:mean(objs.(zindex)), wave:10^loglam, fmean:fmean, fmedian:fmedian, $
         fgeomean:fgeomean, nuse:nuse}

stop
return, stack
end
