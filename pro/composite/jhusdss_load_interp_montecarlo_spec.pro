;+
; NAME:
;   jhusdss_load_interp_spec
; PURPOSE:
;   Given a list of SDSS objects with redshift, load the spectra and project onto the same wavelength grid
; CALLING SEQUENCE:
;   jhusdss_load_interp_spec, objs
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
;   moremask - if given, we chuck more pixels than andmask given 
;              by the original mask based on in ivar.
;              see the code, likely need to be revisited
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
;   17-Nov-2011  Guangtun Zhu, JHU
;-
pro jhusdss_load_interp_montecarlo_spec, nmfver, objs, minwave=minwave, maxwave=maxwave, $
    normminwave=normminwave, normmaxwave=normmaxwave, normoption=normoption, $
    ztag=ztag, platetag=platetag, fibertag=fibertag, mjdtag=mjdtag, $
    path=path, loglam=loglam, moremask=moremask, $
    allflux=allflux, allivar=allivar, allmask=allmask, boss=boss, zuse=zuse

if (n_elements(nmfver) eq 0) then begin
   doc_library, 'jhusdss_load_interp_montecarlo_spec'
   return
endif

if (n_elements(objs) eq 0) then begin
   doc_library, 'jhusdss_load_interp_montecarlo_spec'
   return
endif

nobj = n_elements(objs)

;; set default parameters
if (n_elements(normoption) eq 0) then normoption = 0
if (normoption eq 1) then begin
   if (n_elements(normminwave) eq 0 or n_elements(normmaxwave) eq 0) then $
       message, "You have to provide the wavelength range for normalization!"
endif

;; tags
;if (n_elements(ztag) eq 0) then ztag = jhusdss_get_tags(/ztag)
;zindex = tag_indx(objs[0], ztag)
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

; path
if (n_elements(path) eq 0) then begin
    if (keyword_set(boss)) then begin
        path = jhusdss_get_path(/boss_spectra)
    endif else begin
        path = jhusdss_get_path(/spectra)
    endelse
endif
finalpath = path
mc_path = jhusdss_get_path(/nmfqso)+'/MC_'+string(nmfver, format='(I3.3)')+'/Flux'

;; initialization wavelength grid, in log space

if (n_elements(loglam) eq 0) then begin
   if (n_elements(minwave) eq 0) then minwave = 3000.D
   if (n_elements(maxwave) eq 0) then maxwave = 8000.D
   loglam = jhusdss_get_loglam(minwave, maxwave, nwave=nwave)
endif
nwave = n_elements(loglam)

;; save all the projected flux into these vectors
allflux = fltarr(nobj, nwave)
allivar = fltarr(nobj, nwave)
;; good one are 0, bad ones are 1! -- as opposed to mask/andmask/ormask
allmask = bytarr(nobj, nwave)+1b

for i=0L, nobj-1L do begin

    if (i mod 40 eq 0) then print, '  ' + strtrim(string(i+1),2) + '/' +strtrim(string(nobj),2)

    ;; read out the spectrum
    if (keyword_set(boss)) then finalpath = path+'/'+string(objs[i].(pindex), format='(i4.4)')
;   stop
    readspec, objs[i].(pindex), objs[i].(findex), mjd=objs[i].(mindex), path=finalpath, $
       wave=wave, flux=flux, invvar=ivar, andmask=andmask, ormask=ormask

    ;; get montecarlo flux
    mc_file = mc_path+'/MC_flux_'+string(objs[i].(pindex), format='(i4.4)')+'_'+string(objs[i].(findex), format='(i3.3)')+'_'+string(objs[i].(mindex), format='(i5.5)')+'.fits'
    mc = mrdfits(mc_file, 1)
    flux = mc.flux

    ;; blueshift the spectrum back to rest frame
    if (n_elements(zuse) gt 0) then begin
       wave = wave/(1.+zuse[i])
       flux = flux*(1.+zuse[i])
       ivar = ivar/(1.+zuse[i])^2
    endif else begin
       wave = wave/(1.+objs[i].(zindex))
       flux = flux*(1.+objs[i].(zindex))
       ivar = ivar/(1.+objs[i].(zindex))^2
    endelse

    ;; I don't think this actually works. Need to revisit this
    ;; more pixels masked than andmask?
    usemask = lonarr(n_elements(wave))
    ;; 01/20/2012, kill andmasks, don't think we need it
    ;; e.g., andmasks bit 23 mean bright sky, it's not useful
    if (keyword_set(applymask)) then begin 
        mask = lonarr(n_elements(wave))
        if (keyword_set(moremask)) then begin
            masktmp = double((ivar lt 1.e-3) or (ivar ge 1.d+30)) ;;hardcoded
            masktmp = smooth(masktmp, 5)
            mask = (masktmp ge 1.e-4)
        endif
        usemask = (andmask or mask)

        badmask = where(usemask ne 0, nbad)
        if (nbad gt 0) then ivar[badmask] = 0.
    endif

    goodmask = where(ivar gt 1.e-3)
    if (n_elements(goodmask) lt 5) then begin
       splog, "I found less than 5 good pixels. I am not using this spectrum!"
       continue
    endif

    ;; normalization
    donorm = 1
    case normoption of
       0: begin
          flux_norm = 1.
          end
       ;; won't work if there are too many masked pixels
       1: begin
          inorm = where((wave ge normminwave) and (wave le normmaxwave) and (ivar gt 1.e-3), nnorm)
          if (nnorm le 5) then begin
             splog, "Normalization pixels less than 5, not using this spectrum!"
             donorm = 0
          endif else begin
             flux_norm = median(flux[inorm])
          endelse
          end
       ;; won't work if there are too many masked pixels
       2: begin
          flux_norm = median(flux[goodmask])
          end
       3: begin
          flux_norm = mean(flux[goodmask])
          end
       else: begin
          flux_norm = 1.
          end
    endcase
    if (donorm eq 0) then continue
    flux = flux/flux_norm
    ivar = ivar*flux_norm^2

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
    allmask[i,*] = (allivar[i,*] eq 0.)
endfor

return 
end
