;+
; load and interpolate *ALL* quasar spectra
; All in observer Frame
;; see jhusdss_decompose_allinone.pro
;-
pro jhusdss_write_allskyspec, boss=boss, overwrite=overwrite

;; qsopath
parent_path=jhusdss_get_parent_path()

skypath = jhusdss_get_path(/sky)
outfile = skypath+'/AllInOne/'+jhusdss_allskyspec_filename(boss=boss)

flux_outfile = repstr(outfile, '.fits', '_flux.fits')
;cont_outfile = repstr(outfile, '.fits', '_continuum.fits')
;normresi_outfile = repstr(outfile, '.fits', '_normalized_residual.fits')
;subtresi_outfile = repstr(outfile, '.fits', '_subtracted_residual.fits')

if (file_test(flux_outfile) and ~keyword_set(overwrite)) then begin
   splog, 'File exists. Use Overwrite to overwrite!'
   return
endif else begin
   splog, 'Will write into these files:'
   print, flux_outfile 
;  print, cont_outfile
;  print, normresi_outfile
;  print, subtresi_outfile
endelse

sky = jhusdss_sky_readin(boss=boss)
;stats = jhusdss_qsostats_readin(nmfver, boss=boss)
nsky = n_elements(sky)
strtmp = replicate({z:0.}, nsky)
sky = struct_addtags(sky, strtmp)

;zmin = [0.0, 0.4, 1.8, 2.8]
;zmax = [0.4, 1.8, 2.8, 4.8]
;; ratio = normflux[3020AA:3100AA]/normflux
;;ratio = [1./0.546970, 1., 0.685592, 0.560925*0.685592]

minmaxwave = jhusdss_sdsswave_minmax(boss=boss)
loglam = jhusdss_get_loglam(minwave=minmaxwave[0], maxwave=minmaxwave[1], nwave=nwave)
wave = 10.^loglam

flux_outstr = {wave:wave, flux:fltarr(nsky,nwave), sky:fltarr(nsky, nwave), $
               ivar:fltarr(nsky,nwave), ra:dblarr(nsky), dec:dblarr(nsky), $
               plate:lonarr(nsky), fiber:lonarr(nsky), mjd:lonarr(nsky), $
               spec_snr_median:fltarr(nsky)}

flux_outstr.ra = sky.ra & flux_outstr.dec = sky.dec & flux_outstr.plate = sky.plate & flux_outstr.fiber = sky.fiber & flux_outstr.mjd = sky.mjd 

jhusdss_load_interp_spec, sky, loglam=loglam, $
   zuse=fltarr(nsky), boss=boss, allflux=allflux, allivar=allivar
flux_outstr.flux[*,*] = allflux;/normflux
flux_outstr.ivar[*,*] = allivar;*normflux^2/ratio[i]^2

;for i=0L, 0L do begin
;   iz = where(qso.z ge zmin[i] and qso.z lt zmax[i], nz)
;   print, nz
;    iz = lindgen(nsky)
;    nz = n_elements(iz)

;   normwave = jhusdss_normwave_minmax(option=i+1)

;   iwavemin = value_locate(wave, normwave[0]*(1.+qso[iz].z))
;   iwavemax = value_locate(wave, normwave[1]*(1.+qso[iz].z))
;   normflux = fltarr(nz)
;   for j=0L, nz-1L do normflux[j] = median(allflux[j, iwavemin[j]:iwavemax[j]])
    
;   flux_outstr.norm_ratio[iz] = normflux*ratio[i]
;   subtresi_outstr.subtracted_ivar[iz,*] = allivar

;   for ispec=0L, nz-1L do begin
;       counter, ispec+1, nz
;       ;; load
;       tmpspec = jhusdss_decompose_loadspec(qso[iz[ispec]].plate, qso[iz[ispec]].fiber, nmfver, boss=boss, error=error)
;       if (error) then continue
;       tmpwave = tmpspec.wave*(1.+tmpspec.z)
;       tmpflux = tmpspec.flux
;       tmpivar = tmpspec.ivar
;       tmpivar_residual = tmpspec.ivar*tmpspec.nmf_continuum^2*tmpspec.med_continuum^2
;       tmp_nmf_continuum = tmpspec.nmf_continuum
;       tmp_med_continuum = tmpspec.med_continuum
;       tmp_residual= tmpspec.residual
;       ;; interpolate
;       tmpmask = (tmpivar le 0.)
;       curr_loglam = alog10(tmpwave)

;       ;; get the normalized flux
;       combine1fiber, curr_loglam, tmpflux, tmpivar, newloglam=loglam, $
;                  newflux=finterp, newivar=iinterp, maxiter=0, $
;                  finalmask=tmpmask, andmask=maskinterp
;       tmpratio = median(flux_outstr.flux[iz[ispec],*]/finterp)

;       combine1fiber, curr_loglam, tmp_nmf_continuum, tmpivar, newloglam=loglam, $
;                  newflux=finterp, newivar=iinterp, maxiter=0, $
;                  finalmask=tmpmask, andmask=maskinterp
;       cont_outstr.nmf_continuum[iz[ispec],*] = finterp*tmpratio

;       combine1fiber, curr_loglam, tmp_med_continuum, tmpivar_residual, newloglam=loglam, $
;                  newflux=finterp, newivar=iinterp, maxiter=0, $
;                  finalmask=tmpmask, andmask=maskinterp
;       cont_outstr.med_continuum[iz[ispec],*] = finterp

;       combine1fiber, curr_loglam, tmp_residual, tmpivar_residual, newloglam=loglam, $
;                  newflux=finterp, newivar=iinterp, maxiter=0, $
;                  finalmask=tmpmask, andmask=maskinterp
;       normresi_outstr.residual[iz[ispec],*] = finterp
;       normresi_outstr.residual_ivar[iz[ispec],*] = iinterp
;   endfor
;endfor

;cont_outstr.continuum = cont_outstr.nmf_continuum*cont_outstr.med_continuum
;subtresi_outstr.subtracted_residual = flux_outstr.flux-cont_outstr.continuum

mwrfits, flux_outstr, flux_outfile, /create
;mwrfits, cont_outstr, cont_outfile, /create
;mwrfits, normresi_outstr, normresi_outfile, /create
;mwrfits, subtresi_outstr, subtresi_outfile, /create

end
