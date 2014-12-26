;+
; load and interpolate *ALL* quasar spectra
; All in observer Frame
;-
pro jhusdss_write_alllrgspec, lrgver, boss=boss, overwrite=overwrite

if (n_elements(lrgver) eq 0) then message, 'lrgver required'

;; lrgpath
lrgpath=jhusdss_get_path(/fitlrg)+'/'+string(lrgver, format='(I3.3)')+'/'
;outfile = lrgpath+'/'+jhusdss_alllrgspec_filename(lrgver, boss=boss)
flux_outfile = lrgpath+'/Flux_'+jhusdss_alllrgspec_filename(lrgver, boss=boss)
continuum_outfile = lrgpath+'/Continuum_'+jhusdss_alllrgspec_filename(lrgver, boss=boss)
residual_outfile = lrgpath+'/Residual_'+jhusdss_alllrgspec_filename(lrgver, boss=boss)

if (file_test(outfile) and ~keyword_set(overwrite)) then begin
   splog, 'File exists. Use Overwrite to overwrite!'
   return
endif else begin
   splog, 'Will write into this file:'
   splog, flux_outfile
   splog, continuum_outfile
   splog, residual_outfile
endelse

lrg = jhusdss_lrg_readin(boss=boss)
nlrg = n_elements(lrg)

zmin = [0.01]
zmax = [0.65]
normwave = [5300., 5400.]

;; ratio = normflux[3020AA:3100AA]/normflux
;; ratio = [1./0.546970, 1., 0.685592, 0.560925*0.685592]

minmaxwave = jhusdss_sdsswave_minmax(boss=boss)
loglam = jhusdss_get_loglam(minwave=minmaxwave[0], maxwave=minmaxwave[1], nwave=nwave)
wave = 10.^loglam

;outstr = {wave:wave, flux:fltarr(nlrg,nwave), ivar:fltarr(nlrg,nwave), $
;          nmf_continuum:fltarr(nlrg,nwave), med_continuum:fltarr(nlrg,nwave), $
;          residual:fltarr(nlrg,nwave), norm_ratio: fltarr(nlrg)}

flux_outstr = {wave:wave, flux:fltarr(nlrg, nwave), ivar:fltarr(nlrg, nwave), $
               norm_ratio:fltarr(nlrg)}
continuum_outstr = {wave:wave, nmf_continuum:fltarr(nlrg,nwave), med_continuum:fltarr(nlrg,nwave)}
residual_outstr = {wave:wave, residual:fltarr(nlrg,nwave)}

for i=0L, 0L do begin
    iz = where(lrg.z ge zmin[i] and lrg.z lt zmax[i], nz)
;   iz = iz[3001:3080]
;   nz = n_elements(iz)
    print, nz

    jhusdss_load_interp_spec, lrg[iz], loglam=loglam, $
            zuse=fltarr(nz), boss=boss, allflux=allflux, allivar=allivar

    iwavemin = value_locate(wave, normwave[0]*(1.+lrg[iz].z))
    iwavemax = value_locate(wave, normwave[1]*(1.+lrg[iz].z))
    normflux = fltarr(nz)
    for j=0L, nz-1L do normflux[j] = median(allflux[j, iwavemin[j]:iwavemax[j]])
    
    flux_outstr.flux[iz,*] = allflux;/normflux
    flux_outstr.ivar[iz,*] = allivar;*normflux^2/ratio[i]^2
    flux_outstr.norm_ratio[iz] = normflux;*ratio[i]

    for ispec=0L, nz-1L do begin
        counter, ispec+1, nz
        ;; load
        tmpspec = jhusdss_lrg_decompose_loadspec(lrg[iz[ispec]].plate, lrg[iz[ispec]].fiber, $
                          lrg[iz[ispec]].mjd, lrgver, boss=boss, error=error)
        if (error) then continue
        tmpwave = tmpspec.wave*(1.+tmpspec.z)
        tmpflux = tmpspec.flux
        tmpivar = tmpspec.ivar
        tmpivar_residual = tmpspec.ivar*tmpspec.nmf_continuum^2*tmpspec.med_continuum^2
        tmp_nmf_continuum = tmpspec.nmf_continuum
        tmp_med_continuum = tmpspec.med_continuum
        tmp_residual= tmpspec.residual
        ;; interpolate
        tmpmask = (tmpivar le 0.)
        curr_loglam = alog10(tmpwave)

        ;; get the normalized flux
        combine1fiber, curr_loglam, tmpflux, tmpivar, newloglam=loglam, $
                   newflux=finterp, newivar=iinterp, maxiter=0, $
                   finalmask=tmpmask, andmask=maskinterp
        tmpratio = median(flux_outstr.flux[iz[ispec],*]/finterp)

        combine1fiber, curr_loglam, tmp_nmf_continuum, tmpivar, newloglam=loglam, $
                   newflux=finterp, newivar=iinterp, maxiter=0, $
                   finalmask=tmpmask, andmask=maskinterp
        continuum_outstr.nmf_continuum[iz[ispec],*] = finterp*tmpratio

        combine1fiber, curr_loglam, tmp_med_continuum, tmpivar_residual, newloglam=loglam, $
                   newflux=finterp, newivar=iinterp, maxiter=0, $
                   finalmask=tmpmask, andmask=maskinterp
        continuum_outstr.med_continuum[iz[ispec],*] = finterp

        combine1fiber, curr_loglam, tmp_residual, tmpivar_residual, newloglam=loglam, $
                   newflux=finterp, newivar=iinterp, maxiter=0, $
                   finalmask=tmpmask, andmask=maskinterp
        residual_outstr.residual[iz[ispec],*] = finterp
    endfor
endfor

;mwrfits, outstr, outfile, /create
mwrfits, flux_outstr, flux_outfile, /create
mwrfits, continuum_outstr, continuum_outfile, /create
mwrfits, residual_outstr, residual_outfile, /create

end
