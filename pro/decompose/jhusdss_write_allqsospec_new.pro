;+
; load and interpolate *ALL* quasar spectra
; All in observer Frame
;; see jhusdss_decompose_allinone.pro
;-
pro jhusdss_write_allqsospec_new, nmfver, boss=boss, overwrite=overwrite

if (n_elements(nmfver) eq 0) then message, 'nmfver required'

;; qsopath
parent_path=jhusdss_get_parent_path()

qsopath = parent_path+'/SDSS/QSO/NMF/'+string(nmfver, format='(I3.3)')+'/'
outfile = qsopath+'/AllInOne/'+jhusdss_allqsospec_filename(nmfver, boss=boss)

flux_outfile = repstr(outfile, '.fits', '_flux.fits')
cont_outfile = repstr(outfile, '.fits', '_continuum.fits')
normresi_outfile = repstr(outfile, '.fits', '_normalized_residual.fits')
subtresi_outfile = repstr(outfile, '.fits', '_subtracted_residual.fits')

if (file_test(flux_outfile) and ~keyword_set(overwrite)) then begin
   splog, 'File exists. Use Overwrite to overwrite!'
   return
endif else begin
   splog, 'Will write into these files:'
   print, flux_outfile 
   print, cont_outfile
   print, normresi_outfile
   print, subtresi_outfile
endelse

qso = jhusdss_qso_readin(boss=boss)
stats = jhusdss_qsostats_readin(nmfver, boss=boss)
nqso = n_elements(qso)

zmin = [0.0, 0.4, 1.8, 2.8]
zmax = [0.4, 1.8, 2.8, 4.8]
;; ratio = normflux[3020AA:3100AA]/normflux
ratio = [1./0.546970, 1., 0.685592, 0.560925*0.685592]

minmaxwave = jhusdss_sdsswave_minmax(boss=boss)
loglam = jhusdss_get_loglam(minwave=minmaxwave[0], maxwave=minmaxwave[1], nwave=nwave)
wave = 10.^loglam

flux_outstr = {wave:wave, flux:fltarr(nqso,nwave), ivar:fltarr(nqso,nwave), $
               norm_ratio:fltarr(nqso), $
               ra:dblarr(nqso), dec:dblarr(nqso), $
               plate:lonarr(nqso), fiber:lonarr(nqso), mjd:lonarr(nqso), $
               zqso:fltarr(nqso), err_zqso:fltarr(nqso), $
               spec_snr_median:fltarr(nqso), isitdecomposed:bytarr(nqso), $
               med_sdeviation_red:fltarr(nqso), med_sdeviation_blue:fltarr(nqso)} 

cont_outstr = {wave:wave, continuum:fltarr(nqso,nwave), $
               nmf_continuum:fltarr(nqso,nwave), med_continuum:fltarr(nqso, nwave), $
               ra:dblarr(nqso), dec:dblarr(nqso), $
               plate:lonarr(nqso), fiber:lonarr(nqso), mjd:lonarr(nqso), $
               zqso:fltarr(nqso), err_zqso:fltarr(nqso), $
               spec_snr_median:fltarr(nqso), isitdecomposed:bytarr(nqso), $
               med_sdeviation_red:fltarr(nqso), med_sdeviation_blue:fltarr(nqso)} 

normresi_outstr = {wave:wave, residual:fltarr(nqso,nwave), $
               residual_ivar:fltarr(nqso,nwave), $
               ra:dblarr(nqso), dec:dblarr(nqso), $
               plate:lonarr(nqso), fiber:lonarr(nqso), mjd:lonarr(nqso), $
               zqso:fltarr(nqso), err_zqso:fltarr(nqso), $
               spec_snr_median:fltarr(nqso), isitdecomposed:bytarr(nqso), $
               med_sdeviation_red:fltarr(nqso), med_sdeviation_blue:fltarr(nqso)} 

subtresi_outstr = {wave:wave, subtracted_residual:fltarr(nqso,nwave), $
               subtracted_ivar:fltarr(nqso,nwave), $
               ra:dblarr(nqso), dec:dblarr(nqso), $
               plate:lonarr(nqso), fiber:lonarr(nqso), mjd:lonarr(nqso), $
               zqso:fltarr(nqso), err_zqso:fltarr(nqso), $
               spec_snr_median:fltarr(nqso), isitdecomposed:bytarr(nqso), $
               med_sdeviation_red:fltarr(nqso), med_sdeviation_blue:fltarr(nqso)} 

flux_outstr.ra = stats.ra & flux_outstr.dec = stats.dec & flux_outstr.plate = stats.plate & flux_outstr.fiber = stats.fiber & flux_outstr.mjd = stats.mjd & flux_outstr.zqso = stats.zqso & flux_outstr.err_zqso = stats.err_zqso & flux_outstr.spec_snr_median = stats.spec_snr_median & flux_outstr.isitdecomposed = stats.isitdecomposed & flux_outstr.med_sdeviation_red = stats.med_sdeviation_red & flux_outstr.med_sdeviation_blue = stats.med_sdeviation_blue

cont_outstr.ra = stats.ra & cont_outstr.dec = stats.dec & cont_outstr.plate = stats.plate & cont_outstr.fiber = stats.fiber & cont_outstr.mjd = stats.mjd & cont_outstr.zqso = stats.zqso & cont_outstr.err_zqso = stats.err_zqso & cont_outstr.spec_snr_median = stats.spec_snr_median & cont_outstr.isitdecomposed = stats.isitdecomposed & cont_outstr.med_sdeviation_red = stats.med_sdeviation_red & cont_outstr.med_sdeviation_blue = stats.med_sdeviation_blue

subtresi_outstr.ra = stats.ra & subtresi_outstr.dec = stats.dec & subtresi_outstr.plate = stats.plate & subtresi_outstr.fiber = stats.fiber & subtresi_outstr.mjd = stats.mjd & subtresi_outstr.zqso = stats.zqso & subtresi_outstr.err_zqso = stats.err_zqso & subtresi_outstr.spec_snr_median = stats.spec_snr_median & subtresi_outstr.isitdecomposed = stats.isitdecomposed & subtresi_outstr.med_sdeviation_red = stats.med_sdeviation_red & subtresi_outstr.med_sdeviation_blue = stats.med_sdeviation_blue

normresi_outstr.ra = stats.ra & normresi_outstr.dec = stats.dec & normresi_outstr.plate = stats.plate & normresi_outstr.fiber = stats.fiber & normresi_outstr.mjd = stats.mjd & normresi_outstr.zqso = stats.zqso & normresi_outstr.err_zqso = stats.err_zqso & normresi_outstr.spec_snr_median = stats.spec_snr_median & normresi_outstr.isitdecomposed = stats.isitdecomposed & normresi_outstr.med_sdeviation_red = stats.med_sdeviation_red & normresi_outstr.med_sdeviation_blue = stats.med_sdeviation_blue

for i=0L, 3L do begin
    iz = where(qso.z ge zmin[i] and qso.z lt zmax[i], nz)
    print, nz

    normwave = jhusdss_normwave_minmax(option=i+1)

    jhusdss_load_interp_spec, qso[iz], loglam=loglam, $
            zuse=fltarr(nz), boss=boss, allflux=allflux, allivar=allivar

    iwavemin = value_locate(wave, normwave[0]*(1.+qso[iz].z))
    iwavemax = value_locate(wave, normwave[1]*(1.+qso[iz].z))
    normflux = fltarr(nz)
    for j=0L, nz-1L do normflux[j] = median(allflux[j, iwavemin[j]:iwavemax[j]])
    
    flux_outstr.flux[iz,*] = allflux;/normflux
    flux_outstr.ivar[iz,*] = allivar;*normflux^2/ratio[i]^2
    flux_outstr.norm_ratio[iz] = normflux*ratio[i]
    subtresi_outstr.subtracted_ivar[iz,*] = allivar

    for ispec=0L, nz-1L do begin
        counter, ispec+1, nz
        ;; load
        tmpspec = jhusdss_decompose_loadspec(qso[iz[ispec]].plate, qso[iz[ispec]].fiber, nmfver, boss=boss, error=error)
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
        cont_outstr.nmf_continuum[iz[ispec],*] = finterp*tmpratio

        combine1fiber, curr_loglam, tmp_med_continuum, tmpivar_residual, newloglam=loglam, $
                   newflux=finterp, newivar=iinterp, maxiter=0, $
                   finalmask=tmpmask, andmask=maskinterp
        cont_outstr.med_continuum[iz[ispec],*] = finterp

        combine1fiber, curr_loglam, tmp_residual, tmpivar_residual, newloglam=loglam, $
                   newflux=finterp, newivar=iinterp, maxiter=0, $
                   finalmask=tmpmask, andmask=maskinterp
        normresi_outstr.residual[iz[ispec],*] = finterp
        normresi_outstr.residual_ivar[iz[ispec],*] = iinterp
    endfor
endfor

cont_outstr.continuum = cont_outstr.nmf_continuum*cont_outstr.med_continuum
subtresi_outstr.subtracted_residual = flux_outstr.flux-cont_outstr.continuum

mwrfits, flux_outstr, flux_outfile, /create
mwrfits, cont_outstr, cont_outfile, /create
mwrfits, normresi_outstr, normresi_outfile, /create
mwrfits, subtresi_outstr, subtresi_outfile, /create

end
