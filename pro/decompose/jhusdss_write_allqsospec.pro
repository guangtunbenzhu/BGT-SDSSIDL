;; Obsolete
;+
; load and interpolate *ALL* quasar spectra
; All in observer Frame
;; see jhusdss_decompose_allinone.pro
;-
pro jhusdss_write_allqsospec, nmfver, boss=boss, overwrite=overwrite

if (n_elements(nmfver) eq 0) then message, 'nmfver required'

;; qsopath
qsopath=jhusdss_get_path(/nmfqso)+'/'+string(nmfver, format='(I3.3)')+'/'
outfile = qsopath+'/'+jhusdss_allqsospec_filename(nmfver, boss=boss)

if (file_test(outfile) and ~keyword_set(overwrite)) then begin
   splog, 'File exists. Use Overwrite to overwrite!'
   return
endif else begin
   splog, 'Will write into this file:'
   splog, outfile
endelse

qso = jhusdss_qso_readin(boss=boss)
nqso = n_elements(qso)

zmin = [0.0, 0.4, 1.8, 2.8]
zmax = [0.4, 1.8, 2.8, 4.8]
;; ratio = normflux[3020AA:3100AA]/normflux
ratio = [1./0.546970, 1., 0.685592, 0.560925*0.685592]

minmaxwave = jhusdss_sdsswave_minmax(boss=boss)
loglam = jhusdss_get_loglam(minwave=minmaxwave[0], maxwave=minmaxwave[1], nwave=nwave)
wave = 10.^loglam

outstr = {wave:wave, flux:fltarr(nqso,nwave), ivar:fltarr(nqso,nwave), $
          nmf_continuum:fltarr(nqso,nwave), med_continuum:fltarr(nqso,nwave), $
          residual:fltarr(nqso,nwave), norm_ratio: fltarr(nqso)}

for i=0L, 3L do begin
    iz = where(qso.z ge zmin[i] and qso.z lt zmax[i], nz)
;   iz = iz[3001:3080]
;   nz = n_elements(iz)
    print, nz

    normwave = jhusdss_normwave_minmax(option=i+1)

    jhusdss_load_interp_spec, qso[iz], loglam=loglam, $
            zuse=fltarr(nz), boss=boss, allflux=allflux, allivar=allivar

    iwavemin = value_locate(wave, normwave[0]*(1.+qso[iz].z))
    iwavemax = value_locate(wave, normwave[1]*(1.+qso[iz].z))
    normflux = fltarr(nz)
    for j=0L, nz-1L do normflux[j] = median(allflux[j, iwavemin[j]:iwavemax[j]])
    
    outstr.flux[iz,*] = allflux;/normflux
    outstr.ivar[iz,*] = allivar;*normflux^2/ratio[i]^2
    outstr.norm_ratio[iz] = normflux*ratio[i]

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
        tmpivar_2 = fltarr(n_elements(curr_loglam))+1.
        combine1fiber, curr_loglam, tmpflux, tmpivar, newloglam=loglam, $
                   newflux=finterp, newivar=iinterp, maxiter=0, $
                   finalmask=tmpmask, andmask=maskinterp
        tmpratio = median(outstr.flux[iz[ispec],*]/finterp)

        combine1fiber, curr_loglam, tmp_nmf_continuum, tmpivar, newloglam=loglam, $
                   newflux=finterp, newivar=iinterp, maxiter=0, $
                   finalmask=tmpmask, andmask=maskinterp
        outstr.nmf_continuum[iz[ispec],*] = finterp*tmpratio

        combine1fiber, curr_loglam, tmp_med_continuum, tmpivar_2, newloglam=loglam, $
                   newflux=finterp, newivar=iinterp, maxiter=0, $
                   finalmask=tmpmask, andmask=maskinterp
        outstr.med_continuum[iz[ispec],*] = finterp

        combine1fiber, curr_loglam, tmp_residual, tmpivar_2, newloglam=loglam, $
                   newflux=finterp, newivar=iinterp, maxiter=0, $
                   finalmask=tmpmask, andmask=maskinterp
        outstr.residual[iz[ispec],*] = finterp
    endfor
endfor

mwrfits, outstr, outfile, /create

end
