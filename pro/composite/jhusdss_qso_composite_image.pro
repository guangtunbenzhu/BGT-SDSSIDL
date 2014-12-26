pro jhusdss_qso_composite_image, nmfver, boss=boss, overwrite=overwrite, nsam=nsam

if (n_elements(nmfver) eq 0) then message, 'nmfver required'
if (n_elements(nsam) eq 0) then nsam=25L

qsopath=jhusdss_get_path(/nmfqso)+'/'+string(nmfver, format='(I3.3)')+'/'
fluxfile = qsopath+'/'+jhusdss_qso_composite_image_filename(nmfver, boss=boss)
resifile = qsopath+'/Residual_'+jhusdss_qso_composite_image_filename(nmfver, boss=boss)
if (file_test(fluxfile) and ~keyword_set(overwrite)) then begin
   splog, 'File exists. Use Overwrite to overwrite!'
   return
endif else begin
   splog, 'Will write into these files:'
   splog, fluxfile
   splog, resifile
endelse

qso = jhusdss_qso_readin(boss=boss)
stats = jhusdss_qsostats_readin(nmfver, boss=boss)
allspec = jhusdss_read_allqsospec(nmfver, boss=boss)

; itrim = jhusdss_absorber_trim(stats, /qso)
  itrim = where(stats.spec_snr_median gt 1. $
             and stats.med_sdeviation_red gt 0.00 $
             and stats.med_sdeviation_red le 0.07 $
             and stats.zqso ge 0.1 $
             and stats.zqso lt 4.8)

isort = bsort(qso[itrim].z)
;iuse = where(stats[itrim[isort]].isitdecomposed, nuse)
;index = itrim[isort[iuse]]
index = itrim[isort]
nuse = n_elements(index)

zmin = min(qso[index].z)
zmax = max(qso[index].z)
nbin = nuse/nsam+1L

nwave = n_elements(allspec.wave)
fluximage = {z:fltarr(nbin), wave:allspec.wave, $
             image_mean:fltarr(nwave,nbin), image_median:fltarr(nwave,nbin), $
             image_sigma:fltarr(nwave,nbin), index:lonarr(nsam,nbin)}
resiimage = {z:fltarr(nbin), wave:allspec.wave, $
             image_mean:fltarr(nwave, nbin), image_median:fltarr(nwave, nbin), $
             image_sigma:fltarr(nwave, nbin), index:lonarr(nsam,nbin)}
;outmean= fltarr(nwave, nbin)

for i=0L, nbin-1L do begin
    counter, i+1, nbin
    ibegin = i*nsam
    iend = (i+1L)*nsam-1L
    if (i eq nbin-1L) then iend = nuse-1
    tmpflux = allspec.flux[index[ibegin:iend], *]
    tmpresi = allspec.residual[index[ibegin:iend], *]
    tmpivar = allspec.ivar[index[ibegin:iend], *]

    for j=ibegin, iend do begin
        tmpflux[j-ibegin,*] = allspec.flux[index[j], *]/allspec.norm_ratio[index[j]]
        tmpresi[j-ibegin,*] = allspec.residual[index[j], *]
        tmpivar[j-ibegin,*] = allspec.ivar[index[j], *]*allspec.norm_ratio[index[j]]^2
    endfor

    jhusdss_composite_engine, tmpflux, tmpivar, fmean=fluxmean, fmedian=fluxmedian, fsigma=fluxsigma
    fluximage.z[i] = median(qso[index[ibegin:iend]].z)
    fluximage.image_median[*,i] = fluxmedian
    fluximage.image_mean[*,i] = fluxmean
    fluximage.image_sigma[*,i] = fluxsigma
    fluximage.index[0:iend-ibegin,i] = index[ibegin:iend]

;    zmedian = median(qso[index[ibegin:iend]].z)
;   i_lya_wave = where(allspec.wave le (1216.-15.)*(1.+zmedian), n_lya_wave)
;   if (n_lya_wave gt 21)  then begin
;   for j=ibegin, iend do begin
;       tmpsort = bsort(tmpresi[j-ibegin, i_lya_wave])
;       tmpnorm = median(tmpresi[j-ibegin, i_lya_wave[tmpsort[n_lya_wave-3.*n_lya_wave/4.:n_lya_wave-1]]])
;       tmpresi[j-ibegin, i_lya_wave] = tmpresi[j-ibegin, i_lya_wave]/tmpnorm
;   endfor
;   endif

    jhusdss_composite_engine, tmpresi, tmpivar, fmean=resimean, fmedian=resimedian, fsigma=resisigma
    resiimage.z[i] = median(qso[index[ibegin:iend]].z)
    resiimage.image_median[*,i] = resimedian
    resiimage.image_mean[*,i] = resimean
    resiimage.image_sigma[*,i] = resisigma
    resiimage.index[0:iend-ibegin,i] = index[ibegin:iend]
;   outmean[*,i] = fluxmean
endfor

mwrfits, fluximage, fluxfile, /create
mwrfits, resiimage, resifile, /create
;jhusdss_qso_composite_engine, inloglam, influx, inivar, z, outloglam, ivarweight=ivarweight

end
