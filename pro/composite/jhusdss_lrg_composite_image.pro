pro jhusdss_lrg_composite_image, lrgver, boss=boss, overwrite=overwrite, nsam=nsam

if (n_elements(lrgver) eq 0) then message, 'lrgver required'
if (n_elements(nsam) eq 0) then nsam=50L

lrgpath=jhusdss_get_path(/fitlrg)+'/'+string(lrgver, format='(I3.3)')+'/'

fluxfile = lrgpath+'/'+jhusdss_lrg_composite_image_filename(lrgver, boss=boss)
resifile = lrgpath+'/Residual_'+jhusdss_lrg_composite_image_filename(lrgver, boss=boss)

if (file_test(fluxfile) and ~keyword_set(overwrite)) then begin
   splog, 'File exists. Use Overwrite to overwrite!'
   return
endif else begin
   splog, 'Will write into these files:'
   splog, fluxfile
   splog, resifile
endelse

lrg = jhusdss_lrg_readin(boss=boss)
allspec = jhusdss_read_alllrgspec(lrgver, boss=boss)
allspec_residual = jhusdss_read_alllrgspec(lrgver, /residual, boss=boss)

itrim = jhusdss_lrg_trim(lrg)

;where(lrg.sn_median gt 1. $
;             and lrg.z ge 0.15 $
;             and lrg.z lt 0.55)

isort = bsort(lrg[itrim].z)
;iuse = where(stats[itrim[isort]].isitdecomposed, nuse)
;index = itrim[isort[iuse]]
index = itrim[isort]
nuse = n_elements(index)

zmin = min(lrg[index].z)
zmax = max(lrg[index].z)
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
;for i=4000L, nbin-1L do begin
    counter, i+1, nbin
    ibegin = i*nsam
    iend = (i+1L)*nsam-1L
    if (i eq nbin-1L) then iend = nuse-1
    tmpflux = allspec.flux[index[ibegin:iend], *]
    tmpresi = allspec_residual.residual[index[ibegin:iend], *]
    tmpivar = allspec.ivar[index[ibegin:iend], *]

    for j=ibegin, iend do begin
        tmpflux[j-ibegin,*] = allspec.flux[index[j], *]/allspec.norm_ratio[index[j]]
        tmpresi[j-ibegin,*] = allspec_residual.residual[index[j], *]
        tmpivar[j-ibegin,*] = allspec.ivar[index[j], *]*allspec.norm_ratio[index[j]]^2
    endfor

    jhusdss_composite_engine, tmpflux, tmpivar, fmean=fluxmean, fmedian=fluxmedian, fsigma=fluxsigma
    fluximage.z[i] = median(lrg[index[ibegin:iend]].z)
    fluximage.image_median[*,i] = fluxmedian
    fluximage.image_mean[*,i] = fluxmean
    fluximage.image_sigma[*,i] = fluxsigma
    fluximage.index[0:iend-ibegin,i] = index[ibegin:iend]

;    zmedian = median(lrg[index[ibegin:iend]].z)
;   i_lya_wave = where(allspec.wave le (1216.-15.)*(1.+zmedian), n_lya_wave)
;   if (n_lya_wave gt 21)  then begin
;   for j=ibegin, iend do begin
;       tmpsort = bsort(tmpresi[j-ibegin, i_lya_wave])
;       tmpnorm = median(tmpresi[j-ibegin, i_lya_wave[tmpsort[n_lya_wave-3.*n_lya_wave/4.:n_lya_wave-1]]])
;       tmpresi[j-ibegin, i_lya_wave] = tmpresi[j-ibegin, i_lya_wave]/tmpnorm
;   endfor
;   endif

    jhusdss_composite_engine, tmpresi, tmpivar, fmean=resimean, fmedian=resimedian, fsigma=resisigma
    resiimage.z[i] = median(lrg[index[ibegin:iend]].z)
    resiimage.image_median[*,i] = resimedian
    resiimage.image_mean[*,i] = resimean
    resiimage.image_sigma[*,i] = resisigma
    resiimage.index[0:iend-ibegin,i] = index[ibegin:iend]
;   outmean[*,i] = fluxmean
;   djs_plot, allspec.wave, resimedian, xra=[3800., 9200.], xst=1, yra=[0.5, 3.], yst=1
;   stop
endfor

mwrfits, fluximage, fluxfile, /create
mwrfits, resiimage, resifile, /create
;jhusdss_lrg_composite_engine, inloglam, influx, inivar, z, outloglam, ivarweight=ivarweight

end
