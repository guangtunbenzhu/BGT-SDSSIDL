pro jhusdss_civ_absorber_composite_image, nmfver, boss=boss, overwrite=overwrite, nsam=nsam

if (n_elements(nmfver) eq 0) then message, 'nmfver required'
if (n_elements(nsam) eq 0) then nsam=25L

;; output
if (keyword_set(boss)) then begin
   path=jhusdss_get_path(/nmfqso)+'/'+string(nmfver, format='(I3.3)')+'/Absorbers_BOSS'
endif else begin
   path=jhusdss_get_path(/nmfqso)+'/'+string(nmfver, format='(I3.3)')+'/Absorbers'
endelse

resifile = path+'/CIV_Cooksey_'+jhusdss_absorber_composite_image_filename(nmfver, boss=boss)
if (file_test(resifile) and ~keyword_set(overwrite)) then begin
   splog, 'File exists. Use Overwrite to overwrite!'
   return
endif else begin
   splog, 'Will write into this file:'
   splog, resifile
endelse

qso = jhusdss_qso_readin(boss=boss)
stats = jhusdss_qsostats_readin(nmfver, boss=boss)
;absorbers = jhusdss_absorber_readin(nmfver, boss=boss)
absorbers = mrdfits('~/SDATA/SDSS/CIV/sdss_civ_cookseyetal13_withqsoindex_withzabs.fits', 1)
allspec = jhusdss_read_allqsospec(nmfver, boss=boss)

index = bsort(absorbers.zabs)
nuse = n_elements(index)

zmin = min(absorbers[index].zabs)
zmax = max(absorbers[index].zabs)
nbin = nuse/nsam+1L

nwave = n_elements(allspec.wave)
resiimage = {z:fltarr(nbin), image_mean:fltarr(nwave, nbin), image_median:fltarr(nwave, nbin)}

for i=0L, nbin-1L do begin
    counter, i+1, nbin
    ibegin = i*nsam
    iend = (i+1L)*nsam-1L
    if (i eq nbin-1L) then iend = nuse-1

    tmpzabs = absorbers[index[ibegin:iend]].zabs
    tmpqsoindex = absorbers[index[ibegin:iend]].index_qso
    tmpflux = allspec.flux[tmpqsoindex,*]
    tmpnmfcontinuum = allspec.nmf_continuum[tmpqsoindex,*]
    tmpresi = allspec.residual[tmpqsoindex,*]
    tmpivar = allspec.ivar[tmpqsoindex,*]

    ;; recalculate med_medium
    ;; jhusdss_redecompose_spec.pro
    for j=0L,iend-ibegin do begin
        thiszabs = tmpzabs[j]
        thisflux = reform(tmpflux[j,*])
        thisnmfcontinuum  = reform(tmpnmfcontinuum[j,*])
        thisivar = reform(tmpivar[j, *])
        thisspec = {z:0., wave:allspec.wave, ivar:thisivar, flux:thisflux, nmf_continuum:thisnmfcontinuum}
        tmpresi[j,*] = jhusdss_redecompose_spec(thisspec, thiszabs, new_med_continuum=new_med_continuum)
;       tmpresi[j,*] = thisflux/thisnmfcontinuum/new_med_continuum
    endfor

    jhusdss_composite_engine, tmpresi, tmpivar, fmean=fluxmean, fmedian=fluxmedian
    resiimage.z[i] = median(tmpzabs)
    resiimage.image_median[*,i] = fluxmedian
    resiimage.image_mean[*,i] = fluxmean

endfor

mwrfits, resiimage, resifile, /create
;jhusdss_qso_composite_engine, inloglam, influx, inivar, z, outloglam, ivarweight=ivarweight

end
