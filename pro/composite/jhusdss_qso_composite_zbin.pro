;; pro jhusdss_qso_composite_zbin
;; to-do: de-redden spectra?


nmfver=106
ivarweigth = 1b
overwrite = 1b
snrcut = 2.
sdevcut = 0.10
sigma_cut = 10.

lambda=[3551., 4686., 6165., 7481., 8931.]

qsopath=jhusdss_get_path(/nmfqso)+'/'+string(nmfver, format='(I3.3)')+'/'
fluxfile = qsopath+'/'+jhusdss_qso_composite_image_filename(nmfver, boss=boss)
outfile = repstr(fluxfile, '.fits', '_zbin.fits')
outphotofile = repstr(fluxfile, '.fits', '_photo_zbin.fits')
if (file_test(outfile) and ~keyword_set(overwrite)) then begin
   splog, 'File exists. Use Overwrite to overwrite!'
   return
endif else begin
   splog, 'Will write into these files:'
   splog, outfile
endelse

choice_load_data = 0
read,'load data? [1=yes, 0=no]: ',choice_load_data

if choice_load_data eq 1 then begin
   qso = jhusdss_qso_readin(boss=boss)
   stats = jhusdss_qsostats_readin(nmfver, boss=boss)
   allspec = jhusdss_read_allqsospec(nmfver, /flux, boss=boss)
   qsophoto_infile = getenv('RAW_DATA')+'/SDSS/QSO/HW_dr7qso_newz_photo_more.fits'
   qso_photo = mrdfits(qsophoto_infile, 1)
endif

;; zbins
z_min = 0.10
z_max = 4.8
dz=0.005
delta_z = 0.0005
z_bin = jhusdss_make_bins(z_min, z_max, delta_z, nbin=nbin)
z_bin[nbin-1].max = z_max

;; structure
nwave = n_elements(allspec.wave)
outstr_i_191 = {imag_max:19.1, imag_min:15.0, z:z_bin.mean, wave:allspec.wave, $
             fluxmean:fltarr(nbin, nwave), fluxmedian:fltarr(nbin, nwave), $
             npairs:lonarr(nbin)}
outstr_i_203 = {imag_max:20.3, imag_min:15.0, z:z_bin.mean, wave:allspec.wave, $
             fluxmean:fltarr(nbin, nwave), fluxmedian:fltarr(nbin, nwave), $
             npairs:lonarr(nbin)}
outstr_i_191_203 = {imag_max:20.3, image_min:19.1, z:z_bin.mean, wave:allspec.wave, $
             fluxmean:fltarr(nbin, nwave), fluxmedian:fltarr(nbin, nwave), $
             npairs:lonarr(nbin)}

phoout_i_191 = {imag_max:19.1, imag_min:15.0, z:z_bin.mean, $
             ugriz_mean:fltarr(nbin,5), ugriz_median:fltarr(nbin,5), $
             ugriz_dered_mean:fltarr(nbin,5), ugriz_dered_median:fltarr(nbin,5), $
             npairs:lonarr(nbin)}
phoout_i_203 = {imag_max:20.3, imag_min:15.0, z:z_bin.mean, $
             ugriz_mean:fltarr(nbin,5), ugriz_median:fltarr(nbin,5), $
             ugriz_dered_mean:fltarr(nbin,5), ugriz_dered_median:fltarr(nbin,5), $
             npairs:lonarr(nbin)}
phoout_i_191_203 = {imag_max:20.3, imag_min:19.1, z:z_bin.mean, $
             ugriz_mean:fltarr(nbin,5), ugriz_median:fltarr(nbin,5), $
             ugriz_dered_mean:fltarr(nbin,5), ugriz_dered_median:fltarr(nbin,5), $
             npairs:lonarr(nbin)}

for iz=0L, nbin-1L do begin

   ; QUALITY CONTROL - make sure it is consistent through all studies, shouldn't matter if one uses median
   ithis= where(qso.z gt z_bin[iz].mean-dz and qso.z lt z_bin[iz].mean+dz $
             and qso_photo.bal_flag eq 0L $
             and stats.spec_snr_median gt snrcut $
             and stats.med_sdeviation_red gt 0. $
             and stats.med_sdeviation_red le sdevcut, nthis)

   ithis_i_191 = where(qso_photo[ithis].ugriz_dered[3] lt 19.1, nthis_i_191)
   ithis_i_203 = where(qso_photo[ithis].ugriz_dered[3] lt 20.3, nthis_i_203)
   ithis_i_191_203 = where(qso_photo[ithis].ugriz_dered[3] lt 20.3 and qso_photo[ithis].ugriz_dered[3] ge 19.1, nthis_i_191_203)

   i_191 = ithis[ithis_i_191]
   i_203 = ithis[ithis_i_203]
   i_191_203 = ithis[ithis_i_191_203]

   print, "z=", z_bin[iz].mean, nthis_i_191, nthis_i_203, nthis_i_191_203
   ;; PHOTOMETRY

if (nthis_i_191 gt 1) then begin
   in_photo_191 = transpose(qso_photo[i_191].ugriz)
   in_photo_dered_191 = transpose(qso_photo[i_191].ugriz_dered)
   in_photoivar_191 = transpose(1./qso_photo[i_191].ugriz_err^2)
   weight = fltarr(nthis_i_191)+1.
   jhusdss_composite_engine, in_photo_191, in_photoivar_191, fmean=fmean, fmedian=fmedian, $
       nobjuse=nobjuse, ivarweight=ivarweight, weight=weight, sigma_cut=sigma_cut
   phoout_i_191.ugriz_mean[iz, *] = fmean
   phoout_i_191.ugriz_median[iz, *] = fmedian
   jhusdss_composite_engine, in_photo_dered_191, in_photoivar_191, fmean=fmean, fmedian=fmedian, $
       nobjuse=nobjuse, ivarweight=ivarweight, weight=weight, sigma_cut=sigma_cut
   phoout_i_191.ugriz_dered_mean[iz, *] = fmean
   phoout_i_191.ugriz_dered_median[iz, *] = fmedian
endif

if (nthis_i_203 gt 1) then begin
   in_photo_203 = transpose(qso_photo[i_203].ugriz)
   in_photo_dered_203 = transpose(qso_photo[i_203].ugriz_dered)
   in_photoivar_203 = transpose(1./qso_photo[i_203].ugriz_err^2)
   weight = fltarr(nthis_i_203)+1.
   jhusdss_composite_engine, in_photo_203, in_photoivar_203, fmean=fmean, fmedian=fmedian, $
       nobjuse=nobjuse, ivarweight=ivarweight, weight=weight, sigma_cut=sigma_cut
   phoout_i_203.ugriz_mean[iz, *] = fmean
   phoout_i_203.ugriz_median[iz, *] = fmedian
   jhusdss_composite_engine, in_photo_dered_203, in_photoivar_203, fmean=fmean, fmedian=fmedian, $
       nobjuse=nobjuse, ivarweight=ivarweight, weight=weight, sigma_cut=sigma_cut
   phoout_i_203.ugriz_dered_mean[iz, *] = fmean
   phoout_i_203.ugriz_dered_median[iz, *] = fmedian
endif

if (nthis_i_191_203 gt 1) then begin
   in_photo_191_203 = transpose(qso_photo[i_191_203].ugriz)
   in_photo_dered_191_203 = transpose(qso_photo[i_191_203].ugriz_dered)
   in_photoivar_191_203 = transpose(1./qso_photo[i_191_203].ugriz_err^2)
   weight = fltarr(nthis_i_191_203)+1.
   jhusdss_composite_engine, in_photo_191_203, in_photoivar_191_203, fmean=fmean, fmedian=fmedian, $
       nobjuse=nobjuse, ivarweight=ivarweight, weight=weight, sigma_cut=sigma_cut
   phoout_i_191_203.ugriz_mean[iz, *] = fmean
   phoout_i_191_203.ugriz_median[iz, *] = fmedian
   jhusdss_composite_engine, in_photo_dered_191_203, in_photoivar_191_203, fmean=fmean, fmedian=fmedian, $
       nobjuse=nobjuse, ivarweight=ivarweight, weight=weight, sigma_cut=sigma_cut
   phoout_i_191_203.ugriz_dered_mean[iz, *] = fmean
   phoout_i_191_203.ugriz_dered_median[iz, *] = fmedian
endif

   ;; SPECTROSCOPY

   ; calculate wavelength shift
   ; dz=0.01 is big, let's calculate wavelength shift so that blurring only happens within 1 pixel, ignoring the first and last 50 pixels (losing 100 pixels)
   i6000 = value_locate(allspec.wave, 6000.)

if (nthis_i_191 gt 1) then begin
   wave_shift_191_init = ((value_locate(allspec.wave, 6000.*(1.+qso[i_191].z)/(1.+z_bin[iz].mean)) - i6000) < 50L)
   wave_shift_191 = (wave_shift_191_init > (-50L))
   in_flux_191 = fltarr(nthis_i_191, nwave)
   in_ivar_191 = fltarr(nthis_i_191, nwave)+0.
   print, 'i<19.1 Making composite ... '
   for i=0L+50L, nwave-1L-50L do begin
       in_flux_191[*,i] = allspec.flux[i_191, wave_shift_191+i]
       in_ivar_191[*,i] = allspec.ivar[i_191, wave_shift_191+i]
   endfor
   weight = fltarr(nthis_i_191)+1.
   jhusdss_composite_engine, in_flux_191, in_ivar_191, fmean=fmean, fmedian=fmedian, $
      nobjuse=nobjuse, ivarweight=ivarweight, weight=weight, sigma_cut=sigma_cut
   outstr_i_191.fluxmean[iz, *] = fmean
   outstr_i_191.fluxmedian[iz, *] = fmedian
endif

if (nthis_i_203 gt 1) then begin
   wave_shift_203_init = ((value_locate(allspec.wave, 6000.*(1.+qso[i_203].z)/(1.+z_bin[iz].mean)) - i6000) < 50L)
   wave_shift_203 = (wave_shift_203_init > (-50L))
   in_flux_203 = fltarr(nthis_i_203, nwave)
   in_ivar_203 = fltarr(nthis_i_203, nwave)+0.
   print, 'i<20.3 Making composite ... '
   for i=0L+50L, nwave-1L-50L do begin
       in_flux_203[*,i] = allspec.flux[i_203, wave_shift_203+i]
       in_ivar_203[*,i] = allspec.ivar[i_203, wave_shift_203+i]
   endfor
   weight = fltarr(nthis_i_203)+1.
   jhusdss_composite_engine, in_flux_203, in_ivar_203, fmean=fmean, fmedian=fmedian, $
      nobjuse=nobjuse, ivarweight=ivarweight, weight=weight, sigma_cut=sigma_cut
   outstr_i_203.fluxmean[iz, *] = fmean
   outstr_i_203.fluxmedian[iz, *] = fmedian
endif

if (nthis_i_191_203 gt 1) then begin
   wave_shift_191_203_init = ((value_locate(allspec.wave, 6000.*(1.+qso[i_191_203].z)/(1.+z_bin[iz].mean)) - i6000) < 50L)
   wave_shift_191_203 = (wave_shift_191_203_init > (-50L))
   in_flux_191_203 = fltarr(nthis_i_191_203, nwave)
   in_ivar_191_203 = fltarr(nthis_i_191_203, nwave)+0.
   print, '19.1<i<20.3 Making composite ... '
   for i=0L+50L, nwave-1L-50L do begin
       in_flux_191_203[*,i] = allspec.flux[i_191_203, wave_shift_191_203+i]
       in_ivar_191_203[*,i] = allspec.ivar[i_191_203, wave_shift_191_203+i]
   endfor
   weight = fltarr(nthis_i_191_203)+1.
   jhusdss_composite_engine, in_flux_191_203, in_ivar_191_203, fmean=fmean, fmedian=fmedian, $
      nobjuse=nobjuse, ivarweight=ivarweight, weight=weight, sigma_cut=sigma_cut
   outstr_i_191_203.fluxmean[iz, *] = fmean
   outstr_i_191_203.fluxmedian[iz, *] = fmedian
endif

   ; quicklook
   load_dp, /b

   djs_plot, allspec.wave, smooth(outstr_i_191.fluxmean[iz,*], 5), xra=xra, xst=1, yra=yra, yst=1
   djs_oplot, lambda, 10^(-0.4*(phoout_i_191.ugriz_median[1000,*]-8.9))/3.34E4/lambda^2./1.E-17, psym=4, symsize=8, color='red'
   djs_oplot, allspec.wave, smooth(outstr_i_203.fluxmean[iz,*],5)
   djs_oplot, allspec.wave, smooth(outstr_i_191_203.fluxmean[iz,*],5)
   djs_oplot, allspec.wave, smooth(outstr_i_191.fluxmedian[iz,*],5), color='green'
   djs_oplot, allspec.wave, smooth(outstr_i_203.fluxmedian[iz,*],5), color='green'
   djs_oplot, allspec.wave, smooth(outstr_i_191_203.fluxmedian[iz,*],5), color='green'
endfor

stop
mwrfits, outstr_i_191, repstr(outfile, '.fits', '_191.fits'), /create
mwrfits, outstr_i_203, repstr(outfile, '.fits', '_203.fits'), /create
mwrfits, outstr_i_191_203, repstr(outfile, '.fits', '_191_203.fits'), /create
mwrfits, phoout_i_191, repstr(outphotofile, '.fits', '_191.fits'), /create
mwrfits, phoout_i_203, repstr(outphotofile, '.fits', '_203.fits'), /create
mwrfits, phoout_i_191_203, repstr(outphotofile, '.fits', '_191_203.fits'), /create

end
