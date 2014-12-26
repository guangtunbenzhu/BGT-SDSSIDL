;; pro jhusdss_qso_composite_zbin
;; to-do: de-redden spectra?

nmfver=106
ivarweigth = 1b
overwrite = 1b
snrcut = 2.
sdevcut = 0.10
sigma_cut = 2.

lambda=[3551., 4686., 6165., 7481., 8931.]

parent_path=jhusdss_get_parent_path()

qsopath = parent_path+'/SDSS/QSO/NMF/'+string(nmfver, format='(I3.3)')+'/'
outfile = qsopath+'/AllInOne/'+jhusdss_allqsospec_filename(nmfver, boss=boss)
outfile = repstr(outfile, '.fits', '_imag.fits')
outphotofile = repstr(outfile, '.fits', '_photo.fits')

if (file_test(repstr(outphotofile, '.fits', '_203.fits')) and ~keyword_set(overwrite)) then begin
   splog, 'File exists. Use Overwrite to overwrite!'
   return
endif else begin
   splog, 'Will write into these files:'
   splog, outfile
   splog, outphotofile
endelse

choice_load_data = 0
read,'load data? [1=yes, 0=no]: ',choice_load_data

if choice_load_data eq 1 then begin
   qso = jhusdss_qso_readin(boss=boss)
   stats = jhusdss_qsostats_readin(nmfver, boss=boss)
   allspec = jhusdss_read_allqsospec(nmfver, /flux, boss=boss)
   qsophoto_infile = getenv('RAW_DATA')+'/SDSS/QSO/HW_dr7qso_newz_photo_more.fits'
   qso_photo = mrdfits(qsophoto_infile, 1)
   nqso = n_elements(qso)

   ;; read in composite spectra/photometry @ z
   inqsopath=jhusdss_get_path(/nmfqso)+'/'+string(nmfver, format='(I3.3)')+'/'
   fluxfile = inqsopath+'/'+jhusdss_qso_composite_image_filename(nmfver, boss=boss)
   infile = repstr(fluxfile, '.fits', '_zbin.fits')
   inphotofile = repstr(fluxfile, '.fits', '_photo_zbin.fits')

   outstr_i_191 = mrdfits(repstr(infile, '.fits', '_191.fits'), 1)
   outstr_i_203 = mrdfits(repstr(infile, '.fits', '_203.fits'), 1)
   outstr_i_191_203 = mrdfits(repstr(infile, '.fits', '_191_203.fits'), 1)
   phoout_i_191 = mrdfits(repstr(inphotofile, '.fits', '_191.fits'), 1)
   phoout_i_203 = mrdfits(repstr(inphotofile, '.fits', '_203.fits'), 1)
   phoout_i_191_203 = mrdfits(repstr(inphotofile, '.fits', '_191_203.fits'), 1)
endif

;; output
delta_norm_outstr_i_191 = {imag_max:19.1, imag_min:15.0, wave:allspec.wave, delta:fltarr(nqso,nwave), $
               delta_ivar:fltarr(nqso,nwave), $
               ra:qso.ra, dec:qso.dec, $
               plate:qso.plate, fiber:qso.fiber, mjd:qso.mjd, $
               zqso:qso.z, err_zqso:qso.zerr, $
               spec_snr_median:stats.spec_snr_median, isitdecomposed:stats.isitdecomposed, $
               med_sdeviation_red:stats.med_sdeviation_red, med_sdeviation_blue:stats.med_sdeviation_blue}

delta_sub_outstr_i_191 = {imag_max:19.1, imag_min:15.0, wave:allspec.wave, delta_subtracted:fltarr(nqso,nwave), $
               delta_subtracted_ivar:fltarr(nqso,nwave), $
               ra:qso.ra, dec:qso.dec, $
               plate:qso.plate, fiber:qso.fiber, mjd:qso.mjd, $
               zqso:qso.z, err_zqso:qso.zerr, $
               spec_snr_median:stats.spec_snr_median, isitdecomposed:stats.isitdecomposed, $
               med_sdeviation_red:stats.med_sdeviation_red, med_sdeviation_blue:stats.med_sdeviation_blue}

delta_photo_outstr_i_191 = {imag_max:19.1, imag_min:15.0, wave:lambda, delta_ugriz:fltarr(nqso,5), delta_ugriz_dered:fltarr(nqso,5), $
               delta_ugriz_err:fltarr(nqso,5), $
               ra:qso.ra, dec:qso.dec, $
               plate:qso.plate, fiber:qso.fiber, mjd:qso.mjd, $
               zqso:qso.z, err_zqso:qso.zerr, $
               spec_snr_median:stats.spec_snr_median, isitdecomposed:stats.isitdecomposed, $
               med_sdeviation_red:stats.med_sdeviation_red, med_sdeviation_blue:stats.med_sdeviation_blue}

delta_norm_outstr_i_203 = {imag_max:20.3, imag_min:15.0, wave:allspec.wave, delta:fltarr(nqso,nwave), $
               delta_ivar:fltarr(nqso,nwave), $
               ra:qso.ra, dec:qso.dec, $
               plate:qso.plate, fiber:qso.fiber, mjd:qso.mjd, $
               zqso:qso.z, err_zqso:qso.zerr, $
               spec_snr_median:stats.spec_snr_median, isitdecomposed:stats.isitdecomposed, $
               med_sdeviation_red:stats.med_sdeviation_red, med_sdeviation_blue:stats.med_sdeviation_blue}

delta_sub_outstr_i_203 = {imag_max:20.3, imag_min:15.0, wave:allspec.wave, delta_subtracted:fltarr(nqso,nwave), $
               delta_subtracted_ivar:fltarr(nqso,nwave), $
               ra:qso.ra, dec:qso.dec, $
               plate:qso.plate, fiber:qso.fiber, mjd:qso.mjd, $
               zqso:qso.z, err_zqso:qso.zerr, $
               spec_snr_median:stats.spec_snr_median, isitdecomposed:stats.isitdecomposed, $
               med_sdeviation_red:stats.med_sdeviation_red, med_sdeviation_blue:stats.med_sdeviation_blue}

delta_photo_outstr_i_203 = {imag_max:20.3, imag_min:15.0, wave:lambda, delta_ugriz:fltarr(nqso,5), delta_ugriz_dered:fltarr(nqso,5), $
               delta_ugriz_err:fltarr(nqso,5), $
               ra:qso.ra, dec:qso.dec, $
               plate:qso.plate, fiber:qso.fiber, mjd:qso.mjd, $
               zqso:qso.z, err_zqso:qso.zerr, $
               spec_snr_median:stats.spec_snr_median, isitdecomposed:stats.isitdecomposed, $
               med_sdeviation_red:stats.med_sdeviation_red, med_sdeviation_blue:stats.med_sdeviation_blue}

delta_norm_outstr_i_191_203 = {imag_max:20.3, imag_min:19.1, wave:allspec.wave, delta:fltarr(nqso,nwave), $
               delta_ivar:fltarr(nqso,nwave), $
               ra:qso.ra, dec:qso.dec, $
               plate:qso.plate, fiber:qso.fiber, mjd:qso.mjd, $
               zqso:qso.z, err_zqso:qso.zerr, $
               spec_snr_median:stats.spec_snr_median, isitdecomposed:stats.isitdecomposed, $
               med_sdeviation_red:stats.med_sdeviation_red, med_sdeviation_blue:stats.med_sdeviation_blue}

delta_sub_outstr_i_191_203 = {imag_max:20.3, imag_min:19.1, wave:allspec.wave, delta_subtracted:fltarr(nqso,nwave), $
               delta_subtracted_ivar:fltarr(nqso,nwave), $
               ra:qso.ra, dec:qso.dec, $
               plate:qso.plate, fiber:qso.fiber, mjd:qso.mjd, $
               zqso:qso.z, err_zqso:qso.zerr, $
               spec_snr_median:stats.spec_snr_median, isitdecomposed:stats.isitdecomposed, $
               med_sdeviation_red:stats.med_sdeviation_red, med_sdeviation_blue:stats.med_sdeviation_blue}

delta_photo_outstr_i_191_203 = {imag_max:20.3, imag_min:19.1, wave:lambda, delta_ugriz:fltarr(nqso,5), delta_ugriz_dered:fltarr(nqso,5), $
               delta_ugriz_err:qso_photo.ugriz_err, $
               ra:qso.ra, dec:qso.dec, $
               plate:qso.plate, fiber:qso.fiber, mjd:qso.mjd, $
               zqso:qso.z, err_zqso:qso.zerr, $
               spec_snr_median:stats.spec_snr_median, isitdecomposed:stats.isitdecomposed, $
               med_sdeviation_red:stats.med_sdeviation_red, med_sdeviation_blue:stats.med_sdeviation_blue}


;; get z from zbin
iwhichz = lonarr(nqso)
for iqso=0L, nqso-1L do begin
    tmp = min(abs(qso[iqso].z-outstr_i_203.z), imin)
    iwhichz[iqso] = imin
endfor

;; spectroscopy [parallel programming]
delta_norm_outstr_i_191.delta = allspec.flux/outstr_i_191.fluxmedian[iwhichz, *]
delta_norm_outstr_i_203.delta = allspec.flux/outstr_i_203.fluxmedian[iwhichz, *]
delta_norm_outstr_i_191_203.delta = allspec.flux/outstr_i_191_203.fluxmedian[iwhichz, *]

delta_sub_outstr_i_191.delta_subtracted = allspec.flux-outstr_i_191.fluxmedian[iwhichz, *]
delta_sub_outstr_i_203.delta_subtracted = allspec.flux-outstr_i_203.fluxmedian[iwhichz, *]
delta_sub_outstr_i_191_203.delta_subtracted = allspec.flux-outstr_i_191_203.fluxmedian[iwhichz, *]

;; photometry
delta_photo_outstr_i_191.delta_ugriz = transpose(qso_photo.ugriz)-phoout_i_191.ugriz_median[iwhichz, *]
delta_photo_outstr_i_203.delta_ugriz = transpose(qso_photo.ugriz)-phoout_i_203.ugriz_median[iwhichz, *]
delta_photo_outstr_i_191_203.delta_ugriz = transpose(qso_photo.ugriz)-phoout_i_191_203.ugriz_median[iwhichz, *]
delta_photo_outstr_i_191.delta_ugriz_dered = transpose(qso_photo.ugriz_dered)-phoout_i_191.ugriz_dered_median[iwhichz, *]
delta_photo_outstr_i_203.delta_ugriz_dered = transpose(qso_photo.ugriz_dered)-phoout_i_203.ugriz_dered_median[iwhichz, *]
delta_photo_outstr_i_191_203.delta_ugriz_dered = transpose(qso_photo.ugriz_dered)-phoout_i_191_203.ugriz_dered_median[iwhichz, *]

; quicklook
load_dp, /b

;  djs_plot, allspec.wave, smooth(outstr_i_191.fluxmean[iz,*], 5), xra=xra, xst=1, yra=yra, yst=1
;  djs_oplot, lambda, 10^(-0.4*(phoout_i_191.ugriz_median[1000,*]-8.9))/3.34E4/lambda^2./1.E-17, psym=4, symsize=8, color='red'
;  djs_oplot, allspec.wave, smooth(outstr_i_203.fluxmean[iz,*],5)
;  djs_oplot, allspec.wave, smooth(outstr_i_191_203.fluxmean[iz,*],5)
;  djs_oplot, allspec.wave, smooth(outstr_i_191.fluxmedian[iz,*],5), color='green'
;  djs_oplot, allspec.wave, smooth(outstr_i_203.fluxmedian[iz,*],5), color='green'
;  djs_oplot, allspec.wave, smooth(outstr_i_191_203.fluxmedian[iz,*],5), color='green'

stop
mwrfits, delta_norm_outstr_i_191, repstr(outfile, '.fits', '_normalized_191.fits'), /create
mwrfits, delta_norm_outstr_i_203, repstr(outfile, '.fits', '_normalized_203.fits'), /create
mwrfits, delta_norm_outstr_i_191_203, repstr(outfile, '.fits', '_normalized_191_203.fits'), /create

mwrfits, delta_sub_outstr_i_191, repstr(outfile, '.fits', '_subtracted_191.fits'), /create
mwrfits, delta_sub_outstr_i_203, repstr(outfile, '.fits', '_subtracted_203.fits'), /create
mwrfits, delta_sub_outstr_i_191_203, repstr(outfile, '.fits', '_subtracted_191_203.fits'), /create

mwrfits, delta_photo_outstr_i_191, repstr(outphotofile, '.fits', '_191.fits'), /create
mwrfits, delta_photo_outstr_i_203, repstr(outphotofile, '.fits', '_203.fits'), /create
mwrfits, delta_photo_outstr_i_191_203, repstr(outphotofile, '.fits', '_191_203.fits'), /create

end
