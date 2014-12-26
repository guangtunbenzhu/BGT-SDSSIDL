;pro jhusdss_garching_galqso_stack, nmfver, boss=boss, ivarweight=ivarweight

;read,'BOSS? [1=yes, 0=no]: ',BOSS
Coarse = 0b
DoWeight = 0b
read,'Coarse? [1=yes, 0=no]: ', Coarse
read,'DoWeight? [1=yes, 0=no]: ', DoWeight

nmfver = 106
ivarweigth = 1b
overwrite = 1b
savespec = 0b
savemean = 0b
snrcut = 2.
sdevcut = 0.10
sigma_cut = 5.
whatthewhat = 1

lambda=[3551., 4686., 6165., 7481., 8931.]

;; ca II
ca_line_wave = [3934.79, 3969.59]
;;ca_line_wave = [3935.00, 3969.00]
ca_xra = [3800., 4300.]
;; Na I
na_line_wave = [5890.00, 5896.00]
na_xra = [5700., 6100.]

ha_line_wave = [6563.00, 6584.00]
ha_xra = [6400., 6800.]

line_wave = ca_line_wave
xra = ca_xra

if (~Coarse) then begin
   minrad_tmp = 10.^(alog10(0.020)+findgen(23)*(0.5*alog10(1.5)))
   maxrad_tmp = minrad_tmp*1.5
   minrad = [[0.003, 0.010], minrad_tmp]
   maxrad = [[0.010, 0.020], maxrad_tmp]
   i_indep = [[0,1], lindgen(12)*2+2]
endif else begin
   minrad_tmp = 10.^(alog10(0.003)+findgen(13)*(0.5*alog10(3.0)))
   maxrad_tmp = minrad_tmp*3.0
   minrad = minrad_tmp
   maxrad = maxrad_tmp
;  minrad = [[0.003, 0.010], minrad_tmp]
;  maxrad = [[0.010, 0.020], maxrad_tmp]
;  i_indep = [[0,1], lindgen(11)*2+2]
   i_indep = lindgen(11)
endelse


fix_separation = long((line_wave[1]-line_wave[0])/(line_wave[1]+line_wave[0])/alog(10.)*2E+4)

;; galaxy redshift range
 zgalmin = 0.030
 zgalmax = 0.400
 wave_exclude = 4360.

;; quasar redshift range
 zqsomax=4.

 rpmean = fltarr(n_elements(minrad))
 wmean = fltarr(n_elements(minrad))

; maxrad = minrad*3.
 print, 'radius (kpc): ', minrad*1E3
 print, 'radius (kpc): ', maxrad*1E3

choice_load_data = 0
read,'load data? [1=yes, 0=no]: ',choice_load_data

stackpath = jhusdss_get_path(/nmfqso)+'/'+string(nmfver, format='(I3.3)')+'/Stack/'
outfile = stackpath + jhusdss_garching_galqso_stack_filename(nmfver, boss=boss)
gi_outfile = repstr(outfile, '.fits', '_photo.fits')

if (Coarse) then begin
   gi_outfile=repstr(gi_outfile,'.fits','_coarse.fits')
endif
if (DoWeight) then begin
   gi_outfile=repstr(gi_outfile,'.fits','_weight.fits')
endif

if (file_test(gi_outfile) and ~overwrite) then begin
   splog, 'File already exists, use /overwrite to overwrite' 
   return
endif else begin
   splog, 'Will write into this file: '
   print, gi_outfile
endelse

if choice_load_data eq 1 then begin

   ;; foreground galaxies
   garching_path = jhusdss_get_path(/garching)
   garching_file = garching_path+'/'+'gal_info_dr7_v5_2.fit.gz'
   gal = mrdfits(garching_file, 1)
   uniq_file = garching_path+'/'+'gal_uniq_dr7_v5_2.fits'
   galuniq = mrdfits(uniq_file, 1)

   sfr_file = garching_path+'/'+'gal_totsfr_dr7_v5_2.fits.gz'
   sfr = mrdfits(sfr_file, 1)
   mass_file = garching_path+'/'+'totlgm_dr7_v5_2.fit.gz'
   mass = mrdfits(mass_file, 1)
   ssfr_file = garching_path+'/'+'gal_totspecsfr_dr7_v5_2.fits.gz'
   ssfr = mrdfits(ssfr_file, 1)

   ;; background quasars
   print, 'Load SDSS DR7'
   qso = jhusdss_qso_readin(boss=boss)
   gi_infile = getenv('RAW_DATA')+'/SDSS/QSO/HW_dr7qso_newz_photo.fits'
   qso_gi = mrdfits(gi_infile, 1)
   stat = jhusdss_qsostats_readin(nmfver, boss=boss)
   qsophoto_infile = getenv('RAW_DATA')+'/SDSS/QSO/HW_dr7qso_newz_photo_more.fits'
   qso_photo = mrdfits(qsophoto_infile, 1)
   nqso = n_elements(qso)

   qsopath=jhusdss_get_path(/nmfqso)+'/'+string(nmfver, format='(I3.3)')+'/'
   infile = qsopath+'/AllInOne/'+jhusdss_allqsospec_filename(nmfver, boss=boss)
   infile = repstr(infile, '.fits', '_imag.fits')
   inphotofile = repstr(infile, '.fits', '_photo.fits')

   DoBright = 0b
   read,'DoBright? [1=yes, 0=no]: ', DoBright
   if (DoBright eq 1) then delta_photo = mrdfits(repstr(inphotofile, '.fits', '_191.fits'), 1) $
   else delta_photo = mrdfits(repstr(inphotofile, '.fits', '_203.fits'), 1)

   ;; match
   match0 = jhusdss_galqso_match_readin(boss=boss)
endif

;; double precision
outstr = replicate({delta_gi_mean:0.D0, delta_gi_median:0.D0, delta_uothers_mean:dblarr(4), delta_uothers_median:dblarr(4), delta_gothers_mean:dblarr(3), delta_gothers_median:dblarr(3), $
                    delta_rothers_mean:dblarr(2.), delta_rothers_median:dblarr(2.), delta_iz_mean:0.D0, delta_iz_median:0.D0, $
                    photo_wave:lambda, delta_ugriz_mean:dblarr(5), delta_ugriz_median:dblarr(5), delta_ugriz_dered_mean:dblarr(5), delta_ugriz_dered_median:dblarr(5), $
                    npairs:0L, rp:0., rp_min:0., rp_max:0.}, n_elements(minrad))

nrp = n_elements(minrad)
for irad=0L, nrp-1L do begin

    ;; sdss dr7 
    print, minrad[irad]*1E3, maxrad[irad]*1E3
    isub = where(match0.rp_mpc gt minrad[irad] $
             and match0.rp_mpc le maxrad[irad], nmatch)
    if nmatch eq 0 then message, "Can't find any pair within the annulus"
    match = match0[isub]
    nmatch = n_elements(match)

    ;; quality control
    sdev_red_tmp = stat[match.index_qso].med_sdeviation_red
    sdev_blue_tmp = stat[match.index_qso].med_sdeviation_blue
    snr_tmp = stat[match.index_qso].SPEC_SNR_MEDIAN
    ssfr_tmp = ssfr[match.index_gal].avg
    sfr_tmp = sfr[match.index_gal].avg
    mass_tmp = mass[match.index_gal].avg
    zgal = gal[match.index_gal].z
    zuniq = galuniq[match.index_gal].choose
    zqso = qso[match.index_qso].z
    rp_tmp = match.rp_mpc*1E3

    ;; Now choose
    ;; zgalmin < zgal < zgalmax
    ;; S/N(QSO) > 3
    ;; if CIV(QSO) < CaII(gal) < MgII(QSO) then sdev_red_tmp < 0.1
    ;; if Lya(QSO) < CaII(gal) < CIV(QSO) then sdev_blue_tmp < 0.1

    iall = where((zgal gt zgalmin) and (zgal lt zgalmax) $
;            and ((zgal gt (wave_exclude+15.)/(line_wave[0]-0.)-1.) or (zgal lt (wave_exclude-15.)/(line_wave[1]+0.)-1.)) $
;            and (snr_tmp gt snrcut) $
             and (mass_tmp gt 7. and mass_tmp lt 10. and sfr_tmp gt -3. and sfr_tmp lt 3.) $
             and zuniq $
             and zqso-zgal gt 0.1 and zqso lt zqsomax, $
             nall)

    match = match[iall]
    nmatch = n_elements(match)

    spec_snr = stat[match.index_qso].spec_snr_median
    newzgal = gal[match.index_gal].z
    newzqso = qso[match.index_qso].z
    newdelta_gi = qso_gi[match.index_qso].delta_g_i
    newimag = qso_photo[match.index_qso].ugriz_dered[3]
    newdelta_imag = delta_photo.delta_ugriz[match.index_qso,3]

    ;; ignore -99 choose between [-1, 1]
    iuse = where(newdelta_gi gt -5. and newdelta_gi lt 5. and newimag le 19.1 and newdelta_imag gt -5 and newdelta_imag lt 5)

    ; g-i
    outstr[irad].delta_gi_mean = mean(newdelta_gi[iuse], /double)
    outstr[irad].delta_gi_median = median(newdelta_gi[iuse], /double)

    ; delta ugriz
    outstr[irad].delta_ugriz_mean = mean(delta_photo.delta_ugriz[match[iuse].index_qso,*], /double, dimension=1)
    outstr[irad].delta_ugriz_median = median(delta_photo.delta_ugriz[match[iuse].index_qso,*], /double, dimension=1)
    outstr[irad].delta_ugriz_dered_mean = mean(delta_photo.delta_ugriz_dered[match[iuse].index_qso,*], /double, dimension=1)
    outstr[irad].delta_ugriz_dered_median = median(delta_photo.delta_ugriz_dered[match[iuse].index_qso,*], /double, dimension=1)

    print, 'npairs = ', nmatch
    print, '<rp> = ', median(rp_tmp)
    rpmean[irad] = median(match.rp_mpc*1E3)

    outstr[irad].rp = rpmean[irad]
    outstr[irad].rp_min = minrad[irad]
    outstr[irad].rp_max = maxrad[irad]
    outstr[irad].npairs = nmatch

endfor

load_dp, /b
;djs_plot, outstr.rp, outstr.delta_gi_mean, /xlog, /ylog, yra=[1E-2, 5E-1], yst=1, psym=4, thick=3
djs_plot, outstr.rp, outstr.delta_gi_mean-0.0459, /xlog, /ylog, yra=[1E-5, 5E-1], yst=1, psym=4, thick=3

icolor=['blue', 'green', 'dark green', 'magenta', 'red'] 
zeropoint = [0.055, 0.009, -0.017, -0.030, -0.033]
k_print, filename='delta_tmp.ps', axis_char_scale=1.3, xsize=8, ysize=6
  djs_plot, outstr.rp, outstr.delta_gi_mean-0.0459, /xlog, /ylog, yra=[1E-4, 2E-1], yst=1, psym=4, symsize=1.5, thick=6, $
      xthick=6, ythick=6, xtitle='r_p (kpc)', ytitle='<\Delta g-i>', charsize=1.5, charthick=3, color='red', xra=[20, 3000], xst=1, $
      ytickformat='jhusdss_tick_exponent'
  djs_plot, outstr.rp, outstr.delta_ugriz_dered_mean[0]-zeropoint[0], /xlog, /ylog, yra=[1E-4, 1E0], yst=1, psym=4, symsize=1.5, thick=6, $
          xthick=6, ythick=6, xtitle='r_p (kpc)', ytitle='<ugriz>', charsize=1.5, charthick=3, color='blue', xra=[20, 3000], xst=1, $
          ytickformat='jhusdss_tick_exponent'
  for jc=1L, 4L do  djs_oplot, outstr.rp, outstr.delta_ugriz_dered_mean[jc]-zeropoint[jc], psym=4, symsize=1.5, thick=6, color=icolor[jc]
k_end_print

djs_plot, outstr.rp, (outstr.delta_ugriz_dered_mean[1])-(outstr.delta_ugriz_dered_mean[3])-0.0451, /xlog, /ylog, yra=[1E-4, 1E-1], yst=1, psym=4, symsize=1.5, thick=6, $
          xthick=6, ythick=6, xtitle='r_p (kpc)', ytitle='<ugriz>', charsize=1.5, charthick=3, color='blue', xra=[10, 3000], xst=1, $
          ytickformat='jhusdss_tick_exponent'
djs_oplot, outstr.rp, outstr.delta_gi_mean-0.0418, psym=4, thick=3, color='red'

end

