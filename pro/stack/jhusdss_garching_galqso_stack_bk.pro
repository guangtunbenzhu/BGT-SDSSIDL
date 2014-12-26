;pro jhusdss_garching_galqso_stack, nmfver, boss=boss, ivarweight=ivarweight

read,'BOSS? [1=yes, 0=no]: ',BOSS

nmfver = 106
ivarweigth = 1b
overwrite = 1b
saveall = 1b
saveall = 1b
snrcut = 5.
sdevcut = 0.07

;; ca II
ca_line_wave = [3934.79, 3969.59]
ca_xra = [3800., 4300.]
;; Na I
na_line_wave = [5890.00, 5896.00]
na_xra = [5700., 6100.]

line_wave = ca_line_wave
xra = ca_xra

;; radius bins
minrad = [0.050, 0.060, 0.070, 0.080, 0.090, 0.100, 0.120, 0.150]
maxrad = [0.100, 0.120, 0.140, 0.160, 0.180, 0.200, 0.240, 0.300]

;; galaxy redshift range
 zgalmin = 0.04
 zgalmax = 0.70

 rpmean = fltarr(n_elements(minrad))
 wmean = fltarr(n_elements(minrad))

; maxrad = minrad*3.
 print, 'radius (kpc): ', minrad*1E3
 print, 'radius (kpc): ', maxrad*1E3

choice_load_data = 0
read,'load data? [1=yes, 0=no]: ',choice_load_data

stackpath = jhusdss_get_path(/nmfqso)+'/'+string(nmfver, format='(I3.3)')+'/Stack/'
outfile = stackpath + jhusdss_garching_galqso_stack_filename(nmfver, boss=boss)

if (file_test(outfile) and ~overwrite) then begin
   splog, 'File already exists, use /overwrite to overwrite' 
   return
endif else begin
   splog, 'Will write into this file: '
   print, outfile
endelse

if choice_load_data eq 1 then begin

   ;; foreground galaxies
   garching_path = jhusdss_get_path(/garching)
   garching_file = garching_path+'/'+'gal_info_dr7_v5_2.fit.gz'
   gal = mrdfits(garching_file, 1)
   uniq_file = garching_path+'/'+'gal_uniq_dr7_v5_2.fits'
   galuniq = mrdfits(uniq_file, 1)

;  sfr_file = garching_path+'/'+'gal_totsfr_dr7_v5_2.fit.gz'
;  sfr = mrdfits(sfr_file, 1)
;  mass_file = garching_path+'/'+'totlgm_dr7_v5_2.fits.gz'
;  mass = mrdfits(mass_file, 1)
;  ssfr_file = garching_path+'/'+'gal_totspecsfr_dr7_v5_2.fits.gz'
;  ssfr = mrdfits(ssfr_file, 1)

   ;; background quasars
   qso = jhusdss_qso_readin(boss=boss)
   stat = jhusdss_qsostats_readin(nmfver, boss=boss)
   allspec = jhusdss_read_allqsospec(nmfver, /normresi, boss=boss)

   ;; match
   match0 = jhusdss_galqso_match_readin(boss=boss)

endif

;; only use caii \pm 200 AA
nwave = 600L
wave_shift = 200L
outwave = dblarr(nwave)
outstr = replicate({wave:outwave, fluxmean:fltarr(nwave), fluxmedian:fltarr(nwave), fluxgeomean:fltarr(nwave), $
                   npairs:0L, rp:0., ew_caii_3934:0., err_ew_caii_3934:0.}, n_elements(minrad))

for irad=0L, n_elements(minrad)-1L do begin

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
;   ssfr_tmp = ssfr[match.index_gal].avg
;   sfr_tmp = sfr[match.index_gal].avg
;   mass_tmp = mass[match.index_gal].avg
    zgal = gal[match.index_gal].z
    zuniq = galuniq[match.index_gal].choose
    zqso = qso[match.index_qso].z
    rp_tmp = match.rp_mpc*1E3

    ;; Now choose
    ;; 0.025 < zgal < 0.25
    ;; S/N(QSO) > 3
    ;; if CIV(QSO) < CaII(gal) < MgII(QSO) then sdev_red_tmp < 0.1
    ;; if Lya(QSO) < CaII(gal) < CIV(QSO) then sdev_blue_tmp < 0.1

    iall = where((zgal gt zgalmin) and (zgal lt zgalmax) $
             and (snr_tmp gt snrcut) $
             and zuniq $
;            and ((1600.*(zqso+1.) lt 3930.*(1.+zgal) and 2800.*(zqso+1.) gt 3980.*(1.+zgal) and sdev_red_tmp lt 0.1) $
             and ((1600.*(zqso+1.) lt (line_wave[0]-5.)*(1.+zgal) and sdev_red_tmp lt 0.1) $
              or  (1250.*(zqso+1.) lt (line_wave[0]-5.)*(1.+zgal) and 1500.*(zqso+1.) gt (line_wave[1]+5.)*(1.+zgal) and sdev_blue_tmp lt 0.1)), $
             nall)

    match = match[iall]
    nmatch = n_elements(match)

    print, 'npairs = ', nmatch
    print, '<rp> = ', median(rp_tmp)
    rpmean[irad] = median(rp_tmp)

    newzgal = gal[match.index_gal].z
    newzqso = qso[match.index_qso].z

    print, 'Get wavelength ...'

    icaii = value_locate(allspec.wave, 3934.78*(1.+newzgal))
    for i=0L, nwave-1L do outwave[i] = median(allspec.wave[icaii-wave_shift+i]/(1.+newzgal))
    outstr[irad].wave = outwave

    inresidual = fltarr(nmatch, nwave)
    inivar = fltarr(nmatch, nwave)+1.

    print, 'Get residuals ...'
    for i=0L, nwave-1L do inresidual[*,i] = allspec.residual[match.index_qso, icaii-wave_shift+i]

    print, 'Get ivar...'
    for i=0L, nwave-1L do inivar[*,i] = allspec.residual_ivar[match.index_qso, icaii-wave_shift+i]

    print, 'Make Composite ...'
    jhusdss_composite_engine, inresidual, inivar, fmean=fmean, fmedian=fmedian, $
       fgeomean=fgeomean, nobjuse=nobjuse, ivarweight=ivarweight

    tmp = moment(fmedian, sdev=sdev)

    outstr[irad].fluxmean = fmean
    outstr[irad].fluxmedian = fmedian
    outstr[irad].fluxgeomean = fgeomean

    fmedian=fmean

    ;; make the qa plot
    y = fmedian
    iuse = where(outwave gt 3800. and outwave lt 4300., nuse)
    terror = sqrt((moment(y[iuse]))[1])

    ;; fitting
    in_quadra = 0.
    in_slope = 0.
    in_intercept = median(1.-y[iuse])
    in_center = 3934.78
    in_separation = 3969.59-3934.78
    in_lflux = 0.2
    in_ratio = 0.5
    in_sigma = 2.0

    jhusdss_lowz_doublet_fit2, outwave[iuse], 1.-y[iuse], fltarr(nuse)+1./terror^2, $
       in_quadra, in_slope, in_intercept, in_center, in_separation, in_lflux, in_ratio, in_sigma, $
       quadra=quadra, slope=slope, intercept=intercept, center=center, separation=separation, $
       lflux=lflux, ratio=ratio, sigma=sigma, $
       err_quadra=err_quadra, err_slope=err_slope, err_intercept=err_intercept, err_center=err_center, $
       err_separation=err_separation, err_lflux=err_lflux, err_ratio=err_ratio, err_sigma=err_sigma, $
       maxwidth=maxwidth

    p = [quadra, slope, intercept, center, separation, lflux, ratio, sigma]
    yfit = jhusdss_lowz_doublet_func2(outwave, p)
    clevel = quadra*(outwave-center)^2+slope*(outwave-center)+intercept

    jhusdss_singlet_smooth, 1.-clevel-y, fltarr(n_elements(y))+1./terror^2, $
       outflux=single_newy, outivar=single_newyivar, sigma=2
    jhusdss_doublet_smooth, 1.-clevel-y, fltarr(n_elements(y))+1./terror^2, $
       outflux=newy, outivar=newyivar, sigma=2

    print, lflux, lflux/err_lflux
    wmean[irad] = lflux

    outstr[irad].rp = rpmean[irad]
    outstr[irad].npairs = nmatch
    outstr[irad].ew_caii_3934 = lflux
    outstr[irad].err_ew_caii_3934 = err_lflux

    load_dp, /b
    djs_plot, outwave, 1.-single_newy-clevel+clevel, psym=10, $
        xra=[3800., 4300], xst=1, yra=[1.-4.*terror, 1.+4.*terror]
    djs_oplot, 3934.*[1,1.], !y.crange[0]+[0.5,0.9]*(!y.crange[1]-!y.crange[0])
    djs_oplot, 3970.*[1,1.], !y.crange[0]+[0.5,0.9]*(!y.crange[1]-!y.crange[0])
    ii = where(outwave gt 3930 and outwave lt 3940)
    djs_oplot, outwave[ii], 1.-single_newy[ii]-clevel+clevel, color='red'
    ii = where(outwave gt 3965 and outwave lt 3975)
    djs_oplot, outwave[ii], 1.-single_newy[ii]-clevel+clevel, color='red'

;djs_oplot, outwave, 1.-yfit+clevel, psym=10, color='red'
;djs_oplot, outwave, 1.-newy-clevel+clevel, thick=thick, color='green'

;   djs_oplot, outwave, 1.-sky_fmedian+clevel, color='red'
;djs_oplot, outwave, smooth(sky_fmedian,5)-median(sky_fmedian)+1., color='red'
;djs_oplot, 3934.*[1,1.], !y.crange
;djs_oplot, 3970.*[1,1.], !y.crange
;djs_oplot, sky.wave/(1.+median(newzgal)), 1.-(1.-smooth(sky.fluxmedian, 5))/10., color='red'

;wait, 5
;djs_plot, outwave, smooth(fgeomean,8)-median(fgeomean)+1., xra=[3800., 4300], xst=1, yra=[1.-4.*sdev, 1.+4.*sdev]
;djs_oplot, 3934.*[1,1.], !y.crange
;djs_oplot, 3970.*[1,1.], !y.crange
;mwrfits, match, outfile, /create
;mwrfits, outstr, outfile
endfor

if saveall then mwrfits, outstr, outfile, /create

for i=0,2 do BEEP

end

