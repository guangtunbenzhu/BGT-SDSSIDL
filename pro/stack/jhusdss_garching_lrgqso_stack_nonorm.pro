;pro jhusdss_garching_galqso_stack, nmfver, boss=boss, ivarweight=ivarweight

;read,'BOSS? [1=yes, 0=no]: ',BOSS
Coarse = 0b
DoWeight = 0b
DoNaI = 0b
DoScale = 0b
DoLargeScale = 0b
Scale_factor = 0.2 ;; M* \propto M(halo)^1 \propto r(halo)^3 and r(halo) \propto M*^(1/3); Scale to 10^10.1 Msun
fid_mass = 10.3 ;; M* \propto M(halo)^1 \propto r(halo)^3 and r(halo) \propto M*^(1/3); Scale to 10^10.1 Msun
read,'Coarse? [1=yes, 0=no]: ', Coarse
read,'DoWeight? [1=yes, 0=no]: ', DoWeight
read,'DoNaI? [1=yes, 0=no]: ', DoNaI
read,'DoScale? [1=yes, 0=no]: ', DoScale
read,'DoLargeScale? [1=yes, 0=no]: ', DoLargeScale

nmfver = 106
ivarweigth = 1b
overwrite = 1b
savespec = 1b
savemean = 1b
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

mg_line_wave = [2796.35, 2803.53]
mg_xra = [2500., 3100.]

Fe_line_wave = [2586.65, 2600.17]
Fe_xra = [2500., 2900.]

line_wave = mg_line_wave
xra = mg_xra

if (DoNaI) then begin
   line_wave = na_line_wave
   xra = na_xra
end

if (~Coarse) then begin
   minrad_tmp = 10.^(alog10(0.030)+findgen(31)*(0.5*alog10(1.5)))
   maxrad_tmp = minrad_tmp*1.5
   minrad = [[0.005], minrad_tmp]
   maxrad = [[0.030], maxrad_tmp]
   maxrad[n_elements(maxrad)-1] = 18.
   if (DoLargeScale) then begin
      minrad = [18.]
      maxrad = [25.]
   endif
endif else begin
   minrad_tmp = 10.^(alog10(0.050)+findgen(9)*(0.5*alog10(3.)))
   maxrad_tmp = minrad_tmp*3.
   minrad = [[0.020], minrad_tmp]
   maxrad = [[0.050], maxrad_tmp]
;  minrad = [[0.003, 0.010], minrad_tmp]
;  maxrad = [[0.010, 0.020], maxrad_tmp]
;  i_indep = [[0,1], lindgen(11)*2+2]
;  i_indep = lindgen(11)
endelse

fix_separation = long((line_wave[1]-line_wave[0])/(line_wave[1]+line_wave[0])/alog(10.)*2E+4)

;; galaxy redshift range
 zgalmin = 0.40
 zgalmax = 0.750
 wave_exclude = 4360.

 rpmean = fltarr(n_elements(minrad))
 wmean = fltarr(n_elements(minrad))

; maxrad = minrad*3.
 print, 'radius (kpc): ', minrad*1E3
 print, 'radius (kpc): ', maxrad*1E3

choice_load_data = 0
read,'load data? [1=yes, 0=no]: ',choice_load_data

stackpath = jhusdss_get_path(/nmfqso)+'/'+string(nmfver, format='(I3.3)')+'/Stack/'
outfile = stackpath + jhusdss_garching_lrgqso_stack_filename(nmfver, boss=boss)
spec_outfile = repstr(outfile, '.fits', '_spec_nonorm.fits')
outfile = repstr(outfile, '.fits', '_nonorm.fits')

if (Coarse) then begin
   outfile=repstr(outfile,'.fits','_coarse.fits')
   spec_outfile=repstr(spec_outfile,'.fits','_coarse.fits')
endif
if (DoWeight) then begin
   outfile=repstr(outfile,'.fits','_weight.fits')
   spec_outfile=repstr(spec_outfile,'.fits','_weight.fits')
endif
if (DoNaI) then begin
   outfile=repstr(outfile,'.fits','_NaI.fits')
   spec_outfile=repstr(spec_outfile,'.fits','_NaI.fits')
endif
if (DoScale) then begin
   outfile=repstr(outfile,'.fits','_Scale.fits')
   spec_outfile=repstr(spec_outfile,'.fits','_Scale.fits')
endif
if (DoLargeScale) then begin
   outfile=repstr(outfile,'.fits','_LargeScale.fits')
   spec_outfile=repstr(spec_outfile,'.fits','_LargeScale.fits')
endif


;; both SDSS and BOSS
;; outfile = repstr(outfile, '.fits', '_SDSS.fits')

if (file_test(outfile) and ~overwrite) then begin
   splog, 'File already exists, use /overwrite to overwrite' 
   return
endif else begin
   splog, 'Will write into this file: '
   print, outfile
   print, spec_outfile
endelse

if choice_load_data eq 1 then begin

   ;; foreground galaxies
   gal = jhusdss_lrg_readin(/boss)

   bosslrg_path = jhusdss_get_path(/bosslrg)
;  garching_file = garching_path+'/'+'gal_info_dr7_v5_2.fit.gz'
;  gal = mrdfits(garching_file, 1)
   uniq_file = bosslrg_path+'/'+'uniq_wisconsin_pca_bc03-v5_6_0.fits.gz'
   galuniq = mrdfits(uniq_file, 1)

;; calculate fraction of young stars

;  sfr_file = garching_path+'/'+'gal_totsfr_dr7_v5_2.fits.gz'
;  sfr = mrdfits(sfr_file, 1)
;  mass_file = garching_path+'/'+'totlgm_dr7_v5_2.fit.gz'
;  mass = mrdfits(mass_file, 1)
;  ssfr_file = garching_path+'/'+'gal_totspecsfr_dr7_v5_2.fits.gz'
;  ssfr = mrdfits(ssfr_file, 1)

   ;; background quasars
   print, 'Load SDSS DR7'
   qso = jhusdss_qso_readin(boss=boss)
   stat = jhusdss_qsostats_readin(nmfver, boss=boss)
   allspec = jhusdss_read_allqsospec(nmfver, /flux, boss=boss)
;  allspec = jhusdss_read_allqsospec(nmfver, /subtresi, boss=boss)
   qsophoto_infile = getenv('RAW_DATA')+'/SDSS/QSO/HW_dr7qso_newz_photo_more.fits'
   qso_photo = mrdfits(qsophoto_infile, 1)
   nqso = n_elements(qso)

   qsopath=jhusdss_get_path(/nmfqso)+'/'+string(nmfver, format='(I3.3)')+'/'
   infile = qsopath+'/AllInOne/'+jhusdss_allqsospec_filename(nmfver, boss=boss)
   infile = repstr(infile, '.fits', '_imag.fits')
   inphotofile = repstr(infile, '.fits', '_photo.fits')
   delta_norm_outstr_i_203 = mrdfits(repstr(infile, '.fits', '_normalized_203.fits'), 1)
   delta_photo_outstr_i_203 = mrdfits(repstr(inphotofile, '.fits', '_203.fits'), 1)
   delta_norm_outstr_i_191 = mrdfits(repstr(infile, '.fits', '_normalized_191.fits'), 1)
   delta_photo_outstr_i_191 = mrdfits(repstr(inphotofile, '.fits', '_191.fits'), 1)

   fluxfile = qsopath+'/'+jhusdss_qso_composite_image_filename(nmfver, boss=boss)
   zbininfile = repstr(fluxfile, '.fits', '_zbin.fits')
   zbininphotofile = repstr(fluxfile, '.fits', '_photo_zbin.fits')
   outstr_i_191 = mrdfits(repstr(zbininfile, '.fits', '_191.fits'), 1)
   outstr_i_203 = mrdfits(repstr(zbininfile, '.fits', '_203.fits'), 1)
   outstr_i_191_203 = mrdfits(repstr(zbininfile, '.fits', '_191_203.fits'), 1)
   phoout_i_191 = mrdfits(repstr(zbininphotofile, '.fits', '_191.fits'), 1)
   phoout_i_203 = mrdfits(repstr(zbininphotofile, '.fits', '_203.fits'), 1)
   phoout_i_191_203 = mrdfits(repstr(zbininphotofile, '.fits', '_191_203.fits'), 1)

   ;; match
   match0 = jhusdss_lrgqso_match_readin(boss=boss)
endif
rp_scale = match0.rp_mpc*10.^(Scale_factor*(fid_mass-gal[match0.index_lrg].mstellar_median))
stop

;; only use caii \pm 200 AA
nwave = 3200L
wave_shift = 400L
if (DoLargeScale) then begin
    nwave = 400L 
    wave_shift = 200L
endif
outwave = dblarr(nwave)
outstr = replicate({wave:outwave, fluxmean:fltarr(nwave), fluxmedian:fltarr(nwave), fluxgeomean:fltarr(nwave), $
                    photo_wave:lambda, delta_ugriz_mean:fltarr(5), delta_ugriz_median:fltarr(5), delta_ugriz_dered_mean:fltarr(5), delta_ugriz_dered_median:fltarr(5), $
                   npairs:0L, rp:0., rp_min:0., rp_max:0., ew_caii_3934:0., err_ew_caii_3934:0.}, n_elements(minrad))

nrp = n_elements(minrad)
;for irad=0L, nrp-1L do begin
;for irad=0L, nrp-1L do begin
for irad=0L, 6L do begin

    ;; sdss dr7 
    print, minrad[irad]*1E3, maxrad[irad]*1E3
    isub = where(match0.rp_mpc gt minrad[irad] $
             and match0.rp_mpc le maxrad[irad], nmatch)
    if (DoScale) then begin
       isub = where(rp_scale gt minrad[irad] $
                and rp_scale le maxrad[irad], nmatch)
    endif
;   stop
    if nmatch eq 0 then message, "Can't find any pair within the annulus"

    match = match0[isub]
    nmatch = n_elements(match)

    ;; quality control
    sdev_red_tmp = stat[match.index_qso].med_sdeviation_red
    sdev_blue_tmp = stat[match.index_qso].med_sdeviation_blue
    snr_tmp = stat[match.index_qso].SPEC_SNR_MEDIAN
;   ssfr_tmp = ssfr[match.index_lrg].avg
;   sfr_tmp = sfr[match.index_lrg].avg
;   mass_tmp = mass[match.index_lrg].avg
    mass_tmp = gal[match.index_lrg].mstellar_median
    zgal = gal[match.index_lrg].z
    zuniq = galuniq[match.index_lrg].choose
    zqso = qso[match.index_qso].z
    rp_tmp = match.rp_mpc*1E3
    iqso = qso_photo[match.index_qso].ugriz_dered[3]
    if (DoScale) then rp_tmp = rp_scale[isub]*1E3

    ;; Now choose
    ;; 0.025 < zgal < 0.25
    ;; S/N(QSO) > 3
    ;; if CIV(QSO) < CaII(gal) < MgII(QSO) then sdev_red_tmp < 0.1
    ;; if Lya(QSO) < CaII(gal) < CIV(QSO) then sdev_blue_tmp < 0.1

    iall = where((zgal gt zgalmin) and (zgal lt zgalmax) $
;            and ((zgal gt (wave_exclude+15.)/(line_wave[0]-0.)-1.) or (zgal lt (wave_exclude-15.)/(line_wave[1]+0.)-1.)) $
             and (snr_tmp gt snrcut) $
             and (mass_tmp gt 7. and mass_tmp lt 13.) $
    ;        and (mass_tmp gt 7. and mass_tmp lt 13. and sfr_tmp gt -3. and sfr_tmp lt 3.) $
             and zuniq $
             and iqso lt 19.1 and iqso gt 15. $
             and zqso-zgal gt 0.1 and zqso lt 2.3, $
;            and ((1600.*(zqso+1.) lt 3930.*(1.+zgal) and 2800.*(zqso+1.) gt 3980.*(1.+zgal) and sdev_red_tmp lt 0.1) $
    ;        and ((2830.*(zqso+1.) lt (line_wave[0]-5.)*(1.+zgal) and 3920.*(zqso+1.) gt (line_wave[1]+5.)*(1.+zgal) and sdev_red_tmp lt 0.1) $
    ;        and ((1930.*(zqso+1.) lt (line_wave[0]-5.)*(1.+zgal) and 2770.*(zqso+1.) gt (line_wave[1]+5.)*(1.+zgal) and sdev_red_tmp lt 0.1) $
    ;         or  (1580.*(zqso+1.) lt (line_wave[0]-5.)*(1.+zgal) and 1870.*(zqso+1.) gt (line_wave[1]+5.)*(1.+zgal) and sdev_red_tmp lt 0.1) $
    ;         or  (1250.*(zqso+1.) lt (line_wave[0]-5.)*(1.+zgal) and 1520.*(zqso+1.) gt (line_wave[1]+5.)*(1.+zgal) and sdev_blue_tmp lt 0.1)), $
             nall)

    match = match[iall]
    nmatch = n_elements(match)

    spec_snr = stat[match.index_qso].spec_snr_median
    newzgal = gal[match.index_lrg].z
    newzqso = qso[match.index_qso].z

;   plothist, newzqso, peak=1, bin=0.01, color=djs_icolor('red')
;   stop
;   continue

    icaii0 = value_locate(delta_norm_outstr_i_191.wave, line_wave[0]*(1.+newzgal))
    jcaii0 = icaii0+1
    icaii = icaii0
    for k=0L, nmatch-1L do icaii[k] = (abs(delta_norm_outstr_i_191.wave[icaii0[k]]-line_wave[0]*(1.+newzgal[k])) > abs(delta_norm_outstr_i_191.wave[jcaii0[k]]-line_wave[0]*(1.+newzgal[k]))) ? jcaii0[k] : icaii0[k]

    print, 'npairs = ', nmatch
    print, '<rp> = ', median(rp_tmp[iall])
    print, '<z> = ', median(newzgal)

    ;; get weight
;   hist = histogram(newzgal, bin=0.01)
;   weight = 1.+(newzgal-zgalmin)^2*20.
    weight = fltarr(nmatch)+1.
;   weight = newzgal^2*100.
    if (DoWeight) then weight=1./rp_tmp[iall]*median(rp_tmp[iall])
 
;   rpmean[irad] = median(match.rp_mpc*1E3)
;   rpmean[irad] = median(rp_tmp[iall])
    rpmean[irad] = total(rp_tmp[iall]*weight)/total(weight)
;   if irad eq 5 then stop
    print, '<rp_weighted> = ', rpmean[irad]
    print, 'Get wavelength ...'
    for i=0L, nwave-1L do outwave[i] = median(delta_norm_outstr_i_191.wave[icaii-wave_shift+i]/(1.+newzgal))
    outstr[irad].wave = outwave

    inresidual = fltarr(nmatch, nwave)
    inivar = fltarr(nmatch, nwave)+1.

    print, 'Get residuals ...'
    for i=0L, nwave-1L do inresidual[*,i] = delta_norm_outstr_i_191.delta[match.index_qso, icaii-wave_shift+i]
;   for i=0L, nwave-1L do inresidual[*,i] = allspec.subtracted_residual[match.index_qso, icaii-wave_shift+i]
    print, 'Get ivar...'
    for i=0L, nwave-1L do inivar[*,i] = allspec.ivar[match.index_qso, icaii-wave_shift+i]
;   for i=0L, nwave-1L do inivar[*,i] = allspec.subtracted_ivar[match.index_qso, icaii-wave_shift+i]

;   stop
    print, 'Make Composite ...'
    jhusdss_composite_engine, inresidual, inivar, fmean=fmean, fmedian=fmedian, $
       fgeomean=fgeomean, nobjuse=nobjuse, ivarweight=ivarweight, weight=weight, sigma_cut=sigma_cut
    outstr[irad].fluxmean = fmean
    outstr[irad].fluxmedian = fmedian
    outstr[irad].fluxgeomean = fgeomean

    ;; make the qa plot
;   y = fmedian
    y = fgeomean
    med_continuum = median(y, 71, /even)
    y = (y/med_continuum)

;   photo
    inphoto_resi = delta_photo_outstr_i_191.DELTA_UGRIZ[match.index_qso, *]
    inphoto_ivar = transpose(1./qso_photo[match.index_qso].ugriz_err^2)
    jhusdss_composite_engine, inphoto_resi, inphoto_ivar, fmean=fmean, fmedian=fmedian, $
       fgeomean=fgeomean, nobjuse=nobjuse, ivarweight=ivarweight, weight=weight, sigma_cut=sigma_cut
    outstr[irad].delta_ugriz_mean = fmean
    outstr[irad].delta_ugriz_median = fmedian

    inphoto_resi = delta_photo_outstr_i_191.DELTA_UGRIZ_DERED[match.index_qso, *]
    inphoto_ivar = transpose(1./qso_photo[match.index_qso].ugriz_err^2)
    jhusdss_composite_engine, inphoto_resi, inphoto_ivar, fmean=fmean, fmedian=fmedian, $
       fgeomean=fgeomean, nobjuse=nobjuse, ivarweight=ivarweight, weight=weight, sigma_cut=sigma_cut
    outstr[irad].delta_ugriz_dered_mean = fmean
    outstr[irad].delta_ugriz_dered_median = fmedian

    iuse = where(outwave gt xra[0] and outwave lt xra[1], nuse)
    terror = sqrt((moment(y[iuse]))[1])

    ;; fitting
    in_quadra = 0.
    in_slope = 0.
    in_intercept = median(1.-y[iuse])
    in_center = line_wave[0]
    in_separation = line_wave[1]-line_wave[0]
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

    newy = y+clevel
    ;; in pixel
    jhusdss_singlet_smooth, 1.-newy, fltarr(n_elements(y))+1./terror^2, $
       outflux=single_newy, outivar=single_newyivar, /normalize, sigma=2.0, factor_norm=factor_norm1
;   single_newy = 1.-clevel-y
;   single_newyivar = fltarr(n_elements(y))+1./terror^2
    jhusdss_doublet_smooth, 1.-newy, fltarr(n_elements(y))+1./terror^2, $
       outflux=double_newy, outivar=double_newyivar, /normalize, sigma=2.0, separation=fix_separation, factor_norm=factor_norm2

    print, lflux, lflux/err_lflux
    wmean[irad] = lflux

;   outstr[irad].fluxmean = fmean
;   outstr[irad].fluxmedian = fmedian
;   outstr[irad].fluxgeomean = fgeomean

    outstr[irad].rp = rpmean[irad]
    outstr[irad].rp_min = minrad[irad]
    outstr[irad].rp_max = maxrad[irad]
    outstr[irad].npairs = nmatch
    outstr[irad].ew_caii_3934 = lflux
    outstr[irad].err_ew_caii_3934 = err_lflux

    load_dp, /b

; djs_plot, outwave, fmean
; djs_plot, outwave, y    
; djs_plot, outwave, fmean, xra=[3700, 4000], xst=1  

; i++ & djs_plot, outwave, inresidual[i,*], xra=[3700, 4000], xst=1
; djs_oplot, 3727.*[1,1], !y.crange, color='red'
; djs_oplot, line_wave[0]*[1,1.], !y.crange[0]+[0.5,0.9]*(!y.crange[1]-!y.crange[0]), color='red'
; djs_oplot, line_wave[1]*[1,1.], !y.crange[0]+[0.5,0.9]*(!y.crange[1]-!y.crange[0]), color='red'

; if (keyword_set(whatthewhat)) then begin
    djs_plot, outwave, 1.-single_newy, psym=10, $
        xra=xra, xst=1, yra=[1.-3.*terror, 1.+3.*terror], position=[0.1, 0.5, 0.9, 0.9]
    djs_oplot, line_wave[0]*[1,1.], !y.crange[0]+[0.5,0.9]*(!y.crange[1]-!y.crange[0])
    djs_oplot, line_wave[1]*[1,1.], !y.crange[0]+[0.5,0.9]*(!y.crange[1]-!y.crange[0])
;   djs_oplot, 3727.*[1,1.], !y.crange[0]+[0.5,0.9]*(!y.crange[1]-!y.crange[0])
    ii = where(outwave gt (line_wave[0]-5.) and outwave lt (line_wave[0]+5.))
    djs_oplot, outwave[ii], 1.-single_newy[ii], color='red'
    ii = where(outwave gt (line_wave[1]-5.) and outwave lt (line_wave[1]+5.))
    djs_oplot, outwave[ii], 1.-single_newy[ii], color='red'
    djs_oplot, !x.crange, [1,1], linestyle=2

    djs_plot, outwave, 1.-double_newy, psym=10, $
        xra=xra, xst=1, yra=[1.-3.*terror, 1.+3.*terror], position=[0.1, 0.1, 0.9, 0.5], /noerase
    djs_oplot, line_wave[0]*[1,1.], !y.crange[0]+[0.5,0.9]*(!y.crange[1]-!y.crange[0])
;   djs_oplot, line_wave[1]*[1,1.], !y.crange[0]+[0.5,0.9]*(!y.crange[1]-!y.crange[0])
    ii = where(outwave gt (line_wave[0]-5.) and outwave lt (line_wave[0]+5.))
    djs_oplot, outwave[ii], 1.-double_newy[ii], color='green'
    ii = where(outwave gt (line_wave[1]-5.) and outwave lt (line_wave[1]+5.))
    djs_oplot, outwave[ii], 1.-double_newy[ii], color='green'
    djs_oplot, !x.crange, [1,1], linestyle=2
    djs_xyouts, line_wave[0], !y.crange[0]+0.3*(!y.crange[1]-!y.crange[0]), string(rpmean[irad], format='(f6.2)')

    if savespec then begin
;      spec_outstr = {nrp:nrp, rp:rpmean[irad], rp_min:minrad[irad], rp_max:maxrad[irad], npairs:nmatch, $
;          wave:outwave, residual:inresidual, ivar:inivar, spec_snr:spec_snr, zgal:newzgal, zqso:newzqso}
       spec_outstr = {nrp:nrp, rp:rpmean[irad], rp_min:minrad[irad], rp_max:maxrad[irad], npairs:nmatch, $
           index_qso:match.index_qso, index_lrg:match.index_lrg, zgal:newzgal, zqso:newzqso}
       mwrfits, spec_outstr, spec_outfile, create=(irad eq 0)
    endif
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
;endif

;stop
endfor

if savemean then mwrfits, outstr, outfile, /create


for i=0,2 do BEEP
djs_plot, outstr.rp, outstr.ew_caii_3934, /xlog, /ylog, yra=[1E-5, 2E0], yst=1, psym=4, thick=3

stop
;!p.multi = [0, 1, 6]
;; composite spectra of background quasars

load_lrg_data = 0
read,'load LRG data? [1=yes, 0=no]: ',load_lrg_data
if (load_lrg_data eq 1) then begin
   lrg_composite_zbin = mrdfits('/home/gz323/DATA/SDSS/LRG/101/LRG_composite_image_101.fits', 1)

    filterlist = primus_sdss_filterlist()
    weff = primus_lambda_eff(filterlist=filterlist)
    bandnames = ['u','g','r','i','z']
    filter_colors = ['blue', 'green', 'yellow', 'orange', 'red']
    k_load_filters, filterlist, filter_nlambda, filter_lambda, filter_pass
    fboost=[45., 20., 20., 20., 50.]
endif

thick=8
xthick=8
ythick=8
charsize=1.4
charthick=3

qapath=jhusdss_get_path(/nmfqso)+'/'+string(nmfver, format='(I3.3)')
qapath = qapath+'/QAplots'

zcolor = [djs_icolor('black'), djs_icolor('pink'), djs_icolor('navy'), $
            djs_icolor('orange'), djs_icolor('blue'), djs_icolor('magenta red'), $
            djs_icolor('dark green'), djs_icolor('magenta'), djs_icolor('green'),  $
            djs_icolor('brown'), djs_icolor('dark cyan'), djs_icolor('red')]

psfile=qapath+'/Magnification_Extinction_Absorption_FG.LRG_BG.QSO.ps'
k_print, filename=psfile, xsize=12, ysize=14, axis_char_scale=1.3

nmax = 6
loadct, 4
for i=0L, nmax-1L do zcolor[i] = (i+1)/float(nmax)*200.-1.

xra = [3200, 10000]
yra = [0.95, 1.2]
ytitle = 'f_\lambda/mean(f_\lambda)'
djs_plot, outstr[0].wave*1.52, outstr[0].fluxmedian, xra=xra, yra=yra, xst=1, yst=1, ytitle=ytitle, $
    thick=thick, xthick=xthick, ythick=ythick, /nodata, position=[0.15, 0.28, 0.95, 0.73], xtickformat='(A1)', charsize=charsize-0.2, charthick=charthick
for i=0L, nmax-1L do djs_oplot, outstr[i].wave*1.52, smooth(outstr[i].fluxmedian, 10), thick=thick, color=zcolor[i+1]
for i=0L, nmax-1L do djs_oplot, outstr[i].photo_wave, 10.^(-0.4*(outstr[i].delta_ugriz_median)), psym=4, symsize=3, color=zcolor[i+1], thick=thick
for i=0L, nmax-1L do djs_oplot, outstr[i].photo_wave, 10.^(-0.4*(outstr[i].delta_ugriz_median)), psym=4, symsize=3, color=zcolor[i+1], thick=thick
djs_oplot, !x.crange, [1,1], thick=thick, color='grey', linestyle=1

yra = [0,18]
xtitle = '\lambda_{obs} (\AA)'
ytitle = 'f_\lambda'
tmp = min(abs(1.8-outstr_i_191.z), izqso)
djs_plot, outstr_i_191.wave, outstr_i_191.fluxmedian[izqso,*], xra=xra, yra=yra, xst=1, yst=1, ytitle=ytitle, $
    thick=thick, xthick=xthick, ythick=ythick, position=[0.15, 0.74, 0.95, 0.9], /noerase, xtickformat='(A1)', charsize=charsize, charthick=charthick, color='blue'
djs_axis, xaxis=1, xtitle=xtitle, charsize=charsize, charthick=charthick

yra = [0,13]
tmp = min(abs(0.52-lrg_composite_zbin.z), izlrg)
djs_plot, lrg_composite_zbin.wave, smooth(lrg_composite_zbin.IMAGE_MEDIAN[*,izlrg],10)*10., xra=xra, yra=yra, xst=1, yst=1, xtitle=xtitle, ytitle=ytitle, $
    thick=thick, xthick=xthick, ythick=ythick, position=[0.15, 0.1, 0.95, 0.27], charsize=charsize, charthick=charthick, /noerase, color='red'

for i=0L, 4L do djs_oplot, filter_lambda[*,i], filter_pass[*,i]*fboost[i], thick=thick-2, color=filter_colors[i]

k_end_print

end


