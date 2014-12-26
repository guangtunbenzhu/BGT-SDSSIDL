;pro jhusdss_garching_gallrg_stack, nmfver, boss=boss, ivarweight=ivarweight

read,'BOSS? [1=yes, 0=no]: ',BOSS

lrgver = 101
ivarweigth = 1b
overwrite = 1b
saveall = 1b
snrcut = 1.
sdevcut = 0.10
;; ca II
ca_line_wave = [3934.79, 3969.59]
;;ca_line_wave = [3935.00, 3969.00]
ca_xra = [3800., 4300.]
;; Na I
na_line_wave = [5890.00, 5896.00]
na_xra = [5700., 6100.]

line_wave = ca_line_wave
xra = ca_xra

minrad_tmp = 10.^(alog10(0.020)+findgen(6)*(1.0*alog10(1.5)))
maxrad_tmp = minrad_tmp*3.0
minrad = [[0.003, 0.010], minrad_tmp]
maxrad = [[0.010, 0.020], maxrad_tmp]

;; galaxy redshift range
 zgalmin = 0.030
 zgalmax = 0.400

 rpmean = fltarr(n_elements(minrad))
 wmean = fltarr(n_elements(minrad))

; maxrad = minrad*3.
 print, 'radius (kpc): ', minrad*1E3
 print, 'radius (kpc): ', maxrad*1E3

choice_load_data = 0
read,'load data? [1=yes, 0=no]: ',choice_load_data

stackpath = jhusdss_get_path(/fitlrg)+'/'+string(lrgver, format='(I3.3)')+'/Stack/'
outfile = stackpath + jhusdss_garching_gallrg_stack_filename(lrgver, boss=boss)
spec_outfile = repstr(outfile, '.fits', '_spec.fits')

;; both SDSS and BOSS
;; outfile = repstr(outfile, '.fits', '_SDSS.fits')

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

   ;; background LRGs
   print, 'Load SDSS LRGs DR7'
;  stat = jhusdss_lrgstats_readin(lrgver, boss=boss)
   stat = jhusdss_read_alllrgspec(lrgver, boss=boss, /index)
   allwave = jhusdss_read_alllrgspec(lrgver, boss=boss, /wave)
   allspec = jhusdss_read_alllrgspec(lrgver, boss=boss, /normresi)

   ;; match
   match0 = jhusdss_gallrg_match_readin(boss=boss)
endif

;; only use caii \pm 200 AA
nwave = 600L
wave_shift = 200L
outwave = dblarr(nwave)
outstr = replicate({wave:outwave, fluxmean:fltarr(nwave), fluxmedian:fltarr(nwave), fluxgeomean:fltarr(nwave), $
                   npairs:0L, rp:0., rp_min:0., rp_max:0., ew_caii_3934:0., err_ew_caii_3934:0.}, n_elements(minrad))

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
    sdev_red_tmp = stat[match.index_lrg].med_sdeviation_red
    sdev_blue_tmp = stat[match.index_lrg].med_sdeviation_blue
    snr_tmp = stat[match.index_lrg].SPEC_SNR_MEDIAN
;   ssfr_tmp = ssfr[match.index_gal].avg
;   sfr_tmp = sfr[match.index_gal].avg
;   mass_tmp = mass[match.index_gal].avg
    zgal = gal[match.index_gal].z
    zuniq = galuniq[match.index_gal].choose
    zlrg = stat[match.index_lrg].zlrg
    rp_tmp = match.rp_mpc*1E3

    ;; Now choose
    ;; 0.025 < zgal < 0.25
    ;; S/N(QSO) > 3
    ;; if CIV(QSO) < CaII(gal) < MgII(QSO) then sdev_red_tmp < 0.1
    ;; if Lya(QSO) < CaII(gal) < CIV(QSO) then sdev_blue_tmp < 0.1

    iall = where((zgal gt zgalmin) and (zgal lt zgalmax) $
             and (snr_tmp gt snrcut) $
             and zuniq $
             and zlrg gt zgal+0.1  $
             and ((2400.*(zlrg+1.) lt (line_wave[0]-5.)*(1.+zgal) and (2700.*(zlrg+1.) gt (line_wave[1]-5.)*(1.+zgal) and sdev_blue_tmp lt sdevcut)) $
              or  (2900.*(zlrg+1.) lt (line_wave[0]-5.)*(1.+zgal) and (3650.*(zlrg+1.) gt (line_wave[1]+5.)*(1.+zgal) and sdev_blue_tmp lt sdevcut))), $
             nall)

    match = match[iall]
    nmatch = n_elements(match)

    newzgal = gal[match.index_gal].z
    newzlrg = stat[match.index_lrg].zlrg

    icaii0 = value_locate(allwave.wave, line_wave[0]*(1.+newzgal))
    jcaii0 = icaii0+1
    icaii = icaii0
    for k=0L, nmatch-1L do icaii[k] = (abs(allwave.wave[icaii0[k]]-line_wave[0]*(1.+newzgal[k])) > abs(allwave.wave[jcaii0[k]]-line_wave[0]*(1.+newzgal[k]))) ? jcaii0[k] : icaii0[k]

    print, 'npairs = ', nmatch
    print, '<rp> = ', median(rp_tmp)

    rpmean[irad] = median(rp_tmp)
    print, 'Get wavelength ...'
    for i=0L, nwave-1L do outwave[i] = median(allwave.wave[icaii-wave_shift+i]/(1.+newzgal))
    outstr[irad].wave = outwave

    inresidual = fltarr(nmatch, nwave)
    inivar = fltarr(nmatch, nwave)+1.

    print, 'Get residuals ...'
    for i=0L, nwave-1L do begin
        for j=0L, nmatch-1L do begin
            inresidual[j,i] = allspec[match[j].index_lrg].residual[icaii[j]-wave_shift+i]
        endfor
    endfor

    print, 'Get ivar...'
    for i=0L, nwave-1L do begin
        for j=0L, nmatch-1L do begin
            inivar[j,i] = allspec[match[j].index_lrg].residual_ivar[icaii[j]-wave_shift+i]
        endfor
    endfor

    ;; get weight
;   hist = histogram(newzgal, bin=0.01)
;   weight = 1.+(newzgal-zgalmin)^2*20.
    weight = fltarr(nmatch)+1.
    
    print, 'Make Composite ...'
    jhusdss_composite_engine, inresidual, inivar, fmean=fmean, fmedian=fmedian, $
       fgeomean=fgeomean, nobjuse=nobjuse, ivarweight=ivarweight, weight=weight

    tmp = moment(fmedian, sdev=sdev)

    outstr[irad].fluxmean = fmean
    outstr[irad].fluxmedian = fmedian
    outstr[irad].fluxgeomean = fgeomean

    ;; make the qa plot
;   y = fmedian
    y = fmean
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

    jhusdss_singlet_smooth, 1.-clevel-y, fltarr(n_elements(y))+1./terror^2, $
       outflux=single_newy, outivar=single_newyivar, sigma=1.0
;   single_newy = 1.-clevel-y
;   single_newyivar = fltarr(n_elements(y))+1./terror^2
    fix_separation = long((line_wave[1]-line_wave[0])/(line_wave[1]+line_wave[0])/alog(10.)*2E+4)
    jhusdss_doublet_smooth, 1.-clevel-y, fltarr(n_elements(y))+1./terror^2, $
       outflux=newy, outivar=newyivar, sigma=1.5, separation=fix_separation

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
    djs_plot, outwave, 1.-single_newy-clevel+clevel, psym=10, $
        xra=xra, xst=1, yra=[1.-4.*terror, 1.+4.*terror], position=[0.1, 0.5, 0.9, 0.9]
    djs_oplot, line_wave[0]*[1,1.], !y.crange[0]+[0.5,0.9]*(!y.crange[1]-!y.crange[0])
    djs_oplot, line_wave[1]*[1,1.], !y.crange[0]+[0.5,0.9]*(!y.crange[1]-!y.crange[0])
    ii = where(outwave gt (line_wave[0]-5.) and outwave lt (line_wave[0]+5.))
    djs_oplot, outwave[ii], 1.-single_newy[ii]-clevel+clevel, color='red'
    ii = where(outwave gt (line_wave[1]-5.) and outwave lt (line_wave[1]+5.))
    djs_oplot, outwave[ii], 1.-single_newy[ii]-clevel+clevel, color='red'
    djs_oplot, !x.crange, [1,1], linestyle=2

    djs_plot, outwave, 1.-newy-clevel+clevel, psym=10, $
        xra=xra, xst=1, yra=[1.-4.*terror, 1.+4.*terror], position=[0.1, 0.1, 0.9, 0.5], /noerase
    djs_oplot, line_wave[0]*[1,1.], !y.crange[0]+[0.5,0.9]*(!y.crange[1]-!y.crange[0])
;   djs_oplot, line_wave[1]*[1,1.], !y.crange[0]+[0.5,0.9]*(!y.crange[1]-!y.crange[0])
    ii = where(outwave gt (line_wave[0]-5.) and outwave lt (line_wave[0]+5.))
    djs_oplot, outwave[ii], 1.-newy[ii]-clevel+clevel, color='green'
    ii = where(outwave gt (line_wave[1]-5.) and outwave lt (line_wave[1]+5.))
    djs_oplot, outwave[ii], 1.-newy[ii]-clevel+clevel, color='green'
    djs_oplot, !x.crange, [1,1], linestyle=2
    djs_xyouts, line_wave[0], !y.crange[0]+0.3*(!y.crange[1]-!y.crange[0]), string(rpmean[irad], format='(f6.2)')

    if saveall then begin
       spec_outstr = {nrp:nrp, rp:rpmean[irad], rp_min:minrad[irad], rp_max:maxrad[irad], npairs:nmatch, $
           wave:outwave, residual:inresidual, ivar:inivar}
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
endfor

if saveall then mwrfits, outstr, outfile, /create

for i=0,2 do BEEP

end

