;; This should be just a test version
;pro jhusdss_garching_galqso_match_spec_byrad_docom, nmfver, boss=boss, overwrite=overwrite, $
;       minrad=minrad, maxrad=maxrad

 boss = 0b
 nmfver = 106
 minrad = [0.003, 0.010, 0.020, 0.040, 0.080];, 0.160, 0.320, 0.640];, 1.280, 2.560];, 0.270]
 maxrad = [0.010, 0.020, 0.040, 0.080, 0.160];, 0.320, 0.640, 1.280];, 2.560, 5.120];, 0.810]
;minrad = [0.003, 0.010, 0.015, 0.0225, 0.03375, 0.050625, 0.0759375, 0.113906, 0.170859];, 1.280, 2.560];, 0.270]
;maxrad = [0.010, 0.015, 0.0225, 0.03375, 0.050625, 0.0759375, 0.113906, 0.170859, 0.256289];, 2.560, 5.120];, 0.810]
 zgalmin = 0.025
 zgalmax = 0.25
 rpmean = fltarr(n_elements(minrad))
 wmean = fltarr(n_elements(minrad))

 saveall = 1b

; maxrad = minrad*3.
print, 'radius (kpc): ', minrad*1E3, maxrad*1E3

ivarweight = 1

choice_load_data = 0
read,'load data? [1=yes, 0=no]: ',choice_load_data

;if (n_elements(nmfver) eq 0) then message, 'nmfver required'
;if (n_elements(minrad) eq 0) then minrad = 0.
;if (n_elements(maxrad) eq 0) then maxrad = 0.1
;if (minrad ge maxrad) then message, "minrad can't be larger than maxrad"

if choice_load_data eq 1 then begin

if (keyword_set(boss)) then begin
   lowzpath=jhusdss_get_path(/nmfqso)+'/'+string(nmfver, format='(I3.3)')+'/Lowz_BOSS'
endif else begin
   lowzpath=jhusdss_get_path(/nmfqso)+'/'+string(nmfver, format='(I3.3)')+'/Lowz'
endelse


;outfile = lowzpath+'/Rad'+string(minrad,format='(f4.2)')+'_'+string(maxrad,format='(f4.2)')+'_' $
;        + jhusdss_lowz_spec_filename(nmfver)

;if (file_test(outfile) and ~keyword_set(overwrite)) then begin
;   splog, 'File already exists. Use /overwrite if you want to overwrite it.'
;   return
;endif else begin
;   splog, 'Will write the composite into this file: '
;   print, outfile
;endelse

garching_path = jhusdss_get_path(/garching)
garching_file = garching_path+'/'+'gal_info_dr7_v5_2.fit.gz'
gal = mrdfits(garching_file, 1)
uniq_file = garching_path+'/'+'gal_uniq_dr7_v5_2.fits'
galuniq = mrdfits(uniq_file, 1)
sfr_file = garching_path+'/'+'gal_totsfr_dr7_v5_2.fits.gz'
sfr = mrdfits(sfr_file, 1)
mass_file = garching_path+'/'+'totlgm_dr7_v5_2.fits.gz'
mass = mrdfits(mass_file, 1)
ssfr_file = garching_path+'/'+'gal_totspecsfr_dr7_v5_2.fits.gz'
ssfr = mrdfits(ssfr_file, 1)

;; read in SDSS qso catalog, need ra, dec
qso = jhusdss_qso_readin(boss=boss)
stat = jhusdss_qsostats_readin(nmfver, boss=boss)
allspec = jhusdss_read_allqsospec(nmfver, boss=boss)

sky_file = jhusdss_get_path(/nmfqso)+'/'+ string(nmfver, format='(I3.3)')+$
                       '/Composite/Residual_Composite_'+ string(nmfver, format='(I3.3)')+'.fits'
sky = mrdfits(sky_file, 1)
combine1fiber, alog10(sky.wave), sky.fluxmedian, fltarr(n_elements(sky.wave))+1., newloglam=alog10(allspec.wave), $
    newflux=skyflux, newivar=skyivar, maxiter=0

match0 = jhusdss_galqso_match_readin(boss=boss)

endif

nwave = 600L
wave_shift = 200L
outwave = dblarr(nwave)
outstr = replicate({wave:outwave, fluxmean:fltarr(nwave), fluxmedian:fltarr(nwave), $
                   npairs:0L, rp:0., ew_caii_3934:0., err_ew_caii_3934:0.}, n_elements(minrad))

for irad=0L, n_elements(minrad)-1L do begin

print, minrad[irad]*1E3, maxrad[irad]*1E3
isub = where(match0.rp_mpc gt minrad[irad] $
         and match0.rp_mpc le maxrad[irad], nmatch)
if nmatch eq 0 then message, "Can't find any pair within the annulus"
match = match0[isub]
nmatch = n_elements(match)

    sdev_tmp = stat[match.index_qso].med_sdeviation_red
    snr_tmp = stat[match.index_qso].SPEC_SNR_MEDIAN
;   ssfr_tmp = ssfr[match.index_gal].avg
;   sfr_tmp = sfr[match.index_gal].avg
;   mass_tmp = mass[match.index_gal].avg
    zgal = gal[match.index_gal].z
    zuniq = galuniq[match.index_gal].choose
    zqso = qso[match.index_qso].z
    rp_tmp = match.rp_mpc*1E3
;   flag_tmp = 

    iall = where(sdev_tmp gt 0. and sdev_tmp le 0.07 $
             and snr_tmp gt 3.0 $
             and zgal gt zgalmin $
             and zgal lt zgalmax $
             and zgal lt zqso-0.1 $
;            and zqso gt 0.43 $
;            and zqso lt 2.8 $
;            and (zqso gt 1.02 or zqso lt 0.95) $
             and 1600.*(zqso+1.) lt 5890.*(1.+zgal) $
;            and (~flag_tmp) $
             and zuniq, nall)

match = match[iall]
nmatch = n_elements(match)

print, 'npairs = ', nmatch
print, '<rp> = ', median(rp_tmp)
rpmean[irad] = median(rp_tmp)

;minwave = 3800D0
;maxwave = 4400D0
;loglam = jhusdss_get_loglam(minwave=minwave, maxwave=maxwave)
;nwave = n_elements(loglam)

;outstr = {wave:10.^loglam, allflux:fltarr(nmatch, nwave), allivar:fltarr(nmatch, nwave)}

;nallwave = n_elements(allspec.wave)

newzgal = gal[match.index_gal].z
newzqso = qso[match.index_qso].z

hist = histogram(newzgal, bin=0.002)
hist = float(hist)/max(hist)
jndex = floor((newzgal-zgalmin)/0.002)
weight = 1./hist[jndex]


;newzgal = randomu(seed, nmatch)*0.1+0.05

print, 'Get wavelength ...'

icaii = value_locate(allspec.wave, 5889.95*(1.+newzgal))
for i=0L, nwave-1L do outwave[i] = median(allspec.wave[icaii-wave_shift+i]/(1.+newzgal))
outstr[irad].wave = outwave

inresidual = fltarr(nmatch, nwave)
inivar = fltarr(nmatch, nwave)+1.
insky = fltarr(nmatch, nwave)
inskyivar = fltarr(nmatch, nwave)+1.
ishuffle = ((floor(randomu(seed,nmatch)*nmatch) > 0 ) < (nmatch-1))

sky_jndex = floor((newzgal[ishuffle]-zgalmin)/0.002)
sky_weight = 1./hist[sky_jndex]

print, 'Get residuals ...'
;for i=0L, nwave-1L do inresidual[*,i] = allspec.residual[match.index_qso, icaii-wave_shift+i]-skyflux[icaii-wave_shift+i]
for i=0L, nwave-1L do inresidual[*,i] = allspec.residual[match.index_qso, icaii-wave_shift+i]
;for i=0L, nwave-1L do insky[*,i] = skyflux[icaii-wave_shift+i]
for i=0L, nwave-1L do insky[*,i] = allspec.residual[match[ishuffle].index_qso, icaii-wave_shift+i]

print, 'Get ivar...'
for i=0L, nwave-1L do inivar[*,i] = allspec.ivar[match.index_qso, icaii-wave_shift+i]*allspec.nmf_continuum[match.index_qso, icaii-wave_shift+i]^2*allspec.med_continuum[match.index_qso, icaii-wave_shift+i]^2
for i=0L, nwave-1L do inskyivar[*,i] = allspec.ivar[match[ishuffle].index_qso, icaii-wave_shift+i]*allspec.nmf_continuum[match[ishuffle].index_qso, icaii-wave_shift+i]^2*allspec.med_continuum[match[ishuffle].index_qso, icaii-wave_shift+i]^2

;iout = where(inresidual gt 1.2 or inresidual lt 0.8, nout)
;inivar[iout] = 0.

;outmedian = median(inresidual, dimension=1)

print, 'Make Composite ...'
jhusdss_composite_engine, inresidual, inivar, fmean=fmean, fmedian=fmedian, $
   fgeomean=fgeomean, nobjuse=nobjuse, ivarweight=ivarweight;, weight=weight
jhusdss_composite_engine, insky, inskyivar, fmean=sky_fmean, fmedian=sky_fmedian, $
   fgeomean=sky_fgeomean, nobjuse=nobjuse, ivarweight=ivarweight;, weight=sky_weight

tmp = moment(fmedian, sdev=sdev)
;djs_plot, outwave, smooth(fmedian,5)-median(fmedian)+1., xra=[3800., 4300], xst=1, yra=[1.-4.*sdev, 1.+4.*sdev]
;djs_oplot, 3934.*[1,1.], !y.crange
;djs_oplot, 3970.*[1,1.], !y.crange
;wait, 5

outstr[irad].fluxmean = fmean
outstr[irad].fluxmedian = fmedian

fmedian=fmean
sky_fmedian=sky_fmean

y = fmedian
iuse = where(outwave gt 5700. and outwave lt 6300., nuse)
terror = sqrt((moment(y[iuse]))[1])

   in_quadra = 0.
   in_slope = 0.
   in_intercept = median(1.-y[iuse])
   in_center = 5889.95
   in_separation = 5895.92-5889.95
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
    xra=[5700., 6300.], xst=1, yra=[1.-4.*terror, 1.+4.*terror]
djs_oplot, 5890.*[1,1.], !y.crange[0]+[0.7,0.9]*(!y.crange[1]-!y.crange[0])
djs_oplot, 5896.*[1,1.], !y.crange[0]+[0.7,0.9]*(!y.crange[1]-!y.crange[0])

;djs_oplot, outwave, 1.-yfit+clevel, psym=10, color='red'
;djs_oplot, outwave, 1.-newy-clevel+clevel, thick=thick, color='green'

djs_oplot, outwave, 1.-sky_fmedian+clevel, color='red'
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

if saveall then begin
   outfile = lowzpath + '/' + 'Lowz_composite_spec_rp_qso_nad.fits'
   mwrfits, outstr, outfile, /create
endif

for i=0,2 do BEEP

end

