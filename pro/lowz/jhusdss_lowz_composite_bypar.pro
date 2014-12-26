;; This should be just a test version
pro jhusdss_lowz_composite_bypar, nmfver, boss=boss, overwrite=overwrite

if (n_elements(nmfver) eq 0) then message, 'nmfver required'

if (n_elements(minrad) eq 0) then minrad=0.0
if (n_elements(maxrad) eq 0) then maxrad=0.1

if (keyword_set(boss)) then begin
   lowzpath=jhusdss_get_path(/nmfqso)+'/'+string(nmfver, format='(I3.3)')+'/Lowz_BOSS'
endif else begin
   lowzpath=jhusdss_get_path(/nmfqso)+'/'+string(nmfver, format='(I3.3)')+'/Lowz'
endelse

outfile = lowzpath+'/'+jhusdss_lowz_composite_filename(nmfver)
if (file_test(outfile) and ~keyword_set(overwrite)) then begin
   splog, 'File already exists. Use /overwrite if you want to overwrite it.'
   return
endif else begin
   splog, 'Will write the composite into this file: '
   print, outfile
endelse

garching_path = jhusdss_get_path(/garching)
garching_file = garching_path+'/'+'gal_info_dr7_v5_2.fit.gz'
gal = mrdfits(garching_file, 1)
sfr_file = garching_path+'/'+'gal_totsfr_dr7_v5_2.fits.gz'
sfr = mrdfits(sfr_file, 1)

;; read in SDSS qso catalog, need ra, dec
qso_path = jhusdss_get_path(/qso)
if (keyword_set(boss)) then begin
   qso_file = qso_path+'/'+jhusdss_boss_qsofile()
endif else begin
   qso_file = qso_path+'/'+jhusdss_dr7_qsofile()
endelse
qso = mrdfits(qso_file, 1)

if (keyword_set(boss)) then begin
    stat_path=jhusdss_get_path(/nmfqso)+'/'+string(nmfver, format='(I3.3)')+'/Decompose_BOSS'
endif else begin
    stat_path=jhusdss_get_path(/nmfqso)+'/'+string(nmfver, format='(I3.3)')+'/Decompose'
endelse

stat_file = stat_path+'/'+jhusdss_stat_filename(nmfver, boss=boss)
stat = mrdfits(stat_file, 1)

spec = jhusdss_galqso_match_spec_readin(nmfver, minrad, maxrad, boss=boss, match=match)

;strtmp = {index:0L, zabs:0., plate:0L, fiber:0L, mjd:0L}

;; redshift bins
delta_rp = 25. ; kpc
rp_min = 0.0 ; kpc
rp_max = 100 ; kpc
rp_bin0 = jhusdss_make_bins(rp_min, rp_max, delta_rp, nbin=rp_nbin)

delta_logrp = alog10(2.)
logrp_min = alog10(3.)
logrp_max = alog10(100.)
logrp_bin = jhusdss_make_bins(logrp_min, logrp_max, delta_logrp, nbin=rp_nbin)
rp_bin = replicate(rp_bin0[0], rp_nbin)
rp_bin.min = 10.^logrp_bin.min
rp_bin.max = 10.^logrp_bin.max
rp_min = 10.^logrp_min
rp_max = 10.^logrp_max

rp_bin[0].min=3.00

;rp_bin[0].min = 1.0

minwave = 3850D0
maxwave = 4050D0
loglam = jhusdss_get_loglam(minwave=minwave, maxwave=maxwave)
nwave = n_elements(loglam)

for i=0L, rp_nbin-1L do begin
    print, i+1, rp_bin[i].min, rp_bin[i].max, rp_nbin
    igal = where(match.rp_mpc*1E3 gt rp_bin[i].min and match.rp_mpc*1E3 le rp_bin[i].max, ngal)

;   zstr = replicate(strtmp, ngal)
;   zstr.index = match[igal].index_qso
;   zstr.plate = qso[match[igal].index_qso].plate
;   zstr.fiber = qso[match[igal].index_qso].fiber
;   zstr.mjd = qso[match[igal].index_qso].mjd
    sdev_tmp = stat[match[igal].index_qso].med_sdeviation_red
    sfr_tmp = sfr[match[igal].index_gal].avg
    zgal = gal[match[igal].index_gal].z
    zqso = qso[match[igal].index_qso].z

    iall = where(sdev_tmp gt 0. and sdev_tmp le 0.07 and zgal gt 0.025 $
             and zgal lt zqso-0.025 $
             and 1255.*(zqso+1.) lt 2800.*(1+zgal), nall)

    if (nall le 1) then message, 'I have less than 2 pairs'
    print, 'N(Pairs) = ', nall
    allflux = spec.allflux[igal[iall], *]
    allivar = spec.allivar[igal[iall], *]

;   if (keyword_set(nonsense)) then begin
    a = 'a'
    for j=0L, nall-1L do begin
        djs_plot, 10.^loglam, allflux[j,*], xra=[3850., 4050.], xst=1, yra=[0.0,2.0]
        djs_oplot, 3934.78*[1.,1.], !y.crange
        djs_oplot, 3969.59*[1.,1.], !y.crange
        print, qso[match[igal[iall[j]]].index_qso].plate, qso[match[igal[iall[j]]].index_qso].fiber, qso[match[igal[iall[j]]].index_qso].mjd, qso[match[igal[iall[j]]].index_qso].ra, qso[match[igal[iall[j]]].index_qso].dec
        print, gal[match[igal[iall[j]]].index_gal].plateid, gal[match[igal[iall[j]]].index_gal].fiberid, gal[match[igal[iall[j]]].index_gal].mjd, gal[match[igal[iall[j]]].index_gal].ra, gal[match[igal[iall[j]]].index_gal].dec
        read, a
    endfor
;   endif

    jhusdss_composite_engine, allflux, allivar, fmean=fmean, fmedian=fmedian, fgeomean=fgeomean, /ivarweight

    djs_plot, 10.^loglam, smooth(fmean,5), yra=[0.96, 1.02]
    djs_oplot, 3934.78*[1.,1.], !y.crange
    djs_oplot, 3969.59*[1.,1.], !y.crange
;   if (n_elements(composite) eq 0) then begin
;      composite = jhusdss_lowz_composite_blank(nwave, rp_nbin)
;      composite.wave = composite_tmp.wave
;   endif
;   composite.nabs[i] = composite_tmp.nabs
;   composite.zabs[i] = composite_tmp.zabs
;   composite.rp_min[i] = rp_bin[i].min
;   composite.rp_max[i] = rp_bin[i].max
;   composite.rp[i] = median(match[igal].rp_mpc*1E3)
;   composite.fluxmean[*, i] = composite_tmp.fluxmean
;   composite.fluxmedian[*, i] = composite_tmp.fluxmedian
;   composite.fluxgeomean[*, i] = composite_tmp.fluxgeomean
;   composite.nobjuse[*, i] = composite_tmp.nobjuse
    stop
endfor

;mwrfits, composite, outfile, /create
end
