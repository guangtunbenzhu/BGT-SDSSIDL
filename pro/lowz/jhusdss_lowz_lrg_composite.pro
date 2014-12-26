;; This should be just a test version
pro jhusdss_lowz_lrg_composite, nmfver, boss=boss, overwrite=overwrite

if (n_elements(nmfver) eq 0) then message, 'nmfver required'
if (n_elements(minrad) eq 0) then minrad=0.0
if (n_elements(maxrad) eq 0) then maxrad=0.4
if (n_elements(nssfr) eq 0) then nssfr=2
if (n_elements(nmass) eq 0) then nmass=2

if (keyword_set(boss)) then begin
   lowzpath=jhusdss_get_path(/fitlrg)+'/'+string(nmfver, format='(I3.3)')+'/Lowz_BOSS'
endif else begin
   lowzpath=jhusdss_get_path(/fitlrg)+'/'+string(nmfver, format='(I3.3)')+'/Lowz'
endelse

outfile = lowzpath+'/'+jhusdss_lowz_lrg_composite_filename(nmfver)
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
ssfr_file = garching_path+'/'+'gal_totspecsfr_dr7_v5_2.fits.gz'
ssfr = mrdfits(ssfr_file, 1)
mass_file = garching_path+'/'+'totlgm_dr7_v5_2.fit.gz'
mass = mrdfits(mass_file, 1)
uniq_file = garching_path+'/'+'gal_uniq_dr7_v5_2.fits'
galuniq = mrdfits(uniq_file, 1)

;; read in SDSS qso catalog, need ra, dec
lrg_path = jhusdss_get_path(/garching)
if (keyword_set(boss)) then begin
   lrg_file = lrg_path+'/'+jhusdss_boss_lrgfile()
endif else begin
   lrg_file = lrg_path+'/'+jhusdss_dr7_lrgfile()
endelse
lrg = mrdfits(lrg_file, 1)

;if (keyword_set(boss)) then begin
;    stat_path=jhusdss_get_path(/nmfqso)+'/'+string(nmfver, format='(I3.3)')+'/Decompose_BOSS'
;;endif else begin
;    stat_path=jhusdss_get_path(/nmfqso)+'/'+string(nmfver, format='(I3.3)')+'/Decompose'
;endelse

;stat_file = stat_path+'/'+jhusdss_stat_filename(nmfver, boss=boss)
;stat = mrdfits(stat_file, 1)

spec = jhusdss_gallrg_match_spec_readin(nmfver, minrad, maxrad, boss=boss, match=match);, flag=flag)
;rp_all = 10.^(alog10(match.rp_mpc*1E3) - 0.19*(mass[match.index_gal].avg-10.0))
rp_all = match.rp_mpc*1E3

;strtmp = {index:0L, zabs:0., plate:0L, fiber:0L, mjd:0L}

;; redshift bins
delta_rp = 25. ; kpc
rp_min = 0.0 ; kpc
rp_max = 100 ; kpc
rp_bin0 = jhusdss_make_bins(rp_min, rp_max, delta_rp, nbin=rp_nbin)

delta_logrp = alog10(2.)
logrp_min = alog10(5.)
logrp_max = alog10(401.)
logrp_bin = jhusdss_make_bins(logrp_min, logrp_max, delta_logrp, nbin=rp_nbin)
rp_bin = replicate(rp_bin0[0], rp_nbin)
rp_bin.min = 10.^logrp_bin.min
rp_bin.max = 10.^logrp_bin.max
rp_min = 10.^logrp_min
rp_max = 10.^logrp_max

rp_bin[0].min=3.00

minwave = 3850D0
maxwave = 4050D0
loglam = jhusdss_get_loglam(minwave=minwave, maxwave=maxwave)
nwave = n_elements(loglam)

str_tmp = {wave:10.^loglam, nabs:lonarr(rp_nbin), zabs:fltarr(rp_nbin), $
   rp_min:rp_bin.min, rp_max:rp_bin.max, rp:fltarr(rp_nbin), $
   ssfr:fltarr(rp_nbin), sfr:fltarr(rp_nbin), mass:fltarr(rp_nbin), $
   fluxmean:fltarr(rp_nbin, nwave), fluxmedian:fltarr(rp_nbin, nwave), $
   fluxgeomean:fltarr(rp_nbin, nwave), nobjuse:fltarr(rp_nbin, nwave)}

composite = replicate(str_tmp, (nssfr+1), (nmass+1))

for i=0L, rp_nbin-1L do begin
    print, i+1, rp_bin[i].min, rp_bin[i].max, rp_nbin
    igal = where(rp_all gt rp_bin[i].min and rp_all le rp_bin[i].max, ngal)

;   sdev_tmp = stat[match[igal].index_lrg].med_sdeviation_red
    ssfr_tmp = ssfr[match[igal].index_gal].avg
    sfr_tmp = sfr[match[igal].index_gal].avg
    mass_tmp = mass[match[igal].index_gal].avg
    zgal = gal[match[igal].index_gal].z
    zuniq = galuniq[match[igal].index_gal].choose
    zlrg = lrg[match[igal].index_lrg].z
    rp_tmp = match[igal].rp_mpc*1E3
;   flag_tmp = flag[igal].flag

;   iall = where(sdev_tmp gt 0. and sdev_tmp le 0.10 and zgal gt 0.015 $
    iall = where(zgal gt 0.035 $
             and zgal lt zlrg-0.035 $
             and zlrg gt 0.08 $
             and zuniq, nall)

    if (nall le 1) then message, 'I have less than 2 pairs'
    print, 'N(Pairs) = ', nall

    all_ssfr = ssfr_tmp[iall]
    all_sfr = sfr_tmp[iall]
    all_mass = mass_tmp[iall]
    all_zgal = zgal[iall]
    all_zlrg = zlrg[iall]
    all_rp = rp_tmp[iall]
    allflux = spec.allflux[igal[iall], *]
    allivar = spec.allivar[igal[iall], *]

    jhusdss_lowz_deadalive, all_mass, all_ssfr, index_dead=index_dead, index_alive=index_alive, $
            index_use=index_use

    for issfr=0L, nssfr do begin

        case issfr of
             0: thisindex = index_use
             1: thisindex = index_dead
             2: thisindex = index_alive
        endcase

        ;; all mass first
        mass_divider_tmp = jhusdss_sample_divider(all_mass[thisindex], nmass, binsize=0.0001)
        mass_divider = [7., mass_divider_tmp, 13.]
        for imass=0L, nmass do begin
            if (imass eq 0L) then begin
               useindex = thisindex
            endif else begin
               useindex_tmp = where(all_mass[thisindex] gt mass_divider[imass-1] $
                                and all_mass[thisindex] le mass_divider[imass], ntmp)
               if (ntmp le 0) then continue
               useindex = thisindex[useindex_tmp]
            endelse

            if (n_elements(useindex) le 1) then begin
               splog, "Less than 2 pairs"
               continue
            endif

            print, 'N(Pairs)=', n_elements(useindex)
            allflux_tmp = allflux[useindex, *]
            allivar_tmp = allivar[useindex, *]

    if (keyword_set(nonsense)) then begin
    if (issfr eq 1) then begin
    a = 'a'
    for j=0L, nall-1L do begin
        djs_plot, 10.^loglam, allflux_tmp[j,*], xra=[3850., 4050.], xst=1, yra=[0.0,2.0]
        djs_oplot, 3934.78*[1.,1.], !y.crange
        djs_oplot, 3969.59*[1.,1.], !y.crange
;       print, qso[match[igal[iall[j]]].index_qso].plate, qso[match[igal[iall[j]]].index_qso].fiber, qso[match[igal[iall[j]]].index_qso].mjd, qso[match[igal[iall[j]]].index_qso].ra, qso[match[igal[iall[j]]].index_qso].dec
;       print, gal[match[igal[iall[j]]].index_gal].plateid, gal[match[igal[iall[j]]].index_gal].fiberid, gal[match[igal[iall[j]]].index_gal].mjd, gal[match[igal[iall[j]]].index_gal].ra, gal[match[igal[iall[j]]].index_gal].dec
        read, a
    endfor
    endif
    endif


            jhusdss_composite_engine, allflux_tmp, allivar_tmp, fmean=fmean, fmedian=fmedian, $
               fgeomean=fgeomean, nobjuse=nobjuse, /ivarweight
            composite[issfr, imass].nabs[i] = n_elements(useindex)
            composite[issfr, imass].zabs[i] = median(all_zgal[useindex])
            composite[issfr, imass].rp[i] = median(all_rp[useindex])
            composite[issfr, imass].ssfr[i] = median(all_ssfr[useindex])
            composite[issfr, imass].sfr[i] = median(all_sfr[useindex])
            composite[issfr, imass].mass[i] = median(all_mass[useindex])
            composite[issfr, imass].fluxmean[i, *] = fmean
            composite[issfr, imass].fluxmedian[i, *] = fmedian
            composite[issfr, imass].fluxgeomean[i, *] = fgeomean
            composite[issfr, imass].nobjuse[i, *] = nobjuse
            djs_plot, 10.^loglam, smooth(fmean,5), yra=[0.99, 1.01]
            djs_oplot, 3934.78*[1.,1.], !y.crange
            djs_oplot, 3969.59*[1.,1.], !y.crange
;           a = 'a'
;           read, a
        endfor
    endfor

;   djs_plot, 10.^loglam, smooth(fmean,5), yra=[0.96, 1.02]
;   djs_oplot, 3934.78*[1.,1.], !y.crange
;   djs_oplot, 3969.59*[1.,1.], !y.crange
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
;   a = 'a'
;   read, a
endfor

mwrfits, composite, outfile, /create
end
