;; Version 1: for Nino/X
;; Version 2: Final
;; for all, jhusdss_abs2qso, 107, /keepasso, /addall

pro jhusdss_abs2qso, nmfver, boss=boss, dr12=dr12, feii=feii, overwrite=overwrite, wlim=wlim, snrlimit=snrlimit, $
                     keepasso=keepasso, keeplowzlowsnr=keeplowzlowsnr, addall=addall

    if (n_elements(nmfver) eq 0) then message, 'NMFVER required'
    if (n_elements(wlim) eq 0) then wlim = 0.
    if (n_elements(snrlimit) eq 0) then snrlimit = 0.

    if (keyword_set(dr12)) then begin
       path=jhusdss_get_path(/nmfqso)+'/'+string(nmfver, format='(I3.3)')+'/Absorbers_DR12'
    endif else begin
       if (keyword_set(boss)) then begin
          path=jhusdss_get_path(/nmfqso)+'/'+string(nmfver, format='(I3.3)')+'/Absorbers_BOSS'
       endif else begin
          path=jhusdss_get_path(/nmfqso)+'/'+string(nmfver, format='(I3.3)')+'/Absorbers'
       endelse
    endelse

    outfile = path+'/ALLQSO_Trimmed_QSO_'+jhusdss_absorbers_filename(nmfver, boss=boss, dr12=dr12, /mgii)
    if (keyword_set(addall)) then  $
       outfile = repstr(outfile, '.fits', '_ALL.fits')

    if (file_test(outfile) and ~keyword_set(overwrite)) then begin
       splog, 'File already exists. Use /overwrite if you want to overwrite it.'
       return
    endif else begin
       splog, 'Will write the absorber catalog into this file: '
       print, outfile
    endelse

    qsos = jhusdss_qso_readin(boss=boss, dr12=dr12)
    stats = jhusdss_qsostats_readin(nmfver, boss=boss, dr12=dr12)

    absorbers = jhusdss_absorber_readin(nmfver, boss=boss, dr12=dr12, keepasso=keepasso, keeplowzlowsnr=keeplowzlowsnr, addall=addall, feii=feii)

    strtmp = jhusdss_abs2qso_blank()

    outstr = replicate(strtmp, n_elements(qsos))

    outstr.ra = stats.ra
    outstr.dec = stats.dec
    outstr.plate  = stats.plate
    outstr.fiber = stats.fiber
    outstr.mjd = stats.mjd
    outstr.zqso = stats.zqso
    outstr.err_zqso = stats.err_zqso
    outstr.index_qso = lindgen(n_elements(qsos))
    outstr.spec_snr_median = stats.spec_snr_median
    outstr.med_sdeviation_red = stats.med_sdeviation_red
    outstr.med_sdeviation_blue = stats.med_sdeviation_blue

    for i=0L, n_elements(absorbers)-1L do begin
        if (absorbers[i].rew_mgii_2796 lt wlim) then continue
;       if ((absorbers[i].zabs gt 0.46) and absorbers[i].criterion_mgii_feii eq 0b) then continue

        index = absorbers[i].index_qso
        jabs = outstr[index].nabs
        outstr[index].nabs = outstr[index].nabs+1L

        outstr[index].criterion_mgii[jabs] = absorbers[i].criterion_mgii
        outstr[index].criterion_mgii_feii[jabs] = absorbers[i].criterion_mgii_feii
        outstr[index].criterion_feii[jabs] = absorbers[i].criterion_feii

        outstr[index].signal_mgii_2803[jabs] = absorbers[i].signal_mgii_2803
        outstr[index].signal_mgii_2796[jabs] = absorbers[i].signal_mgii_2796
        outstr[index].snr_mgii_2803[jabs] = absorbers[i].snr_mgii_2803
        outstr[index].snr_mgii_2796[jabs] = absorbers[i].snr_mgii_2796

        outstr[index].zabs[jabs] = absorbers[i].zabs
        outstr[index].err_zabs[jabs] = absorbers[i].err_zabs
        outstr[index].vdisp[jabs] = absorbers[i].vdisp
        outstr[index].err_vdisp[jabs] = absorbers[i].err_vdisp

        outstr[index].rew_mgi_2853[jabs] = absorbers[i].rew_mgi_2853
        outstr[index].err_rew_mgi_2853[jabs] = absorbers[i].err_rew_mgi_2853
        outstr[index].vdisp_mgi_2853[jabs] = absorbers[i].vdisp_mgi_2853
        outstr[index].err_vdisp_mgi_2853[jabs] = absorbers[i].err_vdisp_mgi_2853

        outstr[index].rew_mgii_2803[jabs] = absorbers[i].rew_mgii_2803
        outstr[index].err_rew_mgii_2803[jabs] = absorbers[i].err_rew_mgii_2803
        outstr[index].vdisp_mgii_2803[jabs] = absorbers[i].vdisp_mgii_2803
        outstr[index].err_vdisp_mgii_2803[jabs] = absorbers[i].err_vdisp_mgii_2803

        outstr[index].rew_mgii_2796[jabs] = absorbers[i].rew_mgii_2796
        outstr[index].err_rew_mgii_2796[jabs] = absorbers[i].err_rew_mgii_2796
        outstr[index].vdisp_mgii_2796[jabs] = absorbers[i].vdisp_mgii_2796
        outstr[index].err_vdisp_mgii_2796[jabs] = absorbers[i].err_vdisp_mgii_2796

        outstr[index].rew_feii_2600[jabs] = absorbers[i].rew_feii_2600
        outstr[index].err_rew_feii_2600[jabs] = absorbers[i].err_rew_feii_2600
        outstr[index].vdisp_feii_2600[jabs] = absorbers[i].vdisp_feii_2600
        outstr[index].err_vdisp_feii_2600[jabs] = absorbers[i].err_vdisp_feii_2600

        outstr[index].rew_feii_2586[jabs] = absorbers[i].rew_feii_2586
        outstr[index].err_rew_feii_2586[jabs] = absorbers[i].err_rew_feii_2586
        outstr[index].vdisp_feii_2586[jabs] = absorbers[i].vdisp_feii_2586
        outstr[index].err_vdisp_feii_2586[jabs] = absorbers[i].err_vdisp_feii_2586

        outstr[index].rew_feii_2383[jabs] = absorbers[i].rew_feii_2383
        outstr[index].err_rew_feii_2383[jabs] = absorbers[i].err_rew_feii_2383
        outstr[index].vdisp_feii_2383[jabs] = absorbers[i].vdisp_feii_2383
        outstr[index].err_vdisp_feii_2383[jabs] = absorbers[i].err_vdisp_feii_2383

        outstr[index].rew_feii_2374[jabs] = absorbers[i].rew_feii_2374
        outstr[index].err_rew_feii_2374[jabs] = absorbers[i].err_rew_feii_2374
        outstr[index].vdisp_feii_2374[jabs] = absorbers[i].vdisp_feii_2374
        outstr[index].err_vdisp_feii_2374[jabs] = absorbers[i].err_vdisp_feii_2374

        outstr[index].rew_feii_2344[jabs] = absorbers[i].rew_feii_2344
        outstr[index].err_rew_feii_2344[jabs] = absorbers[i].err_rew_feii_2344
        outstr[index].vdisp_feii_2344[jabs] = absorbers[i].vdisp_feii_2344
        outstr[index].err_vdisp_feii_2344[jabs] = absorbers[i].err_vdisp_feii_2344

    endfor
    
    if (keyword_set(feii)) then begin
    for i=0L, n_elements(feii_absorbers)-1L do begin
        if (feii_absorbers[i].rew_mgii_2796 lt wlim) then continue

        index = feii_absorbers[i].index_qso
        jabs = outstr[index].nabs
        outstr[index].nabs = outstr[index].nabs+1L

        outstr[index].criterion_mgii[jabs] = feii_absorbers[i].criterion_mgii
        outstr[index].criterion_mgii_feii[jabs] = feii_absorbers[i].criterion_mgii_feii
        outstr[index].criterion_feii[jabs] = feii_absorbers[i].criterion_feii

        outstr[index].signal_mgii_2803[jabs] = feii_absorbers[i].signal_mgii_2803
        outstr[index].signal_mgii_2796[jabs] = feii_absorbers[i].signal_mgii_2796
        outstr[index].snr_mgii_2803[jabs] = feii_absorbers[i].snr_mgii_2803
        outstr[index].snr_mgii_2796[jabs] = feii_absorbers[i].snr_mgii_2796

        outstr[index].zabs[jabs] = feii_absorbers[i].zabs
        outstr[index].err_zabs[jabs] = feii_absorbers[i].err_zabs
        outstr[index].vdisp[jabs] = feii_absorbers[i].vdisp
        outstr[index].err_vdisp[jabs] = feii_absorbers[i].err_vdisp

        outstr[index].rew_mgi_2853[jabs] = feii_absorbers[i].rew_mgi_2853
        outstr[index].err_rew_mgi_2853[jabs] = feii_absorbers[i].err_rew_mgi_2853
        outstr[index].vdisp_mgi_2853[jabs] = feii_absorbers[i].vdisp_mgi_2853
        outstr[index].err_vdisp_mgi_2853[jabs] = feii_absorbers[i].err_vdisp_mgi_2853

        outstr[index].rew_mgii_2803[jabs] = feii_absorbers[i].rew_mgii_2803
        outstr[index].err_rew_mgii_2803[jabs] = feii_absorbers[i].err_rew_mgii_2803
        outstr[index].vdisp_mgii_2803[jabs] = feii_absorbers[i].vdisp_mgii_2803
        outstr[index].err_vdisp_mgii_2803[jabs] = feii_absorbers[i].err_vdisp_mgii_2803

        outstr[index].rew_mgii_2796[jabs] = feii_absorbers[i].rew_mgii_2796
        outstr[index].err_rew_mgii_2796[jabs] = feii_absorbers[i].err_rew_mgii_2796
        outstr[index].vdisp_mgii_2796[jabs] = feii_absorbers[i].vdisp_mgii_2796
        outstr[index].err_vdisp_mgii_2796[jabs] = feii_absorbers[i].err_vdisp_mgii_2796

        outstr[index].rew_feii_2600[jabs] = feii_absorbers[i].rew_feii_2600
        outstr[index].err_rew_feii_2600[jabs] = feii_absorbers[i].err_rew_feii_2600
        outstr[index].vdisp_feii_2600[jabs] = feii_absorbers[i].vdisp_feii_2600
        outstr[index].err_vdisp_feii_2600[jabs] = feii_absorbers[i].err_vdisp_feii_2600

        outstr[index].rew_feii_2586[jabs] = feii_absorbers[i].rew_feii_2586
        outstr[index].err_rew_feii_2586[jabs] = feii_absorbers[i].err_rew_feii_2586
        outstr[index].vdisp_feii_2586[jabs] = feii_absorbers[i].vdisp_feii_2586
        outstr[index].err_vdisp_feii_2586[jabs] = feii_absorbers[i].err_vdisp_feii_2586

        outstr[index].rew_feii_2383[jabs] = feii_absorbers[i].rew_feii_2383
        outstr[index].err_rew_feii_2383[jabs] = feii_absorbers[i].err_rew_feii_2383
        outstr[index].vdisp_feii_2383[jabs] = feii_absorbers[i].vdisp_feii_2383
        outstr[index].err_vdisp_feii_2383[jabs] = feii_absorbers[i].err_vdisp_feii_2383

        outstr[index].rew_feii_2374[jabs] = feii_absorbers[i].rew_feii_2374
        outstr[index].err_rew_feii_2374[jabs] = feii_absorbers[i].err_rew_feii_2374
        outstr[index].vdisp_feii_2374[jabs] = feii_absorbers[i].vdisp_feii_2374
        outstr[index].err_vdisp_feii_2374[jabs] = feii_absorbers[i].err_vdisp_feii_2374

        outstr[index].rew_feii_2344[jabs] = feii_absorbers[i].rew_feii_2344
        outstr[index].err_rew_feii_2344[jabs] = feii_absorbers[i].err_rew_feii_2344
        outstr[index].vdisp_feii_2344[jabs] = feii_absorbers[i].vdisp_feii_2344
        outstr[index].err_vdisp_feii_2344[jabs] = feii_absorbers[i].err_vdisp_feii_2344

    endfor
    endif
 
;   ii = where(outstr.nabs gt 0 and outstr.spec_snr_median gt snrlimit, nn)
   
;  sdev_limit = 0.07
   zqso_min = 0.40
   zqso_max = 4.7

;  ii = where(outstr.spec_snr_median gt snrlimit $
;            and outstr.med_sdeviation_red gt 0.00 $
;            and outstr.med_sdeviation_red le sdev_limit $
;            and outstr.zqso ge zqso_min $
;            and outstr.zqso lt zqso_max)
   ii = where(outstr.nabs gt 0)

    mwrfits, outstr[ii], outfile, /create

end
