function jhusdss_absorber_cat_qso2abs, stat0, qso0

;; trim
ihaveabs = where(qso0.nabs gt 0, nhaveabs)
qso = qso0[ihaveabs]
stat = stat0[ihaveabs]
nqso = n_elements(qso)

nabs = long(total(qso.nabs))
abstmp = jhusdss_absorber_finalcat_blank()
absorbers = replicate(abstmp, nabs)

lines = strtrim(qso[0].lines, 2)
index_mgi_2853 =where(strcmp(lines, 'mgi_2853', /fold_case) eq 1, nindex)
index_mgii_2803 =where(strcmp(lines, 'mgii_2803', /fold_case) eq 1, nindex)
index_mgii_2796 =where(strcmp(lines, 'mgii_2796', /fold_case) eq 1, nindex)
index_feii_2600 =where(strcmp(lines, 'feii_2600', /fold_case) eq 1, nindex)
index_feii_2586 =where(strcmp(lines, 'feii_2586', /fold_case) eq 1, nindex)
index_feii_2383 =where(strcmp(lines, 'feii_2383', /fold_case) eq 1, nindex)
index_feii_2374 =where(strcmp(lines, 'feii_2374', /fold_case) eq 1, nindex)
index_feii_2344 =where(strcmp(lines, 'feii_2344', /fold_case) eq 1, nindex)

iabs = 0L
for i=0L, nqso-1L do begin
    for j=0L, qso[i].nabs-1L do begin
        absorbers[iabs].ra = stat[i].ra
        absorbers[iabs].dec = stat[i].dec
        absorbers[iabs].plate = stat[i].plate
        absorbers[iabs].fiber = stat[i].fiber
        absorbers[iabs].mjd = stat[i].mjd
        absorbers[iabs].zqso = stat[i].zqso
        absorbers[iabs].err_zqso = stat[i].err_zqso
        absorbers[iabs].index_qso = ihaveabs[i]

        absorbers[iabs].spec_snr_median = stat[i].spec_snr_median
        absorbers[iabs].med_sdeviation_red = stat[i].med_sdeviation_red
        absorbers[iabs].med_sdeviation_blue = stat[i].med_sdeviation_blue

        absorbers[iabs].zabs = qso[i].zabs[j]
        absorbers[iabs].err_zabs = qso[i].err_zabs[j]

        absorbers[iabs].vdisp = qso[i].vdisp[j]
        absorbers[iabs].err_vdisp = qso[i].err_vdisp[j]
        absorbers[iabs].nabs = qso[i].nabs

        ;; might need to be changed
        absorbers[iabs].criterion_mgii = qso[i].criterion_mgii[j]
        absorbers[iabs].criterion_mgii_feii = qso[i].criterion_mgii_feii[j]
        absorbers[iabs].criterion_feii = qso[i].criterion_feii[j]

        absorbers[iabs].signal_mgii_2803 = qso[i].signal[0,j]
        absorbers[iabs].signal_mgii_2796 = qso[i].signal[1,j]
        absorbers[iabs].snr_mgii_2803 = qso[i].snr[0,j]
        absorbers[iabs].snr_mgii_2796 = qso[i].snr[1,j]

        if (index_mgi_2853 ne -1) then begin
           absorbers[iabs].rew_mgi_2853 = qso[i].ew[index_mgi_2853, j]
           absorbers[iabs].err_rew_mgi_2853 = qso[i].err_ew[index_mgi_2853, j]
           absorbers[iabs].vdisp_mgi_2853 = qso[i].vdisp_all[index_mgi_2853, j]
           absorbers[iabs].err_vdisp_mgi_2853 = qso[i].err_vdisp_all[index_mgi_2853, j]
        endif

        if (index_mgii_2803 ne -1) then begin
           absorbers[iabs].rew_mgii_2803 = qso[i].ew[index_mgii_2803, j]
           absorbers[iabs].err_rew_mgii_2803 = qso[i].err_ew[index_mgii_2803, j]
           absorbers[iabs].vdisp_mgii_2803 = qso[i].vdisp_all[index_mgii_2803, j]
           absorbers[iabs].err_vdisp_mgii_2803 = qso[i].err_vdisp_all[index_mgii_2803, j]
        endif

        if (index_mgii_2796 ne -1) then begin
           absorbers[iabs].rew_mgii_2796 = qso[i].ew[index_mgii_2796, j]
           absorbers[iabs].err_rew_mgii_2796 = qso[i].err_ew[index_mgii_2796, j]
           absorbers[iabs].vdisp_mgii_2796 = qso[i].vdisp_all[index_mgii_2796, j]
           absorbers[iabs].err_vdisp_mgii_2796 = qso[i].err_vdisp_all[index_mgii_2796, j]
        endif

        if (index_feii_2600 ne -1) then begin
           absorbers[iabs].rew_feii_2600 = qso[i].ew[index_feii_2600, j]
           absorbers[iabs].err_rew_feii_2600 = qso[i].err_ew[index_feii_2600, j]
           absorbers[iabs].vdisp_feii_2600 = qso[i].vdisp_all[index_feii_2600, j]
           absorbers[iabs].err_vdisp_feii_2600 = qso[i].err_vdisp_all[index_feii_2600, j]
        endif

        if (index_feii_2586 ne -1) then begin
           absorbers[iabs].rew_feii_2586 = qso[i].ew[index_feii_2586, j]
           absorbers[iabs].err_rew_feii_2586 = qso[i].err_ew[index_feii_2586, j]
           absorbers[iabs].vdisp_feii_2586 = qso[i].vdisp_all[index_feii_2586, j]
           absorbers[iabs].err_vdisp_feii_2586 = qso[i].err_vdisp_all[index_feii_2586, j]
        endif

        if (index_feii_2383 ne -1) then begin
           absorbers[iabs].rew_feii_2383 = qso[i].ew[index_feii_2383, j]
           absorbers[iabs].err_rew_feii_2383 = qso[i].err_ew[index_feii_2383, j]
           absorbers[iabs].vdisp_feii_2383 = qso[i].vdisp_all[index_feii_2383, j]
           absorbers[iabs].err_vdisp_feii_2383 = qso[i].err_vdisp_all[index_feii_2383, j]
        endif

        if (index_feii_2374 ne -1) then begin
           absorbers[iabs].rew_feii_2374 = qso[i].ew[index_feii_2374, j]
           absorbers[iabs].err_rew_feii_2374 = qso[i].err_ew[index_feii_2374, j]
           absorbers[iabs].vdisp_feii_2374 = qso[i].vdisp_all[index_feii_2374, j]
           absorbers[iabs].err_vdisp_feii_2374 = qso[i].err_vdisp_all[index_feii_2374, j]
        endif

        if (index_feii_2344 ne -1) then begin
           absorbers[iabs].rew_feii_2344 = qso[i].ew[index_feii_2344, j]
           absorbers[iabs].err_rew_feii_2344 = qso[i].err_ew[index_feii_2344, j]
           absorbers[iabs].vdisp_feii_2344 = qso[i].vdisp_all[index_feii_2344, j]
           absorbers[iabs].err_vdisp_feii_2344 = qso[i].err_vdisp_all[index_feii_2344, j]
        endif

        iabs++
    endfor
endfor

return, absorbers
end
