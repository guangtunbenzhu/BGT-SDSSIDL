;; see jhusdss_absorber_cat_qso2abs.pro, this is a trimmed version
;; updated with stats info
pro jhusdss_pitts_qso2abs, nmfver

qsofile = jhusdss_get_path(/absorber)+'/MgII/Master_Pitts_Catalog.fits'
qso0 = mrdfits(qsofile, 1)

qsopath=jhusdss_get_path(/nmfqso)+'/'+string(nmfver, format='(I3.3)')+'/Decompose'
statfile = qsopath+'/Pitts_QSO_decompose_NMF_'+string(nmfver, format='(I3.3)')+'stats.fits'
stat0 = mrdfits(statfile, 1)

outpath= jhusdss_get_path(/nmfqso)+'/'+string(nmfver, format='(I3.3)')+'/Absorbers'
outfile = outpath+'/'+'Master_Pitts_Catalog_Absorbers_'+string(nmfver, format='(I3.3)')+'.fits'

ihaveabs = where(qso0.nabs gt 0, nhaveabs)
qso = qso0[ihaveabs]
stat = stat0[ihaveabs]
nqso = n_elements(qso)

nabs = fix(total(qso.nabs))
abstmp = jhusdss_absorber_finalcat_blank()
absorbers = replicate(abstmp, nabs)

lines = strtrim(qso[0].lines, 2)
index_mgi_2853 =where(strcmp(lines, 'mgi_2853', /fold_case) eq 1, nindex)
index_mgii_2803 =where(strcmp(lines, 'mgii_2803', /fold_case) eq 1, nindex)
index_mgii_2796 =where(strcmp(lines, 'mgii_2796', /fold_case) eq 1, nindex)
index_feii_2600 =where(strcmp(lines, 'feii_2600', /fold_case) eq 1, nindex)
index_feii_2586 =where(strcmp(lines, 'feii_2586', /fold_case) eq 1, nindex)

iabs = 0L
for i=0L, nqso-1L do begin
    for j=0L, qso[i].nabs-1L do begin
        absorbers[iabs].ra = qso[i].ra
        absorbers[iabs].dec = qso[i].dec
        absorbers[iabs].plate = qso[i].plate
        absorbers[iabs].fiber = qso[i].fiber
        absorbers[iabs].mjd = qso[i].mjd
        absorbers[iabs].zqso = stat[i].zqso
        absorbers[iabs].err_zqso = qso[i].err_zqso

        absorbers[iabs].spec_snr_median = stat[i].spec_snr_median
        absorbers[iabs].med_sdeviation_red = stat[i].med_sdeviation_red
        absorbers[iabs].med_sdeviation_blue = stat[i].med_sdeviation_blue

        absorbers[iabs].nabs = qso[i].nabs

        absorbers[iabs].zabs = qso[i].zabs[j]
;       absorbers[iabs].err_zabs = qso[i].err_zabs[j]

        if (index_mgi_2853 ne -1) then begin
           absorbers[iabs].rew_mgi_2853 = qso[i].ew[index_mgi_2853, j]
           absorbers[iabs].err_rew_mgi_2853 = qso[i].err_ew[index_mgi_2853, j]
        endif

        if (index_mgii_2803 ne -1) then begin
           absorbers[iabs].rew_mgii_2803 = qso[i].ew[index_mgii_2803, j]
           absorbers[iabs].err_rew_mgii_2803 = qso[i].err_ew[index_mgii_2803, j]
        endif

        if (index_mgii_2796 ne -1) then begin
           absorbers[iabs].rew_mgii_2796 = qso[i].ew[index_mgii_2796, j]
           absorbers[iabs].err_rew_mgii_2796 = qso[i].err_ew[index_mgii_2796, j]
        endif

        if (index_feii_2600 ne -1) then begin
           absorbers[iabs].rew_feii_2600 = qso[i].ew[index_feii_2600, j]
           absorbers[iabs].err_rew_feii_2600 = qso[i].err_ew[index_feii_2600, j]
        endif

        if (index_feii_2586 ne -1) then begin
           absorbers[iabs].rew_feii_2586 = qso[i].ew[index_feii_2586, j]
           absorbers[iabs].err_rew_feii_2586 = qso[i].err_ew[index_feii_2586, j]
        endif
        iabs++
    endfor
endfor

mwrfits, absorbers, outfile, /create
end
