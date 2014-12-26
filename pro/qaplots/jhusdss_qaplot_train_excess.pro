pro jhusdss_qaplot_train_excess, nmfver

path=jhusdss_get_path(/nmfqso)+'/'+string(nmfver, format='(I3.3)')
mfile = path+'/Absorbers/'+'JHU_Pitts_matched_'+string(nmfver, format='(I3.3)')+'.fits'
nomfile = path+'/Absorbers/'+'Pitts_nomatched_'+string(nmfver, format='(I3.3)')+'.fits'
jhu_nomfile = path+'/Absorbers/'+'JHU_nomatched_'+string(nmfver, format='(I3.3)')+'.fits'

outfile = path+'/Absorbers/'+'JHU_nomatched_composite_'+string(nmfver, format='(I3.3)')+'.fits'

pitts0 = mrdfits(mfile, 1)
jhu0 = mrdfits(mfile, 2)
pitts_nomatch0 = mrdfits(nomfile, 1)
jhu_nomatch0 = mrdfits(jhu_nomfile, 1)

;;##########
;; trim
;; This is not accurate -- Guangtun 04-15-2013
;; Let's use ijhu_all instead
inwindow_all = bytarr(n_elements(pitts0))
pitts0.criterion_mgii = 1b
pitts0.criterion_mgii_feii = 1b
pitts0.criterion_feii = 1b
pitts0.snr_mgii_2796 = jhu0.snr_mgii_2796
ipitts_all = jhusdss_absorber_trim(pitts0)
inwindow_all[ipitts_all] = 1b

;; real matches
ireal_matches = bytarr(n_elements(pitts0))
ijhu_all = jhusdss_absorber_trim(jhu0)
pitts_all = pitts0[ijhu_all]
jhu_all = jhu0[ijhu_all]
ireal_matches[ijhu_all] = 1b

;; in window but not real matches
;; Not accurate -- Guangtun 04-15-2013
;; ipitts_all_nomatch = where(inwindow_all and ~ireal_matches)

;; only Mg II
;; pitts_all includes Fe II, pitts_mgii only Mg II
;; Add comp=ipitts_all_nomatch -- Guangtun 04-15-2013
ipitts_mgii = where(jhu_all.criterion_mgii eq 1b, comp=ipitts_all_nomatch)
pitts_mgii = pitts_all[ipitts_mgii]
jhu_mgii = jhu_all[ipitts_mgii]

;; no match in window
pitts_nomatch0.criterion_mgii = 1b
pitts_nomatch0.criterion_mgii_feii = 1b
pitts_nomatch0.criterion_feii = 1b
pitts_nomatch0.snr_mgii_2796 = 0. ;; Low-z shouldn't be a reference, Guangtun 04-15-2013
ipitts_nomatch = jhusdss_absorber_trim(pitts_nomatch0)
;; Change to pitts_all[ipitts_all_nomatch] -- Guangtun 04-15-2013
;; pitts_nomatch = [pitts_nomatch0[ipitts_nomatch], pitts0[ipitts_all_nomatch]]
pitts_nomatch = [pitts_nomatch0[ipitts_nomatch], pitts_all[ipitts_all_nomatch]]

;; no match in window
ijhu_nomatch = jhusdss_absorber_trim(jhu_nomatch0)
jhu_nomatch = jhu_nomatch0[ijhu_nomatch]

;; only Mg II
ijhu_nomatch_mgii = where(jhu_nomatch.criterion_mgii eq 1b)
jhu_nomatch_mgii = jhu_nomatch[ijhu_nomatch_mgii]
;;##########

jhu_nomatch = jhu_nomatch_mgii

ihigh = where(jhu_nomatch.rew_mgii_2796 ge 1.5 and jhu_nomatch.rew_mgii_2796 lt 4.0, nhigh)
imedium = where(jhu_nomatch.rew_mgii_2796 ge 1.1 and jhu_nomatch.rew_mgii_2796 lt 1.5, nmedium)
ilow = where(jhu_nomatch.rew_mgii_2796 ge 0.8 and jhu_nomatch.rew_mgii_2796 lt 1.1, nlow)
iverylow = where(jhu_nomatch.rew_mgii_2796 ge 0.2 and jhu_nomatch.rew_mgii_2796 lt 0.8, nverylow)

i1 = where(jhu_nomatch.rew_mgii_2796 ge 0.65 and jhu_nomatch.rew_mgii_2796 lt 0.8, n1)
i2 = where(jhu_nomatch.rew_mgii_2796 ge 0.5 and jhu_nomatch.rew_mgii_2796 lt 0.65, n2)
i3 = where(jhu_nomatch.rew_mgii_2796 ge 0.4 and jhu_nomatch.rew_mgii_2796 lt 0.5, n3)
i4 = where(jhu_nomatch.rew_mgii_2796 ge 0.2 and jhu_nomatch.rew_mgii_2796 lt 0.3, n4)

print, median(jhu_nomatch[ihigh].rew_mgii_2796), median(jhu_nomatch[imedium].rew_mgii_2796), median(jhu_nomatch[ilow].rew_mgii_2796), median(jhu_nomatch[iverylow].rew_mgii_2796), median(jhu_nomatch[i1].rew_mgii_2796), median(jhu_nomatch[i2].rew_mgii_2796), median(jhu_nomatch[i3].rew_mgii_2796), median(jhu_nomatch[i4].rew_mgii_2796)
stop

comp_high = jhusdss_absorbers_composite_engine(jhu_nomatch[ihigh], nmfver=nmfver)
comp_medium = jhusdss_absorbers_composite_engine(jhu_nomatch[imedium], nmfver=nmfver)
comp_low = jhusdss_absorbers_composite_engine(jhu_nomatch[ilow], nmfver=nmfver)
comp_verylow = jhusdss_absorbers_composite_engine(jhu_nomatch[iverylow], nmfver=nmfver)
comp_1 = jhusdss_absorbers_composite_engine(jhu_nomatch[i1], nmfver=nmfver)
comp_2 = jhusdss_absorbers_composite_engine(jhu_nomatch[i2], nmfver=nmfver)
comp_3 = jhusdss_absorbers_composite_engine(jhu_nomatch[i3], nmfver=nmfver)
comp_4 = jhusdss_absorbers_composite_engine(jhu_nomatch[i4], nmfver=nmfver)

mwrfits, comp_high, outfile, /create
mwrfits, comp_medium, outfile 
mwrfits, comp_low, outfile
mwrfits, comp_verylow, outfile
mwrfits, comp_1, outfile
mwrfits, comp_2, outfile
mwrfits, comp_3, outfile
mwrfits, comp_4, outfile

end
