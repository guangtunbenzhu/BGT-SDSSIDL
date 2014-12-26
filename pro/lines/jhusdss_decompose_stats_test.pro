;+
;-
;; To create a qso parameter file. We will trim those with 'bad' statistics: 
;; weird eigenvalues, large skewness, large width, low s/n 
pro jhusdss_decompose_stats_test, nmfver, overwrite=overwrite, boss=boss

if (n_elements(nmfver) eq 0) then nmfver = jhusdss_get_nmf_version()

;; output
path=jhusdss_get_path(/nmfqso)+'/'+string(nmfver, format='(I3.3)')+'/Decompose'
if (jhusdss_direxist(path) eq 0) then message, "Can't find the directory."
filename = path+'/QSO_decompose_NMF_'+string(nmfver, format='(I3.3)')+'stats.fits'

absorber0 = jhusdss_absorber_readin(nmfver)
stats0 = mrdfits(filename, 1)
ii = where(stats0.isitstated eq 1b, nn)
stats0 = stats0[ii]
absorber0 = absorber0[ii]

stop
ibad = where(stats0.isitstated eq 1b and stats0.med_sdeviation_red gt 0.20 and absorber0.nabs gt 0, nbad); and stats0.sdeviation_red lt 4.0, nbad)
;ibad = where(stats0.isitstated eq 1b and stats0.sdeviation_red gt 4.0, nbad); and stats0.sdeviation_red lt 4.0, nbad)
;ibad = where(stats0.isitstated eq 1b and stats0.sdeviation_red gt 3.0 and $
;             abs(stats0.continuous_diff_red2) gt 0.01, nbad)
;            abs(stats0.continuous_diff_red1)+abs(stats0.continuous_diff_red2)+abs(stats0.continuous_diff_red3) gt 0.04, nbad)
;ibad = where(stats0.isitstated eq 1b and stats0.skewness_red gt 2., nbad)

stats = stats0[ibad]
absorber = absorber0[ibad]
for i=0L, nbad-1L do begin

;  counter, i+1, nbad
   print, i+1, nbad

   ;; load convolved spectra
   spec = jhusdss_decompose_loadspec(stats[i].plate, stats[i].fiber, nmfver, error=error)

   ii = where(spec.ivar gt 0., nn)
;  window, 0
   djs_plot, spec.wave[ii]*(1.+spec.z), spec.flux[ii], xra=[3700., 9200], pos=[0.1, 0.50, 0.9, 0.95], $
             yra=[0,6], xst=1, yst=1
   djs_oplot, spec.wave[ii]*(1.+spec.z), spec.nmf_continuum[ii], color='red'
   djs_oplot, [2800., 2800.]*(1+absorber[i].zabs[0]), !y.crange, color='green'
   djs_plot, spec.wave[ii], spec.med_continuum[ii], xra=[3700., 9200]/(1.+spec.z), pos=[0.1, 0.05, 0.9, 0.45], $
             yra=[0.3, 1.8], xst=1, yst=1, /noerase
;  print, stats[i].sdeviation_red, stats[i].skewness_red
   print, stats[i].med_sdeviation_red, absorber[i].nabs
   print, absorber[i].zabs[0:absorber[i].nabs-1]
   print, 2800.*(1.+absorber[i].zabs[0:absorber[i].nabs-1])
;  print, stats[i].continuous_diff_red1, stats[i].continuous_diff_red2, stats[i].continuous_diff_red3
;  window, 1

   a = 'a'
   read, a
   if (a eq 's') then stop
endfor

end
