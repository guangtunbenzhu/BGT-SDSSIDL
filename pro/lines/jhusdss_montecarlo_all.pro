pro jhusdss_montecarlo_all, nmfver, boss=boss, overwrite=overwrite, noresidual=noresidual

if (n_elements(nmfver) eq 0) then nmfver = jhusdss_get_nmf_version()

;; how many qso we are doing
;; do all of them, it doesn't take much time anyway
;; Yue Shen's catalog
qsopath = jhusdss_get_path(/qso)
if (keyword_set(boss)) then begin
   infile = jhusdss_boss_qsofile()
endif else begin
   infile =  jhusdss_dr7_qsofile()
endelse

qsofile = qsopath+'/'+infile
qsos = mrdfits(qsofile, 1)

;; pool of parent sample
absorbers = jhusdss_montecarlo_absorbers_pool(nmfver) 

jhusdss_montecarlo_spec, qsos, absorbers, nmfver, overwrite=overwrite, boss=boss, noresidual=noresidual

end
