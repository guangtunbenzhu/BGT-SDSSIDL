pro jhusdss_montecarlo_wmin_all, nmfver, boss=boss, overwrite=overwrite, noresidual=noresidual

if (n_elements(nmfver) eq 0) then $
   message, 'nmfver required'

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

jhusdss_montecarlo_wmin_spec, qsos, nmfver, overwrite=overwrite, boss=boss, noresidual=noresidual

end
