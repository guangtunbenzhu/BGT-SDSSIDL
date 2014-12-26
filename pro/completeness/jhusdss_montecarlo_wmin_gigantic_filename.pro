;+
; Documentation needed!
;-

function jhusdss_montecarlo_wmin_gigantic_filename, nmfver, prefix=prefix

;; be careful this will return new values if these arguments are not set.
if (n_elements(prefix) eq 0) then prefix = 'All_in_One_QSO_spec_MonteCarlo_wmin'

return, prefix+'_'+string(nmfver, format='(I3.3)')+'.fits'

end

