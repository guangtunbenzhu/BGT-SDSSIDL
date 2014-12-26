;+
; Documentation needed!
;-

function jhusdss_decompose_name, plate, fiber, prefix=prefix

;; be careful this will return new values if these arguments are not set.
if (n_elements(prefix) eq 0) then prefix = 'QSO_spec_decomposed'

return, prefix+'_'+string(plate, format='(I4.4)')+'_'+string(fiber, format='(I3.3)') $
  +'.fits'

end

