;+
; Documentation needed!
;-

function jhusdss_lrg_decompose_filename, plate, fiber, mjd, prefix=prefix, boss=boss

;; be careful this will return new values if these arguments are not set.
if (n_elements(prefix) eq 0) then prefix = 'LRG_spec_decomposed'

if (not keyword_set(boss)) then begin
    return, prefix+'_'+string(plate, format='(I4.4)')+'_'+string(fiber, format='(I3.3)') $
  +'_'+string(mjd, format='(I5.5)')+'.fits'
endif else begin
    return, prefix+'_'+string(plate, format='(I4.4)')+'_'+string(fiber, format='(I4.4)') $
  +'_'+string(mjd, format='(I5.5)')+'.fits'
endelse


end

