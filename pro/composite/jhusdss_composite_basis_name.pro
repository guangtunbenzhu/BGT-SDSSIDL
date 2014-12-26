;+
; Documentation Needed!
;-
function jhusdss_composite_basis_name, zbin, magbin, prefix=prefix, normminwave=normminwave

;; be careful this will return new values if these arguments are not set.
if (n_elements(prefix) eq 0) then prefix = 'QSO_composite'
if (n_elements(normminwave) eq 0) then normminwave=4150.d
if (n_elements(zbin) ne 1 or n_elements(magbin) ne 1) then $
    message, "Both zbin and magbin should be scalars."

return, prefix+'_z'+string(fix(zbin.min*100.), format='(I3.3)')+'_'+$
        string(fix(zbin.max*100.), format='(I3.3)')$
        +'_Mi'+string(fix(-magbin.min*100.), format='(I4.4)')+'_'+$
        string(fix(-magbin.max*100.), format='(I4.4)')$
        +'_norm'+string(fix(normminwave), format='(I4.4)')+'.fits'

end

