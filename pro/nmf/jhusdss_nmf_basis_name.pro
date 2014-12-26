;+
; Documentation Needed!
;-
function jhusdss_nmf_basis_name, zmin, zmax, prefix=prefix, normminwave=normminwave

;; be careful this will return new values if these arguments are not set.
if (n_elements(prefix) eq 0) then prefix = 'QSO_NMF_basis'
if (n_elements(normminwave) eq 0) then normminwave=4150.d
if (n_elements(zmin) ne 1 or n_elements(zmax) ne 1) then message, "Both zmin and zmax should be scalars."

return, prefix+'_z'+string(fix(zmin*100.), format='(I3.3)')+'_'+string(fix(zmax*100.), format='(I3.3)')$
  +'_norm'+string(fix(normminwave), format='(I4.4)')+'.fits'

end

