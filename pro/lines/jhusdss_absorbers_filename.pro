;+
; Documentation Needed!
;-
function jhusdss_absorbers_filename, nmfver, mgii=mgii, civ=civ, boss=boss, dr12=dr12

;; be careful this will return new values if these arguments are not set.
if (n_elements(prefix) eq 0) then prefix = 'Absorbers'
if (keyword_set(dr12)) then begin
   prefix = prefix+'_DR12'
endif else begin
   if (keyword_set(boss)) then begin
      prefix = prefix+'_BOSS'
   endif
endelse

if (keyword_set(mgii)) then suffix = 'MgII'
if (keyword_set(civ)) then suffix = 'CIV'
if (n_elements(suffix) eq 0) then suffix = 'MgII'

return, prefix+'_NMF_'+string(nmfver, format='(I3.3)')+'_'+suffix+'.fits'

end

