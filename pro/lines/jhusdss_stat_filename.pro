;+
; Documentation Needed!
;-
function jhusdss_stat_filename, nmfver, boss=boss, dr12=dr12

if (keyword_set(dr12)) then begin
   return, 'QSO_decompose_NMF_DR12_'+string(nmfver, format='(I3.3)')+'stats.fits'
endif else begin
   if (keyword_set(boss)) then begin
      return, 'QSO_decompose_NMF_BOSS_'+string(nmfver, format='(I3.3)')+'stats.fits'
   endif else begin
      return, 'QSO_decompose_NMF_'+string(nmfver, format='(I3.3)')+'stats.fits'
   endelse
endelse

end

