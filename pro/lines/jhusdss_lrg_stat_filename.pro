;+
; Documentation Needed!
;-
function jhusdss_lrg_stat_filename, lrgver, boss=boss

if (keyword_set(boss)) then begin
   return, 'LRG_decompose_NMF_BOSS_'+string(lrgver, format='(I3.3)')+'_stats.fits'
endif else begin
   return, 'LRG_decompose_NMF_'+string(lrgver, format='(I3.3)')+'_stats.fits'
endelse

end

