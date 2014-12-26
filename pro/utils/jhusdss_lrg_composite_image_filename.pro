function jhusdss_lrg_composite_image_filename, lrgver, boss=boss

if (n_elements(lrgver) eq 0) then message, 'lrgver required'

if (keyword_set(boss)) then begin
   filename = 'BOSS_LRG_composite_image_'+string(lrgver, format='(I3.3)')+'.fits'
endif else begin
   filename = 'LRG_composite_image_'+string(lrgver, format='(I3.3)')+'.fits'
endelse

return, filename
end
