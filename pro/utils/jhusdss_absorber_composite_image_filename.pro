function jhusdss_absorber_composite_image_filename, nmfver, boss=boss

if (n_elements(nmfver) eq 0) then message, 'nmfver required'

if (keyword_set(boss)) then begin
   filename = 'BOSS_Absorber_composite_image_'+string(nmfver, format='(I3.3)')+'.fits'
endif else begin
   filename = 'Absorber_composite_image_'+string(nmfver, format='(I3.3)')+'.fits'
endelse

return, filename

end
