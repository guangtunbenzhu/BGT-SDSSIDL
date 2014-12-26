function jhusdss_absorber_composite_image_readin, nmfver, boss=boss
if (n_elements(nmfver) eq 0) then message, 'nmfver required'

if (keyword_set(boss)) then begin
   path=jhusdss_get_path(/nmfqso)+'/'+string(nmfver, format='(I3.3)')+'/Absorbers_BOSS'
endif else begin
   path=jhusdss_get_path(/nmfqso)+'/'+string(nmfver, format='(I3.3)')+'/Absorbers'
endelse

fluxfile = path+'/'+jhusdss_absorber_composite_image_filename(nmfver, boss=boss)
return, mrdfits(fluxfile, 1)

end
