function jhusdss_lrg_composite_image_readin, lrgver, residual=residual, boss=boss
if (n_elements(lrgver) eq 0) then message, 'lrgver required'

lrgpath=jhusdss_get_path(/fitlrg)+'/'+string(lrgver, format='(I3.3)')+'/'

if (keyword_set(residual)) then begin
   resifile = lrgpath+'/Residual_'+jhusdss_lrg_composite_image_filename(lrgver, boss=boss)
   return, mrdfits(resifile, 1)
endif else begin
   fluxfile = lrgpath+'/'+jhusdss_lrg_composite_image_filename(lrgver, boss=boss)
   return, mrdfits(fluxfile, 1)
endelse

end
