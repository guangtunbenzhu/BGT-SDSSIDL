function jhusdss_qso_composite_image_readin, nmfver, residual=residual, boss=boss
if (n_elements(nmfver) eq 0) then message, 'nmfver required'

qsopath=jhusdss_get_path(/nmfqso)+'/'+string(nmfver, format='(I3.3)')+'/'

if (keyword_set(residual)) then begin
   resifile = qsopath+'/Residual_'+jhusdss_qso_composite_image_filename(nmfver, boss=boss)
   return, mrdfits(resifile, 1)
endif else begin
   fluxfile = qsopath+'/'+jhusdss_qso_composite_image_filename(nmfver, boss=boss)
   return, mrdfits(fluxfile, 1)
endelse

end
