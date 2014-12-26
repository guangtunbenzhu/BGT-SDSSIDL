function jhusdss_allqsospec_filename, nmfver, flux=flux, continuum=continuum, normresi=normresi, subtresi=subtresi, boss=boss

if (n_elements(nmfver) eq 0) then message, 'NMFVER required!'

if (keyword_set(boss)) then begin
   filename = 'BOSS_ALLQSO_SPEC_'+string(nmfver, format='(I3.3)')+'.fits'
endif else begin
   filename = 'ALLQSO_SPEC_'+string(nmfver, format='(I3.3)')+'.fits'
endelse
   
if (keyword_set(flux)) then return, repstr(filename, '.fits', '_flux.fits')
if (keyword_set(continuum)) then return, repstr(filename, '.fits', '_continuum.fits')
if (keyword_set(normresi)) then return, repstr(filename, '.fits', '_normalized_residual.fits')
if (keyword_set(subtresi)) then return, repstr(filename, '.fits', '_subtracted_residual.fits')

return, filename
end
