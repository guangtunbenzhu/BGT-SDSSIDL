function jhusdss_decompose_allinone_filename, nmfver, $
         boss=boss, continuum=continuum, flux=flux, residual=residual

   if (n_elements(nmfver) eq 0) then message, 'NMFVER required'
   prefix = ''
   if (keyword_set(continuum)) then prefix = 'Continuum_'
   if (keyword_set(flux)) then prefix = 'Flux_'
   if (keyword_set(residual)) then prefix = 'Residual_'
   prefix = 'All_in_one_'+prefix

   suffix = string(nmfver, format='(i3.3)')
   if (keyword_set(boss)) then suffix = suffix+'_BOSS'

   return, prefix+suffix+'.fits'
end
