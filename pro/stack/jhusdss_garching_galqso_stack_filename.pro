function jhusdss_garching_galqso_stack_filename, nmfver, boss=boss

if (n_elements(nmfver) eq 0) then message, 'NMFVER required!'

if (keyword_set(boss)) then begin
   filename = 'GALQSO_STACK_'+string(nmfver, format='(I3.3)')+'_BOSS.fits'
endif else begin
   filename = 'GALQSO_STACK_'+string(nmfver, format='(I3.3)')+'.fits'
endelse
   
return, filename
end
