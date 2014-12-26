function jhusdss_garching_gallrg_stack_filename, lrgver, boss=boss

if (n_elements(lrgver) eq 0) then message, 'NMFVER required!'

if (keyword_set(boss)) then begin
   filename = 'GALLRG_STACK_'+string(lrgver, format='(I3.3)')+'_BOSS.fits'
endif else begin
   filename = 'GALLRG_STACK_'+string(lrgver, format='(I3.3)')+'.fits'
endelse
   
return, filename
end
