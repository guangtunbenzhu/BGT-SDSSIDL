function jhusdss_allskyspec_filename, flux=flux, boss=boss

if (keyword_set(boss)) then begin
   filename = 'BOSS_ALLSKY_SPEC.fits'
endif else begin
   filename = 'ALLSKY_SPEC.fits'
endelse

if (keyword_set(flux)) then return, repstr(filename, '.fits', '_flux.fits')
   
return, filename
end
