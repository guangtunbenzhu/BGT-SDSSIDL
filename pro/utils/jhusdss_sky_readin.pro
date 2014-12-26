function jhusdss_sky_readin, boss=boss

skypath = jhusdss_get_path(/sky)
if (keyword_set(boss)) then begin
   filename = 'all_sky_boss.fits'
endif else begin
   filename =  'all_sky_dr7.fits'
endelse

infile = skypath+'/'+filename

return, mrdfits(infile, 1)

end
