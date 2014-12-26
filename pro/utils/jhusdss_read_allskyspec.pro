;; obsolete
function jhusdss_read_allskyspec, flux=flux, boss=boss

;; qsopath
skypath=jhusdss_get_path(/sky)+'/AllInOne/'
infile = skypath+'/'+jhusdss_allskyspec_filename(flux=flux, boss=boss)

return, mrdfits(infile, 1)

end
