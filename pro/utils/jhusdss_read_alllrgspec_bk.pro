;; obsolete
function jhusdss_read_alllrgspec, lrgver, flux=flux, continuum=continuum, normresi=normresi, subtresi=subtresi, boss=boss

if (n_elements(lrgver) eq 0) then message, 'lrgver required'
;; lrgpath
lrgpath=jhusdss_get_path(/fitlrg)+'/'+string(lrgver, format='(I3.3)')+'/'
infile = lrgpath+'/AllInOne/'+jhusdss_alllrgspec_filename(lrgver, flux=flux, continuum=continuum, normresi=normresi, subtresi=subtresi, boss=boss)

return, mrdfits(infile, 1)

end
