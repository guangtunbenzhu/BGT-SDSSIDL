;; obsolete
function jhusdss_read_alllrgspec, lrgver, wave=wave, index=index, flux=flux, continuum=continuum, normresi=normresi, subtresi=subtresi, corr=corr, boss=boss

if (n_elements(lrgver) eq 0) then message, 'lrgver required'
;lrg; lrgpath
parentpath = jhusdss_get_parent_path()+'SDSS/LRG/'+string(lrgver, format='(I3.3)')+'/'
;lrgpath=jhusdss_get_path(/fitlrg)+'/'+string(lrgver, format='(I3.3)')+'/'
infile = parentpath+'/AllInOne/'+jhusdss_alllrgspec_filename(lrgver, wave=wave, index=index, flux=flux, continuum=continuum, normresi=normresi, subtresi=subtresi, corr=corr, boss=boss)

return, mrdfits(infile, 1)

end
