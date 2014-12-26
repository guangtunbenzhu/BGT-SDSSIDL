function jhusdss_read_allqsospec, nmfver, flux=flux, continuum=continuum, normresi=normresi, subtresi=subtresi, boss=boss

if (n_elements(nmfver) eq 0) then message, 'nmfver required'
;; qsopath
qsopath=jhusdss_get_path(/nmfqso)+'/'+string(nmfver, format='(I3.3)')+'/'
infile = qsopath+'/AllInOne/'+jhusdss_allqsospec_filename(nmfver, flux=flux, continuum=continuum, normresi=normresi, subtresi=subtresi, boss=boss)

return, mrdfits(infile, 1)

end
