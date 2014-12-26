;; obsolete
function jhusdss_read_alllrgspec_old, lrgver, flux=flux, continuum=continuum, normresi=normresi, subtresi=subtresi, boss=boss

if (n_elements(lrgver) eq 0) then message, 'lrgver required'
;; lrgpath
lrgpath=jhusdss_get_path(/fitlrg)+'/'+string(lrgver, format='(I3.3)')+'/'

if (keyword_set(flux)) then filename = 'Flux_ALLLRG_SPEC_101.fits'
if (keyword_set(continuum)) then filename = 'Continuum_ALLLRG_SPEC_101.fits'
if (keyword_set(normresi)) then filename = 'Residual_ALLLRG_SPEC_101.fits'

infile = lrgpath+'/'+filename

return, mrdfits(infile, 1)

end
