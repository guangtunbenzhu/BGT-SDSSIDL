function jhusdss_alllrgspec_filename, lrgver, wave=wave, index=index, flux=flux, continuum=continuum, normresi=normresi, subtresi=subtresi, corr=corr, boss=boss

if (n_elements(lrgver) eq 0) then message, 'LRGVER required!'

if (keyword_set(boss)) then begin
   filename = 'DR9_LRG_ALL.fit'
endif else begin
   filename = 'DR7_LRG_ALL.fit'
endelse

if (keyword_set(wave)) then return, repstr(filename, '.fit', '_WAVE.fit')
if (keyword_set(index)) then return, repstr(filename, '.fit', '_INDEX.fit')
if (keyword_set(flux)) then return, repstr(filename, '.fit', '_FLUX.fit')
if (keyword_set(continuum)) then return, repstr(filename, '.fit', '_CONTINUUM_'+string(lrgver, format='(I3.3)')+'.fit')
if (keyword_set(normresi)) then return, repstr(filename, '.fit', '_NORMALIZED_RESIDUAL_'+string(lrgver, format='(I3.3)')+'.fit')
if (keyword_set(subtresi)) then return, repstr(filename, '.fit', '_SUBTRACTED_RESIDUAL_'+string(lrgver, format='(I3.3)')+'.fit')
if (keyword_set(corr)) then return, repstr(filename, '.fit', '_CORRECTION_'+string(lrgver, format='(I3.3)')+'.fit')

return, filename
end
