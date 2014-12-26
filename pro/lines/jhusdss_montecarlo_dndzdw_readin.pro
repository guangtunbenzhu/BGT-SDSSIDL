function jhusdss_montecarlo_dndzdw_readin, nmfver, boss=boss

if (keyword_set(boss)) then begin
   path=jhusdss_get_path(/nmfqso)+'/'+$
        string(nmfver, format='(I3.3)')+'/MonteCarlo_BOSS'
endif else begin
   path=jhusdss_get_path(/nmfqso)+'/'+$
        string(nmfver, format='(I3.3)')+'/MonteCarlo'
endelse

filename = jhusdss_montecarlo_dndzdw_filename(nmfver)
infile = path+'/'+filename

return, mrdfits(infile, 1)

end
