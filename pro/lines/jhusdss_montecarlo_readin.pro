function jhusdss_montecarlo_readin, nmfver, boss=boss, nocal=nocal

if (keyword_set(boss)) then begin
   path=jhusdss_get_path(/nmfqso)+'/'+$
        string(nmfver, format='(I3.3)')+'/MonteCarlo_BOSS'
endif else begin
   path=jhusdss_get_path(/nmfqso)+'/'+$
        string(nmfver, format='(I3.3)')+'/MonteCarlo'
endelse

filename = jhusdss_montecarlo_completeness_filename(nmfver)
infile = path+'/'+filename
if keyword_set(nocal) then infile = path+'/Nocal_'+filename

return, mrdfits(infile, 1)

end
