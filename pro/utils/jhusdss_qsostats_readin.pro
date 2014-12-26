function jhusdss_qsostats_readin, nmfver, boss=boss, dr12=dr12

if (n_elements(nmfver) eq 0) then message, 'NMFVER required'

if (keyword_set(dr12)) then begin
    stat_path=jhusdss_get_path(/nmfqso)+'/'+string(nmfver, format='(I3.3)')+'/Decompose_DR12'
endif else begin
   if (keyword_set(boss)) then begin
      stat_path=jhusdss_get_path(/nmfqso)+'/'+string(nmfver, format='(I3.3)')+'/Decompose_BOSS'
   endif else begin
      stat_path=jhusdss_get_path(/nmfqso)+'/'+string(nmfver, format='(I3.3)')+'/Decompose'
   endelse
endelse

stat_file = stat_path+'/'+jhusdss_stat_filename(nmfver, boss=boss, dr12=dr12)

return, mrdfits(stat_file, 1)

end
