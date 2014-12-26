function jhusdss_lrgstats_readin, lrgver, boss=boss

if (n_elements(lrgver) eq 0) then message, 'NMFVER required'

if (keyword_set(boss)) then begin
    stat_path=jhusdss_get_path(/fitlrg)+'/'+string(lrgver, format='(I3.3)')+'/Decompose_BOSS'
endif else begin
    stat_path=jhusdss_get_path(/fitlrg)+'/'+string(lrgver, format='(I3.3)')+'/Decompose'
endelse

stat_file = stat_path+'/'+jhusdss_lrg_stat_filename(lrgver, boss=boss)

return, mrdfits(stat_file, 1)

end
