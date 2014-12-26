function jhusdss_galqso_match_spec_readin, nmfver, minrad, maxrad, boss=boss, match=match, flag=flag

if (n_elements(nmfver) eq 0) then message, 'nmfver required!'
if (n_elements(minrad) eq 0) then message, 'minrad required!'
if (n_elements(maxrad) eq 0) then message, 'maxrad required!'

if (keyword_set(boss)) then begin
   lowzpath=jhusdss_get_path(/nmfqso)+'/'+string(nmfver, format='(I3.3)')+'/Lowz_BOSS'
endif else begin
   lowzpath=jhusdss_get_path(/nmfqso)+'/'+string(nmfver, format='(I3.3)')+'/Lowz'
endelse


if (arg_present(flag)) then begin
   flagfile = lowzpath+'/Flag_Rad'+string(minrad,format='(f4.2)')+'_'+string(maxrad,format='(f4.2)')+'_' $
       + jhusdss_lowz_spec_filename(nmfver)
   flag = mrdfits(flagfile, 1)
endif

infile = lowzpath+'/Rad'+string(minrad,format='(f4.2)')+'_'+string(maxrad,format='(f4.2)')+'_' $
       + jhusdss_lowz_spec_filename(nmfver)

if (arg_present(match)) then match = mrdfits(infile, 1)
return, mrdfits(infile, 2)

end
