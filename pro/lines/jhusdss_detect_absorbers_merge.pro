pro jhusdss_detect_absorbers_merge, nmfver, nprocess=nprocess, overwrite=overwrite, boss=boss, dr12=dr12

if (n_elements(nmfver) eq 0) then message, 'nmfver required'

if (n_elements(nprocess) eq 0) then nprocess = 8L

if (keyword_set(dr12)) then begin
   path=jhusdss_get_path(/nmfqso)+'/'+string(nmfver, format='(I3.3)')+'/Absorbers_DR12'
endif else begin
   if (keyword_set(boss)) then begin
      path=jhusdss_get_path(/nmfqso)+'/'+string(nmfver, format='(I3.3)')+'/Absorbers_BOSS'
   endif else begin
      path=jhusdss_get_path(/nmfqso)+'/'+string(nmfver, format='(I3.3)')+'/Absorbers'
   endelse
endelse

outfile = path+'/'+jhusdss_absorbers_filename(nmfver, /mgii, boss=boss, dr12=dr12)

if (file_test(outfile) and ~keyword_set(overwrite)) then begin
      splog, 'File already exists. Use /overwrite if you want to overwrite it.'
      return
endif

for i=1L, nprocess do begin
    filename = path+'/Pro_'+string(i, format='(i1.1)')+'_'+jhusdss_absorbers_filename(nmfver, /mgii, boss=boss, dr12=dr12)
    tmp = mrdfits(filename,1)
    if (n_elements(absorbers) eq 0) then begin
       absorbers = tmp
    endif else begin
       absorbers = [absorbers, tmp]
    endelse 
endfor

mwrfits, absorbers, outfile, /create

end
