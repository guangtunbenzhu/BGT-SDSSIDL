;+
; Documentation needed!
;-
function jhusdss_decompose_donespec, nmfver, boss=boss

if (n_elements(nmfver) eq 0) then message, "NMF version required."
if (keyword_set(boss)) then begin
    path = jhusdss_get_path(/nmfqso)+'/'+string(nmfver, format='(I3.3)')+'/Decompose_BOSS/*/'
endif else begin
    path = jhusdss_get_path(/nmfqso)+'/'+string(nmfver, format='(I3.3)')+'/Decompose/*/'
endelse

files = file_search(path+'/*.fits')
rootpos = strpos(files[0], '.fits')-8L
plate = long(strmid(files, rootpos, 4L))
fiber = long(strmid(files, rootpos+5, 3L))

strtmp = replicate({plate:plate[0], fiber:fiber[0]}, n_elements(plate))
strtmp.plate = plate
strtmp.fiber = fiber

return, strtmp

end
