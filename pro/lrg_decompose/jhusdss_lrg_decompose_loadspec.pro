;+
; Documentation needed!
; Load the decomposed spectra
;-

function jhusdss_lrg_decompose_loadspec, plate, fiber, mjd, nmfver, boss=boss, error=error

;; path
if (keyword_set(boss)) then begin
    path=jhusdss_get_path(/fitlrg)+'/'+string(nmfver, format='(I3.3)')+'/Decompose_BOSS/'+string(plate, format='(i4.4)')
endif else begin
    path=jhusdss_get_path(/fitlrg)+'/'+string(nmfver, format='(I3.3)')+'/Decompose/'+string(plate, format='(i4.4)')
endelse

if (jhusdss_direxist(path) eq 0) then begin
    splog, "Can't find the directory."
    error = 1b
    return, -1
endif

if (n_elements(plate) ne 1 or n_elements(fiber) ne 1) then $
   message, "I can only load one spectrum for now."

filename = path+'/'+jhusdss_lrg_decompose_filename(plate, fiber, mjd, boss=boss)

error = 0b
if (file_test(filename) eq 0) then begin
    error = 1b
    return, -1
endif

return, mrdfits(filename, 1, /silent)

end
