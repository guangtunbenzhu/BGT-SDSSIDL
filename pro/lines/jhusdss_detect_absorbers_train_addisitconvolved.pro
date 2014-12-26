;+
; Don't need this fragile piece of code
;-
pro jhusdss_detect_absorbers_train_addisitconvolved, nmfver, boss=boss, overwrite=overwrite

if (n_elements(nmfver) eq 0) then message, 'nmfver required'

;; the Pittsburgh master catalog
mgiifile = jhusdss_get_path(/absorber)+'/MgII/Master_Pitts_Catalog.fits'
mgii = mrdfits(mgiifile, 1)
nspec = n_elements(mgii)

str_tmp = {isitconvolved:1b}
outstr = replicate(str_tmp, nspec)

for i=0L, nspec-1L do begin

   counter, i+1, n_elements(mgii)

   ;; load convolved spectra
   spec = jhusdss_convolve_loadspec(mgii[i].plate, mgii[i].fiber, nmfver, boss=boss, error=error)
   if error then begin
      outstr[i].isitconvolved = 0b
      splog, "Can't find the convolved spectrum."
      continue
   endif

endfor

if (tag_exist(mgii[0], 'isitconvolved')) then begin
   out = mgii
   out.isitconvolved = outstr.isitconvolved
endif else  begin
   out = struct_addtags(mgii, outstr)
endelse

mwrfits, out, mgiifile, /create

end
