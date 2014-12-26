;+
; Don't need this fragile piece of code
;-
pro jhusdss_detect_absorbers_train_addmedian, nmfver, overwrite=overwrite

if (n_elements(nmfver) eq 0) then message, 'nmfver required'

;; the Pittsburgh master catalog
mgiifile = jhusdss_get_path(/absorber)+'/MgII/Master_Pitts_Catalog.fits'
mgii = mrdfits(mgiifile, 1)
nspec = n_elements(mgii)

str_tmp = {spec_snr_median:0.}
outstr = replicate(str_tmp, nspec)

for i=0L, nspec-1L do begin

   counter, i+1, n_elements(mgii)

   ;; load decomposed spectra
   orispec = jhusdss_decompose_loadspec(mgii[i].plate, mgii[i].fiber, nmfver, boss=boss, error=error)
   if error then begin
      splog, "Can't find the decomposed spectrum."
      continue
   endif

   ii = where(orispec.ivar gt 0., nn)
   if (nn gt 0) then outstr[i].spec_snr_median = median(orispec.flux[ii]*sqrt(orispec.ivar[ii]))

endfor

if (tag_exist(mgii[0], 'spec_snr_median')) then begin
   out = mgii
   out.spec_snr_median = outstr.spec_snr_median
endif else  begin
   out = struct_addtags(mgii, outstr)
endelse

mwrfits, out, mgiifile, /create

end
