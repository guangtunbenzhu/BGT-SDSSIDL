;+
; Documentation needed!
; This is a dangerous piece of code that will modify inputs
;-

;; sanity check
;; input will be modified on the fly
;; flux[nspec, npix], ivar[nspec, npix]
pro jhusdss_nmf_sanity_check, objs, loglam, flux=flux, ivar=ivar, $
       spec_index=spec_index, pix_index=pix_index

;; if a spectrum is missing all pixels, we don't use it
totalflux = total(flux, 2)
spec_index = where(totalflux gt 0., nspec)
if (nspec eq 0) then message, 'No useful spectra!'
if (nspec ne (size(flux))[1]) then begin
   objs = objs[spec_index]
   flux = flux[spec_index, *]
   ivar = ivar[spec_index, *]
endif

;; if a pixel is missing in every galaxy, we don't use it
totalflux = total(flux, 1)
pix_index = where(totalflux gt 0., npix)
if (npix eq 0) then message, 'No useful spectra!'
if (npix ne (size(flux))[2]) then begin
   loglam = loglam[pix_index]
   flux = flux[*, pix_index]
   ivar = ivar[*, pix_index]
endif

end

