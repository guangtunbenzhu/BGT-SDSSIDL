;+
; Documentation Needed!
; This doesn't work ...
;-

;; data[nspec, npix], ivar[nspec, npix], template[ntemplate, npix]
;; reflux[nreflux, npix], reflux[0, npix] has to be 1.,
;; I recommend reflux[1, npix] = x, reflux[2,npix]=x^2 etc...
;; I am not solving this iteratively at this point
function jhusdss_general_iterfit_simple, data, ivar, template=template, $
   chisquares=chisquares, reflux=reflux, ftweak=ftweak, maxiter=maxiter

if (size(data, /n_dim) ne 2) then $
   message, "Data must be in form of [Nspec, Npix], e.g., [1, 30000]"
if (size(template, /n_dim) ne 2 or (size(template))[1] ne 1) then $
   message, "Template must be in form of [1, Npix]."
if (max(reflux) gt 1.) then $
   message, "Reflux vector cannot be larger than unity."
if (n_elements(ftweak) eq 0) then ftweak = 0.05
if (ftweak le 0. or ftweak ge 1.) then $
   message, "Reflux limit has to be between 0 and 1."

nspec = (size(data))[1]
npix = (size(data))[2]
nreflux = (size(reflux))[1]

coeffs = fltarr(nspec, nreflux)
chisquares = fltarr(nspec)

;; Construct template
reflux_template = fltarr(nreflux, npix)
reflux_template[0,*] = template
for j=1L, nreflux-1L do reflux_template[j,*] = (reflux[j,*]+1./ftweak)*template

;iter=1L
;eps = 1.e-4
;while ((iter le maxiter) and (delta_chi2 le eps)) do begin
coeffs = jhusdss_general_linfit(data, ivar, template=reflux_template, $
                    chisquares=chisquares)
coeffs[*,0] = coeffs[*,0] + total(coeffs[*,1:*], 2)/ftweak

return, coeffs
end
