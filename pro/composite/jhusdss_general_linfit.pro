;+
; Documentation Needed!
; Comments:
;    - covariance not added, it's the inverse of alpha
;    - For non-linear least squares fit, use mpfit instead
;-

;; See numerical recipes 15.4 General Linear Least Squares Fit
;; data[nspec, npix], ivar[nspec, npix], template[ntemplate, npix]
function jhusdss_general_linfit, data, ivar, template=template, chisquares=chisquares

if (size(data, /n_dim) ne 2) then $
   message, "Data must be in form of [Nspec, Npix], e.g., [1, 3000]"
if (size(template, /n_dim) ne 2) then $
   message, "Template must be in form of [Nspec, Npix], e.g., [1, 3000]"

ntemplate = (size(template))[1]
if (ntemplate eq 1) then begin
   splog, "There is only one template, calling jhusdss_general_linfit_simple"
   coeff = jhusdss_general_linfit_simple(data, ivar, template=template, chisquares=chisquares)
   return, coeff
endif

nspec = (size(data))[1]
npix = (size(data))[2]

coeffs = fltarr(nspec, ntemplate)
chisquares = fltarr(nspec)
for i=0L, nspec-1L do begin
   ;; construct AA, AA_T, bb, alpha, bbeta
   ;; AA[ntemplate, npix], AA_T[npix, ntemplate]
   AA = template*sqrt(rebin(ivar[i,*], ntemplate, npix))
   AA_T = transpose(AA)
   ;; bb[1, npix]
   bb = data[i,*]*sqrt(ivar[i,*])
   ;; alpha[ntemplate, ntemplate] = AA_T ## AA
   alpha = AA_T##AA
   ;; bbeta[npix] = AA_T ## bb
   bbeta = reform(AA_T##bb)

   ;; use cholesky to invert alpha
   choldc, alpha, p
   coeffs[i, *] = cholsol(alpha, p, bbeta)
   chisquares[i] = total(ivar[i,*]*(data[i,*]-coeffs[i,*]#template)^2)
endfor

return, coeffs
end
