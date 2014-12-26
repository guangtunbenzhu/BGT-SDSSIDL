;+
; Documentation Needed!
;-

;; one template only
;; data[nspec, npix], ivar[nspec, npix], template[1, npix]
function jhusdss_general_linfit_simple, data, ivar, template=template, chisquares=chisquares, $
         reduced_chi2=reduced_chi2

if (size(data, /n_dim) ne 2) then $
   message, "Data must be in form of [Nspec, Npix], e.g., [1, 30000]"
if (size(template, /n_dim) ne 2 or (size(template))[1] ne 1) then $
   message, "Template must be in form of [1, Npix]."

nspec = (size(data))[1]
npix = (size(data))[2]

coeffs = fltarr(nspec)
chisquares = fltarr(nspec)
reduced_chi2 = fltarr(nspec)

for i=0L, nspec-1L do begin
    coeffs[i] = total(ivar[i,*]*data[i,*]*template)/total(ivar[i,*]*template*template)
    chisquares[i] = total(ivar[i,*]*(data[i,*]-coeffs[i]#template)^2)
    reduced_chi2[i] = chisquares[i]/total(ivar[i,*] ne 0.)
endfor

return, coeffs
end
