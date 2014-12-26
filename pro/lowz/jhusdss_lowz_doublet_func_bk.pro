;+
; See jhusdss_doublegaussian_fit
;-
;; same sigma
function jhusdss_lowz_doublet_func, x, p

clevel=p[0]
center1=p[1]
center2=p[2]
flux1=p[3]
flux2=p[4]
sigma=p[5]

model = fltarr(n_elements(x))
model = clevel+flux1*exp(-0.5*(x-center1)^2/sigma^2)/sqrt(2.*!DPI)/sigma $
      + flux2*exp(-0.5*(x-center2)^2/sigma^2)/sqrt(2.*!DPI)/sigma

return, model
end

end
