;+
; See jhusdss_doublegaussian_fit
;-
;; same sigma
function jhusdss_lowz_doublet_func, x, p

slope=p[0]
intercept=p[1]
center=p[2]
separation=p[3]
flux=p[4]
ratio=p[5]
sigma=p[6]

model = fltarr(n_elements(x))
model = slope*(x-center)+intercept+flux*exp(-0.5*(x-center)^2/sigma^2)/sqrt(2.*!DPI)/sigma $
      + flux*ratio*exp(-0.5*(x-center-separation)^2/sigma^2)/sqrt(2.*!DPI)/sigma

return, model
end
