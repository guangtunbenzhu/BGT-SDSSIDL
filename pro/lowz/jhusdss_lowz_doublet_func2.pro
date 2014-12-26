;+
; See jhusdss_doublegaussian_fit
;-
;; same sigma
function jhusdss_lowz_doublet_func2, x, p

quadra=p[0]
slope=p[1]
intercept=p[2]
center=p[3]
separation=p[4]
flux=p[5]
ratio=p[6]
sigma=p[7]

model = fltarr(n_elements(x))
model = quadra*(x-center)^2+slope*(x-center)+intercept+flux*exp(-0.5*(x-center)^2/sigma^2)/sqrt(2.*!DPI)/sigma $
      + flux*ratio*exp(-0.5*(x-center-separation)^2/sigma^2)/sqrt(2.*!DPI)/sigma

return, model
end
