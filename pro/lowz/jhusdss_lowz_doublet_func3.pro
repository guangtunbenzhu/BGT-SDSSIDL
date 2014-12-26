;+
; See jhusdss_doublegaussian_fit
;-
;; same sigma
function jhusdss_lowz_doublet_func3, x, p

quadra=p[0]
slope=p[1]
intercept=p[2]
center=p[3]
separation=p[4]
flux=p[5]
ratio=p[6]
sigma1=p[7]
sigma2=p[8]

model = fltarr(n_elements(x))
model = quadra*(x-center)^2+slope*(x-center)+intercept+flux*exp(-0.5*(x-center)^2/sigma1^2)/sqrt(2.*!DPI)/sigma1 $
      + flux*ratio*exp(-0.5*(x-center-separation)^2/sigma2^2)/sqrt(2.*!DPI)/sigma2

return, model
end
