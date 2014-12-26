;; same sigma
function doublegaussian_line_func, x, p

center1=p[0]
center2=p[1]
flux1=p[2]
flux2=p[3]
sigma=p[4]

model = fltarr(n_elements(x))
model = flux1*exp(-0.5*(x-center1)^2/sigma^2)/sqrt(2.*!DPI)/sigma $
      + flux2*exp(-0.5*(x-center2)^2/sigma^2)/sqrt(2.*!DPI)/sigma

return, model
end
