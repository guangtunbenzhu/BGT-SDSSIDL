;; same sigma
function singlegaussian_line_func, x, p

center=p[0]
flux=p[1]
sigma=p[2]

model = fltarr(n_elements(x))
model = flux*exp(-0.5*(x-center)^2/sigma^2)/sqrt(2.*!DPI)/sigma

return, model
end

pro singlegaussian_line_fit, center, lflux, sigma, err_center, err_lflux, err_sigma

common com_sfit, use_lambda, use_flux, use_err, initline, initlflux, use_maxwidth, sigma_range

start=dblarr(3)

start[0]=initline
start[1]=initlflux

;; init dV = 207 km/s
;; dlog10(lambda) = 3D-4 (3 pixels)
start[2] = 2.0D-4*median(use_lambda)*alog(10.D0)

parinfo=replicate({limited:bytarr(2), limits:dblarr(2)}, 3)

;; center
parinfo[0].limited=[1b, 1b]
parinfo[0].limits=initline+[-use_maxwidth, use_maxwidth]

;; flux
parinfo[1].limited=[1b, 1b]
parinfo[1].limits=[1D-10, 50.]

;; sigma_range = [69km/s, 500km/s] = [1.D-4, 8.D-4]
parinfo[2].limited=[1b, 1b]
parinfo[2].limits=sigma_range

p=mpfitfun('singlegaussian_line_func', use_lambda, use_flux, use_err, $
           start, ftol=1.d-10, bestnorm=chi2, parinfo=parinfo, perror=perror, /double, /quiet)
;print, chi2

center = p[0]
lflux = p[1]
sigma = p[2]

err_center = perror[0]
err_lflux = perror[1]
err_sigma = perror[2]

end

pro jhusdss_singlegaussian_fit, in_lambda, in_flux, in_ivar, in_center, in_line_flux, $
       center=center, err_center=err_center, lflux=lflux, err_lflux=err_lflux, $
       sigma=sigma, err_sigma=err_sigma, maxwidth=maxwidth

common com_sfit

;cc = 2.99792458D5
;velunit = alog(10.D0)*cc

if n_elements(maxwidth) eq 0 then $
   maxwidth = 5.0D-4*median(in_lambda)*alog(10D0)
sigma_range = [0.5D-4, 20.0D-4]*median(in_lambda)*alog(10.D0)

use_maxwidth = maxwidth

lambda = in_lambda
flux = in_flux
ivar = in_ivar
initline = in_center
initlflux = (in_line_flux > 1D-10)

iuse=where(ivar gt 0., nuse)
if(nuse eq 0) then return

use_lambda = lambda[iuse]
use_flux = flux[iuse]
use_err = 1./sqrt(ivar[iuse])

singlegaussian_line_fit, center, lflux, sigma, err_center, err_lflux, err_sigma

;; Monte Carlo Errors?

end
