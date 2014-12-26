;; moved to an independent file
;; same sigma
;function doublegaussian_line_func, x, p

;center1=p[0]
;center2=p[1]
;flux1=p[2]
;flux2=p[3]
;sigma=p[4]

;model = fltarr(n_elements(x))
;model = flux1*exp(-0.5*(x-center1)^2/sigma^2)/sqrt(2.*!DPI)/sigma $
;      + flux2*exp(-0.5*(x-center2)^2/sigma^2)/sqrt(2.*!DPI)/sigma

;return, model
;end

pro doublegaussian_line_fit, center, lflux, sigma, err_center, err_lflux, err_sigma

common com_sfit, use_lambda, use_flux, use_err, initline, initlflux, use_maxwidth, sigma_range

start=dblarr(5)

start[0]=initline[0]
start[1]=initline[1]

start[2]=initlflux[0]
start[3]=initlflux[1]

;; init dV = 207 km/s
;; dlog10(lambda) = 3D-4 (3 pixels)
start[4] = 2.0D-4*median(use_lambda)*alog(10.D0)

parinfo=replicate({limited:bytarr(2), limits:dblarr(2)}, 5)

;; CIV: 2.6/log(10)/1550 = 7.2D-4 ~7 pixels
;; MgII: 7.0/log(10)/2800 = 10.8D-4 ~11 pixels
;; resolution = 1.D-4
;; use maxwidth = 5D-4

;; center of line1
parinfo[0].limited=[1b, 1b]
parinfo[0].limits=initline[0]+[-use_maxwidth, use_maxwidth]

;; center of line2
parinfo[1].limited=[1b, 1b]
parinfo[1].limits=initline[1]+[-use_maxwidth, use_maxwidth]

;; flux of line1
parinfo[2].limited=[1b, 1b]
parinfo[2].limits=[1.D-10, 50.]

;; flux of line2
parinfo[3].limited=[1b, 1b]
parinfo[3].limits=[1.D-10, 50.]

;; sigma_range = [69km/s, 500km/s] = [1.D-4, 8.D-4]
parinfo[4].limited=[1b, 1b]
parinfo[4].limits=sigma_range

p=mpfitfun('doublegaussian_line_func', use_lambda, use_flux, use_err, $
           start, ftol=1.d-10, bestnorm=chi2, parinfo=parinfo, perror=perror, /double, /quiet)
;print, chi2

;; hack 09-15-2014
;; hope this works ...
center = dblarr(2)-999.
lflux = fltarr(2)-999.
err_center = dblarr(2)-999.
err_lflux = fltarr(2)-999.

if n_elements(perror) gt 0L then begin
   center[0] = p[0]
   center[1] = p[1]
   lflux[0] = p[2]
   lflux[1] = p[3]
   sigma = p[4]

   err_center[0] = perror[0]
   err_center[1] = perror[1]
   err_lflux[0] = perror[2]
   err_lflux[1] = perror[3]
   err_sigma = perror[4]
endif

end

pro jhusdss_doublegaussian_fit, in_lambda, in_flux, in_ivar, in_center, in_line_flux, $
       center=center, err_center=err_center, lflux=lflux, err_lflux=err_lflux, $
       sigma=sigma, err_sigma=err_sigma, maxwidth=maxwidth

common com_sfit

;cc = 2.99792458D5
;velunit = alog(10.D0)*cc
if n_elements(maxwidth) eq 0 then $
   maxwidth = 5.0D-4*median(in_lambda)*alog(10D0)

;; CIV: 2.6/log(10)/1550 = 7.2D-4 ~7 pixels
;; MgII: 7.0/log(10)/2800 = 10.8D-4 ~11 pixels
;; resolution = 1.D-4
;; use maxwidth = 5D-4

use_maxwidth = maxwidth
sigma_range = [0.50D-4, 20.0D-4]*median(in_lambda)*alog(10.D0)
;sigma_range = [1.00D-4, 3.0D-4]*median(in_lambda)*alog(10.D0)

lambda = in_lambda
flux = in_flux
ivar = in_ivar
initline = in_center
initlflux = (in_line_flux > 1.D-10)

iuse=where(ivar gt 0., nuse)
if(nuse eq 0) then return

use_lambda = lambda[iuse]
use_flux = flux[iuse]
use_err = 1./sqrt(ivar[iuse])

doublegaussian_line_fit, center, lflux, sigma, err_center, err_lflux, err_sigma

;; Monte Carlo?

end
