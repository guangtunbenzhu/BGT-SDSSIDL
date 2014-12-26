;+
; See jhusdss_doublegaussian_fit
;-
pro jhusdss_lowz_doublet_line_fit, center, lflux, sigma, clevel, err_center, err_lflux, err_sigma, err_clevel

common com_lowz_fit, use_lambda, use_flux, use_err, initline, initlflux, initclevel, use_maxwidth, sigma_range

start=dblarr(6)

start[0]=initclevel

start[1]=initline[0]
start[2]=initline[1]

start[3]=initlflux[0]
start[4]=initlflux[1]

;; init dV = 207 km/s
;; dlog10(lambda) = 3D-4 (3 pixels)
start[5] = 2.0D-4*median(use_lambda)*alog(10.D0)

parinfo=replicate({limited:bytarr(2), limits:dblarr(2)}, 6)

;; CIV: 2.6/log(10)/1550 = 7.2D-4 ~7 pixels
;; MgII: 7.0/log(10)/2800 = 10.8D-4 ~11 pixels
;; resolution = 1.D-4
;; use maxwidth = 5D-4

;; continuum level 
parinfo[0].limited=[1b, 1b]
parinfo[0].limits=[-0.01, 0.01]

;; center of line1
parinfo[1].limited=[1b, 1b]
parinfo[1].limits=initline[0]+[-use_maxwidth, use_maxwidth]

;; center of line2
parinfo[2].limited=[1b, 1b]
parinfo[2].limits=initline[1]+[-use_maxwidth, use_maxwidth]

;; flux of line1
parinfo[3].limited=[1b, 1b]
parinfo[3].limits=[1.D-10, 50.]

;; flux of line2
parinfo[4].limited=[1b, 1b]
parinfo[4].limits=[1.D-10, 50.]

;; sigma_range = [69km/s, 500km/s] = [1.D-4, 8.D-4]
parinfo[5].limited=[1b, 1b]
parinfo[5].limits=sigma_range

p=mpfitfun('jhusdss_lowz_doublet_func', use_lambda, use_flux, use_err, $
           start, ftol=1.d-10, bestnorm=chi2, parinfo=parinfo, perror=perror, /double, /quiet)
;print, chi2

center = dblarr(2)
lflux = fltarr(2)
err_center = dblarr(2)
err_lflux = fltarr(2)

clevel  = p[0]
center[0] = p[1]
center[1] = p[2]
lflux[0] = p[3]
lflux[1] = p[4]
sigma = p[5]

err_clevel = perror[0]
err_center[0] = perror[1]
err_center[1] = perror[2]
err_lflux[0] = perror[3]
err_lflux[1] = perror[4]
err_sigma = perror[5]

end

pro jhusdss_lowz_doublet_fit, in_lambda, in_flux, in_ivar, in_center, in_line_flux, $
       in_clevel, center=center, err_center=err_center, lflux=lflux, err_lflux=err_lflux, $
       sigma=sigma, err_sigma=err_sigma, maxwidth=maxwidth, clevel=clevel, err_clevel=err_clevel

common com_lowz_fit

;cc = 2.99792458D5
;velunit = alog(10.D0)*cc
if n_elements(maxwidth) eq 0 then $
   maxwidth = 2.0D-4*median(in_lambda)*alog(10D0)

;; CIV: 2.6/log(10)/1550 = 7.2D-4 ~7 pixels
;; MgII: 7.0/log(10)/2800 = 10.8D-4 ~11 pixels
;; resolution = 1.D-4
;; use maxwidth = 5D-4

use_maxwidth = maxwidth
;sigma_range = [0.50D-4, 20.0D-4]*median(in_lambda)*alog(10.D0)
sigma_range = [1.00D-4, 10.0D-4]*median(in_lambda)*alog(10.D0)

lambda = in_lambda
flux = in_flux
ivar = in_ivar
initline = in_center
initlflux = (in_line_flux > 1.D-10)
initclevel = ((in_clevel > (-0.01)) < 0.01)

iuse=where(ivar gt 0., nuse)
if(nuse eq 0) then return

use_lambda = lambda[iuse]
use_flux = flux[iuse]
use_err = 1./sqrt(ivar[iuse])

jhusdss_lowz_doublet_line_fit, center, lflux, sigma, clevel, err_center, err_lflux, err_sigma, err_clevel

;; Monte Carlo?

end
