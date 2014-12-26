;+
; See jhusdss_doublegaussian_fit
;-
;pro jhusdss_lowz_doublet_line_fit, center, lflux, sigma, clevel, err_center, err_lflux, err_sigma, err_clevel
pro jhusdss_lowz_doublet_line_fit, slope, intercept, center, separation, lflux, ratio, sigma, $
        err_slope, err_intercept, err_center, err_separation, err_lflux, err_ratio, err_sigma

common com_lowz_fit, use_lambda, use_flux, use_err, $
         init_slope, init_intercept, init_center, init_separation, init_lflux, init_ratio, init_sigma, $
         use_maxwidth, sigma_range

start=dblarr(7)

start[0]=init_slope
start[1]=init_intercept

start[2]=init_center
start[3]=init_separation

start[4]=init_lflux
start[5]=init_ratio

start[6]=init_sigma

;; init dV = 207 km/s
;; dlog10(lambda) = 3D-4 (3 pixels)
;; start[5] = 2.0D-4*median(use_lambda)*alog(10.D0)

parinfo=replicate({limited:bytarr(2), limits:dblarr(2), fixed:0b}, 7)

;; CIV: 2.6/log(10)/1550 = 7.2D-4 ~7 pixels
;; MgII: 7.0/log(10)/2800 = 10.8D-4 ~11 pixels
;; resolution = 1.D-4
;; use maxwidth = 5D-4

;; continuum level slope
parinfo[0].limited=[1b, 1b]
parinfo[0].limits=[-0.1, 0.1]
;parinfo[0].fixed=1b

;; continuum level intercept
parinfo[1].limited=[1b, 1b]
parinfo[1].limits=[-0.01, 0.01]
;parinfo[1].fixed=1b

;; center of line1
parinfo[2].limited=[1b, 1b]
parinfo[2].limits=init_center+[-use_maxwidth, use_maxwidth]
parinfo[2].fixed=1b

;; separation between line1 and line2
parinfo[3].limited=[1b, 1b]
parinfo[3].limits=init_separation+[-use_maxwidth, use_maxwidth]
parinfo[3].fixed=1b

;; flux of line1
parinfo[4].limited=[1b, 1b]
parinfo[4].limits=[1.D-10, 50.]
;parinfo[4].fixed=1b

;; flux ratio between line2 and line1
;; Ca II 3969/3934 = 1./2.
;; Mg II 2803/2796 = 1./2.
parinfo[5].limited=[1b, 1b]
parinfo[5].limits=[1.D-10, 10.]
;parinfo[5].fixed=1b

;; sigma_range = [69km/s, 500km/s] = [1.D-4, 8.D-4]
parinfo[6].limited=[1b, 1b]
parinfo[6].limits=sigma_range
;parinfo[6].fixed=1b

p=mpfitfun('jhusdss_lowz_doublet_func', use_lambda, use_flux, use_err, $
           start, ftol=1.d-10, bestnorm=chi2, parinfo=parinfo, perror=perror, /double, /quiet)
;print, chi2

slope = p[0]
intercept = p[1]
center = p[2]
separation = p[3]
lflux = p[4]
ratio = p[5]
sigma = p[6]

err_slope = perror[0]
err_intercept = perror[1]
err_center = perror[2]
err_separation = perror[3]
err_lflux = perror[4]
err_ratio = perror[5]
err_sigma = perror[6]

end

pro jhusdss_lowz_doublet_fit, in_lambda, in_flux, in_ivar, $
       in_slope, in_intercept, in_center, in_separation, in_lflux, in_ratio, in_sigma, $
       slope=slope, intercept=intercept, center=center, separation=separation, $
       lflux=lflux, ratio=ratio, sigma=sigma, $
       err_slope=err_slope, err_intercept=err_intercept, err_center=err_center, $
       err_separation=err_separation, err_lflux=err_lflux, err_ratio=err_ratio, err_sigma=err_sigma, $
       maxwidth=maxwidth

common com_lowz_fit

;cc = 2.99792458D5
;velunit = alog(10.D0)*cc
if n_elements(maxwidth) eq 0 then $
   maxwidth = 3.0D-4*median(in_lambda)*alog(10D0)

;; CIV: 2.6/log(10)/1550 = 7.2D-4 ~7 pixels
;; MgII: 7.0/log(10)/2800 = 10.8D-4 ~11 pixels
;; resolution = 1.D-4
;; use maxwidth = 5D-4

use_maxwidth = maxwidth
;sigma_range = [0.50D-4, 20.0D-4]*median(in_lambda)*alog(10.D0)
sigma_range = [1.00D-4, 3.0D-4]*median(in_lambda)*alog(10.D0)

lambda = in_lambda
flux = in_flux
ivar = in_ivar

init_slope = in_slope
init_intercept = ((in_intercept > (-0.01)) < 0.01)
init_center = in_center
init_separation = in_separation
init_lflux = (in_lflux > 1.D-10)
init_ratio = in_ratio
init_sigma = in_sigma

iuse=where(ivar gt 0., nuse)
if(nuse eq 0) then return

use_lambda = lambda[iuse]
use_flux = flux[iuse]
use_err = 1./sqrt(ivar[iuse])

jhusdss_lowz_doublet_line_fit, slope, intercept, center, separation, lflux, ratio, sigma, $
        err_slope, err_intercept, err_center, err_separation, err_lflux, err_ratio, err_sigma

;; Monte Carlo?

end
