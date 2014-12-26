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

pro doublegaussian_line_fit, center, lflux, sigma, err_center, err_lflux, err_sigma

common com_sfit, use_loglam, use_flux, use_err, initline, initlflux, maxwidth, sigma_range

start=dblarr(5)

start[0]=initline[0]
start[1]=initline[1]

start[2]=initlflux[0]
start[3]=initlflux[1]

;; init dV = 207 km/s
;; dlog10(lambda) = 3D-4 (3 pixels)
start[4] = 3.0D-4

parinfo=replicate({limited:bytarr(2), limits:dblarr(2)}, 5)

;; CIV: 2.6/log(10)/1550 = 7.2D-4 ~7 pixels
;; MgII: 7.0/log(10)/2800 = 10.8D-4 ~11 pixels
;; resolution = 1.D-4
;; use maxwidth = 5D-4

;; center of line1
parinfo[0].limited=1
parinfo[0].limits=initline[0]+[-maxwidth, maxwidth]

;; center of line2
parinfo[1].limited=1
parinfo[1].limits=initline[1]+[-maxwidth, maxwidth]

;; flux of line1
;; flux of line2

;; sigma_range = [69km/s, 500km/s] = [1.D-4, 8.D-4]
parinfo[4].limited=1
parinfo[4].limits=sigma_range

p=mpfitfun('doublegaussian_line_func', use_loglam, use_flux, use_err, $
           start, ftol=1.d-10, bestnorm=chi2, parinfo=parinfo, perror=perror, /double, /quiet)

center = dblarr(2)
lflux = fltarr(2)
err_center = dblarr(2)
err_lflux = fltarr(2)

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

end

pro jhusdss_doublegaussian_fit, in_lambda, in_flux, in_ivar, in_line, in_lflux, $
               ew=ew, err_ew=err_ew, center=center, err_center=err_center, $
               lflux=lflux, err_lflux=err_lflux, sigma=sigma, err_sigma=err_sigma

;common com_sfit, use_loglam, use_flux, use_err, initline, initlflux, maxwidth, sigma_range
common com_sfit

maxwidth = 5.0D-4
sigma_range = [1.0D-4, 8.0D-4]


;; CIV: 2.6/log(10)/1550 = 7.2D-4 ~7 pixels
;; MgII: 7.0/log(10)/2800 = 10.8D-4 ~11 pixels
;; resolution = 1.D-4
;; use maxwidth = 5D-4

initline = alog10(in_line)
initlflux = in_lflux
lambda = in_lambda
flux = in_flux
ivar = in_ivar

;; deal with badness
ibad=where(ivar ne ivar OR flux ne flux or lambda ne lambda, nbad)
if(nbad gt 0) then begin
    ivar[ibad]=0.
    flux[ibad]=0.
    lambda[ibad]=0.
endif

iuse=where(ivar gt 0., nuse)
if(nuse eq 0) then return

use_loglam = alog10(lambda[iuse])
use_flux = flux[iuse]
use_err = 1./sqrt(ivar[iuse])

doublegaussian_line_fit, logcenter, lflux, sigma, err_logcenter, err_lflux, err_sigma

center = 10.^logcenter
err_center = center*alog(10)*err_logcenter

ew = fltarr(2)
err_ew = fltarr(2)

dlambda = jhusdss_dwave(lambda)
tmplambda = min(abs(lambda-center[0]), itmp)
ew[0] = lflux[0]*dlambda[itmp]
err_ew[0] = err_lflux[0]*dlambda[itmp]
tmplambda = min(abs(lambda-center[1]), itmp)
ew[1] = lflux[1]*dlambda[itmp]
err_ew[1] = err_lflux[1]*dlambda[itmp]

;; Monte Carlo?

end
