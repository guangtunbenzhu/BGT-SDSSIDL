;+
; See 
;-
pro jhusdss_lowz_errf_do_fit, center, sigma, err_center, err_sigma

common com_errf_fit, init_center, init_sigma, use_rew, use_dr, use_err, use_maxwidth, sigma_range

start=dblarr(2)

start[0]=init_center
start[1]=init_sigma

parinfo=replicate({limited:bytarr(2), limits:dblarr(2), fixed:0b}, 2)

parinfo[0].limited=[1b, 1b]
parinfo[0].limits=init_center+[-use_maxwidth, use_maxwidth]
parinfo[0].fixed=0b

parinfo[1].limited=[1b, 1b]
parinfo[1].limits=sigma_range
parinfo[1].fixed=1b

p=mpfitfun('jhusdss_lowz_errf_func', use_rew, use_dr, use_err, $
           start, ftol=1.d-10, bestnorm=chi2, parinfo=parinfo, perror=perror, /double, /quiet)
;print, chi2

center = p[0]
sigma = p[1]

err_center = perror[0]
err_sigma = perror[1]

end

pro jhusdss_lowz_errf_fit, in_rew, in_dr, in_err, in_center, in_sigma,$
       center=center, sigma=sigma, err_center=err_center, err_sigma=err_sigma, $
       maxwidth=maxwidth

common com_errf_fit

if n_elements(maxwidth) eq 0 then $
   maxwidth = 1.

use_maxwidth = maxwidth
sigma_range = [1.00D-4, 1.D0]

init_center = alog10(in_center)
init_sigma = in_sigma

use_rew = alog10(in_rew)
use_dr = in_dr
use_err = in_err

jhusdss_lowz_errf_do_fit, center, sigma, err_center, err_sigma

end
