function powerlaw_func, x, p
slope = p[0]
intercept = p[1]
;model = fltarr(n_elements(x))
model = intercept*x^slope
return, model
end

pro powerlaw_fit, slope, intercept, err_slope, err_intercept, fixslope=fixslope
common com_powerlaw, init_slope, init_intercept, use_rp, use_ew, use_err

;print, use_rp, use_ew, use_err
start = dblarr(2)
start[0]=init_slope
start[1]=init_intercept

parinfo=replicate({limited:bytarr(2), limits:dblarr(2), fixed:0b}, 2)

;; Slope
parinfo[0].limited=[1b, 1b]
parinfo[0].limits=[-3, 1]
;parinfo[0].limits=[-1.385, -1.375]
if (keyword_set(fixslope)) then parinfo[0].fixed=1b

;; Intercept
parinfo[1].limited=[0b, 0b]
parinfo[1].limits=[-1, 2]

p=mpfitfun('powerlaw_func', use_rp, use_ew, use_err, $
           start, ftol=1.d-10, bestnorm=chi2, parinfo=parinfo, perror=perror, /double, /quiet)

print, 'reduced chi-square: ', total(((use_ew-powerlaw_func(use_rp,p)))^2/use_err^2)/(n_elements(use_rp)-2.-1.)

slope = p[0]
intercept = p[1]

err_slope = perror[0]
err_intercept = perror[1]
;print, p, perror

end

;; y = b*x^a
pro jhusdss_powerlaw_fit, in_rp, in_ew, in_ew_err, in_slope=in_slope, in_intercept=in_intercept, $
            slope=slope, intercept=intercept, err_slope=err_slope, err_intercept=err_intercept, fixslope=fixslope
common com_powerlaw, init_slope, init_intercept, use_rp, use_ew, use_err

init_slope = in_slope
init_intercept = in_intercept
use_rp = in_rp
use_ew = in_ew 
use_err = in_ew_err

powerlaw_fit, slope, intercept, err_slope, err_intercept, fixslope=fixslope

end
