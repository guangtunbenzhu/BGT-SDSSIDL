pro dndzdw_noz_fit, n_star, w_star, $
                n_star_err, w_star_err
common com_noz_dfit, use_xw, use_phi, use_err, init_n_star, init_w_star

start = dblarr(2)
start[0] = init_n_star
start[1] = init_w_star

parinfo = replicate({limited:bytarr(2), limits:dblarr(2)}, 2)
;; strong system
;; this is actually n_star/w_star
parinfo[0].limited = [1b, 1b]
parinfo[0].limits = [0.01, 10.0]

;; this is w_star > 0.2
parinfo[1].limited = [1b, 1b]
parinfo[1].limits = [0.01, 10.0]

p=mpfitfun('jhusdss_dndzdw_noz_func', use_xw, use_phi, use_err, $
           start, ftol=1.d-10, bestnorm=chi2, parinfo=parinfo, perror=perror, /double, /quiet)

n_star = p[0]
w_star = p[1]

n_star_err = perror[0]
w_star_err = perror[1]

end

function jhusdss_montecarlo_dndzdw_noz_fit, inphi, wmin=wmin, wmax=wmax

common com_noz_dfit, use_xw, use_phi, use_err, init_n_star, init_w_star

if (n_elements(wmin) eq 0) then wmin=0.60
if (n_elements(wmax) eq 0) then wmax=5.00

z_nbin = n_elements(inphi.z)
n_star = fltarr(z_nbin)
w_star = fltarr(z_nbin)
n_star_err = fltarr(z_nbin)
w_star_err = fltarr(z_nbin)

for i=0L, z_nbin-1L do begin
;   xw = 10.^inphi.median_log10w[*,i]
    xw = inphi.median_w[*,i]
    phi = inphi.phi[*,i]
    err = inphi.phi_poisson_err[*,i]

    ;; only use those within a certain range
    ii = where(xw gt wmin and xw lt wmax and phi gt 0., nn)
    if (nn gt 0) then begin
       use_xw = xw[ii]
       use_phi = phi[ii]
       use_err = err[ii]

       ;; initial guess using Nestor+2005
       init_n_star = 0.9
       init_w_star = 0.5

       dndzdw_noz_fit, n_star_tmp, w_star_tmp, $
                n_star_err_tmp, w_star_err_tmp
       n_star[i] = n_star_tmp
       w_star[i] = w_star_tmp
       n_star_err[i] = n_star_err_tmp
       w_star_err[i] = w_star_err_tmp
    endif
endfor
outstr = {n_star:n_star*w_star, w_star:w_star, $
         n_star_err:sqrt((n_star*w_star_err)^2+(w_star_err*n_star)^2), $
         w_star_err:w_star_err}

return, outstr
end
