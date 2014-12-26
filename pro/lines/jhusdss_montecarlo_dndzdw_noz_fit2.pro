pro dndzdw_noz_fit2, n_star_strong, w_star_strong, $
                n_star_strong_err, w_star_strong_err, $
                n_star_weak, w_star_weak, $
                n_star_weak_err, w_star_weak_err
common com_noz_dfit, use_xw, use_phi, use_err, init_n_star, init_w_star

start = dblarr(4)
start[0] = init_n_star[0]
start[1] = init_w_star[0]
start[2] = init_n_star[1]
start[3] = init_w_star[1]

parinfo = replicate({limited:bytarr(2), limits:dblarr(2)}, 4)
;; strong system
;; this is actually n_star/w_star
parinfo[0].limited = [1b, 1b]
parinfo[0].limits = [0.01, 10.0]

;; this is w_star > 0.2
parinfo[1].limited = [1b, 1b]
parinfo[1].limits = [0.4, 10.0]

;; weak system
;; this is actually n_star/w_star
parinfo[2].limited = [1b, 1b]
parinfo[2].limits = [0.01, 10.0]

;; this is w_star < 0.1
parinfo[3].limited = [1b, 1b]
parinfo[3].limits = [0.0001, 0.3]

p=mpfitfun('jhusdss_dndzdw_noz_func2', use_xw, use_phi, use_err, $
           start, ftol=1.d-10, bestnorm=chi2, parinfo=parinfo, perror=perror, /double, /quiet)

n_star_strong = p[0]
w_star_strong = p[1]
n_star_weak = p[2]
w_star_weak = p[3]

n_star_strong_err = perror[0]
w_star_strong_err = perror[1]
n_star_weak_err = perror[2]
w_star_weak_err = perror[3]

end

function jhusdss_montecarlo_dndzdw_noz_fit2, inphi
common com_noz_dfit, use_xw, use_phi, use_err, init_n_star, init_w_star

z_nbin = n_elements(inphi.z)
n_star_strong = fltarr(z_nbin)
w_star_strong = fltarr(z_nbin)
n_star_strong_err = fltarr(z_nbin)
w_star_strong_err = fltarr(z_nbin)
n_star_weak = fltarr(z_nbin)
w_star_weak = fltarr(z_nbin)
n_star_weak_err = fltarr(z_nbin)
w_star_weak_err = fltarr(z_nbin)

for i=0L, z_nbin-1L do begin
    xw = inphi.median_w[*,i]
;   xw = 10.^inphi.median_log10w[*,i]
    phi = inphi.phi[*,i]
    err = inphi.phi_poisson_err[*,i]

    ;; only use those within a certain range
    ii = where(xw gt 0.20 and xw lt 6. and phi gt 0., nn)
    if (nn gt 0) then begin
       use_xw = xw[ii]
       use_phi = phi[ii]
       use_err = err[ii]

       ;; initial guess using Nestor+2005
       init_n_star = [0.9, 0.7]
       init_w_star = [0.5, 0.05]

       dndzdw_noz_fit2, n_star_strong_tmp, w_star_strong_tmp, $
                n_star_strong_err_tmp, w_star_strong_err_tmp, $
                n_star_weak_tmp, w_star_weak_tmp, $
                n_star_weak_err_tmp, w_star_weak_err_tmp
       n_star_strong[i] = n_star_strong_tmp
       w_star_strong[i] = w_star_strong_tmp
       n_star_strong_err[i] = n_star_strong_err_tmp
       w_star_strong_err[i] = w_star_strong_err_tmp
       n_star_weak[i] = n_star_weak_tmp
       w_star_weak[i] = w_star_weak_tmp
       n_star_weak_err[i] = n_star_weak_err_tmp
       w_star_weak_err[i] = w_star_weak_err_tmp
    endif
endfor
outstr = {n_star_strong:n_star_strong*w_star_strong, w_star_strong:w_star_strong, $
         n_star_strong_err:sqrt((n_star_strong*w_star_strong_err)^2+(w_star_strong_err*n_star_strong)^2), $
         w_star_strong_err:w_star_strong_err, $

         n_star_weak:n_star_weak*w_star_weak, w_star_weak:w_star_weak, $
         n_star_weak_err:sqrt((n_star_weak*w_star_weak_err)^2+(w_star_weak_err*n_star_weak)^2), $
         w_star_weak_err:w_star_weak_err}

return, outstr
end
