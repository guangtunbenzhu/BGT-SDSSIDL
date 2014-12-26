function dndzdw_func, xw, xz, p
n_star = p[0]
w_star = p[1]
alpha = p[2]
bbeta = p[3]

model = n_star*(1.+xz)^alpha*exp(-xw/w_star/(1.+xz)^bbeta)
return, model
end

pro dndzdw_fit, n_star, w_star, alpha, bbeta, n_star_err, w_star_err, alpha_err, bbeta_err
common com_dfit, use_xw, use_xz, use_phi, use_err, init_n_star, init_w_star, init_alpha, init_bbeta 

start = dblarr(4)
start[0] = init_n_star
start[1] = init_w_star
start[2] = init_alpha
start[3] = init_bbeta

parinfo = replicate({limited:bytarr(2), limits:dblarr(2)}, 4)
;; this is actually n_star/w_star
parinfo[0].limited = [1b, 1b]
parinfo[0].limits = [0.01, 10.0]

;; this is w_star
parinfo[1].limited = [1b, 1b]
parinfo[1].limits = [0.01, 10.0]

;; this is actually alpha-bbeta
parinfo[2].limited = [1b, 1b]
parinfo[2].limits = [-10.0, 10.0]

;; this is bbeta
parinfo[3].limited = [1b, 1b]
parinfo[3].limits = [-10.0, 10.0]

p=mpfit2dfun('dndzdw_func', use_xw, use_xz, use_phi, use_err, $
           start, ftol=1.d-10, bestnorm=chi2, parinfo=parinfo, perror=perror, /double, /quiet)

n_star = p[0]
w_star = p[1]
alpha = p[2]
bbeta = p[3]

n_star_err = perror[0]
w_star_err = perror[1]
alpha_err = perror[2]
bbeta_err = perror[3]

end

function jhusdss_montecarlo_dndzdw_fit, inphi
common com_dfit, use_xw, use_xz, use_phi, use_err, init_n_star, init_w_star, init_alpha, init_bbeta 

   xw = 10.^inphi.median_log10w
   xz = inphi.median_z
   phi = inphi.phi
   err = inphi.phi_poisson_err

   ;; only use those within a certain range
   ii = where(xw gt 1.0 and xw lt 4. and phi gt 0. and xz le 2.0, nn)
   use_xw = xw[ii]
   use_xz = xz[ii]
   use_phi = phi[ii]
   use_err = err[ii]
;  stop

   ;; initial guess using Nestor+2005
   init_n_star = 1.001/0.443
   init_w_star = 0.443
   init_alpha = 0.226-0.634
   init_bbeta = 0.634

   dndzdw_fit, n_star, w_star, alpha, bbeta, n_star_err, w_star_err, alpha_err, bbeta_err
   
   outstr = {n_star:n_star*w_star, w_star:w_star, alpha:alpha+bbeta, bbeta:bbeta, $
             n_star_err:sqrt((n_star*w_star_err)^2+(w_star_err*n_star)^2), $
             w_star_err:w_star_err, $
             alpha_err:sqrt(alpha_err^2+bbeta_err^2), $
             bbeta_err:bbeta_err}

return, outstr
end
