;function dndzdw_func2, xw, xz, p
;n_star_strong = p[0]
;w_star_strong = p[1]
;alpha_strong = p[2]
;bbeta_strong = p[3]
;n_star_weak = p[4]
;w_star_weak = p[5]
;alpha_weak = p[6]
;bbeta_weak = p[7]
;
;model = n_star_strong*(1.+xz)^alpha_strong*exp(-xw/w_star_strong/(1.+xz)^bbeta_strong) $
;      + n_star_weak*(1.+xz)^alpha_weak*exp(-xw/w_star_weak/(1.+xz)^bbeta_weak)
;return, model
;end

pro dndzdw_fit2, n_star_strong, w_star_strong, alpha_strong, bbeta_strong, $
                n_star_strong_err, w_star_strong_err, alpha_strong_err, bbeta_strong_err, $
                n_star_weak, w_star_weak, alpha_weak, bbeta_weak, $
                n_star_weak_err, w_star_weak_err, alpha_weak_err, bbeta_weak_err
common com_dfit, use_xw, use_xz, use_phi, use_err, init_n_star, init_w_star, init_alpha, init_bbeta 

start = dblarr(8)
start[0] = init_n_star[0]
start[1] = init_w_star[0]
start[2] = init_alpha[0]
start[3] = init_bbeta[0]
start[4] = init_n_star[1]
start[5] = init_w_star[1]
start[6] = init_alpha[1]
start[7] = init_bbeta[1]

parinfo = replicate({limited:bytarr(2), limits:dblarr(2)}, 8)
;; strong system
;; this is actually n_star/w_star
parinfo[0].limited = [1b, 1b]
parinfo[0].limits = [0.01, 10.0]

;; this is w_star > 0.5
parinfo[1].limited = [1b, 1b]
parinfo[1].limits = [0.3, 10.0]

;; this is actually alpha-bbeta
parinfo[2].limited = [1b, 1b]
parinfo[2].limits = [-10.0, 10.0]

;; this is bbeta
parinfo[3].limited = [1b, 1b]
parinfo[3].limits = [-10.0, 10.0]

;; weak system
;; this is actually n_star/w_star
parinfo[4].limited = [1b, 1b]
parinfo[4].limits = [0.01, 10.0]

;; this is w_star < 0.1
parinfo[5].limited = [1b, 1b]
parinfo[5].limits = [0.001, 0.2]

;; this is actually alpha-bbeta
parinfo[6].limited = [1b, 1b]
parinfo[6].limits = [-10.0, 10.0]

;; this is bbeta
parinfo[7].limited = [1b, 1b]
parinfo[7].limits = [-10.0, 10.0]


p=mpfit2dfun('jhusdss_dndzdw_func2', use_xw, use_xz, use_phi, use_err, $
           start, ftol=1.d-10, bestnorm=chi2, parinfo=parinfo, perror=perror, /double, /quiet)

n_star_strong = p[0]
w_star_strong = p[1]
alpha_strong = p[2]
bbeta_strong = p[3]
n_star_weak = p[4]
w_star_weak = p[5]
alpha_weak = p[6]
bbeta_weak = p[7]

n_star_strong_err = perror[0]
w_star_strong_err = perror[1]
alpha_strong_err = perror[2]
bbeta_strong_err = perror[3]
n_star_weak_err = perror[4]
w_star_weak_err = perror[5]
alpha_weak_err = perror[6]
bbeta_weak_err = perror[7]

end

function jhusdss_montecarlo_dndzdw_fit2, inphi, zmax=zmax, wmin=wmin, wmax=wmax
common com_dfit, use_xw, use_xz, use_phi, use_err, init_n_star, init_w_star, init_alpha, init_bbeta 

   if (n_elements(zmax) eq 0) then zmin=0.425
   if (n_elements(zmax) eq 0) then zmax=1.65
   if (n_elements(wmin) eq 0) then wmin=0.20
   if (n_elements(wmax) eq 0) then wmax=5.00
;  xw = 10.^inphi.median_log10w
   xw = inphi.median_w
   xz = inphi.median_z
   phi = inphi.phi
   err = inphi.phi_poisson_err

   ;; only use those within a certain range
   ii = where(xw gt wmin and xw lt wmax and phi gt 0. and xz gt zmin and xz le zmax, nn)
   use_xw = xw[ii]
   use_xz = xz[ii]
   use_phi = phi[ii]
   use_err = err[ii]
;  stop

   ;; initial guess using Nestor+2005
   init_n_star = [0.9, 0.7]
   init_w_star = [0.5, 0.05]
   init_alpha = [-0.3, -0.3]
   init_bbeta = [0.5, 0.5]

   dndzdw_fit2, n_star_strong, w_star_strong, alpha_strong, bbeta_strong, $
                n_star_strong_err, w_star_strong_err, alpha_strong_err, bbeta_strong_err, $
                n_star_weak, w_star_weak, alpha_weak, bbeta_weak, $
                n_star_weak_err, w_star_weak_err, alpha_weak_err, bbeta_weak_err
   
   outstr = {n_star_strong:n_star_strong*w_star_strong, w_star_strong:w_star_strong, $
             alpha_strong:alpha_strong+bbeta_strong, bbeta_strong:bbeta_strong, $
             n_star_strong_err:sqrt((n_star_strong*w_star_strong_err)^2+(w_star_strong_err*n_star_strong)^2), $
             w_star_strong_err:w_star_strong_err, $
             alpha_strong_err:sqrt(alpha_strong_err^2+bbeta_strong_err^2), $
             bbeta_strong_err:bbeta_strong_err, $

             n_star_weak:n_star_weak*w_star_weak, w_star_weak:w_star_weak, $
             alpha_weak:alpha_weak+bbeta_weak, bbeta_weak:bbeta_weak, $
             n_star_weak_err:sqrt((n_star_weak*w_star_weak_err)^2+(w_star_weak_err*n_star_weak)^2), $
             w_star_weak_err:w_star_weak_err, $
             alpha_weak_err:sqrt(alpha_weak_err^2+bbeta_weak_err^2), $
             bbeta_weak_err:bbeta_weak_err}

return, outstr
end
