pro dndzdw_fit3, f0, f0_alpha, f0_bbeta, f0_ggamma, $
                 f0_err, f0_alpha_err, f0_bbeta_err, f0_ggamma_err, $
                 w0, w0_alpha, w0_bbeta, w0_ggamma, $
                 w0_err, w0_alpha_err, w0_bbeta_err, w0_ggamma_err
common com_dfit3, use_xw, use_xz, use_phi, use_err, init_f0, init_f0_alpha, init_f0_bbeta, init_f0_ggamma, $
                 init_w0, init_w0_alpha, init_w0_bbeta, init_w0_ggamma

start = dblarr(8)
start[0] = init_f0
start[1] = init_f0_alpha
start[2] = init_f0_bbeta
start[3] = init_f0_ggamma
start[4] = init_w0
start[5] = init_w0_alpha
start[6] = init_w0_bbeta
start[7] = init_w0_ggamma

parinfo = replicate({limited:bytarr(2), limits:dblarr(2)}, 8)
;; strong system
parinfo[0].limited = [1b, 1b]
parinfo[0].limits = [0.001, 100.0]

parinfo[1].limited = [1b, 1b]
parinfo[1].limits = [0.001, 100.0]

parinfo[2].limited = [1b, 1b]
parinfo[2].limits = [0.001, 100.0]

parinfo[3].limited = [1b, 1b]
parinfo[3].limits = [0.001, 100.0]

parinfo[4].limited = [1b, 1b]
parinfo[4].limits = [0.001, 100.0]

parinfo[5].limited = [1b, 1b]
parinfo[5].limits = [0.001, 100.0]

parinfo[6].limited = [1b, 1b]
parinfo[6].limits = [0.001, 100.0]

parinfo[7].limited = [1b, 1b]
parinfo[7].limits = [0.001, 100.0]

p=mpfit2dfun('jhusdss_dndzdw_func3', use_xw, use_xz, use_phi, use_err, $
           start, ftol=1.d-10, bestnorm=chi2, parinfo=parinfo, perror=perror, /double, /quiet)

help, use_xw
print, chi2
f0 = p[0]
f0_alpha = p[1]
f0_bbeta = p[2]
f0_ggamma = p[3]
w0 = p[4]
w0_alpha = p[5]
w0_bbeta = p[6]
w0_ggamma = p[7]

f0_err = perror[0]
f0_alpha_err = perror[1]
f0_bbeta_err = perror[2]
f0_ggamma_err = perror[3]
w0_err = perror[4]
w0_alpha_err = perror[5]
w0_bbeta_err = perror[6]
w0_ggamma_err = perror[7]

end

function jhusdss_montecarlo_dndzdw_fit3, inphi, zmax=zmax, wmin=wmin, wmax=wmax, highz=highz
common com_dfit3, use_xw, use_xz, use_phi, use_err, init_f0, init_f0_alpha, init_f0_bbeta, init_f0_ggamma, $
                 init_w0, init_w0_alpha, init_w0_bbeta, init_w0_ggamma

   if (n_elements(zmin) eq 0) then zmin=0.425
   if (n_elements(zmax) eq 0) then zmax=3.00
   if (n_elements(wmin) eq 0) then wmin=0.6
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

   ;; read in matejek, temporary
;  if (keyword_set(highz)) then begin
     readcol, jhusdss_get_path(/nmfqso)+'/matejek.txt', ms_z, ms_w, ms_phi, ms_err
     use_xw = [use_xw, ms_w]
     use_xz = [use_xz, ms_z]
     use_phi = [use_phi, ms_phi]
     use_err = [use_err, ms_err]
;  endif

   ;; initial guess using Nestor+2005
   init_f0 = 1.
   init_f0_alpha = 1.
   init_f0_bbeta = 1.6
   init_f0_ggamma  = 5.
   init_w0 = 0.3
   init_w0_alpha = 1.
   init_w0_bbeta = 1.6
   init_w0_ggamma  = 5.


   dndzdw_fit3,  f0, f0_alpha, f0_bbeta, f0_ggamma, $
                 f0_err, f0_alpha_err, f0_bbeta_err, f0_ggamma_err, $
                 w0, w0_alpha, w0_bbeta, w0_ggamma, $
                 w0_err, w0_alpha_err, w0_bbeta_err, w0_ggamma_err
   
   outstr = {f0:f0, f0_alpha:f0_alpha, f0_bbeta:f0_bbeta, f0_ggamma:f0_ggamma, $
             f0_err:f0_err, f0_alpha_err:f0_alpha_err, f0_bbeta_err:f0_bbeta_err, f0_ggamma_err:f0_ggamma_err, $
             w0:w0, w0_alpha:w0_alpha, w0_bbeta:w0_bbeta, w0_ggamma:w0_ggamma, $
             w0_err:w0_err, w0_alpha_err:w0_alpha_err, w0_bbeta_err:w0_bbeta_err, w0_ggamma_err:w0_ggamma_err}

return, outstr
end
