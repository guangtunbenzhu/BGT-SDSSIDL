;+
; Calculate: 
; d2N/dlog10Wdz = \bar Delta(N) / \bar Delta(Z)
; \bar Delta(N) = \Sigma (1/f_i) / N(QSO)
; \bar Delta(Z) = \Sigma (Ncovered)*dz(MonteCarlo pixel)/dn(MonteCaro pixel) / N(QSO)
; N(QSO) cancels out
; See jhusdss_absorber_completeness.pro
; 
;-
pro jhusdss_montecarlo_dndzdw, nmfver, boss=boss, overwrite=overwrite

if (n_elements(nmfver) eq 0) then begin
   splog, 'nmfver required.'
   return
endif
;; path
if (keyword_set(boss)) then begin
   path=jhusdss_get_path(/nmfqso)+'/'+$
        string(nmfver, format='(I3.3)')+'/MonteCarlo_BOSS'
endif else begin
   path=jhusdss_get_path(/nmfqso)+'/'+$
        string(nmfver, format='(I3.3)')+'/MonteCarlo'
endelse

filename = jhusdss_montecarlo_dndzdw_filename(nmfver)
outfile = path+'/'+filename
if (file_test(outfile) and (not keyword_set(overwrite))) then begin
    splog, outfile+' file exists, not overwriting ...'
    return
endif

;; dz for a MonteCarlo pixel, the same as in convolving 
dz_pixel_mc = 0.0005
;; number of absorbers for each MC pixel
dn_pixel_mc = 1

;; z bins
delta_z = 0.15
z_min = 0.4
z_max = 2.30
z_bin = jhusdss_make_bins(z_min, z_max, delta_z, nbin=z_nbin)
z_bin[z_nbin-1].max = z_max

;; log10(w) bins
delta_log10w = alog10(1.5)/3.
;; ew(2796) bins
log10w_min = alog10(0.20)
;log10w_min = alog10(0.3)-alog10(1.5)/2.
log10w_max = alog10(15.0)
log10w_bin = jhusdss_make_bins(log10w_min, log10w_max, delta_log10w, nbin=log10w_nbin)

;; read in absorbers
absorbers = jhusdss_absorber_readin(nmfver, /byabs, /trim)

;; read in monte carlo completeness
mc = jhusdss_montecarlo_readin(nmfver)
;; post-processing mc
f_log10w_z = float(mc.ndetected)/float((mc.ncovered+(mc.ncovered eq 0)))
nlog10w = n_elements(mc.log10w)

;; for strong absorbers
istrong = where(mc.log10w_min gt alog10(4.), nstrong)
ftmp = total(mc.ndetected[istrong,*], 1)/total(mc.ncovered[istrong,*], 1)
for i=0L, nstrong-1L do f_log10w_z[istrong[i], *] = ftmp

;; for weak absorbers we extrapolate
;iweak = where(mc.log10w_max le alog10(0.25), nweak)
;ftmp = total(mc.ndetected[iweak,*], 1)/total(mc.ncovered[iweak,*], 1)
;for i=0L, nweak-1L do f_log10w_z[iweak[i], *] = ftmp

;; get f for each absorber
log10w_abs = alog10(absorbers.rew_mgii_2796)
z_abs = absorbers.zabs
;; bilinear interpolation
;; virtual subscript
nabs = n_elements(absorbers)
v_log10w_abs = (log10w_abs-min(mc.log10w_min))/mc.dlog10w
v_z_abs = (z_abs-min(mc.z_min))/mc.dz
f_abs = fltarr(nabs)
for i=0L, nabs-1L do f_abs[i] = bilinear(f_log10w_z, v_log10w_abs[i], v_z_abs[i]) 

;; output
outstr = {log10w:log10w_bin.mean, log10w_min:log10w_bin.min, log10w_max:log10w_bin.max, $
          z:z_bin.mean, z_min:z_bin.min, z_max:z_bin.max, $
          zbin_median:fltarr(z_nbin), log10wbin_median:fltarr(log10w_nbin), $
          zbin_mean:fltarr(z_nbin), log10wbin_mean:fltarr(log10w_nbin), $
          zbin_sdev:fltarr(z_nbin), log10wbin_sdev:fltarr(log10w_nbin), $
          median_log10w:fltarr(log10w_nbin, z_nbin), median_z:fltarr(log10w_nbin, z_nbin), $
          mean_log10w:fltarr(log10w_nbin, z_nbin), mean_z:fltarr(log10w_nbin, z_nbin), $
          sdev_log10w:fltarr(log10w_nbin, z_nbin), sdev_z:fltarr(log10w_nbin, z_nbin), $

          nabs:lonarr(log10w_nbin, z_nbin), $
          phi:fltarr(log10w_nbin, z_nbin), $
          phi_poisson_err:fltarr(log10w_nbin, z_nbin), $
          phi_upper_limit:fltarr(log10w_nbin, z_nbin), $
          phi_lower_limit:fltarr(log10w_nbin, z_nbin), $

          ;; cumulative
          cum_median_log10w:fltarr(log10w_nbin, z_nbin), cum_median_z:fltarr(log10w_nbin, z_nbin), $
          cum_nabs:lonarr(log10w_nbin, z_nbin), $
          cum_phi:fltarr(log10w_nbin, z_nbin), $
          cum_phi_poisson_err:fltarr(log10w_nbin, z_nbin), $
          cum_phi_upper_limit:fltarr(log10w_nbin, z_nbin), $
          cum_phi_lower_limit:fltarr(log10w_nbin, z_nbin), $

          ;; fitting of N*(1+z)^alpha/W*(1+z)^beta*exp(-W/W*(1+z)^beta)
          ;; 0.2<w0<5.0
          ;; fit2
          n_star_strong_all:0., w_star_strong_all:0., alpha_strong_all:0., bbeta_strong_all:0., $
          n_star_strong_all_err:0., w_star_strong_all_err:0., alpha_strong_all_err:0., bbeta_strong_all_err:0., $
          n_star_weak_all:0., w_star_weak_all:0., alpha_weak_all:0., bbeta_weak_all:0., $
          n_star_weak_all_err:0., w_star_weak_all_err:0., alpha_weak_all_err:0., bbeta_weak_all_err:0., $

          ;; noz_fit2
          n_star_strong:fltarr(z_nbin), w_star_strong:fltarr(z_nbin), $
          n_star_strong_err:fltarr(z_nbin), w_star_strong_err:fltarr(z_nbin), $
          n_star_weak:fltarr(z_nbin), w_star_weak:fltarr(z_nbin), $
          n_star_weak_err:fltarr(z_nbin), w_star_weak_err:fltarr(z_nbin),  $

          ;; 0.6<w0<5.0
          ;; noz_fit
          n_star:fltarr(z_nbin), w_star:fltarr(z_nbin), $
          n_star_err:fltarr(z_nbin), w_star_err:fltarr(z_nbin), $

          ;; fit3
          ;; 0.6<w0<5.0, with (1+z)/(1+(z/bbeta)^ggamma) form
          f0:0., f0_alpha:0., f0_bbeta:0., f0_ggamma:0., $
          f0_err:0., f0_alpha_err:0., f0_bbeta_err:0., f0_ggamma_err:0., $
          w0:0., w0_alpha:0., w0_bbeta:0., w0_ggamma:0., $
          w0_err:0., w0_alpha_err:0., w0_bbeta_err:0., w0_ggamma_err:0. $
          }

;; LOOP OVER ALL BINS
for jz=0L, z_nbin-1L do begin

    ;; here, redshift bin had best be an integer times mc.
    jj = where(mc.z_min ge z_bin[jz].min and mc.z_max le z_bin[jz].max, mm)
    if (mm eq 0) then message, 'fix this first!'
    bar_delta_z = total(mc.ncovered_all[jj])*dz_pixel_mc/dn_pixel_mc

    this_z_min = min(mc.z_min[jj])
    this_z_max = max(mc.z_max[jj])
   
    tt = where(z_abs gt this_z_min and z_abs lt this_z_max, ss)
    if (ss gt 0) then begin
       outstr.zbin_median[jz] = median(z_abs[tt])
       outstr.zbin_mean[jz] = mean(z_abs[tt])
       outstr.zbin_sdev[jz] = sqrt((moment(z_abs[tt]))[1])
    endif

    for iw=0L, log10w_nbin-1L do begin
        this_log10w_min = log10w_bin[iw].min
        this_log10w_max = log10w_bin[iw].max

        if (jz eq 0) then begin
            tt = where(log10w_abs gt this_log10w_min and log10w_abs lt this_log10w_max, ss)
            if (ss gt 0) then begin
               outstr.log10wbin_median[iw] = median(log10w_abs[tt])
               outstr.log10wbin_mean[iw] = mean(log10w_abs[tt])
               outstr.log10wbin_sdev[iw] = sqrt((moment(log10w_abs[tt]))[1])
            endif
        endif

        ii = where(log10w_abs gt this_log10w_min $
               and log10w_abs le this_log10w_max $
               and z_abs gt this_z_min $
               and z_abs le this_z_max, nn)
        if (nn gt 0) then begin
           outstr.nabs[iw, jz] = nn
           outstr.median_log10w[iw, jz] = median(log10w_abs[ii])
           outstr.mean_log10w[iw, jz] = mean(log10w_abs[ii])
           outstr.sdev_log10w[iw, jz] = sqrt((moment(log10w_abs[ii]))[1])
           outstr.median_z[iw, jz] = median(z_abs[ii])
           outstr.mean_z[iw, jz] = mean(z_abs[ii])
           outstr.sdev_z[iw, jz] = sqrt((moment(z_abs[ii]))[1])

           phi = total(1./(f_abs[ii] + (f_abs[ii] le 0.)))
           phi_poisson_err2 = total(1./(f_abs[ii] + (f_abs[ii] le 0.))^2)
           weff = phi_poisson_err2/phi
           neff = phi/weff
           gehrels_err = im_poisson_limits(neff, 0.8413)
           phi_lower_limit = gehrels_err[0]*weff
           phi_upper_limit = gehrels_err[1]*weff
           outstr.phi[iw, jz] = phi/bar_delta_z/delta_log10w/alog(10.)/10.^outstr.median_log10w[iw, jz]
           outstr.phi_poisson_err[iw, jz] = sqrt(phi_poisson_err2)/bar_delta_z/delta_log10w/alog(10.)/10.^outstr.median_log10w[iw, jz]
           outstr.phi_upper_limit[iw, jz] = phi_upper_limit/bar_delta_z/delta_log10w/alog(10.)/10.^outstr.median_log10w[iw, jz]
           outstr.phi_lower_limit[iw, jz] = phi_lower_limit/bar_delta_z/delta_log10w/alog(10.)/10.^outstr.median_log10w[iw, jz]
        endif

        ;; cumulative phi
        kk = where(log10w_abs gt this_log10w_min $
               and z_abs gt this_z_min $
               and z_abs le this_z_max, ll)
        if (ll gt 0) then begin
           outstr.cum_nabs[iw, jz] = nn
           outstr.cum_median_log10w[iw, jz] = median(log10w_abs[kk])
           outstr.cum_median_z[iw, jz] = median(z_abs[kk])
           phi = total(1./(f_abs[kk] + (f_abs[kk] le 0.)))
           phi_poisson_err2 = total(1./(f_abs[kk] + (f_abs[kk] le 0.))^2)
           weff = phi_poisson_err2/phi
           neff = phi/weff
           gehrels_err = im_poisson_limits(neff, 0.8413)
           phi_lower_limit = gehrels_err[0]*weff
           phi_upper_limit = gehrels_err[1]*weff
           outstr.cum_phi[iw, jz] = phi/bar_delta_z
           outstr.cum_phi_poisson_err[iw, jz] = sqrt(phi_poisson_err2)/bar_delta_z
           outstr.cum_phi_upper_limit[iw, jz] = phi_upper_limit/bar_delta_z
           outstr.cum_phi_lower_limit[iw, jz] = phi_lower_limit/bar_delta_z
        endif

    endfor
endfor

;; fit
;; 0.2<w0<5.0
strtmp = jhusdss_montecarlo_dndzdw_fit2(outstr)
;; 0.2<w0<5.0
strtmp1 = jhusdss_montecarlo_dndzdw_noz_fit2(outstr)
;; 0.6<w0<5.0
str_noweak = jhusdss_montecarlo_dndzdw_noz_fit(outstr)
;; 0.6<w0<5.0
str_newz = jhusdss_montecarlo_dndzdw_fit3(outstr)

;;
;print, strtmp.n_star_strong, strtmp.w_star_strong, strtmp.alpha_strong, strtmp.bbeta_strong
;print, strtmp.n_star_strong_err, strtmp.w_star_strong_err, strtmp.alpha_strong_err, strtmp.bbeta_strong_err
;print, strtmp.n_star_weak, strtmp.w_star_weak, strtmp.alpha_weak, strtmp.bbeta_weak
;print, strtmp.n_star_weak_err, strtmp.w_star_weak_err, strtmp.alpha_weak_err, strtmp.bbeta_weak_err
;for i=0, z_nbin-1 do begin
;djs_plot, 10.^outstr.median_log10w[*,i], outstr.phi[*,i], psym=4, /ylog, xst=1, yst=1, xra=[0,6], yra=[1e-4, 10]
;ii = where(outstr.phi[*,i] gt 0.)
;xw = 10.^outstr.median_log10w[ii, i]
;oploterror, xw, outstr.phi[ii, i], outstr.phi_poisson_err[ii,i] 
;xz = outstr.median_z[ii, i]
;p = [strtmp.n_star_strong/strtmp.w_star_strong, strtmp.w_star_strong, $
;     strtmp.alpha_strong-strtmp.bbeta_strong, strtmp.bbeta_strong, $
;     strtmp.n_star_weak/strtmp.w_star_weak, strtmp.w_star_weak, $
;     strtmp.alpha_weak-strtmp.bbeta_weak, strtmp.bbeta_weak]
;model = jhusdss_dndzdw_func2(xw, xz, p)
;djs_oplot, xw, model, psym=4, color='red'
;
;p1 = [strtmp1.n_star_strong[i]/strtmp1.w_star_strong[i], strtmp1.w_star_strong[i], $
;     strtmp1.n_star_weak[i]/strtmp1.w_star_weak[i], strtmp1.w_star_weak[i]]
;model1 = jhusdss_dndzdw_noz_func2(xw, p1)
;djs_oplot, xw, model1, color='green'
;stop
;endfor
;
;djs_plot, outstr.zbin_mean, strtmp1.w_star_strong, yra=[0.2, 1]
;oploterror, outstr.zbin_mean, strtmp1.w_star_strong, strtmp1.w_star_strong_err
;djs_oplot, outstr.zbin_mean, strtmp.w_star_strong*(1.+outstr.zbin_mean)^strtmp.bbeta_strong
;
outstr.n_star_strong_all = strtmp.n_star_strong
outstr.w_star_strong_all = strtmp.w_star_strong
outstr.n_star_strong_all_err = strtmp.n_star_strong_err
outstr.w_star_strong_all_err = strtmp.w_star_strong_err
outstr.alpha_strong_all = strtmp.alpha_strong
outstr.bbeta_strong_all = strtmp.bbeta_strong
outstr.alpha_strong_all_err = strtmp.alpha_strong_err
outstr.bbeta_strong_all_err = strtmp.bbeta_strong_err

outstr.n_star_weak_all = strtmp.n_star_weak
outstr.w_star_weak_all = strtmp.w_star_weak
outstr.n_star_weak_all_err = strtmp.n_star_weak_err
outstr.w_star_weak_all_err = strtmp.w_star_weak_err
outstr.alpha_weak_all = strtmp.alpha_weak
outstr.bbeta_weak_all = strtmp.bbeta_weak
outstr.alpha_weak_all_err = strtmp.alpha_weak_err
outstr.bbeta_weak_all_err = strtmp.bbeta_weak_err


outstr.n_star_strong = strtmp1.n_star_strong
outstr.w_star_strong = strtmp1.w_star_strong
outstr.n_star_strong_err = strtmp1.n_star_strong_err
outstr.w_star_strong_err = strtmp1.w_star_strong_err
outstr.n_star_weak = strtmp1.n_star_weak
outstr.w_star_weak = strtmp1.w_star_weak
outstr.n_star_weak_err = strtmp1.n_star_weak_err
outstr.w_star_weak_err = strtmp1.w_star_weak_err

outstr.n_star = str_noweak.n_star
outstr.w_star = str_noweak.w_star
outstr.n_star_err = str_noweak.n_star_err
outstr.w_star_err = str_noweak.w_star_err

outstr.f0 = str_newz.f0
outstr.f0_alpha = str_newz.f0_alpha
outstr.f0_bbeta = str_newz.f0_bbeta
outstr.f0_ggamma = str_newz.f0_ggamma
outstr.f0_err = str_newz.f0_err
outstr.f0_alpha_err= str_newz.f0_alpha_err
outstr.f0_bbeta_err = str_newz.f0_bbeta_err
outstr.f0_ggamma_err = str_newz.f0_ggamma_err
outstr.w0 = str_newz.w0
outstr.w0_alpha = str_newz.w0_alpha
outstr.w0_bbeta = str_newz.w0_bbeta
outstr.w0_ggamma = str_newz.w0_ggamma
outstr.w0_err = str_newz.w0_err
outstr.w0_alpha_err= str_newz.w0_alpha_err
outstr.w0_bbeta_err = str_newz.w0_bbeta_err
outstr.w0_ggamma_err = str_newz.w0_ggamma_err

mwrfits, outstr, outfile, /create

stop
end
