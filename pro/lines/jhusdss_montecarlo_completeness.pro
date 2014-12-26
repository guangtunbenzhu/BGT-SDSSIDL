;+
; Calculate: 
; f(w,z) = N(detected, w, z)/N(covered, w, z)
; N(covered, z)
; This depends how we define the final catalog
; See jhusdss_absorber_trim.pro, jhusdss_detect_absorber_newengine.pro
;
; Be careful the maximum of a long integer is ~2^31 
; nocal, nocalcium -- not useful at all!
;-
pro jhusdss_montecarlo_completeness, nmfver, boss=boss, overwrite=overwrite, nocal=nocal

;; path
if (keyword_set(boss)) then begin
   path=jhusdss_get_path(/nmfqso)+'/'+$
        string(nmfver, format='(I3.3)')+'/MonteCarlo_BOSS'
endif else begin
   path=jhusdss_get_path(/nmfqso)+'/'+$
        string(nmfver, format='(I3.3)')+'/MonteCarlo'
endelse

filename = jhusdss_montecarlo_completeness_filename(nmfver) 
outfile = path+'/'+filename
if (keyword_set(nocal)) then outfile = path+'/Nocal_'+filename

if (file_test(outfile) and (not keyword_set(overwrite))) then begin
    splog, outfile+' file exists, not overwriting ...'
    return
endif

qsopath=jhusdss_get_path(/nmfqso)+'/'+string(nmfver, format='(I3.3)')+'/Decompose'
if (jhusdss_direxist(qsopath) eq 0) then message, "Can't find the directory."
statfile = qsopath+'/QSO_decompose_NMF_'+string(nmfver, format='(I3.3)')+'stats.fits'
stat0 = mrdfits(statfile, 1)

;; QSO trim
iqso = jhusdss_absorber_trim(stat0, /qso)
stat = stat0[iqso]
nspec = n_elements(stat)

;; See jhusdss_absorber_trim.pro & jhusdss_absorber_stat_window_mask.pro
wave_limit = 1550.
wave_min = 3800.
mgii_limit = 0.04
wave_ciii = 1909.
ca_wave = [3934.7750, 3969.5901]

;; ew(2796) bins
all_ewmin = alog10(0.2)
all_ewmax = alog10(10.0)
ew_binsize = alog10(1.5)/4.
ew_bin = jhusdss_make_bins(all_ewmin, all_ewmax, ew_binsize, nbin=ew_nbin)

;; z(redshift) bins
all_zmin = 0.36
all_zmax = 2.30
z_binsize = 0.0025
z_bin = jhusdss_make_bins(all_zmin, all_zmax, z_binsize, nbin=z_nbin)

;; output
outstr = {log10w:ew_bin.mean, log10w_min:ew_bin.min, log10w_max:ew_bin.max, $
          z:z_bin.mean, z_min:z_bin.min, z_max:z_bin.max, $
          dz:z_binsize, dlog10w:ew_binsize, $
          ndetected:lonarr(ew_nbin, z_nbin), $
          ncovered:lonarr(ew_nbin, z_nbin), $
          ncovered_all:lonarr(z_nbin), $
          nnotconvolved:0L, $
          nall:nspec}

;; LOOP OVER ALL QSOS
for i=0L, nspec-1L do begin
    counter, i+1, nspec
    spec = jhusdss_montecarlo_loadspec(stat[i].plate, stat[i].fiber, nmfver, error=error)
    if (error) then begin 
       splog, "Can't find the spectrum. "
       outstr.nnotconvolved = outstr.nnotconvolved+1L
       continue
    endif

    ;; covered and within the red window
    icovered = where(spec.isitcovered $
                 and spec.zgrid gt wave_limit*(1.+stat[i].zqso+0.02)/2796.35-1. $
                 and spec.zgrid lt stat[i].zqso-mgii_limit $
                 and spec.zgrid gt wave_min/2796.35-1. $
                 and (spec.zgrid gt wave_ciii*(1.+stat[i].zqso+0.01)/2796.35-1. $
                   or spec.zgrid lt wave_ciii*(1.+stat[i].zqso-0.02)/2803.53-1.), ncovered)

    if (keyword_set(nocal)) then begin
    icovered = where(spec.isitcovered $
                 and spec.zgrid gt wave_limit*(1.+stat[i].zqso+0.02)/2796.35-1. $
                 and spec.zgrid lt stat[i].zqso-mgii_limit $
                 and spec.zgrid gt wave_min/2796.35-1. $
                 and (spec.zgrid gt wave_ciii*(1.+stat[i].zqso+0.01)/2796.35-1. $
                   or spec.zgrid lt wave_ciii*(1.+stat[i].zqso-0.02)/2803.53-1.) $
                 and (abs(spec.zgrid-(ca_wave[0]/2796.35-1.)) gt 0.005)  $
                 and (abs(spec.zgrid-(ca_wave[1]/2796.35-1.)) gt 0.005), ncovered)
    endif

    if (ncovered gt 0) then begin
       index_z = (0>(floor((spec.zgrid[icovered] - all_zmin)/z_binsize) < (z_nbin-1)))
       index_ew = (0>(floor((alog10(spec.rew_mgii_2796[icovered])-all_ewmin)/ew_binsize) < (ew_nbin-1)))
       for j=0L, ncovered-1L do begin
           outstr.ncovered_all[index_z[j]] = outstr.ncovered_all[index_z[j]] + 1L
           outstr.ncovered[index_ew[j], index_z[j]] = outstr.ncovered[index_ew[j], index_z[j]] + 1L
       endfor

       idetected = where(spec.isitdetected[icovered], ndetected)
       if (ndetected gt 0) then begin
          for j=0L, ndetected-1L do $
              outstr.ndetected[index_ew[idetected[j]], index_z[idetected[j]]] = $
              outstr.ndetected[index_ew[idetected[j]], index_z[idetected[j]]] + 1L
       endif
    endif
endfor

mwrfits, outstr, outfile, /create

end
