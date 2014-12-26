function jhusdss_absorber_stat_window_mask, objects, keepciii=keepciii, keepasso=keepasso, keeplowzlowsnr=keeplowzlowsnr, $
         addall=addall

mask = bytarr(n_elements(objects))

;; Beyond 8700. too noisy
;; wave_max = 8700.
;; ii = where(objects.zabs gt wave_max/2796.35-1. $
;;       and objects.criterion_mgii_feii ne 1, nn)
;; if (nn gt 0) then mask[ii] = 1b

;; mask out Calcium II
ca_wave = [3934.7750, 3969.5901]
for ica=0L, n_elements(ca_wave)-1L do begin
    ii = where(abs(objects.zabs - (ca_wave[ica]/2796.35-1.)) le 0.005, nn)
    if (nn gt 0) then mask[ii] = 1b
;   print, nn
endfor

;; Quasar-associated Mg II
if (~keyword_set(keepasso)) then begin
   mgii_limit = 0.04
   ii = where(objects.zabs ge objects.zqso-mgii_limit, nn)
   if (nn gt 0) then mask[ii] = 1b
;   print, nn
endif

if (~keyword_set(keeplowzlowsnr)) then begin
   lowsnr = 4.5
   ii = where(objects.zabs lt 0.41 and objects.snr_mgii_2796 lt lowsnr, nn)
   if (nn gt 0) then mask[ii] = 1b
;   print, nn
endif

;; QSO's CIV
if (~keyword_set(addall)) then begin
   wave_limit = 1550.
   ii = where(objects.zabs le wave_limit*(1.+objects.zqso+0.02)/2796.35-1., nn)
   if (nn gt 0) then mask[ii] = 1b
;   print, nn
endif

;; This is redundant
wave_lower_limit = 1250.
ii = where(objects.zabs le wave_lower_limit*(1.+objects.zqso+0.02)/2796.35-1., nn)
if (nn gt 0) then mask[ii] = 1b
;print, nn

;; wavelength minimum
wave_min = 3800.
ii = where(objects.zabs le wave_min/2796.35-1., nn)
if (nn gt 0) then mask[ii] = 1b
;print, nn

;; QSO's CIII
if (~keyword_set(keepciii)) then begin
   ciii_wave = 1909.
   ii = where(objects.zabs le ciii_wave*(1.+objects.zqso+0.01)/2796.35-1 $
          and objects.zabs gt ciii_wave*(1.+objects.zqso-0.02)/2803.53-1. $
          and objects.criterion_mgii_feii ne 1, nn) 
   if (nn gt 0) then mask[ii] = 1b
;   print, nn
endif

return, mask
end

