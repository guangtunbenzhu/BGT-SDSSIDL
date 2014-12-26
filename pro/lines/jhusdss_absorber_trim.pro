function jhusdss_absorber_trim, objects, qso=qso, keepasso=keepasso, keepciii=keepciii, keeplowzlowsnr=keeplowzlowsnr, addall=addall
   
;; qso selection
sdev_limit = 0.07
sdev_blue_limit = 0.10
zqso_min = 0.40
zqso_max = 4.7
snr_limit = 0. ;; 2012-Apr-19, change to 0 to keep all qsos since this has little effect.

if (~keyword_set(qso)) then begin

   mask= jhusdss_absorber_stat_window_mask(objects, keepasso=keepasso, keepciii=keepciii, keeplowzlowsnr=keeplowzlowsnr, $
         addall=addall)

   if (~keyword_set(addall)) then begin
   index = where(objects.spec_snr_median gt snr_limit $
             and objects.med_sdeviation_red gt 0.00 $
             and objects.med_sdeviation_red le sdev_limit $
             and objects.zqso ge zqso_min $
             and objects.zqso lt zqso_max $
             and mask eq 0b)
   endif else begin
   ;; need to add Fe II selection, omitted for now
   index = where(objects.spec_snr_median gt snr_limit $
             and ((2796.*(objects.zabs+1.) gt 1550.*(objects.zqso+1.+0.02) and objects.med_sdeviation_red gt 0.00 $
             and objects.med_sdeviation_red le sdev_limit) $
              or (2796.*(objects.zabs+1.) lt 1550.*(objects.zqso+1.-0.02) and objects.med_sdeviation_blue gt 0.00 $
             and objects.med_sdeviation_blue le sdev_blue_limit)) $
             and objects.zqso ge zqso_min $
             and objects.zqso lt zqso_max $
             and mask eq 0b)
   endelse

endif else begin
   index = where(objects.spec_snr_median gt snr_limit $
             and objects.med_sdeviation_red gt 0.00 $
             and objects.med_sdeviation_red le sdev_limit $
             and objects.zqso ge zqso_min $
             and objects.zqso lt zqso_max)
endelse

return, index
end
