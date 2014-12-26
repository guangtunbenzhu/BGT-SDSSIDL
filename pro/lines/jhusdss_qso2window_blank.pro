function jhusdss_qso2window_blank, nabsmax=nabsmax, trainlines=trainlines, finallines=finallines, absorber=absorber

if (n_elements(nabsmax) eq 0) then nabsmax = jhusdss_nabsmax()
if (n_elements(trainlines) eq 0) then trainlines = jhusdss_train_lines()
if (n_elements(finallines) eq 0) then finallines = jhusdss_finalpass_lines()

if (~keyword_set(absorber)) then begin
   tmp_qso = {ra:0D0, dec:0D0, plate:0L, fiber:0L, mjd:0L, zqso:0., err_zqso:0., $
              spec_snr_median:0., med_sdeviation_red:-999, med_sdeviation_blue:-999.}
   return, tmp_qso
endif else begin
   tmp_absorber = {nabs:0L, $
                   ;; first pass
;                  zabs_firstpass:fltarr(nabsmax), $
                   snr:fltarr(n_elements(trainlines), nabsmax), $
                   signal:fltarr(n_elements(trainlines), nabsmax), $
;                  ivar:fltarr(n_elements(trainlines), nabsmax), $
                   criterion_mgii:bytarr(nabsmax), $
                   criterion_mgii_feii:bytarr(nabsmax), $
                   criterion_feii:bytarr(nabsmax), $
;                  ew_firstpass:fltarr(n_elements(trainlines), nabsmax), $
;                  err_ew_firstpass:fltarr(n_elements(trainlines), nabsmax), $
;                  lines_firstpass:lines.name, $
                   ;; final pass
                    zabs:fltarr(nabsmax), $
                    err_zabs:fltarr(nabsmax), $
                    zabs_all:fltarr(n_elements(finallines), nabsmax), $
                    ew:fltarr(n_elements(finallines), nabsmax), $
                    err_ew:fltarr(n_elements(finallines), nabsmax), $
                    vdisp:fltarr(nabsmax), $
                    err_vdisp:fltarr(nabsmax), $
                    vdisp_all:fltarr(n_elements(finallines), nabsmax), $
                    err_vdisp_all:fltarr(n_elements(finallines), nabsmax), $
                    lines:finallines.name}

return, tmp_absorber
endelse

end
