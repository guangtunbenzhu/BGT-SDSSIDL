function jhusdss_absorber_blank_byabs, lines=lines, fit_lines=fit_lines

if (n_elements(lines) eq 0) then lines = jhusdss_train_lines()
if (n_elements(fit_lines) eq 0) then fit_lines = jhusdss_finalpass_lines()

absorbers = {ra:0D0, dec:0D0, plate:0L, fiber:0L, mjd:0L, zqso:0., $
             nabs:0L, $
             zabs:0., $
             snr:fltarr(n_elements(lines)), $
             signal:fltarr(n_elements(lines)), $
             ivar:fltarr(n_elements(lines)), $
             magic:fltarr(4), $
             ew:fltarr(n_elements(lines)), $
             err_ew:fltarr(n_elements(lines)), $
             lines:lines.name, $
             zfit:0., $
             err_zfit:0., $
             zabs_fit:fltarr(n_elements(fit_lines)), $
             ew_fit:fltarr(n_elements(fit_lines)), $
             err_ew_fit:fltarr(n_elements(fit_lines)), $
             sigma_fit_mean:0., $
             err_sigma_fit_mean:0., $
             sigma_fit:fltarr(n_elements(fit_lines)), $
             err_sigma_fit:fltarr(n_elements(fit_lines)), $
             lines_fit:fit_lines.name}

return, absorbers
end
