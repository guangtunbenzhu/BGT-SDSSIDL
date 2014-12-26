function jhusdss_absorber_blank, nabsmax=nabsmax, lines=lines
if (n_elements(nabsmax) eq 0) then nabsmax = jhusdss_nabsmax()
if (n_elements(lines) eq 0) then lines = jhusdss_train_lines()

absorbers = {ra:0D0, dec:0D0, plate:0L, fiber:0L, mjd:0L, zqso:0., err_zqso:0., $
             nabs:0L, zabs_firstpass:fltarr(nabsmax), $
             snr:fltarr(n_elements(lines), nabsmax), $
             signal:fltarr(n_elements(lines), nabsmax), $
             ivar:fltarr(n_elements(lines), nabsmax), $
             magic:bytarr(4, nabsmax), $
             ew_firstpass:fltarr(n_elements(lines), nabsmax), $
             err_ew_firstpass:fltarr(n_elements(lines), nabsmax), $
             lines_firstpass:lines.name}

return, absorbers
end
