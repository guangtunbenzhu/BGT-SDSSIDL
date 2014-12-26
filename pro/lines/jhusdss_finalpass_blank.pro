function jhusdss_finalpass_blank, nabsmax=nabsmax, lines=lines
if (n_elements(nabsmax) eq 0) then nabsmax = jhusdss_nabsmax()
if (n_elements(lines) eq 0) then lines = jhusdss_finalpass_lines()

absorbers = {zabs:fltarr(nabsmax), $
             err_zabs:fltarr(nabsmax), $
             zabs_all:fltarr(n_elements(lines), nabsmax), $
             ew:fltarr(n_elements(lines), nabsmax), $
             err_ew:fltarr(n_elements(lines), nabsmax), $
             vdisp:fltarr(nabsmax), $
             err_vdisp:fltarr(nabsmax), $
             vdisp_all:fltarr(n_elements(lines), nabsmax), $
             err_vdisp_all:fltarr(n_elements(lines), nabsmax), $
             lines:lines.name}

return, absorbers
end
