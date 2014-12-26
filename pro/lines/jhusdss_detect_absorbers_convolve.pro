;+
;function jhusdss_detect_absorbers_convolve, wave, dwave, flux, ivar, z, $
;   lines, ew=ew, errew=errew, newivar=newivar
;-
function jhusdss_detect_absorbers_convolve, wave, dwave, flux, ivar, z, $
   lines, ew=ew, errew=errew, newivar=newivar

;; all inputs required
;if (n_elements(lines) eq 0) then lines = jhusdss_abslines_all(/train)

;; index[nexpand, nlines] -> [4,6], 4 pixels, 6 lines
sig = fltarr(n_elements(lines))
;err = fltarr(n_elements(lines))
newivar = fltarr(n_elements(lines))
ew = fltarr(n_elements(lines))
errew = fltarr(n_elements(lines))
index = jhusdss_detect_absorbers_index(wave, lines, z)
index2 = jhusdss_detect_absorbers_index(wave, lines, z, nexpand=10)
for i=0L, n_elements(lines)-1L do begin
    tmpindex = reform(index[*,i])
    tmpindex2 = reform(index2[*,i])
    ii = where(tmpindex ge 0 and tmpindex lt n_elements(wave), nn)
    ii2 = where(tmpindex2 ge 0 and tmpindex2 lt n_elements(wave), nn2)
    ;; at least 2 pixels
    if nn ge 2 then begin
       jj = where(ivar[tmpindex[ii]] gt 0., mm)
       jj2 = where(ivar[tmpindex2[ii2]] gt 0., mm2)
       if mm ge 2 then begin 
          sig[i] = total(1.-flux[tmpindex[ii[jj]]])
;         err[i] = sqrt(total(1./ivar[tmpindex[ii[jj]]]))
          newivar[i] = 1./(total(1./ivar[tmpindex[ii[jj]]]))
          ew[i] = total((1.-flux[tmpindex2[ii2[jj2]]])*dwave[tmpindex2[ii2[jj2]]])/(1.+z)
          errew[i] = sqrt(total(1./ivar[tmpindex2[ii2[jj2]]]*dwave[tmpindex2[ii2[jj2]]]^2))/(1.+z)
       endif
    endif
endfor

;return, sig/(err + (err eq 0.))
return, sig
end
