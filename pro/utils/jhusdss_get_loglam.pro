;+
; Doccumentation needed!
;-
function jhusdss_get_loglam, minwave=minwave, maxwave=maxwave, nwave=nwave

if (n_elements(minwave) eq 0) then begin
   splog, "You didn't give me a minimum wavelength, use default 800 AA."
   minwave = 800.d0
endif
if (n_elements(maxwave) eq 0) then begin
   splog, "You didn't give me a minimum wavelength, use default 9000 AA."
   maxwave = 9000.d0
endif
if (maxwave le minwave) then $
    message, "maxwave smaller than minwave!"

dloglam = 1.D-4
rloglam = (alog10(maxwave)-alog10(minwave))
nwave = long(rloglam/dloglam)
loglam = alog10(minwave)+(dindgen(nwave)+0.5)*dloglam

return, loglam
end
