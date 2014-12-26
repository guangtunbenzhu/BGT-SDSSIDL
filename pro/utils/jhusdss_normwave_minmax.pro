function jhusdss_normwave_minmax, option=option

if n_elements(option) ne 1 then begin
   splog, "You didn't specific the wavelength choice, return default [4150, 4250]"
   option = 1
endif

case option of
  1: normwave = [4150., 4250.]  ; 0<z<1.0
  2: normwave = [3020., 3100.]  ; 0.4<z<1.8
  3: normwave = [2150., 2250.]  ; 0.8<z<2.8
  4: normwave = [1420., 1500.]  ; 2.0<z<4.8
  else: normwave = [4150., 4250.]
endcase

return, normwave

end
