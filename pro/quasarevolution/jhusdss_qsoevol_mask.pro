;+
; It has to been in the observer frame, in Angstrom
;-
function jhusdss_qsoevol_mask, wave, zabs, dwave=dwave

lines = jhusdss_qsoevol_masklines()
if (n_elements(dwave) eq 0) then dwave = 3.

mask = bytarr(n_elements(wave))+1B
for i=0L, n_elements(lines)-1L do begin
    iwave = where((wave ge lines[i].wave-dwave) and (wave le linew[i].wave+dwave), nwave)
    if nwave gt 0 then mask[iwave] = 0B
endfor

return mask
end
