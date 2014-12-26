function jhusdss_detect_absorbers_index, wave, lines, z, nexpand=nexpand

;; wavelength in increasing order
if (n_elements(nexpand) eq 0) then nexpand = 4
index = lonarr(nexpand, n_elements(lines)) - 1L

tmpindex = value_locate(wave, lines.wave*(1.+z))
iuse = where(tmpindex ge 0 and tmpindex lt n_elements(wave), nuse)

if (nuse eq 0) then return, index

for i=0, nuse-1 do index[*, iuse[i]] = tmpindex[iuse[i]]+indgen(nexpand)-nexpand/2

return, index

end
