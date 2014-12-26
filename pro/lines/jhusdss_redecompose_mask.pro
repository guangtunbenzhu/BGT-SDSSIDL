;; see jhusdss_linemask.pro
function jhusdss_redecompose_mask, wave, lines, z

dang = 5.
nwave = n_elements(wave)

mask = bytarr(nwave)
for i=0L, n_elements(z)-1L do begin
    wavemin = (lines.wave-5.)*(1.+z[i])
    wavemax = (lines.wave+5.)*(1.+z[i])
    indexmin = value_locate(wave, wavemin) > 0
    indexmin = indexmin < (nwave-1)
    indexmax = value_locate(wave, wavemax) > 0
    indexmax = indexmax < (nwave-1)

    for j=0, n_elements(lines)-1 do mask[indexmin[j]:indexmax[j]] = 1b
    mask[0] = 1b
    mask[n_elements(wave)-1] = 1b
endfor

return, mask

end
