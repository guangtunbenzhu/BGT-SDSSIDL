function jhusdss_linemask, wave, lines, z

   velshift = 3000. ; km/s
   c0 = 3.e+5
   nwave = n_elements(wave)
   wavemin = lines.wave*(1.+z)*(1.-velshift/c0)
   wavemax = lines.wave*(1.+z)*(1.+velshift/c0)
   indexmin = value_locate(wave, wavemin) > 0
   indexmin = indexmin < (nwave-1)
   indexmax = value_locate(wave, wavemax) > 0
   indexmax = indexmax < (nwave-1)

   mask = bytarr(n_elements(wave))
   for i=0, n_elements(lines)-1 do mask[indexmin[i]:indexmax[i]] = 1b
   mask[0] = 1b
   mask[n_elements(wave)-1] = 1b

   return, mask
end
