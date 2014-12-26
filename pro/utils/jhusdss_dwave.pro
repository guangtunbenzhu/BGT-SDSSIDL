function jhusdss_dwave, wave
   nwave = n_elements(wave)
   dwave = fltarr(nwave)
   wave1 = wave[1:nwave-1]-wave[0:nwave-2]
   dwave[1:nwave-2] = (wave1[0:nwave-3]+wave1[1:nwave-2])/2.
   dwave[0] = wave1[0]
   dwave[nwave-1] = wave1[nwave-2]

   return, dwave
end
