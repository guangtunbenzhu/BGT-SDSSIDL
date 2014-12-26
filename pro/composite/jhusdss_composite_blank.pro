function jhusdss_composite_blank, nwave

outstr = {nobj:0L, wave:wave, fluxmean:fltarr(nwave), fluxmedian:fltarr(nwave), $
          fluxgeomean:fltarr(nwave), zmin:0., zmax:0., mi_min:0., mi_max:0.}

return, outstr

end
