function jhusdss_highz_composite_blank, nwave, rp_nbin
   if (n_elements(nwave) eq 0) then nwave=1000L
   if (n_elements(rp_nbin) eq 0) then rp_nbin=20L

   stack = {nabs: lonarr(rp_nbin), $
            zabs: fltarr(rp_nbin), $
            rp_min: fltarr(rp_nbin), $
            rp_max: fltarr(rp_nbin), $
            rp: fltarr(rp_nbin), $
            nobjuse:lonarr(nwave, rp_nbin), $
            wave:dblarr(nwave), $
            fluxmean:fltarr(nwave, rp_nbin), $
            fluxmedian:fltarr(nwave, rp_nbin), $
            fluxgeomean:fltarr(nwave, rp_nbin)}

   return, stack
end
