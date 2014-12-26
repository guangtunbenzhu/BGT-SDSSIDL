;; mask out detected lines, 02/29/2012 (Leap Year!) -- Guangtun
function jhusdss_decompose_interp_spec, spec, z, loglam 
    wave = spec.wave*(1.+spec.z)/(1.+z)
    flux = spec.residual
    ivar = spec.ivar*spec.med_continuum^2*spec.nmf_continuum^2

    ;; Don't know if mask works. Give it a try
    mask = (ivar eq 0.)
    curr_loglam = alog10(wave)
    combine1fiber, curr_loglam, flux, ivar, newloglam=loglam, $
                   newflux=finterp, newivar=iinterp, maxiter=0, $
                   finalmask=mask, andmask=maskinterp

    return, {flux:finterp, ivar:iinterp}
end
