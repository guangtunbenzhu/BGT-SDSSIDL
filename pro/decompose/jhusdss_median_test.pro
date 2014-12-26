pro jhusdss_median_test, plate, fiber, nmfver=nmfver

    if (n_elements(nmfver) eq 0) then nmfver=004
    spec = jhusdss_decompose_loadspec(plate, fiber, nmfver)
    nmf_residual = reform(spec.flux/spec.nmf_continuum, 1, n_elements(spec.flux))
    tmpivar = reform(spec.ivar*spec.nmf_continuum^2, 1, n_elements(spec.flux))
    nsigma=1.5

    ;; first iteration
    mask = (reform(spec.ivar, 1, n_elements(spec.flux)) eq 0.)
    filter_sizes=[91, 163];, 245]
;   filter_sizes=[61, 143, 215]
;   filter_sizes=[81, 163]
    jhusdss_median_filter, nmf_residual, tmpivar, mask=mask, $
       continuum=tmp_continuum, residual=tmp_residual, filter_sizes=filter_sizes

    djs_plot, spec.wave, smooth(nmf_residual[0,*],3), xra=[3700., 9500.]/(1.+spec.z), xstyle=1, $
              yra=[0.3, 1.5], ystyle=1, position=[0.1, 0.55, 0.9, 0.95]
    djs_oplot, spec.wave, tmp_continuum[0,*], color='red'
    djs_plot, spec.wave, smooth(nmf_residual[0,*]/tmp_continuum[0,*], 3), xra=[3700., 9500.]/(1.+spec.z), xstyle=1, $
              yra=[0.3, 1.5], ystyle=1, position=[0.1, 0.10, 0.9, 0.5], /noerase
    djs_oplot, !x.crange, [1.0, 1.0], color='red'
    djs_oplot, spec.wave, 1.0+1./sqrt(tmpivar), color='gray'
    djs_oplot, spec.wave, 1.0-1./sqrt(tmpivar), color='gray'
    a = 'a'
    read, a

    filter_sizes=[143, 71];, 235]
;   filter_sizes=[71, 153]
    mask = (reform(spec.ivar, 1, n_elements(spec.flux)) eq 0.) or (abs(tmp_continuum-nmf_residual)*sqrt(tmpivar) gt nsigma)
    jhusdss_median_filter, nmf_residual, tmpivar, mask=mask, $
       continuum=tmp_continuum, residual=tmp_residual, filter_sizes=filter_sizes

    djs_plot, spec.wave, smooth(nmf_residual[0,*],3), xra=[3700., 9500.]/(1.+spec.z), xstyle=1, $
              yra=[0.3, 1.5], ystyle=1, position=[0.1, 0.55, 0.9, 0.95]
    djs_oplot, spec.wave, tmp_continuum[0,*], color='red'
    djs_plot, spec.wave, smooth(nmf_residual[0,*]/tmp_continuum[0,*], 3), xra=[3700., 9500.]/(1.+spec.z), xstyle=1, $
              yra=[0.3, 1.5], ystyle=1, position=[0.1, 0.10, 0.9, 0.5], /noerase
    djs_oplot, !x.crange, [1.0, 1.0], color='red'
    djs_oplot, spec.wave, 1.0+1./sqrt(tmpivar), color='gray'
    djs_oplot, spec.wave, 1.0-1./sqrt(tmpivar), color='gray'
    a = 'a'
    read, a


;   filter_sizes=[81, 153, 235]
    filter_sizes=[143, 71];, 225]
;   filter_sizes=[61, 143]
    mask = (reform(spec.ivar, 1, n_elements(spec.flux)) eq 0.) or (abs(tmp_continuum-nmf_residual)*sqrt(tmpivar) gt nsigma)
    jhusdss_median_filter, nmf_residual, tmpivar, mask=mask, $
       continuum=tmp_continuum, residual=tmp_residual, filter_sizes=filter_sizes

    djs_plot, spec.wave, smooth(nmf_residual[0,*],3), xra=[3700., 9500.]/(1.+spec.z), xstyle=1, $
              yra=[0.3, 1.5], ystyle=1, position=[0.1, 0.55, 0.9, 0.95]
    djs_oplot, spec.wave, tmp_continuum[0,*], color='red'
    djs_plot, spec.wave, smooth(nmf_residual[0,*]/tmp_continuum[0,*], 3), xra=[3700., 9500.]/(1.+spec.z), xstyle=1, $
              yra=[0.3, 1.5], ystyle=1, position=[0.1, 0.10, 0.9, 0.5], /noerase
    djs_oplot, !x.crange, [1.0, 1.0], color='red'

    djs_oplot, spec.wave, 1.0+1./sqrt(tmpivar), color='gray'
    djs_oplot, spec.wave, 1.0-1./sqrt(tmpivar), color='gray'


end
