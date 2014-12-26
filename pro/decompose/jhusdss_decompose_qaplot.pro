;+
; Documentation needed!
;-

pro jhusdss_decompose_qaplot, objs, loglam, flux=flux, ivar=ivar, nmf_continuum=nmf_continuum, $
       med_continuum=med_continuum, residual=residual, eigen_values=eigen_values

wave = 10.^loglam
a = 'y'
for i=0L, n_elements(objs)-1L do begin

    ii = where(ivar[i, *] ne 0.)

    djs_plot, wave[ii], flux[i, ii], position=[0.1, 0.55, 0.9, 0.95]
    djs_oplot, wave[ii], nmf_continuum[i, ii], color='red'

    djs_plot, wave[ii], med_continuum[i, ii]*residual[i, ii], position=[0.1, 0.10, 0.9, 0.5], $
        /noerase, yra=[0., 1.5]
    djs_oplot, wave[ii], med_continuum[i, ii], color='red'
    djs_oplot, wave[ii], residual[i, ii]-0.5, color='green'

    print, "Continue: Y; Pause: P (.con); Else: Exit"
    read, a
    if (a eq 'Y' or a eq 'y') then continue
    if (a eq 'P' or a eq 'p') then stop else break

endfor

end
