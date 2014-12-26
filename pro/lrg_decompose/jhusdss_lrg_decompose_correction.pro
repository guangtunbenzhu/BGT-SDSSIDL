;; obsolete

;pro jhusdss_lrg_decompose_correction, lrgver, boss=boss

    lrgver = 101
    lrg_resi_image =  jhusdss_lrg_composite_image_readin(lrgver, /residual)
    nobj = n_elements(lrg_resi_image.z)
;   lrg =  jhusdss_lrg_readin()

    inwave = lrg_resi_image.wave
    ninwave = n_elements(inwave)
    in_iwave = value_locate(inwave, 5000.*(lrg_resi_image.z+1.))

    zmax = 0.6
    loglam = jhusdss_get_loglam(minwave=3800./(1.+zmax), maxwave=9210.)
    outwave = 10.^loglam
    noutwave = n_elements(outwave)
    out_iwave = value_locate(outwave, 5000.)

    newimage = fltarr(ninwave, nobj)

    influx = fltarr(nobj, noutwave)
    inivar = fltarr(nobj, noutwave)+1.

    for i=0L, nobj-1L do begin
        wave_begin = out_iwave-in_iwave[i]
        wave_end = out_iwave-in_iwave[i]+ninwave-1L
        influx[i,wave_begin:wave_end] = lrg_resi_image.image_median[*,i]
    endfor

    jhusdss_composite_engine, influx, inivar, fmean=fmean, fmedian=fmedian

    djs_plot, outwave, fmedian, xst=1, yst=1

    for i=0L, nobj-1L do begin
        wave_begin = out_iwave-in_iwave[i]
        wave_end = out_iwave-in_iwave[i]+ninwave-1L
        newimage[*,i] = lrg_resi_image.image_median[*,i]-fmedian[wave_begin:wave_end]
    endfor

    wait, 1
    Z = ((newimage>0.8) < 1.1)
    contour,(Z),nlevels=my_nlevels,/fill,$
       xmargin=[6,4], ymargin=[6,3], charsize=my_charsize,xstyle=5,ystyle=5

end
