pro jhusdss_qaplot_composite_image, nmfver

    qso_image =  jhusdss_qso_composite_image_readin(nmfver)
    qso_resi_image =  jhusdss_qso_composite_image_readin(nmfver, /residual)
    qso =  jhusdss_qso_readin()

    abs_image =  jhusdss_absorber_composite_image_readin(nmfver)
    absorber = jhusdss_absorber_readin(106)

    loglam = jhusdss_get_loglam(minwave=3700D0, maxwave=9200D0)
    wave = 10.^loglam
    iwave_min = value_locate(wave, 3800.)

    y_axis = indgen(n_elements(qso_image.z))
    y_abs_axis = indgen(n_elements(abs_image.z))
    k_print, filename='temp.ps', xsize=10, ysize=10

       loadct, 3
;      Z = ((qso_image.image_median>0.8) < 1.1)
       Z = qso_image.image_median
       contour, Z[iwave_min:*, *], 10.^loglam[iwave_min:*], y_axis, $
          nlevels=40, /fill, $
          xtitle=xtitle, ytitle=ytitle, xst=9, yst=5, $
          min=0., max=6
       stop

       ZZ = ((qso_resi_image.image_median>0.8) < 1.1)
       contour, ZZ[iwave_min:*, *], 10.^loglam[iwave_min:*], y_axis, $
          nlevels=20, /fill, $
          xtitle=xtitle, ytitle=ytitle, xst=9, yst=5

       xmgii = 2800.*(1.+absorber.zabs)
       ymgii = value_locate(qso_resi_image.z, absorber.zqso)
       djs_oplot, xmgii, ymgii, psym=4, symsize=0.1, color='cyan'

       stop

       ZZZ = ((abs_image.image_median>0.8) < 1.1)
       contour, ZZZ[iwave_min:*, *], 10.^loglam[iwave_min:*], y_abs_axis, $
          nlevels=20, /fill, $
          xtitle=xtitle, ytitle=ytitle, xst=9, yst=5

    k_end_print
end
