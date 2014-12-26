pro jhusdss_qaplot_line_oplot, zabs, wavemin=wavemin, wavemax=wavemax, linepos=linepos

if (n_elements(wavemin) eq 0) then wavemin=1000.
if (n_elements(wavemax) eq 0) then wavemax=3000.
if (n_elements(linepos) eq 0) then linepos = 0.92

linename=['Na D', 'Ca II', 'Mg I', 'Mg II', 'Fe II', 'Fe II', 'Al III', $
           'Al II', 'C IV', 'Si II', 'Si IV', 'C II', 'O I']
linewave=[5890., 3950., 2853.,  2796.,   2590.,   2360.,   1854.,   $
          1661.,   1560.,  1507.,  1390,  1340.,  1294.]
linewave_all = [5890.,  3969.59, 3934.78, 2852.96, 2803.53, 2796.35, 2600.17, $
                2586.55, 2382.77, 2374.46, 2344.21, $
                1862.79, 1854.72, 1670.79, 1550.78, 1548.20, 1526.71, 1402.77, 1393.76, $
                1334.53, 1304.86, 1302.17]

thick=5
charsize=1.5
charthick=2

    for jabs=0L, n_elements(zabs)-1L do begin
        zabstmp = zabs[jabs]
;       jhusdss_qaplot_oplotlines, zabstmp, thiscolor='dark green'
        for iline=0L, n_elements(linewave_all)-1L do begin
            if (linewave_all[iline] gt wavemin) then $
            djs_oplot, replicate(linewave_all[iline]*(1.+zabstmp),2), $
                !y.crange[0]+[linepos-0.08, linepos-0.02]*(!y.crange[1]-!y.crange[0]), $
                color='dark green', thick=thick, linestyle=0
        endfor
        for iline=0L, n_elements(linewave)-1L do begin
            if (linewave[iline] gt wavemin) then $
            djs_xyouts, linewave[iline]*(1.+zabstmp-0.01), $
                !y.crange[0]+linepos*(!y.crange[1]-!y.crange[0]), $
                linename[iline], color='dark green', charsize=0.9, charthick=1.5
        endfor
    endfor

end

