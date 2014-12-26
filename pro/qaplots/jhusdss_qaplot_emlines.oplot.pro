pro jhusdss_qaplot_emline_oplot, zabs, wavemin=wavemin, wavemax=wavemax

if (n_elements(wavemin) eq 0) then wavemin=1000.
if (n_elements(wavemax) eq 0) then wavemin=3000.

linename=['[O I]', '[O I]', '[O II]']
linewave=[6300, 5577., 3727.]
linewave_all = [6300.,  5577., 3729, 3726.]

thick=5
charsize=1.5
charthick=2

    for jabs=0L, n_elements(zabs)-1L do begin
        zabstmp = zabs[jabs]
;       jhusdss_qaplot_oplotlines, zabstmp, thiscolor='dark green'
        for iline=0L, n_elements(linewave_all)-1L do begin
            if (linewave_all[iline] gt wavemin) then $
            djs_oplot, replicate(linewave_all[iline]*(1.+zabstmp),2), $
                !y.crange[0]+[0.72, 0.80]*(!y.crange[1]-!y.crange[0]), $
                color='dark green', thick=thick, linestyle=0
        endfor
        for iline=0L, n_elements(linewave)-1L do begin
            if (linewave[iline] gt wavemin) then $
            djs_xyouts, linewave[iline]*(1.+zabstmp-0.01), $
                !y.crange[0]+0.85*(!y.crange[1]-!y.crange[0]), $
                linename[iline], color='dark green', charsize=0.9, charthick=1.5
        endfor
    endfor

end

