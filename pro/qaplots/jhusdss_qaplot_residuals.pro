pro jhusdss_qaplot_residuals, nmfver

if (n_elements(nmfver) eq 0) then $
    message, "nmf version required so that you won't overwrite files unintentionally."

;; residual
average_residual_file = jhusdss_get_path(/nmfqso)+'/'+ string(nmfver, format='(I3.3)')+$
                        '/Composite/Residual_Composite_'+ string(nmfver, format='(I3.3)')+'.fits'
average_residual =  mrdfits(average_residual_file, 1)

;; init 
thick=4
charsize=1.1
charthick=2.5

path=jhusdss_get_path(/nmfqso)+'/'+string(nmfver, format='(I3.3)')
qapath = path+'/QAplots'
if (jhusdss_direxist(qapath) eq 0) then spawn, 'mkdir -p '+qapath
psfile = qapath+'/'+'Residuals_'+string(nmfver, format='(I3.3)')+'.ps'

k_print, filename=psfile, axis_char_scale=1.3, xsize=11, ysize=4

   xra = [3750, 9250]
   yra = [0.92, 1.18]
   xtitle = 'Observer-frame \lambda (\AA)'
   ytitle = 'Composite Residual of all Quasars'
   x = average_residual.wave
   y = smooth(average_residual.fluxmedian, 5)
   djs_plot, x, y, xra=xra, yra=yra, $
       xtitle=xtitle, title=ytitle, $
       xthick=6, ythick=6, thick=thick, $
       charsize=charsize, charthick=charthick, color='blue'

   linename = ['Na_{ }D', 'Ca_{ }II']
   linewave = [5850., 3900.]
   linewave_all = [5890.,  3969.59, 3934.78]

   zabstmp=0.
   for iline=0L, n_elements(linewave_all)-1L do begin
;      if (linewave_all[iline] gt wavemin) then $
       djs_oplot, replicate(linewave_all[iline]*(1.+zabstmp),2), $
           !y.crange[0]+[0.35, 0.40]*(!y.crange[1]-!y.crange[0]), $
           color='dark green', thick=thick, linestyle=0
   endfor
   for iline=0L, n_elements(linewave)-1L do begin
;      if (linewave[iline] gt wavemin) then $
       djs_xyouts, linewave[iline]*(1.+zabstmp-0.01), $
           !y.crange[0]+0.42*(!y.crange[1]-!y.crange[0]), $
           linename[iline], color='dark green', charsize=0.9, charthick=1.8
   endfor

   linename = ['[O_{ }I] \lambda5577', '[O_{ }I] \lambda6300']
   linewave = [5367., 6090.]
   linewave_all = [5577.,  6300]

   zabstmp=0.
   for iline=0L, n_elements(linewave_all)-1L do begin
;      if (linewave_all[iline] gt wavemin) then $
       djs_oplot, replicate(linewave_all[iline]*(1.+zabstmp),2), $
           !y.crange[0]+[0.22, 0.27]*(!y.crange[1]-!y.crange[0]), $
           color='dark green', thick=thick, linestyle=0
   endfor
   for iline=0L, n_elements(linewave)-1L do begin
;      if (linewave[iline] gt wavemin) then $
       djs_xyouts, linewave[iline]*(1.+zabstmp-0.01), $
           !y.crange[0]+0.16*(!y.crange[1]-!y.crange[0]), $
           linename[iline], color='dark green', charsize=0.9, charthick=1.8
   endfor

   djs_oplot, [6800., 9200.], $
           !y.crange[0]+[0.25, 0.25]*(!y.crange[1]-!y.crange[0]), $
           color='dark green', thick=thick, linestyle=0
   djs_xyouts, 7800., $
           !y.crange[0]+0.18*(!y.crange[1]-!y.crange[0]), $
           'OH bands', color='dark green', charsize=0.9, charthick=1.8
           
;  jhusdss_qaplot_line_oplot, 0, wavemin=3000, wavemax=10000.


k_end_print



end
