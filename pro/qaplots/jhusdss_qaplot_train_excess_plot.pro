pro jhusdss_qaplot_train_excess_plot, nmfver

path=jhusdss_get_path(/nmfqso)+'/'+string(nmfver, format='(I3.3)')
infile = path+'/Absorbers/'+'JHU_nomatched_composite_'+string(nmfver, format='(I3.3)')+'.fits'

comp_high = mrdfits(infile, 1)
comp_medium = mrdfits(infile, 2)
comp_low = mrdfits(infile, 3)
comp_verylow = mrdfits(infile, 4)
comp_1 = mrdfits(infile, 5)
comp_2 = mrdfits(infile, 6)
comp_3 = mrdfits(infile, 7)
comp_4 = mrdfits(infile, 8)

;; init
thick=5
charsize=1.5
charthick=2.5

psfile = path+'/QAplots/'+'JHU_nomatched_composite_one_'+string(nmfver, format='(I3.3)')+'.ps'
k_print, filename=psfile, axis_char_scale=1.4, xsize=14, ysize=6

   x = comp_verylow.wave
   y = comp_verylow.fluxmedian
   y = smooth(comp_verylow.fluxmedian, 5)
   xra = [1350, 3000]
   yra = [0.01, 1.2]
   title = textoidl('Composite spectra of absorbers in JHU but not in Pittsburgh catalog')
   xtitle = textoidl('Absorber-frame \lambda (\AA)')
   ytitle = textoidl('Median Flux Residuals')
;  ytitle = ''

   colors=[djs_icolor('black'), djs_icolor('blue'), djs_icolor('magenta'), djs_icolor('red')]
;  pos = [0.10, 0.70, 0.95, 0.90]
   djs_plot, x+120, y-0.45, xra=xra, yra=yra, $
       xtitle=xtitle, ytitle=ytitle, title=title, $
       thick=thick, xthick=thick, ythick=thick, $
       charsize=charsize, charthick=charthick, $
       xstyle=1, ystyle=1, color=colors[3]
;  legend = textoidl('N(1.5<W^{\lambda2796}_0<4.0 \AA) = '+strtrim(string(comp_high.nabs), 2))
;  djs_xyouts, !x.crange[0]+0.05*(!x.crange[1]-!x.crange[0]), $
;              !y.crange[0]+0.1*(!y.crange[1]-!y.crange[0]), $
;              legend, charsize=charsize, charthick=charthick

   jhusdss_qaplot_line_oplot, 0., wavemin=1350., linepos=0.93
               
   x = comp_low.wave
   y = comp_low.fluxmedian
   y = smooth(comp_low.fluxmedian, 5)
;  yra = [0.45,1.3]
;  pos = [0.10, 0.50, 0.95, 0.70]
   djs_oplot, x+80, y-0.3, xra=xra, yra=yra, $
       xtickformat='(A1)', ytitle='', $
       thick=thick, xthick=thick, ythick=thick, $
       charsize=charsize, charthick=charthick, $
       xstyle=1, ystyle=1, color=colors[2]
;      position=pos, /noerase
;  legend = textoidl('N(1.1<W^{\lambda2796}_0<1.5 \AA) = '+strtrim(string(comp_medium.nabs), 2))
;  djs_xyouts, !x.crange[0]+0.05*(!x.crange[1]-!x.crange[0]), $
;              !y.crange[0]+0.1*(!y.crange[1]-!y.crange[0]), $
;              legend, charsize=charsize, charthick=charthick

   x = comp_medium.wave
   y = comp_medium.fluxmedian
   y = smooth(comp_medium.fluxmedian, 5)
;  yra = [0.45,1.3]
;  pos = [0.10, 0.30, 0.95, 0.50]
   djs_oplot, x+40, y-0.15, xra=xra, yra=yra, $
       xtickformat='(A1)', ytitle='', $
       thick=thick, xthick=thick, ythick=thick, $
       charsize=charsize, charthick=charthick, $
       xstyle=1, ystyle=1, color=colors[1]
;      position=pos, /noerase
;  legend = textoidl('N(0.8<W^{\lambda2796}_0<1.1 \AA) = '+strtrim(string(comp_low.nabs), 2))
;  djs_xyouts, !x.crange[0]+0.05*(!x.crange[1]-!x.crange[0]), $
;              !y.crange[0]+0.1*(!y.crange[1]-!y.crange[0]), $
;              legend, charsize=charsize, charthick=charthick

;  djs_xyouts, !x.crange[0]-0.08*(!x.crange[1]-!x.crange[0]), $
;              !y.crange[0]+0.4*(!y.crange[1]-!y.crange[0]), $
;              ytitle, orientation=90., charsize=1.7, charthick=charthick

   x = comp_high.wave
   y = comp_high.fluxmedian
   y = smooth(comp_high.fluxmedian, 5)
;  yra = [0.45,1.3]
;  pos = [0.10, 0.10, 0.95, 0.30]
   djs_oplot, x, y, xra=xra, yra=yra, $
       xtitle=xtitle, ytitle='', $
       thick=thick, xthick=thick, ythick=thick, $
       charsize=charsize, charthick=charthick, $
       xstyle=1, ystyle=1, color=colors[0]
;      position=pos, /noerase
;  legend = textoidl('N(0.2<W^{\lambda2796}_0<0.8 \AA) = '+strtrim(string(comp_verylow.nabs), 2))
;  djs_xyouts, !x.crange[0]+0.05*(!x.crange[1]-!x.crange[0]), $
;              !y.crange[0]+0.1*(!y.crange[1]-!y.crange[0]), $
;              legend, charsize=charsize, charthick=charthick

   legend_verylow = textoidl('0.2<W_0<0.8 \AA')
;  legend_verylow = textoidl('N(0.2<W_0<0.8 \AA) = ' $
;                 + strtrim(string(comp_verylow.nabs), 2))
   legend_low = textoidl('0.8<W_0<1.1 \AA')
;  legend_low = textoidl('N(0.8<W_0<1.1 \AA) = ' $
;             + strtrim(string(comp_low.nabs), 2))
   legend_medium = textoidl('1.1<W_0<1.5 \AA')
;  legend_medium = textoidl('N(1.1<W_0<1.5 \AA) = ' $
;                + strtrim(string(comp_medium.nabs), 2))
   legend_high = textoidl('1.5<W_0<4.0 \AA')
;  legend_high = textoidl('N(1.5<W_0<4.0 \AA) = ' $
;              + strtrim(string(comp_high.nabs), 2))
   items = [legend_high, legend_medium, legend_low, legend_verylow]
   legend, items, textcolor=colors, box=0, /right, /bottom, $
        charsize=1.2, charthick=3

k_end_print

psfile = path+'/QAplots/'+'JHU_nomatched_composite_one_weak_'+string(nmfver, format='(I3.3)')+'.ps'
k_print, filename=psfile, axis_char_scale=1.4, xsize=14, ysize=6

   x = comp_4.wave
   y = comp_4.fluxmedian
   y = smooth(comp_4.fluxmedian, 5)
   xra = [1350, 3000]
   yra = [0.01, 1.2]
   title = textoidl('Composite spectra of absorbers in JHU but not in Pittsburgh catalog')
   xtitle = textoidl('Absorber-frame \lambda (\AA)')
   ytitle = textoidl('Median Flux Residuals')
;  ytitle = ''

   colors=[djs_icolor('black'), djs_icolor('blue'), djs_icolor('magenta'), djs_icolor('red')]
;  pos = [0.10, 0.70, 0.95, 0.90]
   djs_plot, x+120, y-0.45, xra=xra, yra=yra, $
       xtitle=xtitle, ytitle=ytitle, title=title, $
       thick=thick, xthick=thick, ythick=thick, $
       charsize=charsize, charthick=charthick, $
       xstyle=1, ystyle=1, color=colors[3]
;  legend = textoidl('N(1.5<W^{\lambda2796}_0<4.0 \AA) = '+strtrim(string(comp_high.nabs), 2))
;  djs_xyouts, !x.crange[0]+0.05*(!x.crange[1]-!x.crange[0]), $
;              !y.crange[0]+0.1*(!y.crange[1]-!y.crange[0]), $
;              legend, charsize=charsize, charthick=charthick

   jhusdss_qaplot_line_oplot, 0., wavemin=1350., linepos=0.93
               
   x = comp_3.wave
   y = comp_3.fluxmedian
   y = smooth(comp_3.fluxmedian, 5)
;  yra = [0.45,1.3]
;  pos = [0.10, 0.50, 0.95, 0.70]
   djs_oplot, x+80, y-0.3, xra=xra, yra=yra, $
       xtickformat='(A1)', ytitle='', $
       thick=thick, xthick=thick, ythick=thick, $
       charsize=charsize, charthick=charthick, $
       xstyle=1, ystyle=1, color=colors[2]
;      position=pos, /noerase
;  legend = textoidl('N(1.1<W^{\lambda2796}_0<1.5 \AA) = '+strtrim(string(comp_medium.nabs), 2))
;  djs_xyouts, !x.crange[0]+0.05*(!x.crange[1]-!x.crange[0]), $
;              !y.crange[0]+0.1*(!y.crange[1]-!y.crange[0]), $
;              legend, charsize=charsize, charthick=charthick

   x = comp_2.wave
   y = comp_2.fluxmedian
   y = smooth(comp_2.fluxmedian, 5)
;  yra = [0.45,1.3]
;  pos = [0.10, 0.30, 0.95, 0.50]
   djs_oplot, x+40, y-0.15, xra=xra, yra=yra, $
       xtickformat='(A1)', ytitle='', $
       thick=thick, xthick=thick, ythick=thick, $
       charsize=charsize, charthick=charthick, $
       xstyle=1, ystyle=1, color=colors[1]
;      position=pos, /noerase
;  legend = textoidl('N(0.8<W^{\lambda2796}_0<1.1 \AA) = '+strtrim(string(comp_low.nabs), 2))
;  djs_xyouts, !x.crange[0]+0.05*(!x.crange[1]-!x.crange[0]), $
;              !y.crange[0]+0.1*(!y.crange[1]-!y.crange[0]), $
;              legend, charsize=charsize, charthick=charthick

;  djs_xyouts, !x.crange[0]-0.08*(!x.crange[1]-!x.crange[0]), $
;              !y.crange[0]+0.4*(!y.crange[1]-!y.crange[0]), $
;              ytitle, orientation=90., charsize=1.7, charthick=charthick

   x = comp_1.wave
   y = comp_1.fluxmedian
   y = smooth(comp_1.fluxmedian, 5)
;  yra = [0.45,1.3]
;  pos = [0.10, 0.10, 0.95, 0.30]
   djs_oplot, x, y, xra=xra, yra=yra, $
       xtitle=xtitle, ytitle='', $
       thick=thick, xthick=thick, ythick=thick, $
       charsize=charsize, charthick=charthick, $
       xstyle=1, ystyle=1, color=colors[0]
;      position=pos, /noerase
;  legend = textoidl('N(0.2<W^{\lambda2796}_0<0.8 \AA) = '+strtrim(string(comp_verylow.nabs), 2))
;  djs_xyouts, !x.crange[0]+0.05*(!x.crange[1]-!x.crange[0]), $
;              !y.crange[0]+0.1*(!y.crange[1]-!y.crange[0]), $
;              legend, charsize=charsize, charthick=charthick

   legend_verylow = textoidl('0.2<W_0<0.4 \AA')
;  legend_verylow = textoidl('N(0.2<W_0<0.8 \AA) = ' $
;                 + strtrim(string(comp_verylow.nabs), 2))
   legend_low = textoidl('0.4<W_0<0.5 \AA')
;  legend_low = textoidl('N(0.8<W_0<1.1 \AA) = ' $
;             + strtrim(string(comp_low.nabs), 2))
   legend_medium = textoidl('0.5<W_0<0.65 \AA')
;  legend_medium = textoidl('N(1.1<W_0<1.5 \AA) = ' $
;                + strtrim(string(comp_medium.nabs), 2))
   legend_high = textoidl('0.65<W_0<0.8 \AA')
;  legend_high = textoidl('N(1.5<W_0<4.0 \AA) = ' $
;              + strtrim(string(comp_high.nabs), 2))
   items = [legend_high, legend_medium, legend_low, legend_verylow]
   legend, items, textcolor=colors, box=0, /right, /bottom, $
        charsize=1.2, charthick=3

k_end_print


if (keyword_set(doSeparate)) then begin
psfile = path+'/QAplots/'+'JHU_nomatched_composite_'+string(nmfver, format='(I3.3)')+'.ps'
k_print, filename=psfile, axis_char_scale=1.4, xsize=12, ysize=10

   x = comp_high.wave
   y = comp_high.fluxmedian
   y = smooth(comp_high.fluxmedian, 5)
   xra = [1350, 2900]
   yra = [0.45,1.3]
   title = textoidl('Composite spectra of absorbers in JHU but not in Pittsburgh catalog')
   xtitle = textoidl('Absorber-frame \lambda (\AA)')
   ytitle = textoidl('Composite Residual')

   pos = [0.10, 0.70, 0.95, 0.90]
   djs_plot, x, y, xra=xra, yra=yra, $
       xtickformat='(A1)', ytitle='', title=title, $
       thick=thick, xthick=thick, ythick=thick, $
       charsize=charsize, charthick=charthick, $
       xstyle=1, ystyle=1, $
       position=pos
   legend = textoidl('N(1.5<W^{\lambda2796}_0<4.0 \AA) = '+strtrim(string(comp_high.nabs), 2))
   djs_xyouts, !x.crange[0]+0.05*(!x.crange[1]-!x.crange[0]), $
               !y.crange[0]+0.1*(!y.crange[1]-!y.crange[0]), $
               legend, charsize=charsize, charthick=charthick

   jhusdss_qaplot_line_oplot, 0., wavemin=1350.
               
   x = comp_medium.wave
   y = comp_medium.fluxmedian
   y = smooth(comp_medium.fluxmedian, 5)
   yra = [0.45,1.3]
   pos = [0.10, 0.50, 0.95, 0.70]
   djs_plot, x, y, xra=xra, yra=yra, $
       xtickformat='(A1)', ytitle='', $
       thick=thick, xthick=thick, ythick=thick, $
       charsize=charsize, charthick=charthick, $
       xstyle=1, ystyle=1, $
       position=pos, /noerase
   legend = textoidl('N(1.1<W^{\lambda2796}_0<1.5 \AA) = '+strtrim(string(comp_medium.nabs), 2))
   djs_xyouts, !x.crange[0]+0.05*(!x.crange[1]-!x.crange[0]), $
               !y.crange[0]+0.1*(!y.crange[1]-!y.crange[0]), $
               legend, charsize=charsize, charthick=charthick

   x = comp_low.wave
   y = comp_low.fluxmedian
   y = smooth(comp_low.fluxmedian, 5)
   yra = [0.45,1.3]
   pos = [0.10, 0.30, 0.95, 0.50]
   djs_plot, x, y, xra=xra, yra=yra, $
       xtickformat='(A1)', ytitle='', $
       thick=thick, xthick=thick, ythick=thick, $
       charsize=charsize, charthick=charthick, $
       xstyle=1, ystyle=1, $
       position=pos, /noerase
   legend = textoidl('N(0.8<W^{\lambda2796}_0<1.1 \AA) = '+strtrim(string(comp_low.nabs), 2))
   djs_xyouts, !x.crange[0]+0.05*(!x.crange[1]-!x.crange[0]), $
               !y.crange[0]+0.1*(!y.crange[1]-!y.crange[0]), $
               legend, charsize=charsize, charthick=charthick

   djs_xyouts, !x.crange[0]-0.08*(!x.crange[1]-!x.crange[0]), $
               !y.crange[0]+0.4*(!y.crange[1]-!y.crange[0]), $
               ytitle, orientation=90., charsize=1.7, charthick=charthick

   x = comp_verylow.wave
   y = comp_verylow.fluxmedian
   y = smooth(comp_verylow.fluxmedian, 5)
   yra = [0.45,1.3]
   pos = [0.10, 0.10, 0.95, 0.30]
   djs_plot, x, y, xra=xra, yra=yra, $
       xtitle=xtitle, ytitle='', $
       thick=thick, xthick=thick, ythick=thick, $
       charsize=charsize, charthick=charthick, $
       xstyle=1, ystyle=1, $
       position=pos, /noerase
   legend = textoidl('N(0.2<W^{\lambda2796}_0<0.8 \AA) = '+strtrim(string(comp_verylow.nabs), 2))
   djs_xyouts, !x.crange[0]+0.05*(!x.crange[1]-!x.crange[0]), $
               !y.crange[0]+0.1*(!y.crange[1]-!y.crange[0]), $
               legend, charsize=charsize, charthick=charthick
k_end_print
endif
stop

end
