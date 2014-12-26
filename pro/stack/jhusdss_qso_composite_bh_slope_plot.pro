
nmfver = 106

choice_load_data = 0
read,'load data? [1=yes, 0=no]: ',choice_load_data

if choice_load_data eq 1 then begin
   restore, filename='composite_bh_alpha.sav'
endif

thick=8
xthick=8 
ythick=8
charsize=1.4
charthick=3

qapath = jhusdss_get_path(/nmfqso)+'/'+string(nmfver, format='(I3.3)')+'/QAplots/'
psfile1 = qapath + 'qso_composite_bh_slope_plot_1.ps'
psfile2 = qapath + 'qso_composite_bh_slope_plot_2.ps'
xra = [1250, 5500]
yra = [3E-1, 5]

charsize = 1.5
charthick = 3.0

xtitle = '\lambda (\AA)'
ytitle = 'Normalized Flux'

ntry = 5
k_print, filename=psfile1, axis_char_scale=1.3, xsize=10, ysize=6
  djs_plot, composite_bh_alpha[2*5+0].wave, composite_bh_alpha[2*5+0].fluxmedian, xtitle=xtitle, ytitle=ytitle, $
      xra=xra, yra=yra, charsize=charsize, charthick=charthick, $
      xst=1, yst=1, thick=6, xthick=6, ythick=6, /nodata, /ylog

  loadct, 5
  for i=0L, ntry-1L do begin
      icolor = (200.-(i+1)/float(ntry)*150.)
      y = composite_bh_alpha[2*5+i].fluxmedian
      ii = where(y gt 0.)
      y = smooth(y, 3)
      y = y[ii]
      x = composite_bh_alpha[0].wave
      x = x[ii]
      djs_oplot, x, y, thick=4, color=icolor
  endfor
  loadct, 0
k_end_print

xra = [4000, 5500]
yra = [3E-1, 1]
k_print, filename=psfile2, axis_char_scale=1.3, xsize=10, ysize=6
  djs_plot, composite_bh_alpha[2*5+0].wave, composite_bh_alpha[2*5+0].fluxmedian, xtitle=xtitle, ytitle=ytitle, $
      xra=xra, yra=yra, charsize=charsize, charthick=charthick, $
      xst=1, yst=1, thick=6, xthick=6, ythick=6, /nodata, /ylog

  loadct, 5
  for i=0L, ntry-1, 2 do begin
      icolor = (200.-(i+1)/float(ntry)*150.)
      y = composite_bh_alpha[i*5+2].fluxmedian
      ii = where(y gt 0.)
      y = smooth(y, 3)
      y = y[ii]
      x = composite_bh_alpha[0].wave
      x = x[ii]
      djs_oplot, x, y+i*0.01, thick=6, color=icolor
  endfor
  loadct, 0
k_end_print


end
