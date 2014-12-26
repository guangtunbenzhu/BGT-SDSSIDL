nmfver=106
stackpath = jhusdss_get_path(/nmfqso)+'/'+string(nmfver, format='(I3.3)')+'/Stack/'
infile = stackpath + jhusdss_garching_galqso_stack_filename(nmfver, boss=boss)
infile = repstr(infile, '.fits', '_fit.fits')
a = mrdfits(infile, 1)

ew_nofit_mc = a.ew_nofit_mc
ew_nofit_2_mc = a.ew_nofit_2_mc
  for i=0,11 do begin
      shuffle_infile = repstr(infile, '_fit.fits', '_shuffle_'+string(i+1, format='(i2.2)')+'_fit.fits')
      ashuffle = mrdfits(shuffle_infile, 1)
      ew_nofit_mc = [ew_nofit_mc, ashuffle.ew_nofit_mc]
      ew_nofit_2_mc = [ew_nofit_2_mc, ashuffle.ew_nofit_2_mc]
  endfor

sdev_nofit_mc = fltarr(n_elements(a))
sdev_nofit_2_mc = fltarr(n_elements(a))
for irp=0,n_elements(a)-1 do begin
    tmp = moment(ew_nofit_mc[*,irp], sdev=sdev)
    sdev_nofit_mc[irp] = sdev
    tmp = moment(ew_nofit_2_mc[*,irp], sdev=sdev)
    sdev_nofit_2_mc[irp] = sdev
endfor

line_wave = [3934.79, 3969.59]

psfile1 = repstr(infile, '.fits', '_null_11.ps')
psfile2 = repstr(infile, '.fits', '_null_22.ps')
psfile3 = repstr(infile, '.fits', '_null_33.ps')
psfile4 = repstr(infile, '.fits', '_null_44.ps')
psfile5 = repstr(infile, '.fits', '_null_55.ps')
psfile6 = repstr(infile, '.fits', '_null_66.ps')

;; manually check, independent bins
i_indep = [[0, 1], lindgen(10)*2+2]
i_indep = i_indep[0:8]
;i_indep = lindgen(25)

;best_fit slope
;ifit = i_indep[0:7]
;xfit = alog10(a[ifit].rp)
;yfit = alog10(a[ifit].ew_nofit_2)
;yerror = 0.5*(alog10(a[ifit].ew_nofit_2+a[ifit].sdev_nofit_2_mc)-yfit)+0.5*(yfit-(alog10((a[ifit].ew_nofit_2-a[ifit].sdev_nofit_2_mc)>1.E-10)))
;yfit = alog10(a[ifit].ew_nofit)
;yerror = 0.5*(alog10(a[ifit].ew_nofit+a[ifit].sdev_nofit_mc)-yfit)+0.5*(yfit-(alog10((a[ifit].ew_nofit-a[ifit].sdev_nofit_mc)>1.E-10)))
;coeffs = linfit(xfit, yfit, measure_error=yerror)
;slope = -1.25
;intercept = 0.52
;slope = coeffs[1]
;intercept = coeffs[0]
slope = -1.3853657
intercept = 0.68639344
print, slope, intercept

;iran = [30, 80, 90, 123, 134, 155, 176, 188, 197]
;iran = ((long(randomu(seed, 10)*200) > 0) < 199)
;isort = bsort(iran)
;iran = iran[isort[uniq(isort)]]
;iran = iran[0:8]

;; find top 3 in each bin, not icaii or ired (3969)
nmc = 250
mc_line_wave = findgen(nmc)*1.5+line_wave[0]-60.
itmp = where((mc_line_wave lt line_wave[0]-5.) $
          or ((mc_line_wave gt line_wave[0]+5.) and (mc_line_wave lt line_wave[1]-5.)) $
          or (mc_line_wave gt line_wave[1]+5.), $
          ntmp)
itop3 = lonarr(3,5)
irp_begin=2
irp_end=irp_begin+4
for irp=irp_begin,irp_end do begin
    isort = reverse(bsort(a[i_indep[irp]].ew_nofit_2_mc[itmp]))
;   print, a[i_indep[irp]].ew_nofit_mc[itmp[isort[[0,2,4]]]]
;   itop3[*, irp-2] = itmp[isort[[0,3,7]]]
    itop3[*, irp-irp_begin] = itmp[isort[[0,1,2]]]
;   if (irp eq 1) then itop3[*, irp-irp_begin] = itmp[isort[[0,2,3]]]
    if (irp eq 2) then itop3[*, irp-irp_begin] = itmp[isort[[0,2,3]]]
    if (irp eq 3) then itop3[*, irp-irp_begin] = itmp[isort[[0,3,6]]]
    if (irp eq 4) then itop3[*, irp-irp_begin] = itmp[isort[[0,2,4]]]
    if (irp eq 5) then itop3[*, irp-irp_begin] = itmp[isort[[0,1,6]]]
    if (irp eq 6) then itop3[*, irp-irp_begin] = itmp[isort[[1,3,5]]]
;   if (irp eq 6) then itop3[*, irp-irp_begin] = itmp[isort[[6,7,8]]]
    print, itop3[*,irp-irp_begin]
;   if (irp eq 6) then itop3[*, irp-2] = itmp[isort[[3,5,6]]]
endfor

thick=6
xthick=6
ythick=6
charsize=1.5
charthick=3

xra=[5,1500]
yra=[-3E0, 6E0]
xtitle='r_p (Kpc)' 
ytitle='<W_0^{\lambda}>/\sigma_{<W>}'
title='Single Gaussian Line Profile Measurement'

xx = 10.^(findgen(1000)*0.0024+1.)


pos = [0.20, 0.15, 0.95, 0.9]
k_print, filename=psfile1, axis_char_scale=1.5, xsize=8, ysize=8
  !p.multi=[0,3,5]
  !x.margin=0
  !y.margin=0
  for i=0,4 do begin
  for j=0,2 do begin
      djs_plot, a[i_indep].rp, a[i_indep].ew_nofit, xtickformat='(A1)', ytickformat='(A1)', $
         xra=xra, yra=yra, xst=1, yst=1, thick=thick, xthick=xthick, ythick=ythick, $
         charsize=charsize, charthick=charthick, /nodata, /xlog;, pos=pos
      ibegin=8
      iend=23 
;     polyfill, [a[ibegin:iend].rp, reverse(a[ibegin:iend].rp)], [a[ibegin:iend].sdev_nofit_mc, -reverse(a[ibegin:iend].sdev_nofit_mc)], color=djs_icolor('light gray')
;     djs_oplot, xx, 10.^(alog10(xx)*slope+intercept), thick=thick
;     oploterror, a[i_indep].rp, a[i_indep].ew_nofit, a[i_indep].sdev_nofit_mc, $
;         psym=8, symsize=1.5, color=djs_icolor('blue'), thick=thick, errcolor=djs_icolor('blue'), errthick=thick
      polyfill, [xra, reverse(xra)], [1,1,-1,-1], color=djs_icolor('light gray')
      djs_oplot, a[i_indep].rp, a[i_indep].ew_nofit_mc[itop3[j,i]]/a[i_indep].sdev_nofit_mc, psym=5, symsize=1.1, thick=thick, color='red'
      plotsym, 0, /fill
      djs_oplot, a[i_indep].rp, a[i_indep].ew_nofit/a[i_indep].sdev_nofit_mc, psym=8, symsize=1.5, thick=thick, color='blue'
      plotsym, 4, /fill
      djs_oplot, a[i_indep[2+i]].rp*[1,1], a[i_indep[2+i]].ew_nofit_mc[itop3[j,i]]/a[i_indep[2+i]].sdev_nofit_mc*[1,1], psym=8, symsize=1.4, thick=thick, color='magenta'
      djs_oplot, xra, [0, 0], color='grey'
      if (j eq 0) then djs_axis, yaxis=0, ytitle=ytitle, charsize=charsize, charthick=charthick
      if (i eq 4) then djs_axis, xaxis=0, xtitle=xtitle, charsize=charsize, charthick=charthick
      if (i eq 0 and j eq 2) then legend, [textoidl('Measurement at \lambda=3934 (Ca II)'), textoidl('Null test at \lambda_{top3}')], $
         psym=[8,4], colors=[djs_icolor('blue'), djs_icolor('red')], thick=8, charsize=1.3, charthick=charthick, $
         /top, /right, box=0, textcolors=[djs_icolor('blue'), djs_icolor('red')], pos=[2600,9.5]
      djs_xyouts, 200, 4.5, $
          '\lambda='+string(a[0].wave_mc[itop3[j,i]], format='(i4.4)')+' \AA', $
          charsize=charsize-0.6, charthick=charthick, color='red'
  endfor
  endfor

k_end_print

xra=[5,400]
yra=[-3E0, 6E0]
xtitle='r_p (Kpc)' 
ytitle='<W_0^{\lambda}>/\sigma_{<W>}'
title='Double Gaussian Line Profile Measurement'

k_print, filename=psfile4, axis_char_scale=1.3, xsize=8, ysize=8
  !p.multi=[0,3,5]
  !x.margin=0
  !y.margin=0
  for i=0,4 do begin
  for j=0,2 do begin
      djs_plot, a[i_indep].rp, a[i_indep].ew_nofit_2, xtickformat='(A1)', ytickformat='(A1)', $
         xra=xra, yra=yra, xst=1, yst=1, thick=thick, xthick=xthick, ythick=ythick, $
         charsize=charsize, charthick=charthick, /nodata, /xlog, xticklen=0.04;, pos=pos
;     ibegin=8
;     iend=23 
;     polyfill, [a[ibegin:iend].rp, reverse(a[ibegin:iend].rp)], [a[ibegin:iend].sdev_nofit_mc, -reverse(a[ibegin:iend].sdev_nofit_mc)], color=djs_icolor('light gray')
;     djs_oplot, xx, 10.^(alog10(xx)*slope+intercept), thick=thick
;     oploterror, a[i_indep].rp, a[i_indep].ew_nofit, a[i_indep].sdev_nofit_mc, $
;         psym=8, symsize=1.5, color=djs_icolor('blue'), thick=thick, errcolor=djs_icolor('blue'), errthick=thick
      polyfill, [xra, reverse(xra)], [1,1,-1,-1], color=djs_icolor('light gray')
      plotsym, 0, /fill
      djs_oplot, a[i_indep].rp, a[i_indep].ew_nofit_2/a[i_indep].sdev_nofit_2_mc, psym=8, symsize=1.2, thick=thick, color='light red'
      djs_oplot, a[i_indep].rp, a[i_indep].ew_nofit_2/a[i_indep].sdev_nofit_2_mc, thick=thick, color='light red'
      djs_oplot, a[i_indep].rp, a[i_indep].ew_nofit_2_mc[itop3[j,i]]/a[i_indep].sdev_nofit_2_mc, psym=symcat(14), symsize=2.0, thick=thick, color='blue'
      djs_oplot, a[i_indep].rp, a[i_indep].ew_nofit_2_mc[itop3[j,i]]/a[i_indep].sdev_nofit_2_mc, thick=thick, color='blue'
;     djs_oplot, a.rp, a.ew_nofit_2_mc[itop3[j,i]]/a.sdev_nofit_2_mc, psym=symcat(14), symsize=2.0, thick=thick, color='blue'
;     djs_oplot, a.rp, a.ew_nofit_2_mc[itop3[j,i]]/a.sdev_nofit_2_mc, thick=thick, color='blue'
      djs_oplot, a[i_indep[irp_begin+i]].rp*[1,1], a[i_indep[irp_begin+i]].ew_nofit_2_mc[itop3[j,i]]/a[i_indep[irp_begin+i]].sdev_nofit_2_mc*[1,1], psym=symcat(14), symsize=2.0, thick=thick, color='dark green'
      djs_oplot, xra, [0, 0], color='grey'
      if (j eq 0) then djs_axis, yaxis=0, ytitle=ytitle, charsize=charsize, charthick=charthick
      if (i eq 4) then djs_axis, xaxis=0, xtitle=xtitle, charsize=charsize, charthick=charthick
      if (i eq 0 and j eq 2) then legend, [textoidl('Measurement at \lambda=3934 (Ca II)'), textoidl('Null test at \lambda_{top3}')], $
         psym=[8,symcat(14)], colors=[djs_icolor('light red'), djs_icolor('blue')], thick=8, charsize=1.3, charthick=charthick, $
         /top, /right, box=0, textcolors=[djs_icolor('light red'), djs_icolor('blue')], pos=[500,9.5]
      plotsym, 0, /fill
      if (i eq 0 and j eq 2) then legend, [textoidl('                                   '), textoidl('                           ')], $
         psym=[8,3], colors=[djs_icolor('light red'), djs_icolor('blue')], thick=8, charsize=1.3, charthick=charthick, $
         /top, /right, box=0, textcolors=[djs_icolor('light red'), djs_icolor('blue')], pos=[500,9.5]
      djs_xyouts, 70, 4.5, $
          '\lambda='+string(a[0].wave_mc[itop3[j, i]], format='(i4.4)')+' \AA', $
          charsize=charsize-0.5, charthick=charthick, color='blue'
;     if (i eq 4) then $
;     djs_xyouts, 26, -0.055, '30', charsize=charsize-0.3, charthick=charthick

  endfor
  endfor
k_end_print

xra=[20,400]
ytitle='<W_0^{\lambda}>'
k_print, filename=psfile5, axis_char_scale=1.6, xsize=8, ysize=8
  !p.multi=[0,3,5]
  !x.margin=0
  !y.margin=0
  !x.omargin=[15, 5.0]
  for i=0,4 do begin
  for j=0,2 do begin
      djs_plot, a[i_indep].rp, a[i_indep].ew_nofit_2, xtickformat='(A1)', ytickformat='(A1)', $
         xra=xra, yra=[-4.0E-2, 5E-2], xst=1, yst=1, thick=thick, xthick=xthick, ythick=ythick, $
         charsize=charsize, charthick=charthick, /nodata, /xlog;, pos=pos
      ibegin=2
      iend=15
      polyfill, [a[ibegin:iend].rp, reverse(a[ibegin:iend].rp)], [sdev_nofit_2_mc[ibegin:iend], -reverse(sdev_nofit_2_mc[ibegin:iend])], color=djs_icolor('light gray')

      djs_oplot, a[i_indep].rp, a[i_indep].ew_nofit_2, psym=8, symsize=1.2, thick=thick, color='light red'
      djs_oplot, a[i_indep].rp, a[i_indep].ew_nofit_2, thick=10, color='light red'

      djs_oplot, a[i_indep].rp, a[i_indep].ew_nofit_2_mc[itop3[j,i]], psym=symcat(14), symsize=2.0, thick=thick, color='blue'
      djs_oplot, a[i_indep].rp, a[i_indep].ew_nofit_2_mc[itop3[j,i]], thick=10, color='blue'
      djs_oplot, a[i_indep[irp_begin+i]].rp*[1,1], a[i_indep[irp_begin+i]].ew_nofit_2_mc[itop3[j,i]]*[1,1], psym=symcat(14), symsize=2.0, thick=thick, color='dark green'
      djs_oplot, xra, [0, 0], color='grey'
      if (j eq 0) then djs_axis, yaxis=0, ytitle=ytitle, charsize=charsize, charthick=charthick
      if (i eq 4) then djs_axis, xaxis=0, xtitle=xtitle, charsize=charsize, charthick=charthick
      if (i eq 0 and j eq 2) then legend, [textoidl('Measurement at \lambda=3934 (Ca II)'), textoidl('Null test at random \lambda')], $
         psym=[8,symcat(14)], colors=[djs_icolor('light red'), djs_icolor('blue')], thick=8, charsize=1.3, charthick=charthick, $
         /top, /right, box=0, textcolors=[djs_icolor('light red'), djs_icolor('blue')], pos=[400,0.085]
      plotsym, 0, /fill
      if (i eq 0 and j eq 2) then legend, [textoidl('                                   '), textoidl('                           ')], $
         psym=[8,3], colors=[djs_icolor('light red'), djs_icolor('blue')], thick=8, charsize=1.3, charthick=charthick, $
         /top, /right, box=0, textcolors=[djs_icolor('light red'), djs_icolor('blue')], pos=[400,0.085]
      djs_xyouts, 100, 0.035, $
          '\lambda='+string(a[0].wave_mc[itop3[j, i]], format='(i4.4)')+' \AA', $
          charsize=charsize-0.5, charthick=charthick, color='blue'
      if (i eq 4) then $
      djs_xyouts, 26, -0.055, '30', charsize=charsize-0.3, charthick=charthick
  endfor
  endfor
k_end_print

end
