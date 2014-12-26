nmfver=106
stackpath = jhusdss_get_path(/nmfqso)+'/'+string(nmfver, format='(I3.3)')+'/Stack/'
infile = stackpath + jhusdss_garching_lrgqso_stack_filename(nmfver, boss=boss)
boot_infile = repstr(infile, '.fits', '_bootstrap_fit_samesigma_moment.fits')
infile = repstr(infile, '.fits', '_fit_samesigma.fits')
a = mrdfits(infile, 1)
a_boot = mrdfits(boot_infile, 1)
;; The last two bins did not get calculated, manually set:
a_boot[28:31].err_ew = a[28:31].sdev_nofit_lr_1_s5_mc*2.0
a_boot[29:31].err_sigma = [2.2, 2.2, 2.2]

ew_mc = a.ew_mc
ew_nofit_mc = a.ew_nofit_mc
ew_nofit_2_mc = a.ew_nofit_2_mc

;; temporary

for i=3,3+12 do begin
;   print, i
    shuffle_infile = repstr(infile, '_fit_samesigma.fits', '_shuffle_'+string(i+1, format='(i2.2)')+'_fit_samesigma.fits')
    ashuffle = mrdfits(shuffle_infile, 1)
    ew_nofit_mc = [ew_nofit_mc, ashuffle.ew_nofit_mc]
    ew_nofit_2_mc = [ew_nofit_2_mc, ashuffle.ew_nofit_2_mc]
endfor

sdev_mc = fltarr(n_elements(a))
sdev_nofit_mc = fltarr(n_elements(a))
sdev_nofit_2_mc = fltarr(n_elements(a))
for irp=0,n_elements(a)-1 do begin
    tmp = moment(ew_mc[*,irp], sdev=sdev)
    sdev_mc[irp] = sdev
    tmp = moment(ew_nofit_mc[*,irp], sdev=sdev)
    sdev_nofit_mc[irp] = sdev
    tmp = moment(ew_nofit_2_mc[*,irp], sdev=sdev)
    sdev_nofit_2_mc[irp] = sdev
endfor

psfile1 = repstr(infile, '.fits', '_shuffle_null_1.ps')
psfile2 = repstr(infile, '.fits', '_shuffle_null_2.ps')
psfile3 = repstr(infile, '.fits', '_shuffle_null_3.ps')
psfile4 = repstr(infile, '.fits', '_shuffle_null_4.ps')
psfile5 = repstr(infile, '.fits', '_shuffle_null_5.ps')
psfile6 = repstr(infile, '.fits', '_shuffle_null_6.ps')

;; manually check, independent bins
i_indep = lindgen(16)*2+1

;best_fit slope
;ifit = i_indep[0:6]
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
;print, slope, intercept
slope = -1.3853657
intercept = 0.68639344

;iran = [30, 80, 90, 123, 134, 155, 176, 188, 197]
;iran = ((long(randomu(seed, 10)*200) > 0) < 199)
;isort = bsort(iran)
;iran = iran[isort[uniq(isort)]]
;iran = iran[0:8]

thick=8
xthick=8
ythick=8
charsize=1.5
charthick=3

xra=[5,1500]
yra=[-3E0, 6E0]
xtitle='r_p (kpc)' 
ytitle='<W_0^{\lambda}>/\sigma_{<W>}'
title='Single Gaussian Line Profile Measurement'

xx = 10.^(findgen(1000)*0.0024+1.)

plotsym, 0, /fill

pos = [0.20, 0.15, 0.95, 0.9]
k_print, filename=psfile1, axis_char_scale=1.5, xsize=8, ysize=8
  !p.multi=[0,3,3]
  !x.margin=0
  !y.margin=0
  for i=0,8 do begin
      shuffle_infile = repstr(infile, '_fit_samesigma.fits', '_shuffle_'+string(i+4, format='(i2.2)')+'_fit_samesigma.fits')
      ashuffle = mrdfits(shuffle_infile, 1)

      djs_plot, a[i_indep].rp, a[i_indep].ew_nofit, xtickformat='(A1)', ytickformat='(A1)', $
         xra=xra, yra=yra, xst=1, yst=1, thick=thick, xthick=xthick, ythick=ythick, $
         charsize=charsize, charthick=charthick, /nodata, /xlog;, pos=pos
      ibegin=8
      iend=23 
;     polyfill, [a[ibegin:iend].rp, reverse(a[ibegin:iend].rp)], [a[ibegin:iend].sdev_nofit_mc, -reverse(a[ibegin:iend].sdev_nofit_mc)], color=djs_icolor('light gray')
;     djs_oplot, xx, 10.^(alog10(xx)*slope+intercept), thick=thick
;     oploterror, a[i_indep].rp, a[i_indep].ew_nofit, a[i_indep].sdev_nofit_mc, $
;         psym=8, symsize=1.5, color=djs_icolor('dark gray'), thick=thick, errcolor=djs_icolor('dark gray'), errthick=thick
      polyfill, [xra, reverse(xra)], [1,1,-1,-1], color=djs_icolor('light gray')
      djs_oplot, a[i_indep].rp, ashuffle[i_indep].ew_nofit/a[i_indep].sdev_nofit_mc, psym=5, symsize=1.2, thick=thick, color='red'
      djs_oplot, a[i_indep].rp, a[i_indep].ew_nofit/a[i_indep].sdev_nofit_mc, psym=8, symsize=1.6, thick=thick, color='dark gray'
;     djs_oplot, a[i_indep].rp, ashuffle[i_indep].ew_nofit/a[i_indep].sdev_nofit_mc, psym=5, symsize=1.2, thick=thick, color='red'
;     djs_oplot, a[i_indep].rp, a[i_indep].ew_nofit/a[i_indep].sdev_nofit_mc, psym=8, symsize=1.6, thick=thick, color='dark gray'
      djs_oplot, xra, [0, 0], color='grey'
      if (i eq 0 or i eq 3 or i eq 6) then djs_axis, yaxis=0, ytitle=ytitle, charsize=charsize, charthick=charthick
      if (i eq 6 or i eq 7 or i eq 8) then djs_axis, xaxis=0, xtitle=xtitle, charsize=charsize, charthick=charthick
      if (i eq 2) then legend, [textoidl('Measurement at \lambda=3934 (Ca II)'), textoidl('Null tests using random quasars')], $
         psym=[8,4], colors=[djs_icolor('dark gray'), djs_icolor('red')], thick=8, charsize=1.3, charthick=charthick, $
         /top, /right, box=0, textcolors=[djs_icolor('dark gray'), djs_icolor('red')], pos=[2600,7.9]

  endfor
k_end_print

if (keyword_set(plotall)) then begin
k_print, filename=psfile2, axis_char_scale=1.3, xsize=8, ysize=6
  djs_plot, a.rp, a.ew_nofit_mc[iran], xtitle=xtitle, ytitle=ytitle, title=title, $
      xra=xra, yra=yra, xst=1, yst=1, thick=thick, xthick=xthick, ythick=ythick, $
      charsize=charsize, charthick=charthick, /nodata, /xlog, pos=pos
  ibegin=1
  iend=23 
  polyfill, [a[ibegin:iend].rp, reverse(a[ibegin:iend].rp)], [a[ibegin:iend].sdev_nofit_mc, -reverse(a[ibegin:iend].sdev_nofit_mc)], color=djs_icolor('light gray')
  djs_oplot, xx, 10.^(alog10(xx)*slope+intercept), thick=thick
  oploterror, a.rp, a.ew_nofit_mc[iran], a.sdev_nofit_mc, $
      psym=8, symsize=2, color=djs_icolor('dark gray'), thick=thick, errcolor=djs_icolor('dark gray'), errthick=thick
k_end_print

k_print, filename=psfile3, axis_char_scale=1.3, xsize=8, ysize=6
  xra=[220, 1500]
  yra=[-4E-3, 4E-3]
  djs_plot, a.rp, a.ew_nofit_mc[iran], xtitle=xtitle, ytitle=ytitle, title=title, $
      xra=xra, yra=yra, xst=1, yst=1, thick=thick, xthick=xthick, ythick=ythick, $
      charsize=charsize, charthick=charthick, /nodata, /xlog, pos=pos;, /ylog
  ibegin=13
  iend=23 
  polyfill, [a[ibegin:iend].rp, reverse(a[ibegin:iend].rp)], [a[ibegin:iend].sdev_nofit_mc, -reverse(a[ibegin:iend].sdev_nofit_mc)], color=djs_icolor('light gray')
  djs_oplot, xx, 10.^(alog10(xx)*slope+intercept), thick=thick
  oploterror, a.rp, a.ew_nofit_mc[iran], a.sdev_nofit_mc, $
      psym=8, symsize=2, color=djs_icolor('dark gray'), thick=thick, errcolor=djs_icolor('dark gray'), errthick=thick
k_end_print
endif

;ichoose=[1, 2, 3, 4, 7, 8, 9, 10, 12]
;ichoose=[1, 2, 3, 4, 5, 6, 7, 8, 9]
ichoose=lindgen(12)+1

xra=[20,30000]/1E3
yra=[-3E0, 7.5E0]
xtitle='r_p' 
ytitle='<W_0^{\lambda}>/\sigma_{<W>}'
title='Double Gaussian Line Profile Measurement'

;; one panel
k_print, filename=psfile2, axis_char_scale=1.6, xsize=8, ysize=6
; !p.multi=[0,3,4]
  !p.multi=[0,1,1]
  !x.margin=0
  !y.margin=0
  !x.omargin=[10, 4]
  !y.omargin=[6, 2]
  djs_plot, a[i_indep].rp, a[i_indep].ew, xtickformat='(A1)', ytitle=ytitle, $
      xra=xra, yra=yra, thick=thick, xthick=xthick, ythick=ythick, $
      charsize=charsize, charthick=charthick, /nodata, /xlog, ystyle=1, xstyle=5

  nmax=12
  zcolor=fltarr(nmax)
  loadct, 0
  for i=0,11 do begin
      zcolor[i] = (i+1)/float(nmax)*200.-1.
      shuffle_infile = repstr(infile, '_fit_samesigma.fits', '_shuffle_'+string(ichoose[i], format='(i2.2)')+'_fit_samesigma.fits')
      ashuffle = mrdfits(shuffle_infile, 1)

      djs_oplot, a[i_indep].rp/1E3, ashuffle[i_indep].ew_nofit_lr_1_s5/a[i_indep].sdev_nofit_lr_1_s5_mc, psym=symcat(14), symsize=1.8, thick=thick, color=zcolor[i]
      djs_oplot, a[i_indep].rp/1E3, ashuffle[i_indep].ew_nofit_lr_1_s5/a[i_indep].sdev_nofit_lr_1_s5_mc, thick=10, color=zcolor[i]
;     djs_oplot, a[i_indep].rp/1E3, ashuffle[i_indep].ew_nofit_lr_1_s4*2/0.94/a_boot[i_indep].err_ew, psym=symcat(14), symsize=1.8, thick=thick, color=zcolor[i]
;     djs_oplot, a[i_indep].rp/1E3, ashuffle[i_indep].ew_nofit_lr_1_s4*2/0.94/a_boot[i_indep].err_ew, thick=10, color=zcolor[i]
;     djs_oplot, a[i_indep].rp/1E3, ashuffle[i_indep].ew*(1./ashuffle[i_indep].line_ratio+1.)/a_boot[i_indep].err_ew, psym=symcat(14), symsize=1.8, thick=thick, color=zcolor[i]
;     djs_oplot, a[i_indep].rp/1E3, ashuffle[i_indep].ew*(1./ashuffle[i_indep].line_ratio+1.)/a_boot[i_indep].err_ew, thick=10, color=zcolor[i]

;     djs_oplot, a[i_indep].rp, ashuffle[i_indep].ew_nofit_2/a[i_indep].sdev_nofit_2_mc, psym=symcat(14), symsize=1.8, thick=thick, color=zcolor[i]
;     djs_oplot, a[i_indep].rp, ashuffle[i_indep].ew_nofit_2/a[i_indep].sdev_nofit_2_mc, thick=thick, color=zcolor[i], linestyle=0
      djs_oplot, xra, [0, 0], color='grey'
;     if (i mod 3 eq 0) then djs_axis, yaxis=0, ytitle=ytitle, charsize=charsize, charthick=charthick
;     if (i ge 9) then djs_axis, xaxis=0, xtitle=xtitle, charsize=charsize, charthick=charthick
      if (i eq 2) then legend, [textoidl('Measurement at \lambda(Mg II)'), textoidl('Null tests using random quasars')], $
         psym=[8,symcat(14)], colors=[djs_icolor('blue'), djs_icolor('dark gray')], thick=8, charsize=1.5, charthick=charthick, $
         /top, /right, box=0, textcolors=[djs_icolor('blue'), djs_icolor('dark gray')], pos=[28,7.2]
      plotsym, 0, /fill
      if (i eq 2) then legend, [textoidl('                                   '), textoidl('                           ')], $
         psym=[8,3], colors=[djs_icolor('blue'), djs_icolor('dark gray')], thick=8, charsize=1.5, charthick=charthick, $
         /top, /right, box=0, textcolors=[djs_icolor('blue'), djs_icolor('dark gray')], pos=[28,7.2]
  endfor
  loadct, 0
; djs_oplot, a[i_indep].rp, a[i_indep].ew_nofit_2/a[i_indep].sdev_nofit_2_mc, psym=8, symsize=3.0, thick=thick, color='blue'
; djs_oplot, a[i_indep].rp, a[i_indep].ew_nofit_2/a[i_indep].sdev_nofit_2_mc, thick=10, color='blue'
  djs_oplot, a[i_indep].rp/1E3, a[i_indep].ew*(1.+1./a[i_indep].line_ratio)/a[i_indep].sdev_nofit_lr_1_s5_mc/2.0, psym=8, symsize=3, thick=thick, color='blue'
  djs_oplot, a[i_indep].rp/1E3, a[i_indep].ew*(1.+1./a[i_indep].line_ratio)/a[i_indep].sdev_nofit_lr_1_s5_mc/2.0, thick=10, color='blue'
; djs_oplot, a[i_indep].rp/1E3, a[i_indep].ew*(1.+1./a[i_indep].line_ratio)/a_boot[i_indep].err_ew, psym=8, symsize=3, thick=thick, color='blue'
; djs_oplot, a[i_indep].rp/1E3, a[i_indep].ew*(1.+1./a[i_indep].line_ratio)/a_boot[i_indep].err_ew, thick=10, color='blue'

  xtickv = [0.1, 1., 10.]
  xticknames = ['100 kpc', '1 Mpc', '10 Mpc']
  djs_axis, xaxis=0, xtitle=xtitle, charsize=charsize, charthick=charthick, xra=xra, xthick=xthick, $
      xtickv=xtickv, xtickname=xticknames;, xticks=n_elements(xtickv)
  djs_axis, xaxis=1, xtitle='', charsize=charsize, charthick=charthick, xra=xra, xthick=xthick, xtickformat='(A1)'

k_end_print


k_print, filename=psfile4, axis_char_scale=1.6, xsize=8, ysize=8
  !p.multi=[0,3,4]
  !x.margin=0
  !y.margin=0
  !x.omargin=[15, 4]
  for i=0,11 do begin
      shuffle_infile = repstr(infile, '_fit_samesigma.fits', '_shuffle_'+string(ichoose[i], format='(i2.2)')+'_fit_samesigma.fits')
      ashuffle = mrdfits(shuffle_infile, 1)

      djs_plot, a[i_indep].rp, a[i_indep].ew_nofit_2, xtickformat='(A1)', ytickformat='(A1)', $
         xra=xra, yra=yra, xst=1, yst=1, thick=thick, xthick=xthick, ythick=ythick, $
         charsize=charsize, charthick=charthick, /nodata, /xlog;, pos=pos
;     ibegin=8
;     iend=23 
;     polyfill, [a[ibegin:iend].rp, reverse(a[ibegin:iend].rp)], [a[ibegin:iend].sdev_nofit_mc, -reverse(a[ibegin:iend].sdev_nofit_mc)], color=djs_icolor('light gray')
;     djs_oplot, xx, 10.^(alog10(xx)*slope+intercept), thick=thick
;     oploterror, a[i_indep].rp, a[i_indep].ew_nofit, a[i_indep].sdev_nofit_mc, $
;         psym=8, symsize=1.5, color=djs_icolor('dark gray'), thick=thick, errcolor=djs_icolor('dark gray'), errthick=thick
      polyfill, [xra, reverse(xra)], [1,1,-1,-1], color=djs_icolor('light gray')
      djs_oplot, a[i_indep].rp, a[i_indep].ew_nofit_2/a[i_indep].sdev_nofit_2_mc, psym=8, symsize=1.4, thick=thick, color='blue'
      djs_oplot, a[i_indep].rp, a[i_indep].ew_nofit_2/a[i_indep].sdev_nofit_2_mc, thick=10, color='blue'
      djs_oplot, a[i_indep].rp, ashuffle[i_indep].ew_nofit_2/a[i_indep].sdev_nofit_2_mc, psym=symcat(14), symsize=1.5, thick=thick, color='dark gray'
      djs_oplot, a[i_indep].rp, ashuffle[i_indep].ew_nofit_2/a[i_indep].sdev_nofit_2_mc, thick=10, color='dark gray'
      djs_oplot, xra, [0, 0], color='grey'
      if (i mod 3 eq 0) then djs_axis, yaxis=0, ytitle=ytitle, charsize=charsize, charthick=charthick
      if (i ge 9) then djs_axis, xaxis=0, xtitle=xtitle, charsize=charsize, charthick=charthick
      if (i eq 2) then legend, [textoidl('Measurement at \lambda=3934 (Ca II)'), textoidl('Null tests using random quasars')], $
         psym=[8,symcat(14)], colors=[djs_icolor('blue'), djs_icolor('dark gray')], thick=8, charsize=1.5, charthick=charthick, $
         /top, /right, box=0, textcolors=[djs_icolor('blue'), djs_icolor('dark gray')], pos=[500,9.0]
      plotsym, 0, /fill
      if (i eq 2) then legend, [textoidl('                                   '), textoidl('                           ')], $
         psym=[8,3], colors=[djs_icolor('blue'), djs_icolor('dark gray')], thick=8, charsize=1.5, charthick=charthick, $
         /top, /right, box=0, textcolors=[djs_icolor('blue'), djs_icolor('dark gray')], pos=[500,9.0]
  endfor
k_end_print

xra=[20,400]
ytitle='<W_0^{\lambda}>'
k_print, filename=psfile5, axis_char_scale=1.6, xsize=8, ysize=8
  !p.multi=[0,3,3]
  !x.margin=0
  !y.margin=0
  !x.omargin=[15, 5.0]
  for i=0,8 do begin
      shuffle_infile = repstr(infile, '_fit_samesigma.fits', '_shuffle_'+string(ichoose[i], format='(i2.2)')+'_fit_samesigma.fits')
      ashuffle = mrdfits(shuffle_infile, 1)

      djs_plot, a[i_indep].rp, a[i_indep].ew_nofit_2, xtickformat='(A1)', ytickformat='(A1)', $
         xra=xra, yra=[-4.0E-2, 5E-2], xst=1, yst=1, thick=thick, xthick=xthick, ythick=ythick, $
         charsize=charsize, charthick=charthick, /nodata, /xlog;, pos=pos
      ibegin=2
      iend=15
      polyfill, [a[ibegin:iend].rp, reverse(a[ibegin:iend].rp)], [sdev_nofit_2_mc[ibegin:iend], -reverse(sdev_nofit_2_mc[ibegin:iend])], color=djs_icolor('light gray')
;     ibegin=8
;     iend=23 
;     polyfill, [a[ibegin:iend].rp, reverse(a[ibegin:iend].rp)], [a[ibegin:iend].sdev_nofit_mc, -reverse(a[ibegin:iend].sdev_nofit_mc)], color=djs_icolor('light gray')
;     djs_oplot, xx, 10.^(alog10(xx)*slope+intercept), thick=thick
;     oploterror, a[i_indep].rp, a[i_indep].ew_nofit, a[i_indep].sdev_nofit_mc, $
;         psym=8, symsize=1.5, color=djs_icolor('dark gray'), thick=thick, errcolor=djs_icolor('dark gray'), errthick=thick
      djs_oplot, a[i_indep].rp, a[i_indep].ew_nofit_2, psym=8, symsize=1.2, thick=thick, color='blue'
      djs_oplot, a[i_indep].rp, a[i_indep].ew_nofit_2, thick=10, color='blue'

      djs_oplot, a[i_indep].rp, ashuffle[i_indep].ew_nofit_2, psym=symcat(14), symsize=2.0, thick=thick, color='dark gray'
      djs_oplot, a[i_indep].rp, ashuffle[i_indep].ew_nofit_2, thick=10, color='dark gray'

      djs_oplot, xra, [0, 0], color='grey'
      if (i eq 0 or i eq 3 or i eq 6) then djs_axis, yaxis=0, ytitle=ytitle, charsize=charsize, charthick=charthick
      if (i eq 6 or i eq 7 or i eq 8) then djs_axis, xaxis=0, xtitle=xtitle, charsize=charsize, charthick=charthick
      if (i eq 2) then legend, [textoidl('Measurement at \lambda=3934 (Ca II)'), textoidl('Null test using random quasars')], $
         psym=[8,symcat(14)], colors=[djs_icolor('blue'), djs_icolor('dark gray')], thick=8, charsize=1.3, charthick=charthick, $
         /top, /right, box=0, textcolors=[djs_icolor('blue'), djs_icolor('dark gray')], pos=[400,0.07]
      plotsym, 0, /fill
      if (i eq 2) then legend, [textoidl('                                   '), textoidl('                           ')], $
         psym=[8,3], colors=[djs_icolor('blue'), djs_icolor('dark gray')], thick=8, charsize=1.3, charthick=charthick, $
         /top, /right, box=0, textcolors=[djs_icolor('blue'), djs_icolor('dark gray')], pos=[400,0.07]
      if (i eq 6 or i eq 7 or i eq 8) then $
      djs_xyouts, 26, -0.049, '30', charsize=charsize-0.3, charthick=charthick
  endfor
k_end_print

if (keyword_set(plotall)) then begin
k_print, filename=psfile5, axis_char_scale=1.3, xsize=8, ysize=6
  djs_plot, a.rp, a.ew_nofit_2_mc[iran], xtitle=xtitle, ytitle=ytitle, title=title, $
      xra=xra, yra=yra, xst=1, yst=1, thick=thick, xthick=xthick, ythick=ythick, $
      charsize=charsize, charthick=charthick, /nodata, /xlog, pos=pos
  ibegin=1
  iend=23 
  polyfill, [a[ibegin:iend].rp, reverse(a[ibegin:iend].rp)], [a[ibegin:iend].sdev_nofit_2_mc, -reverse(a[ibegin:iend].sdev_nofit_2_mc)], color=djs_icolor('light gray')
  djs_oplot, xx, 10.^(alog10(xx)*slope+intercept), thick=thick
  oploterror, a.rp, a.ew_nofit_2_mc[iran], a.sdev_nofit_2_mc, $ psym=8, symsize=2, color=djs_icolor('dark gray'), thick=thick, errcolor=djs_icolor('dark gray'), errthick=thick k_end_print

k_print, filename=psfile6, axis_char_scale=1.3, xsize=8, ysize=6
  xra=[220, 1500]
  yra=[-4E-3, 4E-3]
  djs_plot, a.rp, a.ew_nofit_2_mc[iran], xtitle=xtitle, ytitle=ytitle, title=title, $
      xra=xra, yra=yra, xst=1, yst=1, thick=thick, xthick=xthick, ythick=ythick, $
      charsize=charsize, charthick=charthick, /nodata, /xlog, pos=pos;, /ylog
  ibegin=13
  iend=23 
  polyfill, [a[ibegin:iend].rp, reverse(a[ibegin:iend].rp)], [a[ibegin:iend].sdev_nofit_2_mc, -reverse(a[ibegin:iend].sdev_nofit_2_mc)], color=djs_icolor('light gray')
  djs_oplot, xx, 10.^(alog10(xx)*slope+intercept), thick=thick
  oploterror, a.rp, a.ew_nofit_2_mc[iran], a.sdev_nofit_2_mc, $
      psym=8, symsize=2, color=djs_icolor('dark gray'), thick=thick, errcolor=djs_icolor('dark gray'), errthick=thick
k_end_print
endif

end
