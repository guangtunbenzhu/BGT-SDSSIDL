lrgver=101
stackpath = jhusdss_get_path(/fitlrg)+'/'+string(lrgver, format='(I3.3)')+'/Stack/'
infile = stackpath + jhusdss_garching_gallrg_stack_filename(lrgver, boss=boss)
infile = repstr(infile, '.fits', '_fit.fits')
a = mrdfits(infile, 1)

psfile1 = repstr(infile, '.fits', '_1.ps')
psfile2 = repstr(infile, '.fits', '_2.ps')
psfile3 = repstr(infile, '.fits', '_3.ps')
psfile4 = repstr(infile, '.fits', '_4.ps')
psfile5 = repstr(infile, '.fits', '_5.ps')
psfile6 = repstr(infile, '.fits', '_6.ps')

;; manually check, independent bins
i_indep = [[1], lindgen(12)*2+2]

;best_fit slope
slope = -1.29629
intercept = 0.636813

thick=8
xthick=8
ythick=8
charsize=1.4
charthick=3

xra=[10,500]
yra=[1.E-3, 2.E0]
xtitle='r_p (Kpc)' 
ytitle='W_0^{\lambda 3934} (\AA)'
title='Single Gaussian Line Profile Measurement'

xx = 10.^(findgen(1000)*0.0024+1.)

plotsym, 0, /fill

pos = [0.20, 0.15, 0.95, 0.9]
k_print, filename=psfile1, axis_char_scale=1.3, xsize=8, ysize=6
  djs_plot, a[i_indep].rp, a[i_indep].ew_nofit, xtitle=xtitle, ytitle=ytitle, title=title, $
      xra=xra, yra=yra, xst=1, yst=1, thick=thick, xthick=xthick, ythick=ythick, $
      charsize=charsize, charthick=charthick, /nodata, /xlog, /ylog, pos=pos
  djs_oplot, xx, 10.^(alog10(xx)*slope+intercept), thick=thick
  oploterror, a[i_indep].rp, a[i_indep].ew_nofit, a[i_indep].sdev_nofit_mc, $
      psym=8, symsize=2, color=djs_icolor('blue'), thick=thick, errcolor=djs_icolor('blue'), errthick=thick
k_end_print

k_print, filename=psfile2, axis_char_scale=1.3, xsize=8, ysize=6
  djs_plot, a.rp, a.ew_nofit, xtitle=xtitle, ytitle=ytitle, title=title, $
      xra=xra, yra=yra, xst=1, yst=1, thick=thick, xthick=xthick, ythick=ythick, $
      charsize=charsize, charthick=charthick, /nodata, /xlog, /ylog, pos=pos
  ibegin=1
  iend=23 
  polyfill, [a[ibegin:iend].rp, reverse(a[ibegin:iend].rp)], [a[ibegin:iend].sdev_nofit_mc, -reverse(a[ibegin:iend].sdev_nofit_mc)], color=djs_icolor('light gray')
  djs_oplot, xx, 10.^(alog10(xx)*slope+intercept), thick=thick
  oploterror, a.rp, a.ew_nofit, a.sdev_nofit_mc, $
      psym=8, symsize=2, color=djs_icolor('blue'), thick=thick, errcolor=djs_icolor('blue'), errthick=thick
k_end_print

k_print, filename=psfile3, axis_char_scale=1.3, xsize=8, ysize=6
  xra=[220, 2000]
  yra=[-4E-3, 4E-3]
  djs_plot, a.rp, a.ew_nofit, xtitle=xtitle, ytitle=ytitle, title=title, $
      xra=xra, yra=yra, xst=1, yst=1, thick=thick, xthick=xthick, ythick=ythick, $
      charsize=charsize, charthick=charthick, /nodata, /xlog, pos=pos;, /ylog
  ibegin=13
  iend=23 
  polyfill, [a[ibegin:iend].rp, reverse(a[ibegin:iend].rp)], [a[ibegin:iend].sdev_nofit_mc, -reverse(a[ibegin:iend].sdev_nofit_mc)], color=djs_icolor('light gray')
  djs_oplot, xx, 10.^(alog10(xx)*slope+intercept), thick=thick
  oploterror, a.rp, a.ew_nofit, a.sdev_nofit_mc, $
      psym=8, symsize=2, color=djs_icolor('blue'), thick=thick, errcolor=djs_icolor('blue'), errthick=thick

k_end_print

xra=[10,500]
yra=[1.E-3, 2.E-0]
xtitle='r_p (Kpc)' 
ytitle='W_0^{\lambda 3934} (\AA)'
title='Double Gaussian Line Profile Measurement'

k_print, filename=psfile4, axis_char_scale=1.3, xsize=8, ysize=6
  djs_plot, a[i_indep].rp, a[i_indep].ew_nofit_2, xtitle=xtitle, ytitle=ytitle, title=title, $
      xra=xra, yra=yra, xst=1, yst=1, thick=thick, xthick=xthick, ythick=ythick, $
      charsize=charsize, charthick=charthick, /nodata, /xlog, /ylog, pos=pos
  djs_oplot, xx, 10.^(alog10(xx)*slope+intercept), thick=thick
  oploterror, a[i_indep].rp, a[i_indep].ew_nofit_2, a[i_indep].sdev_nofit_2_mc, $
      psym=8, symsize=2, color=djs_icolor('blue'), thick=thick, errcolor=djs_icolor('blue'), errthick=thick
k_end_print

k_print, filename=psfile5, axis_char_scale=1.3, xsize=8, ysize=6
  djs_plot, a.rp, a.ew_nofit_2, xtitle=xtitle, ytitle=ytitle, title=title, $
      xra=xra, yra=yra, xst=1, yst=1, thick=thick, xthick=xthick, ythick=ythick, $
      charsize=charsize, charthick=charthick, /nodata, /xlog, /ylog, pos=pos
  ibegin=1
  iend=23 
  polyfill, [a[ibegin:iend].rp, reverse(a[ibegin:iend].rp)], [a[ibegin:iend].sdev_nofit_2_mc, -reverse(a[ibegin:iend].sdev_nofit_2_mc)], color=djs_icolor('light gray')
  djs_oplot, xx, 10.^(alog10(xx)*slope+intercept), thick=thick
  oploterror, a.rp, a.ew_nofit_2, a.sdev_nofit_2_mc, $
      psym=8, symsize=2, color=djs_icolor('blue'), thick=thick, errcolor=djs_icolor('blue'), errthick=thick
k_end_print

k_print, filename=psfile6, axis_char_scale=1.3, xsize=8, ysize=6
  xra=[220, 2000]
  yra=[-4E-3, 4E-3]
  djs_plot, a.rp, a.ew_nofit_2, xtitle=xtitle, ytitle=ytitle, title=title, $
      xra=xra, yra=yra, xst=1, yst=1, thick=thick, xthick=xthick, ythick=ythick, $
      charsize=charsize, charthick=charthick, /nodata, /xlog, pos=pos;, /ylog
  ibegin=13
  iend=23 
  polyfill, [a[ibegin:iend].rp, reverse(a[ibegin:iend].rp)], [a[ibegin:iend].sdev_nofit_2_mc, -reverse(a[ibegin:iend].sdev_nofit_2_mc)], color=djs_icolor('light gray')
  djs_oplot, xx, 10.^(alog10(xx)*slope+intercept), thick=thick
  oploterror, a.rp, a.ew_nofit_2, a.sdev_nofit_2_mc, $
      psym=8, symsize=2, color=djs_icolor('blue'), thick=thick, errcolor=djs_icolor('blue'), errthick=thick

k_end_print

end
