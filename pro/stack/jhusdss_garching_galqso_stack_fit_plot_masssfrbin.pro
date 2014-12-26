nmfver=106
stackpath = jhusdss_get_path(/nmfqso)+'/'+string(nmfver, format='(I3.3)')+'/Stack/'
massive = 0b
quiescent = 0b

infile0 = stackpath + jhusdss_garching_galqso_stack_filename(nmfver, boss=boss)
infile = infile0
if (massive) then infile = repstr(infile, '.fits', '_highmass.fits') else infile = repstr(infile, '.fits', '_lowmass.fits')
if (quiescent) then infile = repstr(infile, '.fits', '_quiescent.fits') else infile = repstr(infile, '.fits', '_sf.fits')

infile = repstr(infile, '.fits', '_fit.fits')
a = mrdfits(infile, 1)

mcmax = 200
ew_nofit_mc = a.ew_nofit_mc[0:mcmax]
ew_nofit_2_mc = a.ew_nofit_2_mc[0:mcmax]
  for i=0,11 do begin
      shuffle_infile = repstr(infile0, '.fits', '_shuffle_'+string(i+1, format='(i2.2)')+'_fit.fits')
      ashuffle = mrdfits(shuffle_infile, 1)
      ew_nofit_mc = [ew_nofit_mc, ashuffle[1:24].ew_nofit_mc[0:mcmax]*2.]
      ew_nofit_2_mc = [ew_nofit_2_mc, ashuffle[1:24].ew_nofit_2_mc[0:mcmax]*2.]
  endfor

sdev_nofit_mc = fltarr(n_elements(a))
sdev_nofit_2_mc = fltarr(n_elements(a))
for irp=0,n_elements(a)-1 do begin
    tmp = moment(ew_nofit_mc[*,irp], sdev=sdev)
    sdev_nofit_mc[irp] = sdev
    tmp = moment(ew_nofit_2_mc[*,irp], sdev=sdev)
    sdev_nofit_2_mc[irp] = sdev
endfor 
;print, sdev_nofit_2_mc

psfile1 = repstr(infile, '.fits', '_1.ps')
psfile2 = repstr(infile, '.fits', '_2.ps')
psfile3 = repstr(infile, '.fits', '_3.ps')
psfile4 = repstr(infile, '.fits', '_4.ps')
psfile5 = repstr(infile, '.fits', '_5.ps')
psfile6 = repstr(infile, '.fits', '_6.ps')

;; manually check, independent bins
i_indep = [[0], lindgen(12)*2+1]

;best_fit slope
ifit = i_indep[0:4]
;ifit = i_indep[0:11]

if (keyword_set(firstfit)) then begin
   xfit = alog10(a[ifit].rp/100.)
   ;yfit = alog10(a[ifit].ew_nofit_2)
   ;yerror = 0.5*(alog10(a[ifit].ew_nofit_2+a[ifit].sdev_nofit_2_mc)-yfit)+0.5*(yfit-(alog10((a[ifit].ew_nofit_2-a[ifit].sdev_nofit_2_mc)>1.E-10)))
   yfit = alog10(a[ifit].ew_nofit_2)
   ;yerror = 0.5*(alog10(a[ifit].ew_nofit+a[ifit].sdev_nofit_mc)-yfit)+0.5*(yfit-(alog10((a[ifit].ew_nofit-a[ifit].sdev_nofit_mc)>1.E-10)))
   yerror = (alog10(a[ifit].ew_nofit+a[ifit].sdev_nofit_mc)-yfit);+0.5*(yfit-(alog10((a[ifit].ew_nofit-a[ifit].sdev_nofit_mc)>1.E-10)))
   coeffs = linfit(xfit, yfit, measure_error=yerror)
   ;slope = -1.25
   ;intercept = 0.52
   slope = coeffs[1]
   intercept = coeffs[0]
endif
slope = -1.38
intercept = 0.01

jhusdss_powerlaw_fit, a[ifit].rp/100., a[ifit].ew_nofit_2, sdev_nofit_2_mc[ifit], in_slope=slope, in_intercept=10.^intercept, $
;jhusdss_powerlaw_fit, a[ifit].rp, a[ifit].ew_nofit, sdev_nofit_mc[ifit], in_slope=slope, in_intercept=10.^intercept, $
   slope=out_slope, intercept=out_intercept, err_slope=err_slope, err_intercept=err_intercept
print, slope, intercept
print, out_slope, alog10(out_intercept)
slope = out_slope
intercept = alog10(out_intercept)
print, err_slope, out_intercept, err_intercept
fid_slope = -1.38
fid_intercept = alog10(0.01)
;stop

thick=8
xthick=8
ythick=8
charsize=1.4
charthick=3

;xra=[5,300]
xra=[5,1000]
yra=[1.E-3, 1.E-0]
xtitle='r_p (kpc)' 
ytitle='W_0^{\lambda 3934} (Ca II) [\AA]'
title='Single Gaussian Line Profile Measurement'

xx = 10.^(findgen(1500)*0.0024+0.1)
print, minmax(xx)

plotsym, 0, /fill

pos = [0.17, 0.20, 0.85, 0.90]
k_print, filename=psfile1, axis_char_scale=1.3, xsize=8, ysize=5
  djs_plot, a[i_indep].rp, a[i_indep].ew_nofit_2, xtitle=xtitle, ytitle=ytitle, $
      xra=xra, yra=yra, xst=1, yst=1, thick=thick, xthick=xthick, ythick=ythick, $
      charsize=charsize, charthick=charthick, /nodata, /xlog, /ylog, pos=pos, $
      ytickformat='jhusdss_tick_exponent'
  djs_oplot, xx, 10.^(alog10(xx/100.)*slope+intercept), thick=thick, color='blue'
  oploterror, a[i_indep].rp, a[i_indep].ew_nofit_2, sdev_nofit_2_mc[i_indep], $
      psym=8, symsize=2, color=djs_icolor('blue'), thick=thick, errcolor=djs_icolor('blue'), errthick=thick

   ytitle1 = 'N (Ca II) [cm^{-2}]'
   ;; see Eq. 9.15 in Bruce Drain
   factor = 1.13E12*1E8/0.6480/(3934.78)^2
   djs_axis, yaxis=1, yra=yra*factor, /ylog, ythick=thick, ytitle=ytitle1, $
       charsiz=charsize, charthick=charthick

    djs_xyouts, 6, 0.005, 'Average Ca II Column Density', charsize=charsize, charthick=charthick+1
    djs_xyouts, 6, 0.003, 'Profile in the CGM', charsize=charsize, charthick=charthick+1
    djs_xyouts, 170, 0.000475, '200', charsize=charsize+0.35, charthick=charthick

k_end_print

xra=[5,2000]
yra=[5.E-5, 6.E-1]
k_print, filename=psfile2, axis_char_scale=1.3, xsize=8, ysize=6
  djs_plot, a.rp, a.ew_nofit, xtitle=xtitle, ytitle=ytitle, title=title, $
      xra=xra, yra=yra, xst=1, yst=1, thick=thick, xthick=xthick, ythick=ythick, $
      charsize=charsize, charthick=charthick, /nodata, /xlog, /ylog, pos=pos
  ibegin=1
  iend=23 
; polyfill, [a[ibegin:iend].rp, reverse(a[ibegin:iend].rp)], [a[ibegin:iend].sdev_nofit_mc, -reverse(a[ibegin:iend].sdev_nofit_mc)], color=djs_icolor('light gray')
  djs_oplot, a[i_indep].rp, sdev_nofit_mc[i_indep], linestyle=1, thick=thick, color='gray'
  djs_oplot, xx, 10.^(alog10(xx/100.)*slope+intercept), thick=thick, color='blue'
  oploterror, a[i_indep].rp, a[i_indep].ew_nofit, sdev_nofit_mc[i_indep], $
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
  polyfill, [a[ibegin:iend].rp, reverse(a[ibegin:iend].rp)], [sdev_nofit_mc[ibegin:iend], -reverse(sdev_nofit_mc[ibegin:iend])], color=djs_icolor('light gray')
  djs_oplot, xx, 10.^(alog10(xx/100.)*slope+intercept), thick=thick, color='blue'
  oploterror, a[i_indep].rp, a[i_indep].ew_nofit, sdev_nofit_mc[i_indep], $
      psym=8, symsize=2, color=djs_icolor('blue'), thick=thick, errcolor=djs_icolor('blue'), errthick=thick
; djs_axis, xaxis=0, charsize=charsize, charthick=charthick, thick=thick, $
;     xtickv=[300], xticknames='300'
  djs_xyouts, 270, -0.0046, '300', charsize=charsize+0.4, charthick=charthick

k_end_print

xra=[5,2000]
yra=[5.E-5, 1.E0]
xtitle='r_p (kpc)' 
ytitle='W_0^{\lambda 3934} (Ca II) [\AA]'
title='Double Gaussian Line Profile Measurement'

k_print, filename=psfile4, axis_char_scale=1.3, xsize=10, ysize=6
  djs_plot, a[i_indep].rp, a[i_indep].ew_nofit_2, xtitle=xtitle, ytitle=ytitle, $
      xra=xra, yra=yra, xst=1, yst=1, thick=thick, xthick=xthick, ythick=ythick, $
      charsize=charsize, charthick=charthick, /nodata, /xlog, /ylog, pos=pos, $
      ytickformat='jhusdss_tick_exponent', ystyle=9
  djs_oplot, xx, 10.^(alog10(xx/100.)*slope+intercept), thick=thick, color='red', linestyle=1
  djs_oplot, xx, 10.^(alog10(xx/100.)*fid_slope+fid_intercept), thick=thick+5, color='blue', linestyle=0
  oploterror, a[i_indep].rp*1.08, a[i_indep].ew_nofit, sdev_nofit_mc[i_indep], $
      psym=4, symsize=2, color=djs_icolor('gray'), thick=thick, errcolor=djs_icolor('gray'), errthick=thick
  oploterror, a[i_indep].rp, a[i_indep].ew_nofit_2, sdev_nofit_2_mc[i_indep], $
      psym=8, symsize=2, color=djs_icolor('red'), thick=thick, errcolor=djs_icolor('red'), errthick=thick
  legend, ['Double Gaussian fitting', 'Single Gaussian fitting'], psym=[8,4], $
      colors=[djs_icolor('red'), djs_icolor('gray')], $
      textcolors=[djs_icolor('red'), djs_icolor('gray')], $
      box=0, /bottom, /left, charsize=charsize-0.3, charthick=charthick

   text_relation = 'W_0=('+string(out_intercept*100, format='(f3.1)')+'\pm'+string(err_intercept*100,format='(f3.1)')+')\times10^{-2}\times(r_p/100 kpc)^{'+string(slope, format='(f5.2)')+'\pm'+string(err_slope,format='(f5.2)')+'}'
   djs_xyouts, 6.5, 4.E-4, text_relation, charsize=charsize-0.20, charthick=charthick-0.3, color='red'
   if (massive) then begin
      if (quiescent) then text_massive_quiescent = 'massive quiescent'
      if (~quiescent) then text_massive_quiescent = 'massive star-forming'
   endif
   if (~massive) then begin
      if (quiescent) then text_massive_quiescent = 'lowmass quiescent'
      if (~quiescent) then text_massive_quiescent = 'lowmass star-forming'
   endif
   djs_xyouts, 6.5, 6.E-4, text_massive_quiescent, charsize=charsize-0.20, charthick=charthick-0.3, color='red'

   ytitle1 = 'N (Ca II) [cm^{-2}]'
   ;; see Eq. 9.15 in Bruce Drain
   factor = 1.13E12*1E8/0.6480/(3934.78)^2
   djs_axis, yaxis=1, yra=yra*factor, /ylog, ythick=thick, ytitle=ytitle1, $
       charsiz=charsize, charthick=charthick


  xra=[220,1250]
  yra=[-2.6E-3, 2.6E-3]
  djs_plot, a.rp, a.ew_nofit_2, $
      xra=xra, yra=yra, xst=1, yst=1, thick=thick-3.8, xthick=xthick-3.8, ythick=ythick-3.8, $
      charsize=charsize-0.7, charthick=charthick-0.8, /nodata, /xlog, pos=[0.61, 0.61, 0.81, 0.81], /noerase
;     charsize=charsize-0.6, charthick=charthick-0.5, /nodata, /xlog, pos=[0.30, 0.22, 0.55, 0.47], /noerase
  ibegin=12
  iend=20 
  polyfill, [a[ibegin:iend].rp, reverse(a[ibegin:iend].rp)], [sdev_nofit_2_mc[ibegin:iend], -reverse(sdev_nofit_2_mc[ibegin:iend])], color=djs_icolor('light gray')
  djs_oplot, xx, 10.^(alog10(xx/100.)*slope+intercept), thick=thick-3, color='red'
; djs_oplot, !x.crange, [0,0], color='light gray', thick=thick-3, linestyle=1
  oploterror, a[i_indep].rp*1.08, a[i_indep].ew_nofit, sdev_nofit_mc[i_indep], $
      psym=4, symsize=1.2, color=djs_icolor('gray'), thick=thick-2, errcolor=djs_icolor('gray'), errthick=thick-2
  oploterror, a[i_indep].rp, a[i_indep].ew_nofit_2, sdev_nofit_2_mc[i_indep], $
      psym=8, symsize=1.2, color=djs_icolor('red'), thick=thick-2, errcolor=djs_icolor('red'), errthick=thick-2
  djs_xyouts, 260, -0.0034, '300', charsize=charsize-0.50, charthick=charthick-0.5
  djs_xyouts, 410, -0.0038, xtitle, charsize=charsize-0.50, charthick=charthick-0.5
k_end_print

xra=[5,2000]
yra=[5.E-5, 6.E-1]
k_print, filename=psfile5, axis_char_scale=1.3, xsize=8, ysize=6
  djs_plot, a.rp, a.ew_nofit_2, xtitle=xtitle, ytitle=ytitle, title=title, $
      xra=xra, yra=yra, xst=1, yst=1, thick=thick, xthick=xthick, ythick=ythick, $
      charsize=charsize, charthick=charthick, /nodata, /xlog, /ylog, pos=pos
  ibegin=1
  iend=23 
; polyfill, [a[ibegin:iend].rp, reverse(a[ibegin:iend].rp)], [a[ibegin:iend].sdev_nofit_2_mc, -reverse(a[ibegin:iend].sdev_nofit_2_mc)], color=djs_icolor('light gray')
  djs_oplot, a[i_indep].rp, sdev_nofit_2_mc[i_indep], linestyle=1, thick=thick, color='gray'
  djs_oplot, xx, 10.^(alog10(xx/100.)*slope+intercept), thick=thick, color='blue'
  oploterror, a[i_indep].rp, a[i_indep].ew_nofit_2, sdev_nofit_2_mc[i_indep], $
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
  polyfill, [a[ibegin:iend].rp, reverse(a[ibegin:iend].rp)], [sdev_nofit_2_mc[ibegin:iend], -reverse(sdev_nofit_2_mc[ibegin:iend])], color=djs_icolor('light gray')
  djs_oplot, xx, 10.^(alog10(xx/100.)*slope+intercept), thick=thick
  djs_oplot, !x.crange, [0,0], color='light gray', thick=thick, linestyle=1
  oploterror, a[i_indep].rp, a[i_indep].ew_nofit_2, sdev_nofit_2_mc[i_indep], $
      psym=8, symsize=2, color=djs_icolor('blue'), thick=thick, errcolor=djs_icolor('blue'), errthick=thick
  djs_xyouts, 270, -0.0046, '300', charsize=charsize+0.4, charthick=charthick
k_end_print

end
