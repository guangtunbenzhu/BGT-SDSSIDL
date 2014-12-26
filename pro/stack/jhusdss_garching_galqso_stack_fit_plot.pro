;Coarse = 0b
DoWeight = 0b
;DoNaI = 0b
DoScale = 0b
Scale_factor = 0.2 ;; M* \propto M(halo)^1 \propto r(halo)^3 and r(halo) \propto M*^(1/3); Scale to 10^10.2 Msun
;read,'Coarse? [1=yes, 0=no]: ', Coarse
read,'DoWeight? [1=yes, 0=no]: ', DoWeight
;read,'DoNaI? [1=yes, 0=no]: ', DoNaI
read,'DoScale? [1=yes, 0=no]: ', DoScale
;rvir = 200. ;; kpc for 10^10.3 (10^12)
if (DoScale) then rvir=200. else rvir=200.

;readcol, '/home/gz323/Code/jhu-sdss/pro/stack/caiipairs2.dat', rho, ewcaii, err_ewcaii, nlit, format='(f,f,f, i)'
readcol, '/home/gz323/Code/jhu-sdss/pro/stack/caiipairs2.dat', rho, ewcaii, err_ewcaii, format='(f,f,f)'

nmfver=106
stackpath = jhusdss_get_path(/nmfqso)+'/'+string(nmfver, format='(I3.3)')+'/Stack/'
infile0 = stackpath + jhusdss_garching_galqso_stack_filename(nmfver, boss=boss)
infile = infile0
if (DoWeight) then infile = repstr(infile, '.fits', '_weight.fits')
if (DoScale) then infile = repstr(infile, '.fits', '_Scale.fits')
infile = repstr(infile, '.fits', '_fit.fits')
a = mrdfits(infile, 1)

mcmax = 200
ew_nofit_mc = a.ew_nofit_mc[0:mcmax]
ew_nofit_2_mc = a.ew_nofit_2_mc[0:mcmax]
  for i=0,11 do begin
      shuffle_infile = repstr(infile0, '.fits', '_shuffle_'+string(i+1, format='(i2.2)')+'_fit.fits')
      ashuffle = mrdfits(shuffle_infile, 1)
      ew_nofit_mc = [ew_nofit_mc, ashuffle.ew_nofit_mc[0:mcmax]]
      ew_nofit_2_mc = [ew_nofit_2_mc, ashuffle.ew_nofit_2_mc[0:mcmax]]
  endfor

sdev_nofit_mc = fltarr(n_elements(a))
sdev_nofit_2_mc = fltarr(n_elements(a))

njump = n_elements(ew_nofit_mc[*,0])
;ijump = lindgen(njump/2)*2

for irp=0,n_elements(a)-1 do begin
    tmp = moment(ew_nofit_mc[*,irp], sdev=sdev)
;   tmp = moment(ew_nofit_mc[ijump,irp], sdev=sdev)
    sdev_nofit_mc[irp] = sdev
    tmp = moment(ew_nofit_2_mc[*,irp], sdev=sdev)
;   tmp = moment(ew_nofit_2_mc[ijump,irp], sdev=sdev)
    sdev_nofit_2_mc[irp] = sdev
endfor 
print, sdev_nofit_2_mc

psfile1 = repstr(infile, '.fits', '_1.ps')
psfile2 = repstr(infile, '.fits', '_2.ps')
psfile3 = repstr(infile, '.fits', '_3.ps')
psfile4 = repstr(infile, '.fits', '_4.ps')
psfile5 = repstr(infile, '.fits', '_5.ps')
psfile6 = repstr(infile, '.fits', '_6.ps')

;; manually check, independent bins
i_indep = [[0, 1], lindgen(11)*2+2]

;best_fit slope
ifit = i_indep[0:7]
;ifit = i_indep[0:11]

if (keyword_set(firstfit)) then begin
   xfit = alog10(a[ifit].rp/rvir)
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
slope = -1.4
intercept = 0.8

jhusdss_powerlaw_fit, a[ifit].rp/rvir*2., a[ifit].ew_nofit_2, sdev_nofit_2_mc[ifit], in_slope=slope, in_intercept=10.^intercept, $
;jhusdss_powerlaw_fit, a[ifit].rp, a[ifit].ew_nofit, sdev_nofit_mc[ifit], in_slope=slope, in_intercept=10.^intercept, $
   slope=out_slope, intercept=out_intercept, err_slope=err_slope, err_intercept=err_intercept
print, slope, intercept
print, out_slope, alog10(out_intercept)
slope = out_slope
intercept = alog10(out_intercept)
print, err_slope, out_intercept, err_intercept
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
;ytitle='W_0^{\lambda 3934} (Ca II) [\AA]'
ytitle='<W_0^{K}(Ca II)> (\AA)'
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
  ixx = where(xx lt 250., comp=jxx)
  djs_oplot, xx[ixx], 10.^(alog10(xx[ixx]/rvir*2.)*slope+intercept), thick=thick, color='blue'
  djs_oplot, xx[jxx], 10.^(alog10(xx[jxx]/rvir*2.)*slope+intercept), thick=thick, color='blue', linestyle=1
  oploterror, a[i_indep].rp, a[i_indep].ew_nofit_2, sdev_nofit_2_mc[i_indep], $
      psym=8, symsize=2, color=djs_icolor('blue'), thick=thick, errcolor=djs_icolor('blue'), errthick=thick

   ytitle1 = '<N(Ca II)> (cm^{-2})'
   ;; see Eq. 9.15 in Bruce Drain
   factor = 1.13E12*1E8/0.6267/(3934.78)^2
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
  djs_oplot, xx, 10.^(alog10(xx/rvir*2.)*slope+intercept), thick=thick, color='blue'
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
  djs_oplot, xx, 10.^(alog10(xx/rvir*2.)*slope+intercept), thick=thick, color='blue'
  oploterror, a[i_indep].rp, a[i_indep].ew_nofit, sdev_nofit_mc[i_indep], $
      psym=8, symsize=2, color=djs_icolor('blue'), thick=thick, errcolor=djs_icolor('blue'), errthick=thick
; djs_axis, xaxis=0, charsize=charsize, charthick=charthick, thick=thick, $
;     xtickv=[300], xticknames='300'
  djs_xyouts, 270, -0.0046, '300', charsize=charsize+0.4, charthick=charthick

k_end_print

xra=[3,2000]
yra=[5.E-5, 2.5E0]
xtitle='r_p (kpc)' 
;ytitle='W_0^{\lambda 3934} (Ca II) [\AA]'
ytitle='<W_0^{K}(Ca II)> (\AA)'
title='Double Gaussian Line Profile Measurement'

k_print, filename=psfile4, axis_char_scale=1.3, xsize=10, ysize=6
  djs_plot, a[i_indep].rp, a[i_indep].ew_nofit_2, xtitle=xtitle, ytitle=ytitle, $
      xra=xra, yra=yra, xst=5, yst=1, thick=thick, xthick=xthick, ythick=ythick, $
      charsize=charsize, charthick=charthick, /nodata, /xlog, /ylog, pos=pos, $
      ytickformat='jhusdss_tick_exponent', ystyle=9

  if (DoScale) then begin 
      djs_axis, xaxis=0, xtitle='r_p/r_{vir}', charsize=charsize, charthick=charthick, xra=xra/rvir, xthick=xthick
      djs_axis, xaxis=1, charsize=charsize, charthick=charthick, xra=xra/rvir, xthick=xthick, xtickformat='(A1)'
  endif else begin
      djs_axis, xaxis=0, xtitle=xtitle, charsize=charsize, charthick=charthick, xra=xra, xthick=xthick
      djs_axis, xaxis=1, charsize=charsize, charthick=charthick, xra=xra, xthick=xthick, xtickformat='(A1)'
  endelse

  rdisk = 8.0
; djs_xyouts, rdisk*0.94, 1.05, 'I', charsize=charsize+2, charthick=charthick+2, color='gray'
; djs_xyouts, rdisk*1.10, 1.3, '--', charsize=charsize, charthick=charthick+1, color='blue'
; djs_xyouts, rdisk*1.40, 1.3, '> Halo', charsize=charsize, charthick=charthick+1, color='blue'
; djs_xyouts, rdisk*0.61, 1.3, 'Disk', charsize=charsize, charthick=charthick+1, color='blue'

  djs_oplot, rdisk*[1,1], [6.5E-3,yra[1]], linestyle=2, color='gray', thick=thick+1
  djs_xyouts, rdisk*0.79, 4.0E-3, '<r_{90}>', charsize=charsize-0.3, charthick=charthick, color='gray'
; djs_xyouts, rdisk*1.20, 1.03E-2, '--', charsize=charsize, charthick=charthick+1, color='blue'
; djs_xyouts, rdisk*1.55, 1.03E-2, '>', charsize=charsize, charthick=charthick+1, color='blue'
  djs_xyouts, rdisk*1.25, 6.5E-3, 'Halo', charsize=charsize, charthick=charthick+1, color='gray'
; djs_xyouts, rdisk*0.60, 1.03E-2, '--', charsize=charsize, charthick=charthick+1, color='blue'
; djs_xyouts, rdisk*0.54, 1.03E-2, '<', charsize=charsize, charthick=charthick+1, color='blue'
  djs_xyouts, rdisk*0.52, 6.5E-3, 'Disk', charsize=charsize, charthick=charthick+1, color='gray'

; oploterror, rho, ewcaii, err_ewcaii, psym=5, symsize=1.3, color=djs_icolor('light red'), errcolor=djs_icolor('light red'), thick=thick-2, errthick=errthick

; plotsym, 7
; djs_oplot, rdisk*[1,1], 1.5*[1,1], psym=8, symsize=2, thick=thick
; plotsym, 0, /fill
; one_arrow,212,224,270,'S',charsize=3

  ixx = where(xx lt 250., comp=jxx)
  djs_oplot, xx[ixx], 10.^(alog10(xx[ixx]/rvir*2.)*slope+intercept), thick=thick, color='blue'
  djs_oplot, xx[jxx], 10.^(alog10(xx[jxx]/rvir*2.)*slope+intercept), thick=thick, color='blue', linestyle=1
  oploterror, a[i_indep].rp*1.08, a[i_indep].ew_nofit, sdev_nofit_mc[i_indep], $
      psym=4, symsize=2, color=djs_icolor('gray'), thick=thick, errcolor=djs_icolor('gray'), errthick=thick
  oploterror, a[i_indep].rp, a[i_indep].ew_nofit_2, sdev_nofit_2_mc[i_indep], $
      psym=8, symsize=2, color=djs_icolor('blue'), thick=thick, errcolor=djs_icolor('blue'), errthick=thick
  legend, ['Doublet profile fitting (K & H)', textoidl('Single Gaussian fitting for \lambda3934 (K)')], psym=[8,4], $
      colors=[djs_icolor('blue'), djs_icolor('gray')], $
      textcolors=[djs_icolor('blue'), djs_icolor('gray')], $
      box=0, /bottom, /left, charsize=charsize-0.3, charthick=charthick


;  text_relation = 'W_0=('+string(out_intercept*100, format='(f3.1)')+'\pm'+string(err_intercept*100,format='(f3.1)')+')\times10^{-2}\times(r_p/100 kpc)^{'+string(slope, format='(f5.2)')+'\pm'+string(err_slope,format='(f4.2)')+'}'
   text_relation = '<W_0>='+string(out_intercept, format='(f5.3)')+'\times(r_p/100 kpc)^{'+string(slope, format='(f5.2)')+'}'
;  if (DoScale) then text_relation = 'W_0=('+string(out_intercept*1000, format='(f3.1)')+'\pm'+string(err_intercept*1000,format='(f3.1)')+')\times10^{-3}\times(r_p/r_{vir})^{'+string(slope, format='(f5.2)')+'\pm'+string(err_slope,format='(f4.2)')+'}'
   if (DoScale) then text_relation = '<W_0>='+string(out_intercept, format='(f6.4)')+'\times(r_p/r_{vir})^{'+string(slope, format='(f5.2)')+'}'
   djs_xyouts, 6.0, 4.E-4, text_relation, charsize=charsize-0.20, charthick=charthick-0.3, color='blue'
   djs_oplot, [3.5, 5], 4.E-4*[1.2, 1.2], thick=thick, color='blue'

   ytitle1 = '<N(Ca II)> (cm^{-2})'
   ;; see Eq. 9.15 in Bruce Drain
   factor = 1.13E12*1E8/0.6267/(3934.78)^2
   djs_axis, yaxis=1, yra=yra*factor, /ylog, ythick=thick, ytitle=ytitle1, $
       charsiz=charsize, charthick=charthick

  djs_plot, alog10(xra), alog10(yra), /nodata, /noerase, xtickformat='(A1)', ytickformat='(A1)', charsize=charsize, charthick=charthick, pos=pos, xst=5, yst=5
  one_arrow, alog10(rdisk*1.2), alog10(1.5E-2), 0., ' ', arrowsize=[50,15,40], $
      charsize=charsize, thick=thick, color=djs_icolor('gray'), /data
  one_arrow, alog10(rdisk/1.2), alog10(1.5E-2), 180., ' ', arrowsize=[50,15,40], $
      charsize=charsize, thick=thick, color=djs_icolor('gray'), /data

  xra=[220,1250]
  yra=[-2.6E-3, 2.6E-3]
  ytickv=[textoidl('-2\times10^{-3}'), textoidl('-1\times10^{-3}'), '0', textoidl('1\times10^{-3}'), textoidl('2\times10^{-3}')]
  djs_plot, a.rp, a.ew_nofit_2, $
      xra=xra, yra=yra, xst=1, yst=1, thick=thick-3.8, xthick=xthick-3.8, ythick=ythick-3.8, $
      charsize=charsize-0.7, charthick=charthick-0.8, /nodata, /xlog, pos=[0.63, 0.63, 0.81, 0.81], /noerase, ytickformat='(A1)', yticklen=0.03
; djs_xyouts, 150, -2E-3, textoidl('-2\times10^{-3}'), charsize=charsize-0.7, charthick=charthick-0.8
  djs_xyouts, 133, -1.3E-3, textoidl('-10^{-3}'), charsize=charsize-0.3, charthick=charthick-0.8
  djs_xyouts, 170, -0.3E-3, textoidl('0'), charsize=charsize-0.3, charthick=charthick-0.8
  djs_xyouts, 150, 0.7E-3, textoidl('10^{-3}'), charsize=charsize-0.3, charthick=charthick-0.8
; djs_xyouts, 150, 2E-3, textoidl('2\times10^{-3}'), charsize=charsize-0.7, charthick=charthick-0.8
;     charsize=charsize-0.6, charthick=charthick-0.5, /nodata, /xlog, pos=[0.30, 0.22, 0.55, 0.47], /noerase
  ibegin=13
  iend=21 
; polyfill, [a[ibegin:iend].rp, reverse(a[ibegin:iend].rp)], [sdev_nofit_2_mc[ibegin:iend], -reverse(sdev_nofit_2_mc[ibegin:iend])], color=djs_icolor('light gray')
  djs_oplot, [10, 3000], [1., 1.]*0., linestyle=1, thick=thick-3.8
  djs_oplot, xx, 10.^(alog10(xx/rvir*2.)*slope+intercept), thick=thick-3, color='blue', linestyle=1
; djs_oplot, !x.crange, [0,0], color='light gray', thick=thick-3, linestyle=1
  oploterror, a[i_indep].rp*1.08, a[i_indep].ew_nofit, sdev_nofit_mc[i_indep], $
      psym=4, symsize=1.2, color=djs_icolor('gray'), thick=thick-2, errcolor=djs_icolor('gray'), errthick=thick-2
  oploterror, a[i_indep].rp, a[i_indep].ew_nofit_2, sdev_nofit_2_mc[i_indep], $
      psym=8, symsize=1.2, color=djs_icolor('blue'), thick=thick-2, errcolor=djs_icolor('blue'), errthick=thick-2
  djs_xyouts, 260, -0.0034, '300', charsize=charsize-0.50, charthick=charthick-0.5
  djs_xyouts, 410, -0.0038, xtitle, charsize=charsize-0.50, charthick=charthick-0.5


k_end_print

xra=[3,2000]
yra=[5.E-5, 2.5E0]
xtitle='r_p (kpc)' 
;ytitle='W_0^{\lambda 3934} (Ca II) [\AA]'
ytitle='<W_0^{K}(Ca II)> (\AA)'
title='Double Gaussian Line Profile Measurement'

k_print, filename=repstr(psfile4, '.ps', '_2.ps'), axis_char_scale=1.3, xsize=10, ysize=6
  djs_plot, a[i_indep].rp, a[i_indep].ew_nofit_2, xtitle=xtitle, ytitle=ytitle, $
      xra=xra, yra=yra, xst=5, yst=1, thick=thick, xthick=xthick, ythick=ythick, $
      charsize=charsize, charthick=charthick, /nodata, /xlog, /ylog, pos=pos, $
      ytickformat='jhusdss_tick_exponent', ystyle=9

  if (DoScale) then begin 
      djs_axis, xaxis=0, xtitle='r_p/r_{vir}', charsize=charsize, charthick=charthick, xra=xra/rvir, xthick=xthick
      djs_axis, xaxis=1, charsize=charsize, charthick=charthick, xra=xra/rvir, xthick=xthick, xtickformat='(A1)'
  endif else begin
      djs_axis, xaxis=0, xtitle=xtitle, charsize=charsize, charthick=charthick, xra=xra, xthick=xthick
      djs_axis, xaxis=1, charsize=charsize, charthick=charthick, xra=xra, xthick=xthick, xtickformat='(A1)'
  endelse

  rdisk = 8.0
; djs_xyouts, rdisk*0.94, 1.05, 'I', charsize=charsize+2, charthick=charthick+2, color='gray'
; djs_xyouts, rdisk*1.10, 1.3, '--', charsize=charsize, charthick=charthick+1, color='blue'
; djs_xyouts, rdisk*1.40, 1.3, '> Halo', charsize=charsize, charthick=charthick+1, color='blue'
; djs_xyouts, rdisk*0.61, 1.3, 'Disk', charsize=charsize, charthick=charthick+1, color='blue'

  djs_oplot, rdisk*[1,1], [6.5E-3,yra[1]], linestyle=2, color='gray', thick=thick+1
  djs_xyouts, rdisk*0.79, 4.0E-3, '<r_{90}>', charsize=charsize-0.3, charthick=charthick, color='gray'
; djs_xyouts, rdisk*1.20, 1.03E-2, '--', charsize=charsize, charthick=charthick+1, color='blue'
; djs_xyouts, rdisk*1.55, 1.03E-2, '>', charsize=charsize, charthick=charthick+1, color='blue'
  djs_xyouts, rdisk*1.25, 6.5E-3, 'Halo', charsize=charsize, charthick=charthick+1, color='gray'
; djs_xyouts, rdisk*0.60, 1.03E-2, '--', charsize=charsize, charthick=charthick+1, color='blue'
; djs_xyouts, rdisk*0.54, 1.03E-2, '<', charsize=charsize, charthick=charthick+1, color='blue'
  djs_xyouts, rdisk*0.52, 6.5E-3, 'Disk', charsize=charsize, charthick=charthick+1, color='gray'

; oploterror, rho, ewcaii, err_ewcaii, psym=5, symsize=1.3, color=djs_icolor('light red'), errcolor=djs_icolor('light red'), thick=thick-2, errthick=errthick

; plotsym, 7
; djs_oplot, rdisk*[1,1], 1.5*[1,1], psym=8, symsize=2, thick=thick
; plotsym, 0, /fill
; one_arrow,212,224,270,'S',charsize=3

  ixx = where(xx lt 250., comp=jxx)
  djs_oplot, xx[ixx], 10.^(alog10(xx[ixx]/rvir*2.)*slope+intercept), thick=thick, color='blue'
  djs_oplot, xx[jxx], 10.^(alog10(xx[jxx]/rvir*2.)*slope+intercept), thick=thick, color='blue', linestyle=1
; oploterror, a[i_indep].rp*1.08, a[i_indep].ew_nofit, sdev_nofit_mc[i_indep], $
;     psym=4, symsize=2, color=djs_icolor('gray'), thick=thick, errcolor=djs_icolor('gray'), errthick=thick
  oploterror, a[i_indep].rp, a[i_indep].ew_nofit_2, sdev_nofit_2_mc[i_indep], $
      psym=8, symsize=2, color=djs_icolor('blue'), thick=thick, errcolor=djs_icolor('blue'), errthick=thick
; legend, ['Doublet profile fitting (K & H)', textoidl('Single Gaussian fitting for \lambda3934 (K)')], psym=[8,4], $
;     colors=[djs_icolor('blue'), djs_icolor('gray')], $
;     textcolors=[djs_icolor('blue'), djs_icolor('gray')], $
;     box=0, /bottom, /left, charsize=charsize-0.3, charthick=charthick


;  text_relation = 'W_0=('+string(out_intercept*100, format='(f3.1)')+'\pm'+string(err_intercept*100,format='(f3.1)')+')\times10^{-2}\times(r_p/100 kpc)^{'+string(slope, format='(f5.2)')+'\pm'+string(err_slope,format='(f4.2)')+'}'
   text_relation = '<W_0>='+string(out_intercept, format='(f5.3)')+'\times(r_p/100 kpc)^{'+string(slope, format='(f5.2)')+'}'
;  if (DoScale) then text_relation = 'W_0=('+string(out_intercept*1000, format='(f3.1)')+'\pm'+string(err_intercept*1000,format='(f3.1)')+')\times10^{-3}\times(r_p/r_{vir})^{'+string(slope, format='(f5.2)')+'\pm'+string(err_slope,format='(f4.2)')+'}'
   if (DoScale) then text_relation = '<W_0>='+string(out_intercept, format='(f6.4)')+'\times(r_p/r_{vir})^{'+string(slope, format='(f5.2)')+'}'
   djs_xyouts, 6.0, 4.E-4, text_relation, charsize=charsize-0.20, charthick=charthick-0.3, color='blue'
   djs_oplot, [3.5, 5], 4.E-4*[1.2, 1.2], thick=thick, color='blue'

   ytitle1 = '<N(Ca II)> (cm^{-2})'
   ;; see Eq. 9.15 in Bruce Drain
   factor = 1.13E12*1E8/0.6267/(3934.78)^2
   djs_axis, yaxis=1, yra=yra*factor, /ylog, ythick=thick, ytitle=ytitle1, $
       charsiz=charsize, charthick=charthick

  djs_plot, alog10(xra), alog10(yra), /nodata, /noerase, xtickformat='(A1)', ytickformat='(A1)', charsize=charsize, charthick=charthick, pos=pos, xst=5, yst=5
  one_arrow, alog10(rdisk*1.2), alog10(1.5E-2), 0., ' ', arrowsize=[50,15,40], $
      charsize=charsize, thick=thick, color=djs_icolor('gray'), /data
  one_arrow, alog10(rdisk/1.2), alog10(1.5E-2), 180., ' ', arrowsize=[50,15,40], $
      charsize=charsize, thick=thick, color=djs_icolor('gray'), /data

  xra=[220,1250]
  yra=[-2.6E-3, 2.6E-3]
  ytickv=[textoidl('-2\times10^{-3}'), textoidl('-1\times10^{-3}'), '0', textoidl('1\times10^{-3}'), textoidl('2\times10^{-3}')]
  djs_plot, a.rp, a.ew_nofit_2, $
      xra=xra, yra=yra, xst=1, yst=1, thick=thick-3.8, xthick=xthick-3.8, ythick=ythick-3.8, $
      charsize=charsize-0.7, charthick=charthick-0.8, /nodata, /xlog, pos=[0.63, 0.63, 0.81, 0.81], /noerase, ytickformat='(A1)', yticklen=0.03
; djs_xyouts, 150, -2E-3, textoidl('-2\times10^{-3}'), charsize=charsize-0.7, charthick=charthick-0.8
  djs_xyouts, 133, -1.3E-3, textoidl('-10^{-3}'), charsize=charsize-0.3, charthick=charthick-0.8
  djs_xyouts, 170, -0.3E-3, textoidl('0'), charsize=charsize-0.3, charthick=charthick-0.8
  djs_xyouts, 150, 0.7E-3, textoidl('10^{-3}'), charsize=charsize-0.3, charthick=charthick-0.8
; djs_xyouts, 150, 2E-3, textoidl('2\times10^{-3}'), charsize=charsize-0.7, charthick=charthick-0.8
;     charsize=charsize-0.6, charthick=charthick-0.5, /nodata, /xlog, pos=[0.30, 0.22, 0.55, 0.47], /noerase
  ibegin=13
  iend=21 
; polyfill, [a[ibegin:iend].rp, reverse(a[ibegin:iend].rp)], [sdev_nofit_2_mc[ibegin:iend], -reverse(sdev_nofit_2_mc[ibegin:iend])], color=djs_icolor('light gray')
  djs_oplot, [10, 3000], [1., 1.]*0., linestyle=1, thick=thick-3.8
  djs_oplot, xx, 10.^(alog10(xx/rvir*2.)*slope+intercept), thick=thick-3, color='blue', linestyle=1
; djs_oplot, !x.crange, [0,0], color='light gray', thick=thick-3, linestyle=1
; oploterror, a[i_indep].rp*1.08, a[i_indep].ew_nofit, sdev_nofit_mc[i_indep], $
;     psym=4, symsize=1.2, color=djs_icolor('gray'), thick=thick-2, errcolor=djs_icolor('gray'), errthick=thick-2
  oploterror, a[i_indep].rp, a[i_indep].ew_nofit_2, sdev_nofit_2_mc[i_indep], $
      psym=8, symsize=1.2, color=djs_icolor('blue'), thick=thick-2, errcolor=djs_icolor('blue'), errthick=thick-2
  djs_xyouts, 260, -0.0034, '300', charsize=charsize-0.50, charthick=charthick-0.5
  djs_xyouts, 410, -0.0038, xtitle, charsize=charsize-0.50, charthick=charthick-0.5


k_end_print


xra=[2.5, 250]
yra=[5.E-3, 5.0E0]
xtitle='r_p (kpc)' 
ytitle='<W_0^{K}(Ca II)> (\AA)'
title='Double Gaussian Line Profile Measurement'

k_print, filename=psfile5, axis_char_scale=1.3, xsize=10, ysize=6
  djs_plot, a[i_indep].rp, a[i_indep].ew_nofit_2, xtitle=xtitle, ytitle=ytitle, $
      xra=xra, yra=yra, xst=1, yst=1, thick=thick, xthick=xthick, ythick=ythick, $
      charsize=charsize, charthick=charthick+1, /nodata, /xlog, /ylog, $
      pos=[0.22, 0.17, 0.9, 0.85], $
      ytickformat='jhusdss_tick_exponent', title='Individual Ca II absorbers from the literature'

; rdisk = 8.0
; djs_xyouts, rdisk*0.94, 1.05, 'I', charsize=charsize+2, charthick=charthick+2, color='gray'
; djs_xyouts, rdisk*1.10, 1.3, '--', charsize=charsize, charthick=charthick+1, color='blue'
; djs_xyouts, rdisk*1.40, 1.3, '> Halo', charsize=charsize, charthick=charthick+1, color='blue'
; djs_xyouts, rdisk*0.61, 1.3, 'Disk', charsize=charsize, charthick=charthick+1, color='blue'

; djs_oplot, rdisk*[1,1], [6.5E-3,yra[1]], linestyle=2, color='gray', thick=thick+1
; djs_xyouts, rdisk*1.20, 1.03E-2, '--', charsize=charsize, charthick=charthick+1, color='blue'
; djs_xyouts, rdisk*1.55, 1.03E-2, '>', charsize=charsize, charthick=charthick+1, color='blue'
; djs_xyouts, rdisk*1.25, 6.5E-3, 'Halo', charsize=charsize, charthick=charthick+1, color='gray'
; djs_xyouts, rdisk*0.60, 1.03E-2, '--', charsize=charsize, charthick=charthick+1, color='blue'
; djs_xyouts, rdisk*0.54, 1.03E-2, '<', charsize=charsize, charthick=charthick+1, color='blue'
; djs_xyouts, rdisk*0.52, 6.5E-3, 'Disk', charsize=charsize, charthick=charthick+1, color='gray'

; djs_xyouts, 20, 1E0, 'Individual Ca II Absorbers from the literature', charsize=charsize, charthick=charthick
  ;oploterror, rho[0:5], ewcaii[0:5], err_ewcaii[0:5], psym=4, symsize=1.8, color=djs_icolor('black'), errcolor=djs_icolor('black'), thick=thick, errthick=thick
  ;oploterror, rho[6:9], ewcaii[6:9], err_ewcaii[6:9], psym=6, symsize=1.8, color=djs_icolor('magenta'), errcolor=djs_icolor('magenta'), thick=thick, errthick=thick
  ;djs_oplot, rho[10:16], ewcaii[10:16], psym=5, symsize=1.8, color='dark green', errcolor='dark green', thick=thick, symthick=symthick
  oploterror, rho[0:19], ewcaii[0:19], err_ewcaii[0:19], psym=4, symsize=1.8, color=djs_icolor('dark green'), errcolor=djs_icolor('dark green'), thick=thick, errthick=thick
  djs_oplot, rho[20:26], ewcaii[20:26], psym=4, symsize=1.8, color=djs_icolor('dark green'), thick=thick
; djs_oplot, xx, 10.^(alog10(xx/rvir)*slope+intercept), thick=thick, color='blue'

k_end_print


if (keyword_set(movingon)) then begin
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
  djs_oplot, xx, 10.^(alog10(xx/rvir*2.)*slope+intercept), thick=thick, color='blue'
  oploterror, a[i_indep].rp, a[i_indep].ew_nofit_2, sdev_nofit_2_mc[i_indep], $
      psym=8, symsize=2, color=djs_icolor('blue'), thick=thick, errcolor=djs_icolor('blue'), errthick=thick
k_end_print
endif

k_print, filename=psfile6, axis_char_scale=1.3, xsize=8, ysize=6
  xra=[220, 2000]
  yra=[-4E-3, 4E-3]
  djs_plot, a.rp, a.ew_nofit_2, xtitle=xtitle, ytitle=ytitle, title=title, $
      xra=xra, yra=yra, xst=1, yst=1, thick=thick, xthick=xthick, ythick=ythick, $
      charsize=charsize, charthick=charthick, /nodata, /xlog, pos=pos;, /ylog
  ibegin=13
  iend=23 
  polyfill, [a[ibegin:iend].rp, reverse(a[ibegin:iend].rp)], [sdev_nofit_2_mc[ibegin:iend], -reverse(sdev_nofit_2_mc[ibegin:iend])], color=djs_icolor('light gray')
  djs_oplot, xx, 10.^(alog10(xx/rvir*2.)*slope+intercept), thick=thick
  djs_oplot, !x.crange, [0,0], color='light gray', thick=thick, linestyle=1
  oploterror, a[i_indep].rp, a[i_indep].ew_nofit_2, sdev_nofit_2_mc[i_indep], $
      psym=8, symsize=2, color=djs_icolor('blue'), thick=thick, errcolor=djs_icolor('blue'), errthick=thick
  djs_xyouts, 270, -0.0046, '300', charsize=charsize+0.4, charthick=charthick
k_end_print

;; total Ca II mass between rpmin and rpmax
rpmin=10.
rpmax=200.
kpc = 3.08D21 ;cm
factor = 1.13E12*1E8/0.6485/(3934.78)^2 ; turn REW to column density
Mca40 = 40.*1.67D-27
Msolar = 2.D30
;; Int (factor*out_intercept*(r/rvir)^slope*2.*pi*rdr) = factor*out_intercept*2.*pi*(1/rvir)^slope*r^(2.+slope)
number_caii = factor*out_intercept*2.*!pi*(1./rvir*2.)^slope*kpc^2*(rpmax^(2.+slope)-rpmin^(2.+slope))/(2.+slope)
err_number_caii = factor*err_intercept*2.*!pi*(1./rvir*2.)^slope*kpc^2*(rpmax^(2.+slope)-rpmin^(2.+slope))/(2.+slope)
mass_caii = number_caii*Mca40/Msolar
err_mass_caii = err_number_caii*Mca40/Msolar
print, mass_caii, err_mass_caii
;; total HI mass
;; solar: N(H)/N(Ca) = 10^5.66 -> M(H)/M(Ca) = 10^5.66/40
mass_HI_solar = 10.D0^5.66/40.*mass_caii
print, mass_hi_solar
;; MW: log10[N(CaII)] = +0.22log10[N(HI)]+7.45
;;     log10[N(HI)] = (log10[N(CaII)]-7.45)/0.22
;;     N(HI) = N(CaII)^(1./0.22)*10.^(-7.45/0.22)

end
