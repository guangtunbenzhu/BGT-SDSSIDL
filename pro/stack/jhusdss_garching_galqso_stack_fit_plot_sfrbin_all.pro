minrad_tmp = 10.^(alog10(0.010)+findgen(13)*(0.5*alog10(2.0)))
maxrad_tmp = minrad_tmp*2.0
minrad = minrad_tmp
maxrad = maxrad_tmp
;minrad = [[0.010], minrad_tmp]
;maxrad = [[0.020], maxrad_tmp]

;Coarse = 0b
DoWeight = 0b
;DoNaI = 0b
DoScale = 0b
Scale_factor = 0.2 ;; M* \propto M(halo)^1 \propto r(halo)^3 and r(halo) \propto M*^(1/3); Scale to 10^10.2 Msun
;read,'Coarse? [1=yes, 0=no]: ', Coarse
read,'DoWeight? [1=yes, 0=no]: ', DoWeight
;read,'DoNaI? [1=yes, 0=no]: ', DoNaI
read,'DoScale? [1=yes, 0=no]: ', DoScale
;rvir = 228. ;; kpc for 10^10.3 (10^12)
if (DoScale) then rvir=200. else rvir=200.

nmfver=106
stackpath = jhusdss_get_path(/nmfqso)+'/'+string(nmfver, format='(I3.3)')+'/Stack/'

infile0 = stackpath + jhusdss_garching_galqso_stack_filename(nmfver, boss=boss)

quiescent = 1b 
infile = infile0
if (quiescent) then infile = repstr(infile, '.fits', '_quiescent.fits') else infile = repstr(infile, '.fits', '_sf.fits')
if (DoScale) then infile=repstr(infile,'.fits','_Scale.fits')
infile = repstr(infile, '.fits', '_fit.fits')
a_quiescent= mrdfits(infile, 1)

quiescent = 0b 
infile = infile0
if (quiescent) then infile = repstr(infile, '.fits', '_quiescent.fits') else infile = repstr(infile, '.fits', '_sf.fits')
if (DoScale) then infile=repstr(infile,'.fits','_Scale.fits')
infile = repstr(infile, '.fits', '_fit.fits')
a_sf= mrdfits(infile, 1)

mcmax = 200
ew_nofit_mc = [a_quiescent.ew_nofit_mc[0:mcmax], a_sf.ew_nofit_mc[0:mcmax]]
ew_nofit_2_mc = [a_quiescent.ew_nofit_2_mc[0:mcmax], a_sf.ew_nofit_2_mc[0:mcmax]]

; for i=0,11 do begin
;     shuffle_infile = repstr(infile0, '.fits', '_shuffle_'+string(i+1, format='(i2.2)')+'_fit.fits')
;     ashuffle = mrdfits(shuffle_infile, 1)
;     ew_nofit_mc = [ew_nofit_mc, ashuffle[1:24].ew_nofit_mc[0:mcmax]*sqrt(2.)]
;     ew_nofit_2_mc = [ew_nofit_2_mc, ashuffle[1:24].ew_nofit_2_mc[0:mcmax]*sqrt(2.)]
; endfor

sdev_nofit_mc = fltarr(n_elements(a_quiescent))
sdev_nofit_2_mc = fltarr(n_elements(a_quiescent))
for irp=0,n_elements(a_quiescent)-1 do begin
    tmp = moment(ew_nofit_mc[*,irp], sdev=sdev)
    sdev_nofit_mc[irp] = sdev
    tmp = moment(ew_nofit_2_mc[*,irp], sdev=sdev)
    sdev_nofit_2_mc[irp] = sdev
endfor 
;print, sdev_nofit_2_mc

infile1 = infile0
if (DoScale) then infile1=repstr(infile1,'.fits','_Scale.fits')
psfile1 = repstr(infile1, '.fits', '_fit_sfrbin_1.ps')
psfile2 = repstr(infile1, '.fits', '_fit_sfrbin_2.ps')
psfile3 = repstr(infile1, '.fits', '_fit_sfrbin_3.ps')
psfile4 = repstr(infile1, '.fits', '_fit_sfrbin_4.ps')
psfile5 = repstr(infile1, '.fits', '_fit_sfrbin_5.ps')
psfile6 = repstr(infile1, '.fits', '_fit_sfrbin_6.ps')

;; manually check, independent bins
;; total Ca II mass between rpmin and rpmax
rpmin=10.
rpmax=200.
kpc = 3.08D21 ;cm
factor = 1.13E12*1E8/0.6267/(3934.78)^2 ; turn REW to column density
Mca40 = 40.*1.67D-27
Msolar = 2.D30

;; manually check, independent bins
i_indep = [[0], lindgen(6)*2+1]
;i_indep = lindgen(13)
ifit = i_indep[0:6]

;i_indep = [[0], lindgen(10)*1+1]
;ifit = i_indep[0:8]
;maxrad[ifit[7]] = 0.2
;maxrad[ifit[8]] = 0.2

;i_indep = [[0], lindgen(6)*2+2]
;ifit = i_indep[0:3]
;maxrad[ifit[3]] = 0.2

;; Int (factor*out_intercept*(r/100.)^slope*2.*pi*rdr) = factor*out_intercept*2.*pi*(1/100.)^slope*r^(2.+slope)
annulus_area = factor*kpc^2*((maxrad*1000.D0)^2-(minrad*1000.D0)^2)*!dpi
;; if not interpolating
annulus_area[ifit[0]]=factor*kpc^2*!dpi*[((maxrad[ifit[0]]*1000.D0)^2-(minrad[ifit[0]]*1000.D0)^2) $
                     +((minrad[ifit[1]]*1000.D0)^2-(minrad[ifit[0]]*1000.D0)^2)]
;annulus_area[ifit[7]]=factor*kpc^2*!dpi*[((maxrad[ifit[7]]*1000.D0)^2-(minrad[ifit[7]]*1000.D0)^2) $
;                     +((maxrad[ifit[7]]*1000.D0)^2-(maxrad[ifit[6]]*1000.D0)^2)]
annulus_area = annulus_area/2.

thick=8
xthick=8
ythick=8
charsize=1.4
charthick=3

;xra=[11,120]
xra=[10,300]
yra=[1.E-3, 7.E-1]
;yra=[2.E-3, 8.E-1]
xtitle='r_p (kpc)' 
ytitle='<W_0^{K}(Ca II)> (\AA)'
title='Single Gaussian Line Profile Measurement'

xx = 10.^(findgen(150)*0.024+0.1)
print, minmax(xx)

if (DoScale) then begin
   fid_slope = -1.38
   fid_intercept = alog10(0.0032)
endif else begin
   fid_slope = -1.38
   fid_intercept = alog10(0.011)
endelse

;; get the fitting for SF
jhusdss_powerlaw_fit, a_sf[ifit].rp/rvir*2., a_sf[ifit].ew_nofit_2, sdev_nofit_2_mc[ifit], in_slope=fid_slope, in_intercept=10.^fid_intercept, $
   slope=sf_slope, intercept=sf_intercept, err_slope=err_sf_slope, err_intercept=err_sf_intercept, /fixslope
;; get the fitting for quiescent
jhusdss_powerlaw_fit, a_quiescent[ifit].rp/rvir*2., a_quiescent[ifit].ew_nofit_2, sdev_nofit_2_mc[ifit], in_slope=fid_slope, in_intercept=10.^fid_intercept, $
   slope=quiescent_slope, intercept=quiescent_intercept, err_slope=err_quiescent_slope, err_intercept=err_quiescent_intercept, /fixslope

print, sf_slope, sf_intercept
print, quiescent_slope, quiescent_intercept
print, err_quiescent_intercept, err_sf_intercept
sf_intercept=alog10(sf_intercept)
quiescent_intercept=alog10(quiescent_intercept)

pos = [0.17, 0.20, 0.90, 0.90]
k_print, filename=psfile1, axis_char_scale=1.3, xsize=7, ysize=5
  djs_plot, a_sf[i_indep].rp, a_sf[i_indep].ew_nofit_2, xtitle=xtitle, ytitle=ytitle, $
      xra=xra, yra=yra, xst=5, yst=1, thick=thick, xthick=xthick, ythick=ythick, $
      charsize=charsize, charthick=charthick, /nodata, /xlog, /ylog, pos=pos, $
      ytickformat='jhusdss_tick_exponent'

  if (DoScale) then begin
      djs_axis, xaxis=0, xtitle='r_p/r_{vir}', charsize=charsize, charthick=charthick, xra=xra/rvir, xthick=xthick
      djs_axis, xaxis=1, charsize=charsize, charthick=charthick, xra=xra/rvir, xthick=xthick, xtickformat='(A1)'
  endif else begin
      djs_axis, xaxis=0, xtitle=xtitle, charsize=charsize, charthick=charthick, xra=xra, xthick=xthick
      djs_axis, xaxis=1, charsize=charsize, charthick=charthick, xra=xra, xthick=xthick, xtickformat='(A1)'
  endelse

; ibegin=0
; iend=11 
; polyfill, [a_sf[ibegin:iend].rp, reverse(a_sf[ibegin:iend].rp)], [sdev_nofit_2_mc[ibegin:iend], -reverse(sdev_nofit_2_mc[ibegin:iend])], color=djs_icolor('light gray')
; djs_axis, xaxis=0, xthick=xthick, xtickformat='(A1)'

  djs_oplot, xx, 10.^(alog10(xx/rvir*2.)*sf_slope+sf_intercept), thick=thick+2, color='navy', linestyle=0
  djs_oplot, xx, 10.^(alog10(xx/rvir*2.)*quiescent_slope+quiescent_intercept), thick=thick+2, color='red', linestyle=2

  rpscale=0.95
; djs_oplot, xx, 10.^(alog10(xx/rvir)*fid_slope+fid_intercept), thick=thick+2, color='blue'
  plotsym, 4, thick=2;, /fill
  djs_oplot, a_quiescent[i_indep].rp*rpscale, a_quiescent[i_indep].ew_nofit_2,  $
      psym=8, symsize=2, color=djs_icolor('red'), thick=thick, errcolor=djs_icolor('red'), errthick=thick
; djs_oplot, a_quiescent[i_indep].rp*rpscale, (a_quiescent[i_indep].ew_nofit_2>1.E-5),  $
;     color=djs_icolor('red'), thick=thick, linestyle=2
  oploterror, a_quiescent[i_indep].rp*rpscale, a_quiescent[i_indep].ew_nofit_2,  sdev_nofit_2_mc[i_indep], $
      psym=8, symsize=2, color=djs_icolor('red'), thick=thick, errcolor=djs_icolor('red'), errthick=thick-2, errstyle=0
  plotsym, 0, thick=2, /fill
  djs_oplot, a_sf[i_indep].rp, a_sf[i_indep].ew_nofit_2, $
      psym=8, symsize=2, color=djs_icolor('navy'), thick=thick, errcolor=djs_icolor('navy'), errthick=thick
  oploterror, a_sf[i_indep].rp, a_sf[i_indep].ew_nofit_2, sdev_nofit_2_mc[i_indep], $
      psym=8, symsize=2, color=djs_icolor('navy'), thick=thick, errcolor=djs_icolor('navy'), errthick=thick
; djs_oplot, a_sf[i_indep].rp, (a_sf[i_indep].ew_nofit_2>1.E-5), $
;     color=djs_icolor('navy'), thick=thick


;  ytitle1 = 'n (ca ii) [cm^{-2}]'
   ;; see eq. 9.15 in bruce drain
;  factor = 1.13e12*1e8/0.6267/(3934.78)^2
;  djs_axis, yaxis=1, yra=yra*factor, /ylog, ythick=thick, ytitle=ytitle1, $
;      charsiz=charsize, charthick=charthick

  plotsym, 0, /fill
  djs_xyouts, 90, 4E-1, 'Star-forming', color='navy', charsize=charsize-0.2, charthick=charthick
  djs_oplot, 250*[1,1], 4.5E-1*[1,1], color='navy', psym=8, symsize=1.3, thick=thick
  plotsym, 4, thick=thick
  djs_xyouts, 90, 4E-1*0.79^2, '    Quiescent', color='red', charsize=charsize-0.2, charthick=charthick
  djs_oplot, 250*[1,1], 4.5E-1*0.77^2*[1,1], color='red', psym=8, symsize=1, thick=thick

; if (~DoScale) then djs_xyouts, 17.5, 0.000042, '20', charsize=charsize+0.35, charthick=charthick
k_end_print

k_print, filename=psfile2, axis_char_scale=1.3, xsize=7, ysize=5
  djs_plot, a_sf[i_indep].rp, a_sf[i_indep].ew_nofit_2, xtitle=xtitle, ytitle=ytitle, $
      xra=xra, yra=yra, xst=1, yst=1, thick=thick, xthick=xthick, ythick=ythick, $
      charsize=charsize, charthick=charthick, /nodata, /xlog, /ylog, pos=pos, $
      ytickformat='jhusdss_tick_exponent'
  djs_oplot, xx, 10.^(alog10(xx/rvir)*fid_slope+fid_intercept), thick=thick+2, color='blue'
  rpscale=1.4
  plotsym, 4, thick=2;, /fill
  djs_oplot, a_quiescent[i_indep].rp/rpscale, a_quiescent[i_indep].ew_nofit_2,  $
      psym=8, symsize=2, color=djs_icolor('red'), thick=thick, errcolor=djs_icolor('red'), errthick=thick
  djs_oplot, a_quiescent[i_indep].rp/rpscale, (a_quiescent[i_indep].ew_nofit_2>1.E-5),  $
      color=djs_icolor('red'), thick=thick, linestyle=2
  oploterror, a_quiescent[i_indep].rp/rpscale, a_quiescent[i_indep].ew_nofit_2,  sdev_nofit_2_mc[i_indep], $
      psym=8, symsize=2, color=djs_icolor('red'), thick=thick, errcolor=djs_icolor('red'), errthick=thick-2, errstyle=0
  plotsym, 0, thick=2, /fill
  djs_oplot, a_sf[i_indep].rp, a_sf[i_indep].ew_nofit_2, $
      psym=8, symsize=2, color=djs_icolor('navy'), thick=thick, errcolor=djs_icolor('navy'), errthick=thick
  oploterror, a_sf[i_indep].rp, a_sf[i_indep].ew_nofit_2, sdev_nofit_2_mc[i_indep], $
      psym=8, symsize=2, color=djs_icolor('navy'), thick=thick, errcolor=djs_icolor('navy'), errthick=thick
  djs_oplot, a_sf[i_indep].rp, (a_sf[i_indep].ew_nofit_2>1.E-5), $
      color=djs_icolor('navy'), thick=thick

;  ytitle1 = 'n (ca ii) [cm^{-2}]'
   ;; see eq. 9.15 in bruce drain
;  factor = 1.13e12*1e8/0.6267/(3934.78)^2
;  djs_axis, yaxis=1, yra=yra*factor, /ylog, ythick=thick, ytitle=ytitle1, $
;      charsiz=charsize, charthick=charthick

  plotsym, 0, /fill
  djs_xyouts, 300, 4E-1, 'star-forming', color='navy', charsize=charsize-0.2, charthick=charthick
  djs_oplot, 1200*[1,1], 4.5E-1*[1,1], color='navy', psym=8, symsize=1.3, thick=thick
  plotsym, 4, thick=thick
  djs_xyouts, 300, 4E-1*0.8^2, 'quiescent', color='red', charsize=charsize-0.2, charthick=charthick
  djs_oplot, 1200*[1,1], 4.5E-1*0.8^2*[1,1], color='red', psym=8, symsize=1, thick=thick

  djs_xyouts, 17.5, 0.000042, '20', charsize=charsize+0.4, charthick=charthick
k_end_print


if keyword_set(movingon) then begin
xra=[5,2000]
yra=[5.e-5, 6.e-1]
k_print, filename=psfile2, axis_char_scale=1.3, xsize=8, ysize=6
  djs_plot, a.rp, a.ew_nofit, xtitle=xtitle, ytitle=ytitle, title=title, $
      xra=xra, yra=yra, xst=1, yst=1, thick=thick, xthick=xthick, ythick=ythick, $
      charsize=charsize, charthick=charthick, /nodata, /xlog, /ylog, pos=pos
  ibegin=1
  iend=23 
; polyfill, [a[ibegin:iend].rp, reverse(a[ibegin:iend].rp)], [a[ibegin:iend].sdev_nofit_mc, -reverse(a[ibegin:iend].sdev_nofit_mc)], color=djs_icolor('light gray')
  djs_oplot, a[i_indep].rp, sdev_nofit_mc[i_indep], linestyle=1, thick=thick, color='gray'
  djs_oplot, xx, 10.^(alog10(xx/rvir)*slope+intercept), thick=thick, color='light blue'
  oploterror, a[i_indep].rp, a[i_indep].ew_nofit, sdev_nofit_mc[i_indep], $
      psym=8, symsize=2, color=djs_icolor('light blue'), thick=thick, errcolor=djs_icolor('light blue'), errthick=thick
k_end_print

k_print, filename=psfile3, axis_char_scale=1.3, xsize=8, ysize=6
  xra=[220, 2000]
  yra=[-4e-3, 4e-3]
  djs_plot, a.rp, a.ew_nofit, xtitle=xtitle, ytitle=ytitle, title=title, $
      xra=xra, yra=yra, xst=1, yst=1, thick=thick, xthick=xthick, ythick=ythick, $
      charsize=charsize, charthick=charthick, /nodata, /xlog, pos=pos;, /ylog
  ibegin=13
  iend=23 
  polyfill, [a[ibegin:iend].rp, reverse(a[ibegin:iend].rp)], [sdev_nofit_mc[ibegin:iend], -reverse(sdev_nofit_mc[ibegin:iend])], color=djs_icolor('light gray')
  djs_oplot, xx, 10.^(alog10(xx/rvir)*slope+intercept), thick=thick, color='light blue'
  oploterror, a[i_indep].rp, a[i_indep].ew_nofit, sdev_nofit_mc[i_indep], $
      psym=8, symsize=2, color=djs_icolor('light blue'), thick=thick, errcolor=djs_icolor('light blue'), errthick=thick
; djs_axis, xaxis=0, charsize=charsize, charthick=charthick, thick=thick, $
;     xtickv=[300], xticknames='300'
  djs_xyouts, 270, -0.0046, '300', charsize=charsize+0.4, charthick=charthick

k_end_print

xra=[5,2000]
yra=[5.e-5, 1.e0]
xtitle='r_p (kpc)' 
ytitle='w_0^{\lambda 3934} (ca ii) [\aa]'
title='double gaussian line profile measurement'

k_print, filename=psfile4, axis_char_scale=1.3, xsize=10, ysize=6
  djs_plot, a[i_indep].rp, a[i_indep].ew_nofit_2, xtitle=xtitle, ytitle=ytitle, $
      xra=xra, yra=yra, xst=1, yst=1, thick=thick, xthick=xthick, ythick=ythick, $
      charsize=charsize, charthick=charthick, /nodata, /xlog, /ylog, pos=pos, $
      ytickformat='jhusdss_tick_exponent', ystyle=9
  djs_oplot, xx, 10.^(alog10(xx/rvir)*slope+intercept), thick=thick, color='red', linestyle=1
  djs_oplot, xx, 10.^(alog10(xx/rvir)*fid_slope+fid_intercept), thick=thick+5, color='light blue', linestyle=0
  oploterror, a[i_indep].rp*1.08, a[i_indep].ew_nofit, sdev_nofit_mc[i_indep], $
      psym=4, symsize=2, color=djs_icolor('gray'), thick=thick, errcolor=djs_icolor('gray'), errthick=thick
  oploterror, a[i_indep].rp, a[i_indep].ew_nofit_2, sdev_nofit_2_mc[i_indep], $
      psym=8, symsize=2, color=djs_icolor('red'), thick=thick, errcolor=djs_icolor('red'), errthick=thick
  legend, ['Double Gaussian fitting', 'Single Gaussian fitting'], psym=[8,4], $
      colors=[djs_icolor('red'), djs_icolor('gray')], $
      textcolors=[djs_icolor('red'), djs_icolor('gray')], $
      box=0, /bottom, /left, charsize=charsize-0.3, charthick=charthick

   text_relation = 'W_0=('+string(out_intercept*100, format='(f3.1)')+'\pm'+string(err_intercept*100,format='(f3.1)')+')\times10^{-2}\times(r_p/200 kpc)^{'+string(slope, format='(f5.2)')+'\pm'+string(err_slope,format='(f5.2)')+'}'
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
   factor = 1.13E12*1E8/0.6267/(3934.78)^2
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
  djs_oplot, xx, 10.^(alog10(xx/rvir)*slope+intercept), thick=thick-3, color='red'
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
  djs_oplot, xx, 10.^(alog10(xx/rvir)*slope+intercept), thick=thick, color='light blue'
  oploterror, a[i_indep].rp, a[i_indep].ew_nofit_2, sdev_nofit_2_mc[i_indep], $
      psym=8, symsize=2, color=djs_icolor('light blue'), thick=thick, errcolor=djs_icolor('light blue'), errthick=thick
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
  djs_oplot, xx, 10.^(alog10(xx/rvir)*slope+intercept), thick=thick
  djs_oplot, !x.crange, [0,0], color='light gray', thick=thick, linestyle=1
  oploterror, a[i_indep].rp, a[i_indep].ew_nofit_2, sdev_nofit_2_mc[i_indep], $
      psym=8, symsize=2, color=djs_icolor('light blue'), thick=thick, errcolor=djs_icolor('light blue'), errthick=thick
  djs_xyouts, 270, -0.0046, '300', charsize=charsize+0.4, charthick=charthick
k_end_print
endif

;; total Ca II mass between rpmin and rpmax
rpmin=10.
rpmax=200.
kpc = 3.08D21 ;cm
factor = 1.13E12*1E8/0.6267/(3934.78)^2 ; turn REW to column density
Mca40 = 40.*1.67D-27
Msolar = 2.D30
;; Int (factor*out_intercept*(r/rvir)^slope*2.*pi*rdr) = factor*out_intercept*2.*pi*(1/rvir)^slope*r^(2.+slope)
annulus_area = factor*kpc^2*((maxrad*1000.D0)^2-(minrad*1000.D0)^2)*!dpi
number_caii = total(a_sf[ifit].ew_nofit_2*annulus_area[ifit])
error_number_caii = sqrt(total((sdev_nofit_2_mc[ifit]*annulus_area[ifit])^2))
error_mass_caii = error_number_caii*Mca40/Msolar
mass_caii = number_caii*Mca40/Msolar
print, 'star-forming M_caii: ', mass_caii, error_mass_caii, alog10(mass_caii), $
  (alog10(mass_caii+error_mass_caii)-alog10(mass_caii-error_mass_caii))/2.
;; total HI mass
;; solar: N(H)/N(Ca) = 10^5.69 -> M(H)/M(Ca) = 10^5.69/40
mass_HI_solar = 10.D0^5.69/40.*mass_caii
print, mass_hi_solar

number_caii = total(a_quiescent[ifit].ew_nofit_2*annulus_area[ifit])
mass_caii = number_caii*Mca40/Msolar
print, 'quiescent M_caii: ', mass_caii, error_mass_caii, alog10(mass_caii), $
  (alog10(mass_caii+error_mass_caii)-alog10(mass_caii-error_mass_caii))/2.
;; total HI mass
;; solar: N(H)/N(Ca) = 10^5.69 -> M(H)/M(Ca) = 10^5.69/40
mass_HI_solar = 10.D0^5.69/40.*mass_caii
print, mass_hi_solar


end
