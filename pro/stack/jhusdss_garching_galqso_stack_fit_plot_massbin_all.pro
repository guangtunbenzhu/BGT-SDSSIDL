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
;rvir = 200. ;; kpc for 10^10.3 (10^12)
if (DoScale) then rvir=200. else rvir=200.

nmfver=106
stackpath = jhusdss_get_path(/nmfqso)+'/'+string(nmfver, format='(I3.3)')+'/Stack/'

infile0 = stackpath + jhusdss_garching_galqso_stack_filename(nmfver, boss=boss)

massive = 1b 
infile = infile0
if (massive) then infile = repstr(infile, '.fits', '_highmass.fits') else infile = repstr(infile, '.fits', '_lowmass.fits')
if (DoScale) then infile=repstr(infile,'.fits','_Scale.fits')
infile = repstr(infile, '.fits', '_fit.fits')
a_highmass = mrdfits(infile, 1)

massive = 0b 
infile = infile0
if (massive) then infile = repstr(infile, '.fits', '_highmass.fits') else infile = repstr(infile, '.fits', '_lowmass.fits')
if (DoScale) then infile=repstr(infile,'.fits','_Scale.fits')
infile = repstr(infile, '.fits', '_fit.fits')
a_lowmass= mrdfits(infile, 1)

mcmax = 200
ew_nofit_mc = [a_highmass.ew_nofit_mc[0:mcmax], a_lowmass.ew_nofit_mc[0:mcmax]]
ew_nofit_2_mc = [a_highmass.ew_nofit_2_mc[0:mcmax], a_lowmass.ew_nofit_2_mc[0:mcmax]]

; for i=0,11 do begin
;     shuffle_infile = repstr(infile0, '.fits', '_shuffle_'+string(i+1, format='(i2.2)')+'_fit.fits')
;     ashuffle = mrdfits(shuffle_infile, 1)
;     ew_nofit_mc = [ew_nofit_mc, ashuffle[1:24].ew_nofit_mc[0:mcmax]*sqrt(2.)]
;     ew_nofit_2_mc = [ew_nofit_2_mc, ashuffle[1:24].ew_nofit_2_mc[0:mcmax]*sqrt(2.)]
; endfor

sdev_nofit_mc = fltarr(n_elements(a_highmass))
sdev_nofit_2_mc = fltarr(n_elements(a_highmass))
for irp=0,n_elements(a_highmass)-1 do begin
    tmp = moment(ew_nofit_mc[*,irp], sdev=sdev)
    sdev_nofit_mc[irp] = sdev
    tmp = moment(ew_nofit_2_mc[*,irp], sdev=sdev)
    sdev_nofit_2_mc[irp] = sdev
endfor 
;print, sdev_nofit_2_mc

infile1 = infile0
if (DoScale) then infile1=repstr(infile1,'.fits','_Scale.fits')
psfile1 = repstr(infile1, '.fits', '_fit_massbin_1.ps')
psfile2 = repstr(infile1, '.fits', '_fit_massbin_2.ps')
psfile3 = repstr(infile1, '.fits', '_fit_massbin_3.ps')
psfile4 = repstr(infile1, '.fits', '_fit_massbin_4.ps')
psfile5 = repstr(infile1, '.fits', '_fit_massbin_5.ps')
psfile6 = repstr(infile1, '.fits', '_fit_massbin_6.ps')

;; total Ca II mass between rpmin and rpmax
rpmin=10.
rpmax=200.
kpc = 3.08D21 ;cm
factor = 1.13E12*1E8/0.6485/(3934.78)^2 ; turn REW to column density
Mca40 = 40.*1.67D-27
Msolar = 2.D30

;; manually check, independent bins
i_indep = [[0], lindgen(6)*2+1]
;i_indep = [[0], lindgen(6)*2+1]
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

;xra=[11, 120]
xra=[10, 300]
yra=[1.E-3, 7.E-1]
;yra=[2.E-3, 8.E-1]
xtitle='r_p (kpc)' 
ytitle='<W_0^{K}(Ca II)> (\AA)'
title='Single Gaussian Line Profile Measurement'

xx = 10.^(findgen(150)*0.024+0.10)
print, minmax(xx)

if (DoScale) then begin
   fid_slope = -1.38
   fid_intercept = alog10(0.0032)
endif else begin
   fid_slope = -1.38
   fid_intercept = alog10(0.011)
endelse

;; get the fitting for high mass
jhusdss_powerlaw_fit, a_highmass[ifit].rp/rvir*2., a_highmass[ifit].ew_nofit_2, sdev_nofit_2_mc[ifit], in_slope=fid_slope, in_intercept=10.^fid_intercept, $
   slope=high_slope, intercept=high_intercept, err_slope=err_high_slope, err_intercept=err_high_intercept, /fixslope
;; get the fitting for low mass
jhusdss_powerlaw_fit, a_lowmass[ifit].rp/rvir*2., a_lowmass[ifit].ew_nofit_2, sdev_nofit_2_mc[ifit], in_slope=fid_slope, in_intercept=10.^fid_intercept, $
   slope=low_slope, intercept=low_intercept, err_slope=err_low_slope, err_intercept=err_low_intercept, /fixslope

print, high_slope, high_intercept, err_high_slope, err_high_intercept
print, low_slope, low_intercept, err_low_slope, err_low_intercept
high_intercept=alog10(high_intercept)
low_intercept=alog10(low_intercept)

pos = [0.17, 0.20, 0.90, 0.90]
k_print, filename=psfile1, axis_char_scale=1.3, xsize=7, ysize=5
  djs_plot, a_lowmass[i_indep].rp, a_lowmass[i_indep].ew_nofit_2, xtitle=xtitle, ytitle=ytitle, $
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
; polyfill, [a_lowmass[ibegin:iend].rp, reverse(a_lowmass[ibegin:iend].rp)], [sdev_nofit_2_mc[ibegin:iend], -reverse(sdev_nofit_2_mc[ibegin:iend])], color=djs_icolor('light gray')
; djs_axis, xaxis=0, xthick=xthick, xtickformat='(A1)'

; djs_oplot, xx, 10.^(alog10(xx/rvir)*fid_slope+fid_intercept), thick=thick+2, color='blue', linestyle=0

; ixx = where(xx lt 250., comp=jxx)
; djs_oplot, xx[ixx], 10.^(alog10(xx[ixx]/rvir)*high_slope+high_intercept), thick=thick, color='dark green', linestyle=0
; djs_oplot, xx[jxx], 10.^(alog10(xx[jxx]/rvir)*high_slope+high_intercept), thick=thick, color='blue', linestyle=1
; djs_oplot, xx[ixx], 10.^(alog10(xx[ixx]/rvir)*low_slope+low_intercept), thick=thick, color='magenta', linestyle=2
; djs_oplot, xx[jxx], 10.^(alog10(xx[jxx]/rvir)*low_slope+low_intercept), thick=thick, color='magenta', linestyle=1

  djs_oplot, xx, 10.^(alog10(xx/rvir*2.)*high_slope+high_intercept), thick=thick+2, color='dark green', linestyle=0
  djs_oplot, xx, 10.^(alog10(xx/rvir*2.)*low_slope+low_intercept), thick=thick+2, color='magenta', linestyle=2
  rpscale = 0.95
  plotsym, 4, thick=2;, /fill
  djs_oplot, a_lowmass[i_indep].rp*rpscale, a_lowmass[i_indep].ew_nofit_2,  $
      psym=8, symsize=2, color=djs_icolor('magenta'), thick=thick, errcolor=djs_icolor('magenta'), errthick=thick
; djs_oplot, a_lowmass[i_indep].rp*rpscale, (a_lowmass[i_indep].ew_nofit_2>1.E-5),  $
;     color=djs_icolor('magenta'), thick=thick, linestyle=2
  oploterror, a_lowmass[i_indep].rp*rpscale, a_lowmass[i_indep].ew_nofit_2,  $
      sdev_nofit_2_mc[i_indep], $
      psym=8, symsize=2, color=djs_icolor('magenta'), thick=thick, errcolor=djs_icolor('magenta'), errthick=thick-2, errstyle=0
; oploterror, a_lowmass[i_indep].rp*rpscale, a_lowmass[i_indep].ew_nofit_2,  $
;     (a_lowmass[i_indep].rp-minrad[i_indep]*1000.), sdev_nofit_2_mc[i_indep], /lobar, $
;     psym=8, symsize=2, color=djs_icolor('magenta'), thick=thick, errcolor=djs_icolor('magenta'), errthick=thick-2, errstyle=0
; oploterror, a_lowmass[i_indep].rp*rpscale, a_lowmass[i_indep].ew_nofit_2,  $
;     (maxrad[i_indep]*1000.-a_lowmass[i_indep].rp), sdev_nofit_2_mc[i_indep], /hibar, $
;     psym=8, symsize=2, color=djs_icolor('magenta'), thick=thick, errcolor=djs_icolor('magenta'), errthick=thick-2, errstyle=0
; djs_oplot, a_massive_sf[i_indep].rp, a_massive_sf[i_indep].ew_nofit_2,  $
;     color=djs_icolor('dark green'), thick=thick, errcolor=djs_icolor('dark green'), errthick=thick
; oploterror, a_massive_sf[i_indep].rp, a_massive_sf[i_indep].ew_nofit_2, sdev_nofit_2_mc[i_indep], $
;     psym=8, symsize=2, color=djs_icolor('dark green'), thick=thick, errcolor=djs_icolor('dark green'), errthick=thick
  plotsym, 0, thick=2, /fill
  djs_oplot, a_highmass[i_indep].rp, a_highmass[i_indep].ew_nofit_2, $
      psym=8, symsize=2, color=djs_icolor('dark green'), thick=thick, errcolor=djs_icolor('dark green'), errthick=thick
; djs_oplot, a_highmass[i_indep].rp, (a_highmass[i_indep].ew_nofit_2>1.E-5), $
;     color=djs_icolor('dark green'), thick=thick
  oploterror, a_highmass[i_indep].rp, a_highmass[i_indep].ew_nofit_2, $
      sdev_nofit_2_mc[i_indep], $
      psym=8, symsize=2, color=djs_icolor('dark green'), thick=thick, errcolor=djs_icolor('dark green'), errthick=thick
; oploterror, a_highmass[i_indep].rp, a_highmass[i_indep].ew_nofit_2, $
;     (a_highmass[i_indep].rp-minrad[i_indep]*1000.), sdev_nofit_2_mc[i_indep], /lobar, $
;     psym=8, symsize=2, color=djs_icolor('dark green'), thick=thick, errcolor=djs_icolor('dark green'), errthick=thick
; oploterror, a_highmass[i_indep].rp, a_highmass[i_indep].ew_nofit_2, $
;     (maxrad[i_indep]*1000.-a_highmass[i_indep].rp), sdev_nofit_2_mc[i_indep], /hibar, $
;     psym=8, symsize=2, color=djs_icolor('dark green'), thick=thick, errcolor=djs_icolor('dark green'), errthick=thick

; oploterror, a_lowmass_sf[i_indep].rp*1.05, a_lowmass_sf[i_indep].ew_nofit_2, sdev_nofit_2_mc[i_indep], $
;     psym=8, symsize=2, color=djs_icolor('light blue'), thick=thick, errcolor=djs_icolor('light blue'), errthick=thick
; plotsym, 0, thick=thick
; djs_oplot, a_massive_quiescent[i_indep].rp*1.05^2, a_massive_quiescent[i_indep].ew_nofit_2, $
;     psym=8, symsize=1.5, color=djs_icolor('magenta'), thick=thick, errcolor=djs_icolor('magenta'), errthick=thick, errstyle=2
; djs_oplot, a_massive_quiescent[i_indep].rp*1.05^2, a_massive_quiescent[i_indep].ew_nofit_2, $
;     color=djs_icolor('magenta'), thick=thick, errcolor=djs_icolor('magenta'), errthick=thick, errstyle=2
; oploterror, a_massive_quiescent[i_indep].rp*1.05^2, a_massive_quiescent[i_indep].ew_nofit_2, sdev_nofit_2_mc[i_indep], $
;     psym=8, symsize=1.5, color=djs_icolor('magenta'), thick=thick, errcolor=djs_icolor('magenta'), errthick=thick, errstyle=2
; plotsym, 4, thick=thick
; djs_oplot, a_lowmass_quiescent[i_indep].rp*1.05^3, a_lowmass_quiescent[i_indep].ew_nofit_2, $
;     psym=8, symsize=1.5, color=djs_icolor('magenta'), thick=thick, errcolor=djs_icolor('magenta'), errthick=thick, errstyle=2
; djs_oplot, a_lowmass_quiescent[i_indep].rp*1.05^3, a_lowmass_quiescent[i_indep].ew_nofit_2, $
;     color=djs_icolor('magenta'), thick=thick, errcolor=djs_icolor('magenta'), errthick=thick, errstyle=2
; oploterror, a_lowmass_quiescent[i_indep].rp*1.05^3, a_lowmass_quiescent[i_indep].ew_nofit_2, sdev_nofit_2_mc[i_indep], $
;     psym=8, symsize=1.5, color=djs_icolor('magenta'), thick=thick, errcolor=djs_icolor('magenta'), errthick=thick, errstyle=2

;  ytitle1 = 'n (ca ii) [cm^{-2}]'
   ;; see eq. 9.15 in bruce drain
;  factor = 1.13e12*1e8/0.6267/(3934.78)^2
;  djs_axis, yaxis=1, yra=yra*factor, /ylog, ythick=thick, ytitle=ytitle1, $
;      charsiz=charsize, charthick=charthick

  plotsym, 0, /fill
  djs_xyouts, 110, 4E-1, 'High-mass', color='dark green', charsize=charsize-0.2, charthick=charthick
  djs_oplot, 250*[1,1], 4.5E-1*[1,1], color='dark green', psym=8, symsize=1.3, thick=thick
  plotsym, 4, thick=thick
  djs_xyouts, 110, 4E-1*0.79^2, ' Low-mass', color='magenta', charsize=charsize-0.2, charthick=charthick
  djs_oplot, 250*[1,1], 4.5E-1*0.77^2*[1,1], color='magenta', psym=8, symsize=1, thick=thick

; if (~DoScale) then djs_xyouts, 17.5, 0.00052, '20', charsize=charsize+0.35, charthick=charthick

k_end_print

number_caii = total(a_lowmass[ifit].ew_nofit_2*annulus_area[ifit])
error_number_caii = sqrt(total((sdev_nofit_2_mc[ifit]*annulus_area[ifit])^2))
error_mass_caii = error_number_caii*Mca40/Msolar
mass_caii = number_caii*Mca40/Msolar
print, 'low-mass M_caii: ', mass_caii, error_mass_caii, alog10(mass_caii), $
  (alog10(mass_caii+error_mass_caii)-alog10(mass_caii-error_mass_caii))/2.
;; total HI mass
;; solar: N(H)/N(Ca) = 10^5.66 -> M(H)/M(Ca) = 10^5.66/40
mass_HI_solar = 10.D0^5.66/40.*mass_caii
print, mass_hi_solar

number_caii = total(a_highmass[ifit].ew_nofit_2*annulus_area[ifit])
mass_caii = number_caii*Mca40/Msolar
print, 'high-mass M_caii: ', mass_caii, error_mass_caii, alog10(mass_caii), $
  (alog10(mass_caii+error_mass_caii)-alog10(mass_caii-error_mass_caii))/2.
;; total HI mass
;; solar: N(H)/N(Ca) = 10^5.66 -> M(H)/M(Ca) = 10^5.66/40
mass_HI_solar = 10.D0^5.66/40.*mass_caii
print, mass_hi_solar

end
