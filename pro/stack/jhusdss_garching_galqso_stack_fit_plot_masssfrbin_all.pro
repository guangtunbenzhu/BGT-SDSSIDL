;minrad_tmp = 10.^(alog10(0.020)+findgen(11)*(0.5*alog10(2.0)))
;maxrad_tmp = minrad_tmp*2.0
;minrad = [[0.010], minrad_tmp]
;maxrad = [[0.020], maxrad_tmp]

;minrad_tmp = 10.^(alog10(0.010)+findgen(11)*(0.5*alog10(2.0)))
;maxrad_tmp = minrad_tmp*2.0
;minrad = minrad_tmp
;maxrad = maxrad_tmp

minrad_tmp = 10.^(alog10(0.010)+findgen(10)*alog10(1.35))
maxrad_tmp = minrad_tmp*1.35
minrad = minrad_tmp
maxrad = maxrad_tmp

;; interpolation
;minrad_tmp = 10.^(alog10(0.010)+findgen(11)*(alog10(1.2)))
;maxrad_tmp = minrad_tmp*1.2
;minrad = minrad_tmp
;maxrad = maxrad_tmp
;meanrad = (maxrad+minrad)/2.
;lowmass_sf = interpol(

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


massive = 1b & quiescent = 0b
infile = infile0
if (massive) then infile = repstr(infile, '.fits', '_highmass.fits') else infile = repstr(infile, '.fits', '_lowmass.fits')
if (quiescent) then infile = repstr(infile, '.fits', '_quiescent.fits') else infile = repstr(infile, '.fits', '_sf.fits')
infile = repstr(infile, '.fits', '_fit.fits')
a_massive_sf = mrdfits(infile, 1)

massive = 0b & quiescent = 0b
infile = infile0
if (massive) then infile = repstr(infile, '.fits', '_highmass.fits') else infile = repstr(infile, '.fits', '_lowmass.fits')
if (quiescent) then infile = repstr(infile, '.fits', '_quiescent.fits') else infile = repstr(infile, '.fits', '_sf.fits')
infile = repstr(infile, '.fits', '_fit.fits')
a_lowmass_sf= mrdfits(infile, 1)

massive = 1b & quiescent = 1b
infile = infile0
if (massive) then infile = repstr(infile, '.fits', '_highmass.fits') else infile = repstr(infile, '.fits', '_lowmass.fits')
if (quiescent) then infile = repstr(infile, '.fits', '_quiescent.fits') else infile = repstr(infile, '.fits', '_sf.fits')
infile = repstr(infile, '.fits', '_fit.fits')
a_massive_quiescent= mrdfits(infile, 1)

massive = 0b & quiescent = 1b
infile = infile0
if (massive) then infile = repstr(infile, '.fits', '_highmass.fits') else infile = repstr(infile, '.fits', '_lowmass.fits')
if (quiescent) then infile = repstr(infile, '.fits', '_quiescent.fits') else infile = repstr(infile, '.fits', '_sf.fits')
infile = repstr(infile, '.fits', '_fit.fits')
a_lowmass_quiescent= mrdfits(infile, 1)

mcmax = 200
ew_nofit_mc = [a_lowmass_sf.ew_nofit_mc[0:mcmax], a_massive_sf.ew_nofit_mc[0:mcmax], a_lowmass_quiescent.ew_nofit_mc[0:mcmax], a_massive_quiescent.ew_nofit_mc[0:mcmax]]
ew_nofit_2_mc = [a_lowmass_sf.ew_nofit_2_mc[0:mcmax], a_massive_sf.ew_nofit_2_mc[0:mcmax], a_lowmass_quiescent.ew_nofit_2_mc[0:mcmax], a_massive_quiescent.ew_nofit_2_mc[0:mcmax]]
; for i=0,11 do begin
;     shuffle_infile = repstr(infile0, '.fits', '_shuffle_'+string(i+1, format='(i2.2)')+'_fit.fits')
;     ashuffle = mrdfits(shuffle_infile, 1)
;     ew_nofit_mc = [ew_nofit_mc, ashuffle[1:24].ew_nofit_mc[0:mcmax]*2.]
;     ew_nofit_2_mc = [ew_nofit_2_mc, ashuffle[1:24].ew_nofit_2_mc[0:mcmax]*2.]
; endfor

sdev_nofit_mc = fltarr(n_elements(a_massive_quiescent))
sdev_nofit_2_mc = fltarr(n_elements(a_massive_quiescent))
a_massive_quiescent_sdev_nofit_2_mc = fltarr(n_elements(a_massive_quiescent))
a_lowmass_quiescent_sdev_nofit_2_mc = fltarr(n_elements(a_massive_quiescent))
a_massive_sf_sdev_nofit_2_mc = fltarr(n_elements(a_massive_quiescent))
a_lowmass_sf_sdev_nofit_2_mc = fltarr(n_elements(a_massive_quiescent))
sdev_nofit_2_mc = fltarr(n_elements(a_massive_quiescent))
for irp=0,n_elements(a_massive_quiescent)-1 do begin
    tmp = moment(ew_nofit_mc[*,irp], sdev=sdev)
    sdev_nofit_mc[irp] = sdev
    tmp = moment(ew_nofit_2_mc[*,irp], sdev=sdev)
    sdev_nofit_2_mc[irp] = sdev

    tmp = moment(a_massive_quiescent[irp].ew_nofit_2_mc[0:mcmax], sdev=sdev)
    a_massive_quiescent_sdev_nofit_2_mc[irp] = sdev
    tmp = moment(a_massive_sf[irp].ew_nofit_2_mc[0:mcmax], sdev=sdev)
    a_massive_sf_sdev_nofit_2_mc[irp] = sdev
    tmp = moment(a_lowmass_sf[irp].ew_nofit_2_mc[0:mcmax], sdev=sdev)
    a_lowmass_sf_sdev_nofit_2_mc[irp] = sdev
    tmp = moment(a_lowmass_quiescent[irp].ew_nofit_2_mc[0:mcmax], sdev=sdev)
    a_lowmass_quiescent_sdev_nofit_2_mc[irp] = sdev
endfor 
;print, sdev_nofit_2_mc

infile1 = infile0
if (DoScale) then infile1=repstr(infile1,'.fits','_Scale.fits')
psfile1 = repstr(infile1, '.fits', '_fit_masssfrbin_1.ps')
psfile2 = repstr(infile1, '.fits', '_fit_masssfrbin_2.ps')
psfile3 = repstr(infile1, '.fits', '_fit_masssfrbin_3.ps')
psfile4 = repstr(infile1, '.fits', '_fit_masssfrbin_4.ps')
psfile5 = repstr(infile1, '.fits', '_fit_masssfrbin_5.ps')
psfile6 = repstr(infile1, '.fits', '_fit_masssfrbin_6.ps')

;; total Ca II mass between rpmin and rpmax
rpmin=10.
rpmax=200.
kpc = 3.08D21 ;cm
factor = 1.13E12*1E8/0.6267/(3934.78)^2 ; turn REW to column density
Mca40 = 40.*1.67D-27
Msolar = 2.D30

;; manually check, independent bins
;i_indep = [[0], lindgen(10)*1+1]
;ifit = i_indep[0:8]
;maxrad[ifit[7]] = 0.2
;maxrad[ifit[8]] = 0.2
;i_indep = lindgen(10)
;ifit = i_indep

;; manually check, independent bins
i_indep = [[0], lindgen(6)*2+1]
;i_indep = lindgen(13)
ifit = i_indep[0:6]

;i_indep = [[0], lindgen(6)*2+2]
;ifit = i_indep[0:3]
;maxrad[ifit[3]] = 0.2

;; Int (factor*out_intercept*(r/100.)^slope*2.*pi*rdr) = factor*out_intercept*2.*pi*(1/100.)^slope*r^(2.+slope)
annulus_area = factor*kpc^2*((maxrad*1000.D0)^2-(minrad*1000.D0)^2)*!dpi
;; if not interpolating
;annulus_area[ifit[0]]=factor*kpc^2*!dpi*[((maxrad[ifit[0]]*1000.D0)^2-(minrad[ifit[0]]*1000.D0)^2) $
;                     +((minrad[ifit[1]]*1000.D0)^2-(minrad[ifit[0]]*1000.D0)^2)]
;annulus_area[ifit[7]]=factor*kpc^2*!dpi*[((maxrad[ifit[7]]*1000.D0)^2-(minrad[ifit[7]]*1000.D0)^2) $
;                     +((maxrad[ifit[7]]*1000.D0)^2-(maxrad[ifit[6]]*1000.D0)^2)]
;annulus_area = annulus_area/2.

;number_caii = total(a_lowmass_quiescent[ifit].ew_nofit_2*annulus_area[ifit]/a_lowmass_quiescent_sdev_nofit_2_mc[ifit]^2)/total(1./a_lowmass_quiescent_sdev_nofit_2_mc[ifit]^2)
number_caii = total(a_lowmass_quiescent[ifit].ew_nofit_2*annulus_area[ifit])
error_number_caii = sqrt(total((a_lowmass_quiescent_sdev_nofit_2_mc[ifit]*annulus_area[ifit])^2))
error_mass_caii = error_number_caii*Mca40/Msolar
mass_caii = number_caii*Mca40/Msolar
print, 'low-mass quiescent M_caii: ', mass_caii, error_mass_caii, alog10(mass_caii), $
 (alog10(mass_caii+error_mass_caii)-alog10(mass_caii-error_mass_caii))/2.
;; total HI mass
;; solar: N(H)/N(Ca) = 10^5.69 -> M(H)/M(Ca) = 10^5.69/40
mass_HI_solar = 10.D0^5.69/40.*mass_caii
print, mass_hi_solar

number_caii = total(a_massive_quiescent[ifit].ew_nofit_2*annulus_area[ifit])
error_number_caii = sqrt(total((a_massive_quiescent_sdev_nofit_2_mc[ifit]*annulus_area[ifit])^2))
error_mass_caii = error_number_caii*Mca40/Msolar
mass_caii = number_caii*Mca40/Msolar
print, 'high-mass quiescent M_caii: ', mass_caii, error_mass_caii, alog10(mass_caii), $
  (alog10(mass_caii+error_mass_caii)-alog10(mass_caii-error_mass_caii))/2.
;; total HI mass
;; solar: N(H)/N(Ca) = 10^5.69 -> M(H)/M(Ca) = 10^5.69/40
mass_HI_solar = 10.D0^5.69/40.*mass_caii
print, mass_hi_solar

number_caii = total(a_massive_sf[ifit].ew_nofit_2*annulus_area[ifit])
error_number_caii = sqrt(total((a_massive_sf_sdev_nofit_2_mc[ifit]*annulus_area[ifit])^2))
error_mass_caii = error_number_caii*Mca40/Msolar
mass_caii = number_caii*Mca40/Msolar
print, 'high-mass star-forming M_caii: ', mass_caii, error_mass_caii, alog10(mass_caii), $
  (alog10(mass_caii+error_mass_caii)-alog10(mass_caii-error_mass_caii))/2.
;; total HI mass
;; solar: N(H)/N(Ca) = 10^5.69 -> M(H)/M(Ca) = 10^5.69/40
mass_HI_solar = 10.D0^5.69/40.*mass_caii
print, mass_hi_solar

number_caii = total(a_lowmass_sf[ifit].ew_nofit_2*annulus_area[ifit])
error_number_caii = sqrt(total((a_lowmass_sf_sdev_nofit_2_mc[ifit]*annulus_area[ifit])^2))
error_mass_caii = error_number_caii*Mca40/Msolar
mass_caii = number_caii*Mca40/Msolar
print, 'low-mass star forming M_caii: ', mass_caii, error_mass_caii, alog10(mass_caii), $
  (alog10(mass_caii+error_mass_caii)-alog10(mass_caii-error_mass_caii))/2.
;; total HI mass
;; solar: N(H)/N(Ca) = 10^5.69 -> M(H)/M(Ca) = 10^5.69/40
mass_HI_solar = 10.D0^5.69/40.*mass_caii
print, mass_hi_solar

thick=8
xthick=8
ythick=8
charsize=1.4
charthick=3

;xra=[11,120]
xra=[11,300]
;yra=[5.E-4, 8.E-1]
yra=[5.E-5, 8.E-1]
xtitle='r_p (kpc)' 
ytitle='W_0^{\lambda 3934} (Ca II) [\AA]'
title='Single Gaussian Line Profile Measurement'

xx = 10.^(findgen(1500)*0.0024+0.1)
print, minmax(xx)

if (DoScale) then begin
   fid_slope = -1.38
   fid_intercept = alog10(0.0032)
endif else begin
   fid_slope = -1.38
   fid_intercept = alog10(0.011)
endelse

;; get the fitting for SF
jhusdss_powerlaw_fit, a_lowmass_sf[ifit].rp/rvir*2., a_lowmass_sf[ifit].ew_nofit_2, sdev_nofit_2_mc[ifit], in_slope=fid_slope, in_intercept=10.^fid_intercept, $
   slope=lowmass_sf_slope, intercept=lowmass_sf_intercept, err_slope=err_lowmass_sf_slope, err_intercept=err_lowmass_sf_intercept, /fixslope
;; get the fitting for quiescent
jhusdss_powerlaw_fit, a_lowmass_quiescent[ifit].rp/rvir*2., a_lowmass_quiescent[ifit].ew_nofit_2, sdev_nofit_2_mc[ifit], in_slope=fid_slope, in_intercept=10.^fid_intercept, $
   slope=lowmass_quiescent_slope, intercept=lowmass_quiescent_intercept, err_slope=err_lowmass_quiescent_slope, err_intercept=err_lowmass_quiescent_intercept, /fixslope
;; get the fitting for SF
jhusdss_powerlaw_fit, a_massive_sf[ifit].rp/rvir*2., a_massive_sf[ifit].ew_nofit_2, sdev_nofit_2_mc[ifit], in_slope=fid_slope, in_intercept=10.^fid_intercept, $
   slope=massive_sf_slope, intercept=massive_sf_intercept, err_slope=err_massive_sf_slope, err_intercept=err_massive_sf_intercept, /fixslope
;; get the fitting for quiescent
jhusdss_powerlaw_fit, a_massive_quiescent[ifit].rp/rvir*2., a_massive_quiescent[ifit].ew_nofit_2, sdev_nofit_2_mc[ifit], in_slope=fid_slope, in_intercept=10.^fid_intercept, $
   slope=massive_quiescent_slope, intercept=massive_quiescent_intercept, err_slope=err_massive_quiescent_slope, err_intercept=err_massive_quiescent_intercept, /fixslope

print, massive_sf_slope, massive_sf_intercept, lowmass_sf_slope, lowmass_sf_intercept
print, massive_quiescent_slope, massive_quiescent_intercept, lowmass_quiescent_slope, lowmass_quiescent_intercept
print, err_massive_quiescent_intercept, err_massive_sf_intercept, err_lowmass_quiescent_intercept, err_lowmass_sf_intercept

;; total Ca II mass between rpmin and rpmax

rpmin=10.
rpmax=200.
kpc = 3.08D21 ;cm
factor = 1.13E12*1E8/0.6485/(3934.78)^2 ; turn REW to column density
Mca40 = 40.*1.67D-27
Msolar = 2.D30

;; Int (factor*out_intercept*(r/rvir)^slope*2.*pi*rdr) = factor*out_intercept*2.*pi*(1/rvir)^slope*r^(2.+slope)

rpmin=10.*10.^(0.2*(9.57-10.30))
rpmax=200.*10.^(0.2*(9.57-10.30))
lowmass_sf_number_caii = factor*lowmass_sf_intercept*2.*!pi*(1./rvir*2.)^fid_slope*kpc^2*(rpmax^(2.+fid_slope)-rpmin^(2.+fid_slope))/(2.+fid_slope)
rpmin=10.*10.^(0.2*(10.35-10.30))
rpmax=200.*10.^(0.2*(10.35-10.30))
lowmass_quiescent_number_caii = factor*lowmass_quiescent_intercept*2.*!pi*(1./rvir*2.)^fid_slope*kpc^2*(rpmax^(2.+fid_slope)-rpmin^(2.+fid_slope))/(2.+fid_slope)
rpmin=10.*10.^(0.2*(10.52-10.30))
rpmax=200.*10.^(0.2*(10.52-10.30))
massive_sf_number_caii = factor*massive_sf_intercept*2.*!pi*(1./rvir*2.)^fid_slope*kpc^2*(rpmax^(2.+fid_slope)-rpmin^(2.+fid_slope))/(2.+fid_slope)
rpmin=10.*10.^(0.2*(10.98-10.30))
rpmax=200.*10.^(0.2*(10.98-10.30))
massive_quiescent_number_caii = factor*massive_quiescent_intercept*2.*!pi*(1./rvir*2.)^fid_slope*kpc^2*(rpmax^(2.+fid_slope)-rpmin^(2.+fid_slope))/(2.+fid_slope)
lowmass_sf_mass_caii = lowmass_sf_number_caii*Mca40/Msolar
lowmass_quiescent_mass_caii = lowmass_quiescent_number_caii*Mca40/Msolar
massive_sf_mass_caii = massive_sf_number_caii*Mca40/Msolar
massive_quiescent_mass_caii = massive_quiescent_number_caii*Mca40/Msolar
print, lowmass_sf_mass_caii, lowmass_quiescent_mass_caii, massive_sf_mass_caii, massive_quiescent_mass_caii
err_lowmass_sf_number_caii = factor*err_lowmass_sf_intercept*2.*!pi*(1./rvir*2.)^fid_slope*kpc^2*(rpmax^(2.+fid_slope)-rpmin^(2.+fid_slope))/(2.+fid_slope)
err_lowmass_quiescent_number_caii = factor*err_lowmass_quiescent_intercept*2.*!pi*(1./rvir*2.)^fid_slope*kpc^2*(rpmax^(2.+fid_slope)-rpmin^(2.+fid_slope))/(2.+fid_slope)
err_massive_sf_number_caii = factor*err_massive_sf_intercept*2.*!pi*(1./rvir*2.)^fid_slope*kpc^2*(rpmax^(2.+fid_slope)-rpmin^(2.+fid_slope))/(2.+fid_slope)
err_massive_quiescent_number_caii = factor*err_massive_quiescent_intercept*2.*!pi*(1./rvir*2.)^fid_slope*kpc^2*(rpmax^(2.+fid_slope)-rpmin^(2.+fid_slope))/(2.+fid_slope)
err_lowmass_sf_mass_caii = err_lowmass_sf_number_caii*Mca40/Msolar
err_lowmass_quiescent_mass_caii = err_lowmass_quiescent_number_caii*Mca40/Msolar
err_massive_sf_mass_caii = err_massive_sf_number_caii*Mca40/Msolar
err_massive_quiescent_mass_caii = err_massive_quiescent_number_caii*Mca40/Msolar
print, err_lowmass_sf_mass_caii, err_lowmass_quiescent_mass_caii, err_massive_sf_mass_caii, err_massive_quiescent_mass_caii
;; total HI mass
;; solar: N(H)/N(Ca) = 10^5.69 -> M(H)/M(Ca) = 10^5.69/40
;mass_HI_solar = 10.D0^5.66/40.*mass_caii
;print, mass_hi_solar

stop

;lowmass_sf_intercept=alog10(lowmass_sf_intercept)
;lowmass_quiescent_intercept=alog10(lowmass_quiescent_intercept)
;massive_sf_intercept=alog10(massive_sf_intercept)
;massive_quiescent_intercept=alog10(massive_quiescent_intercept)


pos = [0.17, 0.20, 0.85, 0.90]
k_print, filename=psfile1, axis_char_scale=1.3, xsize=8, ysize=5
  djs_plot, a_massive_sf[i_indep].rp, a_massive_sf[i_indep].ew_nofit_2, xtitle=xtitle, ytitle=ytitle, $
      xra=xra, yra=yra, xst=1, yst=1, thick=thick, xthick=xthick, ythick=ythick, $
      charsize=charsize, charthick=charthick, /nodata, /xlog, /ylog, pos=pos, $
      ytickformat='jhusdss_tick_exponent'

  ibegin=0
  iend=9 
  polyfill, [a_massive_sf[ibegin:iend].rp, reverse(a_massive_sf[ibegin:iend].rp)], [sdev_nofit_2_mc[ibegin:iend], -reverse(sdev_nofit_2_mc[ibegin:iend])], color=djs_icolor('light gray')

  djs_oplot, xx, 10.^(alog10(xx/100.)*fid_slope+fid_intercept), thick=thick+2, color='blue'
  plotsym, 0, /fill
  djs_oplot, a_massive_sf[i_indep].rp, a_massive_sf[i_indep].ew_nofit_2,  $
      psym=8, symsize=2, color=djs_icolor('navy'), thick=thick, errcolor=djs_icolor('navy'), errthick=thick
  djs_oplot, a_massive_sf[i_indep].rp, (a_massive_sf[i_indep].ew_nofit_2>1.E-10),  $
      color=djs_icolor('navy'), thick=thick, errcolor=djs_icolor('navy'), errthick=thick
  oploterror, a_massive_sf[i_indep].rp, a_massive_sf[i_indep].ew_nofit_2, sdev_nofit_2_mc[i_indep], $
      psym=8, symsize=2, color=djs_icolor('navy'), thick=thick, errcolor=djs_icolor('navy'), errthick=thick
  plotsym, 4, /fill
  djs_oplot, a_lowmass_sf[i_indep].rp*1.05, (a_lowmass_sf[i_indep].ew_nofit_2), $
      psym=8, symsize=2, color=djs_icolor('light blue'), thick=thick, errcolor=djs_icolor('light blue'), errthick=thick
  djs_oplot, a_lowmass_sf[i_indep].rp*1.05, (a_lowmass_sf[i_indep].ew_nofit_2>1.E-10), $
      color=djs_icolor('light blue'), thick=thick, errcolor=djs_icolor('light blue'), errthick=thick
  oploterror, a_lowmass_sf[i_indep].rp*1.05, a_lowmass_sf[i_indep].ew_nofit_2, sdev_nofit_2_mc[i_indep], $
      psym=8, symsize=2, color=djs_icolor('light blue'), thick=thick, errcolor=djs_icolor('light blue'), errthick=thick
  plotsym, 0, thick=thick
  djs_oplot, a_massive_quiescent[i_indep].rp*1.10^2, a_massive_quiescent[i_indep].ew_nofit_2, $
      psym=8, symsize=1.5, color=djs_icolor('red'), thick=thick, errcolor=djs_icolor('red'), errthick=thick, errstyle=2
  djs_oplot, a_massive_quiescent[i_indep].rp*1.10^2, (a_massive_quiescent[i_indep].ew_nofit_2>1.E-10), $
      color=djs_icolor('red'), thick=thick, errcolor=djs_icolor('red'), errthick=thick, errstyle=2
  oploterror, a_massive_quiescent[i_indep].rp*1.10^2, a_massive_quiescent[i_indep].ew_nofit_2, sdev_nofit_2_mc[i_indep], $
      psym=8, symsize=1.5, color=djs_icolor('red'), thick=thick, errcolor=djs_icolor('red'), errthick=thick, errstyle=2
  plotsym, 4, thick=thick
  djs_oplot, a_lowmass_quiescent[i_indep].rp*1.15^3, a_lowmass_quiescent[i_indep].ew_nofit_2, $
      psym=8, symsize=1.5, color=djs_icolor('magenta'), thick=thick, errcolor=djs_icolor('magenta'), errthick=thick, errstyle=2
  djs_oplot, a_lowmass_quiescent[i_indep].rp*1.15^3, (a_lowmass_quiescent[i_indep].ew_nofit_2>1.E-10), $
      color=djs_icolor('magenta'), thick=thick, errcolor=djs_icolor('magenta'), errthick=thick, errstyle=2
  oploterror, a_lowmass_quiescent[i_indep].rp*1.15^3, a_lowmass_quiescent[i_indep].ew_nofit_2, sdev_nofit_2_mc[i_indep], $
      psym=8, symsize=1.5, color=djs_icolor('magenta'), thick=thick, errcolor=djs_icolor('magenta'), errthick=thick, errstyle=2

   ytitle1 = 'n (ca ii) [cm^{-2}]'
   ;; see eq. 9.15 in bruce drain
   factor = 1.13e12*1e8/0.6267/(3934.78)^2
   djs_axis, yaxis=1, yra=yra*factor, /ylog, ythick=thick, ytitle=ytitle1, $
       charsiz=charsize, charthick=charthick

  plotsym, 0, /fill
  djs_xyouts, 38, 5E-1, 'high-mass star-forming', color='navy', charsize=charsize-0.4, charthick=charthick
  djs_oplot, 100*[1,1], 5.5E-1*[1,1], color='navy', psym=8, symsize=1.3, thick=thick
  plotsym, 4, /fill
  djs_xyouts, 38, 5E-1*0.7, 'low-mass  star-forming', color='light blue', charsize=charsize-0.4, charthick=charthick
  djs_oplot, 100*[1,1], 5.5E-1*0.7*[1,1], color='light blue', psym=8, symsize=1.3, thick=thick
  plotsym, 0, thick=thick
  djs_xyouts, 38, 5E-1*0.7^2, 'high-mass  quiescent', color='red', charsize=charsize-0.4, charthick=charthick
  djs_oplot, 100*[1,1], 5.5E-1*0.7^2*[1,1], color='red', psym=8, symsize=1, thick=thick
  plotsym, 4, thick=thick
  djs_xyouts, 38, 5E-1*0.7^3, 'low-mass   quiescent', color='magenta', charsize=charsize-0.4, charthick=charthick
  djs_oplot, 100*[1,1], 5.5E-1*0.7^3*[1,1], color='magenta', psym=8, symsize=1, thick=thick

  djs_xyouts, 18.5, 0.00022, '20', charsize=charsize+0.4, charthick=charthick
;   djs_xyouts, 6, 0.005, 'average ca ii column density', charsize=charsize, charthick=charthick+1
;   djs_xyouts, 6, 0.003, 'profile in the cgm', charsize=charsize, charthick=charthick+1
;   djs_xyouts, 170, 0.000475, '200', charsize=charsize+0.35, charthick=charthick

k_end_print

if keyword_set(movingon) then begin
xra=[5,2000]
yra=[5.e-5, 6.e-1]
k_print, filename=psfile2, axis_char_scale=1.3, xsize=8, ysize=6
  djs_plot, a.rp, a.ew_nofit, xtitle=xtitle, ytitle=ytitle, title=title, $
      xra=xra, yra=yra, xst=1, yst=1, thick=thick, xthick=xthick, ythick=ythick, $
      charsize=charsize, charthick=charthick, /nodata, /xlog, /ylog, pos=pos
  ibegin=0
  iend=10 
  polyfill, [a[ibegin:iend].rp, reverse(a[ibegin:iend].rp)], [a[ibegin:iend].sdev_nofit_mc, -reverse(a[ibegin:iend].sdev_nofit_mc)], color=djs_icolor('light gray')
  djs_oplot, a[i_indep].rp, sdev_nofit_mc[i_indep], linestyle=1, thick=thick, color='gray'
  djs_oplot, xx, 10.^(alog10(xx/100.)*slope+intercept), thick=thick, color='light blue'
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
  djs_oplot, xx, 10.^(alog10(xx/100.)*slope+intercept), thick=thick, color='light blue'
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
  djs_oplot, xx, 10.^(alog10(xx/100.)*slope+intercept), thick=thick, color='red', linestyle=1
  djs_oplot, xx, 10.^(alog10(xx/100.)*fid_slope+fid_intercept), thick=thick+5, color='light blue', linestyle=0
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
  djs_oplot, xx, 10.^(alog10(xx/100.)*slope+intercept), thick=thick, color='light blue'
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
  djs_oplot, xx, 10.^(alog10(xx/100.)*slope+intercept), thick=thick
  djs_oplot, !x.crange, [0,0], color='light gray', thick=thick, linestyle=1
  oploterror, a[i_indep].rp, a[i_indep].ew_nofit_2, sdev_nofit_2_mc[i_indep], $
      psym=8, symsize=2, color=djs_icolor('light blue'), thick=thick, errcolor=djs_icolor('light blue'), errthick=thick
  djs_xyouts, 270, -0.0046, '300', charsize=charsize+0.4, charthick=charthick
k_end_print
endif

end
