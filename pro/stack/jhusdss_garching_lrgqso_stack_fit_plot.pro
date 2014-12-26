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
if (DoScale) then rvir=1000. else rvir=1000.

nmfver=106
stackpath = jhusdss_get_path(/nmfqso)+'/'+string(nmfver, format='(I3.3)')+'/Stack/'
infile = stackpath + jhusdss_garching_lrgqso_stack_filename(nmfver, boss=boss)
;infile = repstr(infile, '.fits', '_fit_bk.fits')
comp = mrdfits(infile,1)
if (keyword_set(freedom)) then begin
   boot_infile = repstr(infile, '.fits', '_bootstrap_fit_moment.fits')
   infile = repstr(infile, '.fits', '_fit.fits') 
endif else begin
   boot_infile = repstr(infile, '.fits', '_bootstrap_fit_samesigma_moment.fits')
   infile = repstr(infile, '.fits', '_fit_samesigma.fits')
endelse

a = mrdfits(infile, 1)
a_boot = mrdfits(boot_infile, 1)
;; The last two bins did not get calculated, manually set:
a_boot[29:31].err_ew = a[29:31].sdev_nofit_lr_1_s5_mc*2.0
a_boot[29:31].err_sigma = [2.2, 2.2, 2.2]

;a.ew_nofit_2 = a.ew_nofit_2*2.
;a.sdev_nofit_2_mc= a.sdev_nofit_2_mc*2.
;bootstrap_corr = fltarr(n_elements(a))+1.
;bootstrap_corr[1:18] = [1.41721,1.20982,1.23868,1.26567,0.989393,1.17334,0.947738,0.964353,0.909729,0.746676, 0.847485,0.761055,0.901519,1.16273,0.915482,1.11239,1.47041,0.81654]
;bootstrap_corr[1:6] = [1.41721,1.20982,1.23868,1.26567,1.,1.17334]
;bootstrap_corr[1:18] = [1.63686, 1.36330, 1.41082, 1.43796, 1.12306, 1.30891, 1.02353, 1.06469, 1.04978,0.857614,0.976261,0.917107, 1.02336, 1.22142,0.969464, 1.23497, 1.72883,0.961026]
;bootstrap_corr[1:9] = [1.63686, 1.36330, 1.41082, 1.43796, 1.12306, 1.30891, 1.02353, 1.06469, 1.04978]
;bootstrap_corr[1:16] = [1.41325, 1.18590, 1.22311, 1.23959, 0.958721,1.13099, 0.931704,0.990080,0.922841,0.757122,0.834926,0.830137,0.893792,1.04526, 0.801743,1.05222]

;bootstrap_corr[1:31]  = [1.39546, 1.59855, 1.26812, 1.17987, 0.898434,1.10963, 1.01236, 0.901248,0.980111,0.961529,0.780865,0.760658,1.01504, 0.904289,0.785204,1.26466, $
;                         1.11514, 0.907185,0.773530,0.797137,0.758422,0.762385,0.904971,0.709069,0.706981,0.884198,1.11134, 1.10, 1.10, 1.10, 1.10]
;sigma_bootstrap_corr[1:31]  = [1.39546, 1.59855, 1.26812, 1.17987, 0.898434,1.10963, 1.01236, 0.901248,0.980111,0.961529,0.780865,0.760658,1.01504, 0.904289,0.785204,1.26466, $
;                         1.11514, 0.907185,0.773530,0.797137,0.758422,0.762385,0.904971,0.709069,0.706981,0.884198,1.11134, 1.10, 1.10, 1.10, 1.10]

caii_infile = stackpath + jhusdss_garching_galqso_stack_filename(nmfver, boss=boss)
caii_infile = repstr(caii_infile, '.fits', '_fit.fits')
a_caii = mrdfits(caii_infile, 1)
a_caii_i_indep = [[0, 1], lindgen(11)*2+2]
a_caii_i_fit = a_caii_i_indep[1:7]
 
psfile1 = repstr(infile, '.fits', '_1.ps')
psfile11 = repstr(infile, '.fits', '_11.ps')
psfile2 = repstr(infile, '.fits', '_2.ps')
psfile3 = repstr(infile, '.fits', '_3.ps')
psfile4 = repstr(infile, '.fits', '_4.ps')
psfile5 = repstr(infile, '.fits', '_5.ps')
psfile6 = repstr(infile, '.fits', '_6.ps')
psfile7 = repstr(infile, '.fits', '_7.ps')

;; manually check, independent bins
;i_indep = [[0], lindgen(17)*2+2]
i_indep = lindgen(16)*2+1
;i_indep = lindgen(16)*1+2
;i_indep[15] = i_indep[15]-1
n_indep = n_elements(i_indep)

;best_fit slope
;best_fit slope
ifit = i_indep[0:15]

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
slope = -1.3
intercept = -2.

jhusdss_powerlaw_fit, a[ifit].rp/rvir*2., a[ifit].ew_nofit_2, a[ifit].sdev_nofit_2_mc, in_slope=slope, in_intercept=10.^intercept, $
;jhusdss_powerlaw_fit, a[ifit].rp, a[ifit].ew_nofit, sdev_nofit_mc[ifit], in_slope=slope, in_intercept=10.^intercept, $
   slope=out_slope, intercept=out_intercept, err_slope=err_slope, err_intercept=err_intercept
print, slope, intercept
print, out_slope, alog10(out_intercept)
slope = out_slope
intercept = alog10(out_intercept)
print, err_slope, out_intercept, err_intercept

print, slope, intercept
slope_caii = -1.38
intercept_caii = 0.011

rp_cs = [31., 63, 103.]
ew_cs_lya = [2.01, 1.23, 0.92]
ew_cs_lya_err = [0.15, 0.20, 0.12]

ew_cs_siii1260 = [0.42, 0.41, 0.05];upper limit
ew_cs_siii1260_err = [0.06, 0.09, 0.05];upper limit

ew_cs_cii1334 = [0.90, 0.67, 0.12];upper limit
ew_cs_cii1334_err = [0.08, 0.12, 0.12];upper limit

ew_cs_siiv1393 = [0.39, 0.19, 0.12]
ew_cs_siiv1393_err = [0.08, 0.08, 0.06]

ew_cs_siii1526 = [0.37, 0.15, 0.04];upper limits
ew_cs_siii1526_err = [0.06, 0.15, 0.04];upper limits

ew_cs_civ1549 = [2.13, 1.18, 0.13]
ew_cs_civ1549_err = [0.15, 0.15, 0.05]

ew_cs_alii1670 = [0.40, 0.20, 0.10];upper limits
ew_cs_alii1670_err = [0.08, 0.20, 0.10]

rp_bor = [29., 50.1, 70.9, 102.3, 140.2];
ew_bor = [0.15, 0.13, 0.20, 0.09, 0.11]
ew_bor_err = [0.040, 0.04, 0.02, 0.05, 0.02]

rp_bor = [37., 58.3, 72.9]
ew_bor = [0.20, 0.13, 0.09]
ew_bor_err = [0.07, 0.03, 0.03]

;ew_bor = [0.165, 0.13, 0.18, 0.075]
;ew_bor_err = [0.032, 0.022, 0.014, 0.010]
;ew_bor_err = [0.05, 0.04, 0.02, 0.020]
readcol, '/home/gz323/Code/jhu-sdss/pro/mandelbaum2006_sm6.dat', msm7_rp, msm7_Sigma, msm7_Sigma_err, format='(f,f,f)'
sm_factor=1.
if (keyword_set(read_sm7)) then begin
   readcol, '/home/gz323/Code/jhu-sdss/pro/mandelbaum2006_sm7.dat', msm7_rp, msm7_Sigma, msm7_Sigma_err, format='(f,f,f)'
   sm_factor=0.7
endif
msm7_rp = msm7_rp/0.7*sm_factor
msm7_Sigma = msm7_Sigma*0.7*sm_factor
msm7_Sigma_err = msm7_Sigma_err*0.7*sm_factor

;; turn REW to surface density (minimum) in MSun/pc^2
pccm = 3.08D18 ;cm
;factor = 1.13E20/0.6155/2803.^2
factor = 1.13E20/0.3058/2803.^2
Mmg24 = 24.305*1.67D-27
Msolar = 2.D30
mass_factor = factor*Mmg24/Msolar*pccm^2*10^(12.-7.6)*(9.8/7.4)/2. ;; H => H&He
;; Mg solar 10^7.60
;; ionization*refractory correction
fbaryon = 0.167*0.0018

; print out the table
for i=0L, n_elements(i_indep)-1L do begin
  print, "($"+string(comp[i_indep[i]].rp_min, format='(f6.3)')+","+string(comp[i_indep[i]].rp_max, format='(f6.3)')+"$] & $"+ $
         string(comp[i_indep[i]].rp/1E3, format='(f6.3)')+"$ & $"+ $
         string(comp[i_indep[i]].npairs, format='(i7)')+"$ & $"+ $
         string(a[i_indep[i]].ew*(1.+1./a[i_indep[i]].line_ratio)*1E3, format='(f7.2)')+" \pm "+string(a_boot[i_indep[i]].err_ew*1E3, format='(f7.2)')+"$ \\"
;        string(a[i_indep[i]].ew*(1.+1./a[i_indep[i]].line_ratio)*1E3, format='(f7.2)')+" \pm "+string(a[i_indep[i]].sdev_nofit_lr_1_s5_mc*2.0*bootstrap_corr[i_indep[i]]*1E3, format='(f7.2)')+"$ \\"
;        string(a[i_indep[i]].ew*(1.+1./a[i_indep[i]].line_ratio)/3.*factor*Mmg24/Msolar*pccm^2*1E9, format='(f7.2)')+" \pm "+string(a[i_indep[i]].sdev_nofit_lr_1_s5_mc*2.0*bootstrap_corr[i_indep[i]]/3.*factor*Mmg24/Msolar*pccm^2*1E9, format='(f7.2)')+"$ \\"
endfor
print, " "

for i=0L, n_elements(i_indep)-1L do begin
  print, "($"+string(comp[i_indep[i]].rp_min, format='(f6.3)')+","+string(comp[i_indep[i]].rp_max, format='(f6.3)')+"$] & $"+ $
         string(comp[i_indep[i]].rp/1E3, format='(f6.3)')+"$ & $"+ $
         string(comp[i_indep[i]].npairs, format='(i7)')+"$ & $"+ $
         string(a[i_indep[i]].ew*(1.+1./a[i_indep[i]].line_ratio)*1E3, format='(f7.2)')+" $ & $  "+string(a_boot[i_indep[i]].err_ew*1E3, format='(f7.2)')+"$ \\"
;        string(a[i_indep[i]].ew*(1.+1./a[i_indep[i]].line_ratio)*1E3, format='(f7.2)')+" $ & $  "+string(a[i_indep[i]].sdev_nofit_lr_1_s5_mc*2.0*bootstrap_corr[i_indep[i]]*1E3, format='(f7.2)')+"$ \\"
;        string(a[i_indep[i]].ew*(1.+1./a[i_indep[i]].line_ratio)/3.*factor*Mmg24/Msolar*pccm^2*1E9, format='(f7.2)')+" \pm "+string(a[i_indep[i]].sdev_nofit_lr_1_s5_mc*2.0*bootstrap_corr[i_indep[i]]/3.*factor*Mmg24/Msolar*pccm^2*1E9, format='(f7.2)')+"$ \\"
endfor
print, " "

for i=0L, n_elements(i_indep)-1L do begin
  print, string(comp[i_indep[i]].rp/1E3, format='(f6.3)')+" "+$
         string(comp[i_indep[i]].npairs, format='(i7)')+" "+$
         string(a[i_indep[i]].ew*(1.+1./a[i_indep[i]].line_ratio)*1E3, format='(f7.3)')+" "+$
         string(a_boot[i_indep[i]].err_ew*1E3, format='(f7.3)')+" "+$
;        string(a[i_indep[i]].sdev_nofit_lr_1_s5_mc*2.0*bootstrap_corr[i_indep[i]]*1E3, format='(f7.3)')+" "+$
         string(a[i_indep[i]].ew*(1.+1./a[i_indep[i]].line_ratio)/3.*factor*Mmg24/Msolar*pccm^2*1E9, format='(f7.2)')+" "+$
         string(a_boot[i_indep[i]].err_ew/3.*factor*Mmg24/Msolar*pccm^2*1E9, format='(f7.2)')
;        string(a[i_indep[i]].sdev_nofit_lr_1_s5_mc*2.0*bootstrap_corr[i_indep[i]]/3.*factor*Mmg24/Msolar*pccm^2*1E9, format='(f7.2)')
endfor

sigma_factor = 1./alog(10.)/2800.*1.D+4*69.03
for i=0L, n_elements(i_indep)-1L do begin
  print, string(comp[i_indep[i]].rp/1E3, format='(f6.3)')+" "+$
         string(comp[i_indep[i]].npairs, format='(i7)')+" "+$
         string(sqrt((a[i_indep[i]].sigma*sigma_factor)^2-69.03^2), format='(f6.2)')+" "+$
         string(a_boot[i_indep[i]].err_sigma*sigma_factor*(a[i_indep[i]].sigma*sigma_factor)/sqrt((a[i_indep[i]].sigma*sigma_factor)^2-69.03^2), format='(f6.2)')
;        string(a[i_indep[i]].err_sigma*sigma_factor*(a[i_indep[i]].sigma*sigma_factor)/sqrt((a[i_indep[i]].sigma*sigma_factor)^2-69.03^2), format='(f6.2)')
endfor

;; print out the table
;for i=0L, n_elements(i_indep)-1L do begin
;  print, "($"+string(comp[i_indep[i]].rp_min, format='(f6.3)')+","+string(comp[i_indep[i]].rp_max, format='(f6.3)')+"$] & $"+ $
;         string(comp[i_indep[i]].rp/1E3, format='(f6.3)')+"$ & $"+ $
;         string(comp[i_indep[i]].npairs, format='(i7)')+"$ & $"+ $
;         string(a[i_indep[i]].ew*(1.+1./a[i_indep[i]].line_ratio)*1E3, format='(f7.2)')+" \pm "+string(a[i_indep[i]].sdev_nofit_lr_1_s5_mc*2.0*bootstrap_corr[i_indep[i]]*1E3, format='(f7.2)')+"$ & $"+ $
;         string(a[i_indep[i]].ew*(1.+1./a[i_indep[i]].line_ratio)/3.*factor*Mmg24/Msolar*pccm^2*1E9, format='(f7.2)')+" \pm "+string(a[i_indep[i]].sdev_nofit_lr_1_s5_mc*2.0*bootstrap_corr[i_indep[i]]/3.*factor*Mmg24/Msolar*pccm^2*1E9, format='(f7.2)')+"$ \\"
;endfor
;print, " "
;for i=0L, n_elements(i_indep)-1L do begin
;  print, string(comp[i_indep[i]].rp/1E3, format='(f6.3)')+" "+$
;         string(comp[i_indep[i]].npairs, format='(i7)')+" "+$
;         string(a[i_indep[i]].ew*(1.+1./a[i_indep[i]].line_ratio)*1E3, format='(f7.3)')+" "+$
;         string(a[i_indep[i]].sdev_nofit_lr_1_s5_mc*2.0*bootstrap_corr[i_indep[i]]*1E3, format='(f7.3)')+" "+$
;         string(a[i_indep[i]].ew*(1.+1./a[i_indep[i]].line_ratio)/3.*factor*Mmg24/Msolar*pccm^2*1E9, format='(f7.2)')+" "+$
;         string(a[i_indep[i]].sdev_nofit_lr_1_s5_mc*2.0*bootstrap_corr[i_indep[i]]/3.*factor*Mmg24/Msolar*pccm^2*1E9, format='(f7.2)')
;endfor

;rp_bor = [40., 60., 75.]
;ew_bor = [0.20, 0.13, 0.09]
;ew_bor_err = [0.07, 0.03, 0.03]
;ew_bor_blue = [0.88, 0.32, 0.11]
;ew_bor_blue_err = [0.07, 0.05, 0.04]

thick=8
xthick=8
ythick=8
charsize=1.4
charthick=3

xra=[15,15000]
yra=[3.E-4, 7.E-1]
xtitle='r_p (kpc)' 
ytitle='W_0^{\lambda 2796} (Mg II) [\AA]'
title='Single Gaussian Line Profile Measurement'

xx = 10.^(findgen(1500)*0.0024+1.)

plotsym, 0, /fill

;pos = [0.20, 0.15, 0.95, 0.9]
 pos = [0.17, 0.20, 0.85, 0.90]
k_print, filename=psfile1, axis_char_scale=1.3, xsize=7, ysize=6
  djs_plot, a[i_indep].rp, a[i_indep].ew_nofit_2, xtitle=xtitle, ytitle=ytitle, $
      xra=xra, yra=yra, xst=1, yst=1, thick=thick, xthick=xthick, ythick=ythick, $
      charsize=charsize, charthick=charthick, /nodata, /xlog, /ylog, pos=pos, $
      ytickformat='jhusdss_tick_exponent'
  djs_oplot, xx, 10.^(alog10(xx/rvir*2.)*slope+intercept), thick=thick, linestyle=1, color='blue'
  oploterror, a[i_indep].rp, a[i_indep].ew_nofit, a[i_indep].sdev_nofit_mc, $
      psym=8, symsize=2, color=djs_icolor('blue'), thick=thick, errcolor=djs_icolor('blue'), errthick=thick
; djs_oplot, rp_bor, ew_bor, psym=5, symsize=1.8, thick=thick-1, color='red'
  oploterror, rp_bor, ew_bor, ew_bor_err, psym=5, symsize=1.8, thick=thick-1, color=djs_icolor('gray'), $
     errcolor=djs_icolor('gray'), errthick=thick


;  ytitle1 = 'N (Mg II) [cm^{-2}]'
   ;; see Eq. 9.15 in Bruce Drain
;  factor = 1.13E12*1E8/0.6155/(2976.35)^2
;  djs_axis, yaxis=1, yra=yra*factor, /ylog, ythick=thick, ytitle=ytitle1, $
;      charsiz=charsize, charthick=charthick

    djs_xyouts, 30, 0.0007, 'LRG-Mg II Correlation', charsize=charsize, charthick=charthick+1
;   djs_xyouts, 6, 0.003, 'Profile in the CGM', charsize=charsize, charthick=charthick+1
;   djs_xyouts, 170, 0.000475, '200', charsize=charsize+0.35, charthick=charthick

  djs_xyouts, 300., 0.2, $
             'Bordoloi et al. 2011 (red, group)', color='gray', charsize=charsize, charthick=charthick
  djs_xyouts, 300., 0.3, $
             'Zhu, Menard & BOSS', color='blue', charsize=charsize, charthick=charthick

k_end_print

k_print, filename=psfile11, axis_char_scale=1.3, xsize=7, ysize=6
  djs_plot, a[i_indep].rp, a[i_indep].ew_nofit_2, xtitle=xtitle, ytitle=ytitle, $
      xra=xra, yra=yra, xst=1, yst=1, thick=thick, xthick=xthick, ythick=ythick, $
      charsize=charsize, charthick=charthick, /nodata, /xlog, /ylog, pos=pos, $
      ytickformat='jhusdss_tick_exponent'
  djs_oplot, xx, 10.^(alog10(xx/rvir*2.)*slope+intercept), thick=thick, linestyle=1, color='blue'
  oploterror, a[i_indep].rp, a[i_indep].ew_nofit, a[i_indep].sdev_nofit_mc, $
      psym=8, symsize=2, color=djs_icolor('blue'), thick=thick, errcolor=djs_icolor('blue'), errthick=thick
; djs_oplot, rp_bor, ew_bor, psym=5, symsize=2, thick=thick, color='red'
  oploterror, rp_bor, ew_bor, ew_bor_err, psym=5, symsize=1.8, thick=thick-1, color=djs_icolor('gray'), $
     errcolor=djs_icolor('gray'), errthick=thick

;  ytitle1 = 'N (Mg II) [cm^{-2}]'
   ;; see Eq. 9.15 in Bruce Drain
;  factor = 1.13E12*1E8/0.6155/(2976.35)^2
;  djs_axis, yaxis=1, yra=yra*factor, /ylog, ythick=thick, ytitle=ytitle1, $
;      charsiz=charsize, charthick=charthick

    djs_xyouts, 30, 0.0007, 'LRG-Mg II Correlation', charsize=charsize, charthick=charthick+1
;   djs_xyouts, 6, 0.003, 'Profile in the CGM', charsize=charsize, charthick=charthick+1
;   djs_xyouts, 170, 0.000475, '200', charsize=charsize+0.35, charthick=charthick

; djs_xyouts, 300., 0.2, $
;            'Bordoloi et al. 2011', color='red', charsize=charsize, charthick=charthick
  djs_xyouts, 300., 0.3, $
             'Zhu, Menard & BOSS', color='blue', charsize=charsize, charthick=charthick

k_end_print

;k_print, filename=psfile1, axis_char_scale=1.3, xsize=8, ysize=6
; djs_plot, a[i_indep].rp, a[i_indep].ew_nofit, xtitle=xtitle, ytitle=ytitle, title=title, $
;     xra=xra, yra=yra, xst=1, yst=1, thick=thick, xthick=xthick, ythick=ythick, $
;     charsize=charsize, charthick=charthick, /nodata, /xlog, /ylog, pos=pos
; djs_oplot, xx, 10.^(alog10(xx)*slope+intercept), thick=thick
; djs_oplot, xx, 10.^(alog10(xx)*slope_caii+intercept_caii), thick=thick, linestyle=1
; djs_oplot, rp_bor, ew_bor, psym=5, symsize=2, thick=thick, color='red'
; oploterror, a[i_indep].rp, a[i_indep].ew_nofit, a[i_indep].sdev_nofit_mc, $
;     psym=8, symsize=2, color=djs_icolor('blue'), thick=thick, errcolor=djs_icolor('blue'), errthick=thick
; djs_xyouts, 500., 0.4, $
;            'Bordoloi et al. 2011', color='red', charsize=charsize, charthick=charthick
; djs_xyouts, 500., 0.7, $
;            'Zhu & Menard this work', color='blue', charsize=charsize, charthick=charthick
;_end_print

k_print, filename=psfile2, axis_char_scale=1.3, xsize=8, ysize=6
  djs_plot, a.rp, a.ew_nofit, xtitle=xtitle, ytitle=ytitle, $
      xra=xra, yra=yra, xst=1, yst=1, thick=thick, xthick=xthick, ythick=ythick, $
      charsize=charsize, charthick=charthick, /nodata, /xlog, /ylog, pos=pos
  ibegin=1
  iend=23 
  polyfill, [a[ibegin:iend].rp, reverse(a[ibegin:iend].rp)], [a[ibegin:iend].sdev_nofit_mc, -reverse(a[ibegin:iend].sdev_nofit_mc)], color=djs_icolor('light gray')
  djs_oplot, xx, 10.^(alog10(xx/rvir*2.)*slope+intercept), thick=thick, linestyle=1, color='blue'
  djs_oplot, rp_bor, ew_bor, psym=5, symsize=2, thick=thick, color='red'
  oploterror, a.rp, a.ew_nofit, a.sdev_nofit_mc, $
      psym=8, symsize=2, color=djs_icolor('blue'), thick=thick, errcolor=djs_icolor('blue'), errthick=thick
k_end_print

k_print, filename=psfile3, axis_char_scale=1.3, xsize=8, ysize=6
  xra=[220, 20000]
  yra=[-4E-3, 4E-3]
  djs_plot, a.rp, a.ew_nofit, xtitle=xtitle, ytitle=ytitle, $
      xra=xra, yra=yra, xst=1, yst=1, thick=thick, xthick=xthick, ythick=ythick, $
      charsize=charsize, charthick=charthick, /nodata, /xlog, pos=pos;, /ylog
  ibegin=13
  iend=23 
  polyfill, [a[ibegin:iend].rp, reverse(a[ibegin:iend].rp)], [a[ibegin:iend].sdev_nofit_mc, -reverse(a[ibegin:iend].sdev_nofit_mc)], color=djs_icolor('light gray')
  djs_oplot, xx, 10.^(alog10(xx/rvir*2.)*slope+intercept), thick=thick, linestyle=1, color='blue'
  oploterror, a.rp, a.ew_nofit, a.sdev_nofit_mc, $
      psym=8, symsize=2, color=djs_icolor('blue'), thick=thick, errcolor=djs_icolor('blue'), errthick=thick

k_end_print

xra=[8,33000]/1E3
yra=[1.E-4, 2.E-0]
xtitle='r_p' 
;ytitle='W_0^{\lambda 2796} (Mg II) [\AA]'
ytitle='<W_0> (\AA)'
title='Double Gaussian Line Profile Measurement'

pos = [0.17, 0.20, 0.85, 0.90]
plotsym, 0, /fill
k_print, filename=psfile4, axis_char_scale=1.3, xsize=10, ysize=6
  djs_plot, a[i_indep].rp/1E3, a[i_indep].ew_nofit_2, xtitle=xtitle, ytitle=ytitle, $
      xra=xra, yra=yra, thick=thick, xthick=xthick, ythick=ythick, $
      charsize=charsize, charthick=charthick, /nodata, /xlog, /ylog, pos=pos, $
      ytickformat='jhusdss_tick_exponent', ystyle=9, xstyle=5

  if (DoScale) then begin
      djs_axis, xaxis=0, xtitle='r_p/r_{vir}', charsize=charsize, charthick=charthick, xra=xra/rvir*1E3, xthick=xthick
      djs_axis, xaxis=1, charsize=charsize, charthick=charthick, xra=xra/rvir*1E3, xthick=xthick, xtickformat='(A1)'
  endif else begin
      xtickv = [0.01, 0.1, 1., 10.]
      xticknames = ['10 kpc', '100 kpc', '1 Mpc', '10 Mpc']
      djs_axis, xaxis=0, xtitle=xtitle, charsize=charsize, charthick=charthick, xra=xra, xthick=xthick, $
          xtickv=xtickv, xtickname=xticknames;, xticks=n_elements(xtickv)
      djs_axis, xaxis=1, xtitle='', charsize=charsize, charthick=charthick, xra=xra, xthick=xthick, xtickformat='(A1)'
  endelse
  djs_xyouts, 0.10, 2.35, 'Galaxy-Gas Correlation', charsize=charsize, charthick=charthick

; ixx = where(xx lt 250., comp=jxx)
; djs_oplot, xx[ixx], 10.^(alog10(xx[ixx]/rvir*2.)*slope+intercept), thick=thick, color='blue'
; djs_oplot, xx[jxx], 10.^(alog10(xx[jxx]/rvir*2.)*slope+intercept), thick=thick, color='blue', linestyle=1

; djs_oplot, xx, 10.^(alog10(xx/rvir*2.)*slope+intercept), thick=thick, linestyle=1, color='gray'
; djs_oplot, xx, 10.^(alog10(xx)*slope_caii+intercept_caii), thick=thick, linestyle=1
  plotsym, 4
  oploterror, rp_bor/1E3, ew_bor, ew_bor_err, psym=8, symsize=1., color=djs_icolor('orange'), $
     errcolor=djs_icolor('orange')
  djs_oplot, [0.015, 0.015], [4.5E-4, 4.5E-4], psym=8, symsiz=1., color=djs_icolor('orange')
  plotsym, 3
  oploterror, rp_cs/1E3, ew_cs_siiv1393, ew_cs_siiv1393_err, psym=8, symsize=1., color=djs_icolor('magenta'), $
     errcolor=djs_icolor('magenta')
  djs_oplot, [0.015, 0.015], [7.0E-4, 7.0E-4], psym=8, symsiz=1., color=djs_icolor('magenta')
; plotsym, 4
; oploterror, rp_cs/1E3, ew_cs_civ1549, ew_cs_civ1549_err, psym=8, symsize=1., color=djs_icolor('magenta'), $
;    errcolor=djs_icolor('magenta')
; djs_oplot, [0.015, 0.015], [4.5E-4, 4.5E-4], psym=8, symsiz=1., color=djs_icolor('magenta')
; oploterror, rp_bor, ew_bor_blue, ew_bor_blue_err, psym=8, symsize=1.8, thick=thick-1, color=djs_icolor('gray'), $
;    errcolor=djs_icolor('gray'), errthick=thick
  plotsym, 0, /fill
; oploterror, a[i_indep].rp/1E3, a[i_indep].ew*(1./a[i_indep].line_ratio+1.), a[i_indep].sdev_nofit_lr_1_s5_mc*2.0*bootstrap_corr[i_indep], $
  oploterror, a[i_indep].rp/1E3, a[i_indep].ew*(1./a[i_indep].line_ratio+1.), a_boot[i_indep].err_ew, $
      psym=8, symsize=1.4, color=djs_icolor('blue'), thick=thick, errcolor=djs_icolor('blue'), errthick=thick, hatlength=!d.x_vsize/80.
  djs_oplot, [3.20, 3.20], [0.525, 0.525], psym=8, symsiz=1.2, color=djs_icolor('blue'), thick=thick

  ;; Mg I
  plotsym, 8, /fill
  oploterror, a[i_indep[0]].rp/1E3, [0.0831], [0.0645], $
      psym=8, symsize=1.2, color=djs_icolor('blue'), thick=2, errcolor=djs_icolor('blue'), errthick=2, hatlength=!d.x_vsize/80.
  djs_oplot, [3.20, 3.20], [0.325, 0.325], psym=8, symsiz=1., color=djs_icolor('blue')
; djs_oplot, [0.015, 0.015], [4.5E-4, 4.5E-4], psym=8, symsiz=1.2, color=djs_icolor('blue'), symthick=2

; plotsym, 0
; oploterror, a[i_indep].rp/1E3, a[i_indep].ew_nofit_2_s4*1.5, a[i_indep].sdev_nofit_2_s4_mc*1.5, $
;     psym=8, symsize=1.4, color=djs_icolor('dark green'), thick=thick, errcolor=djs_icolor('dark green'), errthick=thick, hatlength=!d.x_vsize/80.
; plotsym, 0
; oploterror, a[i_indep].rp/1E3, a[i_indep].ew_nofit_lr_1_s5*2.0, a[i_indep].sdev_nofit_lr_1_s5_mc*2.0, $
;     psym=8, symsize=1.4, color=djs_icolor('red'), thick=thick, errcolor=djs_icolor('red'), errthick=thick, hatlength=!d.x_vsize/80.

  plotsym, 0
; plotsym, 0
; oploterror, a[i_indep[0:4]].rp/1E3, a[i_indep[0:4]].ew*(1./a[i_indep[0:4]].line_ratio+1.)*3./2., a[i_indep[0:4]].sdev_nofit_2_s4_mc*1.5, $
;     psym=8, symsize=1.3, color=djs_icolor('blue'), thick=thick, errcolor=djs_icolor('blue'), errthick=thick, hatlength=!d.x_vsize/50.
  oploterror, a_caii[a_caii_i_fit].rp/1E3, a_caii[a_caii_i_fit].ew_nofit_2, a_caii[a_caii_i_fit].sdev_nofit_2_mc, $
     psym=8, symsize=1., color=djs_icolor('light blue'), errcolor=djs_icolor('light blue')
  djs_oplot, [0.015, 0.015], [2.75E-4, 2.75E-4], psym=8, symsiz=1., color=djs_icolor('light blue')
; plotsym, 0
; oploterror, msm7_rp/1E3, msm7_Sigma/mass_factor*fbaryon, msm7_Sigma_err/mass_factor*fbaryon, $
;     psym=8, symsize=2, color=djs_icolor('gray'), thick=thick, errcolor=djs_icolor('gray'), errthick=thick
; djs_xyouts, 1.25, 0.475, $
  djs_xyouts, 0.018, 6.00E-4, $
             'Si IV around LBGs at z\sim2.2, (Steidel et al. 2010)', color='magenta', charsize=charsize-0.3, charthick=charthick
  djs_xyouts, 0.018, 4.00E-4, $
             'Mg II around red galaxies at z\sim0.7 (Bordoloi et al. 2011)', color='orange', charsize=charsize-0.3, charthick=charthick
; djs_xyouts, 1.25, 0.850, $
;            'Zhu et al. (2013) [Mg II]', color='blue', charsize=charsize-0.2, charthick=charthick+0.3
; djs_xyouts, 3.00, 0.850, $
  djs_xyouts, 1.5, 0.475, $
             'Mg II    Mg I around Luminous Red Galaxies at \sim0.5 Zhu et al. (2013)', color='blue', charsize=charsize-0.2, charthick=charthick+0.3
; djs_xyouts, 4.0, 0.475, $
;            'Mg II', color='blue', charsize=charsize-0.2, charthick=charthick+0.3
; djs_xyouts, 4.0, 0.275, $
;            'Mg I', color='blue', charsize=charsize-0.2, charthick=charthick+0.3
; djs_xyouts, 0.018, 4.00E-4, $
;            'Zhu et al. (2013) [Mg I]', color='blue', charsize=charsize-0.2, charthick=charthick+0.3
  djs_xyouts, 0.018, 2.5E-4, $
             ANSI_value('Ca II around all galaxies at z\sim0.1 (Zhu & M!3'+string("351B)+'!Xnard 2013)'), color='light blue', charsize=charsize-0.3, charthick=charthick
;            ANSI_value('Zhu & M!3'+string("351B)+'!Xnard (2013,')+' Ca II)', color='gray', charsize=charsize-0.3, charthick=charthick

; djs_xyouts, 300./1E3, 0.4, $
;            "Mandelbaum & SDSS (Dark Matter)", color='gray', charsize=charsize-0.3, charthick=charthick
; djs_xyouts, 300./1E3, 0.7, $
;            "Zhu, Menard & SDSS ('Dark' Baryons)", color='blue', charsize=charsize-0.3, charthick=charthick

; text_relation = '<W_0>='+string(out_intercept, format='(f5.3)')+'\times(r_p/500 kpc)^{'+string(slope, format='(f5.2)')+'}'
; djs_oplot, [15,25], [1,1]*4.5E-4, thick=thick, linestyle=1, color='gray'
; djs_xyouts, 30, 4.E-4, text_relation, charsize=charsize-0.00, charthick=charthick-0.0, color='black'

   ytitle1 = '<\Sigma_{Mg II}> (M_\odot pc^{-2})'
   ;; see Eq. 9.15 in Bruce Drain
   djs_axis, yaxis=1, yra=yra/3.*factor*Mmg24/Msolar*pccm^2, /ylog, ythick=thick, ytitle=ytitle1, $
       charsiz=charsize, charthick=charthick

k_end_print

xra=[8,48000]/1E3
yra=[6.E-4, 5.E-0]
pos = [0.20, 0.20, 0.90, 0.90]
plotsym, 0, /fill
k_print, filename=repstr(psfile4, '.ps', '_2.ps'), axis_char_scale=1.3, xsize=10, ysize=6
  djs_plot, a[i_indep].rp/1E3, a[i_indep].ew_nofit_2, xtitle=xtitle, ytitle=ytitle, $
      xra=xra, yra=yra, thick=thick, xthick=xthick, ythick=ythick, $
      charsize=charsize, charthick=charthick, /nodata, /xlog, /ylog, pos=pos, $
      ytickformat='jhusdss_tick_exponent', ystyle=1, xstyle=5

  if (DoScale) then begin
      djs_axis, xaxis=0, xtitle='r_p/r_{vir}', charsize=charsize, charthick=charthick, xra=xra/rvir*1E3, xthick=xthick
      djs_axis, xaxis=1, charsize=charsize, charthick=charthick, xra=xra/rvir*1E3, xthick=xthick, xtickformat='(A1)'
  endif else begin
      xtickv = [0.01, 0.1, 1., 10.]
      xticknames = ['10 kpc', '100 kpc', '1 Mpc', '10 Mpc']
      djs_axis, xaxis=0, xtitle=xtitle, charsize=charsize, charthick=charthick, xra=xra, xthick=xthick, $
          xtickv=xtickv, xtickname=xticknames;, xticks=n_elements(xtickv)
      djs_axis, xaxis=1, xtitle='', charsize=charsize, charthick=charthick, xra=xra, xthick=xthick, xtickformat='(A1)'
  endelse
  djs_xyouts, 0.055, 5.7, 'The Galaxy-Gas (Metal) Correlation', charsize=charsize, charthick=charthick

  plotsym, 3
  oploterror, rp_bor/1E3, ew_bor, ew_bor_err, psym=8, symsize=1., color=djs_icolor('orange'), $
     errcolor=djs_icolor('orange')
; djs_oplot, [0.011, 0.011], [2.2E-4, 2.2E-4], psym=8, symsiz=1., color=djs_icolor('orange')
  djs_oplot, [0.10, 0.10], [5.3E-1, 5.3E-1]*1.5, psym=8, symsiz=1., color=djs_icolor('orange')

  plotsym, 4
  oploterror, rp_cs/1E3, ew_cs_siiv1393, ew_cs_siiv1393_err, psym=8, symsize=1., color=djs_icolor('magenta'), $
     errcolor=djs_icolor('magenta')
; djs_oplot, [0.011, 0.011], [3.5E-4, 3.5E-4], psym=8, symsiz=1., color=djs_icolor('magenta')
  djs_oplot, [0.10, 0.10], [5.3E-1, 5.3E-1]*1.5^2, psym=8, symsiz=1., color=djs_icolor('magenta')

  plotsym, 0
  oploterror, a_caii[a_caii_i_fit].rp/1E3, a_caii[a_caii_i_fit].ew_nofit_2, a_caii[a_caii_i_fit].sdev_nofit_2_mc, $
     psym=8, symsize=1., color=djs_icolor('light blue'), errcolor=djs_icolor('light blue')
; djs_oplot, [0.011, 0.011], [1.35E-4, 1.35E-4], psym=8, symsiz=1., color=djs_icolor('light blue')
  djs_oplot, [0.10, 0.10], [5.3E-1, 5.3E-1], psym=8, symsiz=1., color=djs_icolor('light blue')

  plotsym, 0, /fill
; oploterror, a[i_indep].rp/1E3, a[i_indep].ew*(1./a[i_indep].line_ratio+1.), a[i_indep].sdev_nofit_lr_1_s5_mc*2.0*bootstrap_corr[i_indep], $
  oploterror, a[i_indep].rp/1E3, a[i_indep].ew*(1./a[i_indep].line_ratio+1.), a_boot[i_indep].err_ew, $
      psym=8, symsize=1.6, color=djs_icolor('blue'), thick=thick, errcolor=djs_icolor('blue'), errthick=thick, hatlength=!d.x_vsize/80.
  djs_oplot, [0.033, 0.033], [2.56, 2.56], psym=8, symsize=1.3, color=djs_icolor('blue'), thick=thick

  ;; Mg I
; plotsym, 8, /fill
; oploterror, a[i_indep[0]].rp/1E3, [0.0831], [0.0645], $
;     psym=8, symsize=1.2, color=djs_icolor('blue'), thick=2, errcolor=djs_icolor('blue'), errthick=2, hatlength=!d.x_vsize/80.
; djs_oplot, [0.25, 0.25], [1.6, 1.6], psym=8, symsiz=1., color=djs_icolor('blue')

; djs_oplot, [0.35, 0.35], [0.9, 0.9], psym=8, symsiz=1., color=djs_icolor('blue')
; djs_oplot, [0.015, 0.015], [4.5E-4, 4.5E-4], psym=8, symsiz=1.2, color=djs_icolor('blue'), symthick=2

; plotsym, 0
; oploterror, a[i_indep].rp/1E3, a[i_indep].ew_nofit_2_s4*1.5, a[i_indep].sdev_nofit_2_s4_mc*1.5, $
;     psym=8, symsize=1.4, color=djs_icolor('dark green'), thick=thick, errcolor=djs_icolor('dark green'), errthick=thick, hatlength=!d.x_vsize/80.
; plotsym, 0
; oploterror, a[i_indep].rp/1E3, a[i_indep].ew_nofit_lr_1_s5*2.0, a[i_indep].sdev_nofit_lr_1_s5_mc*2.0, $
;     psym=8, symsize=1.4, color=djs_icolor('red'), thick=thick, errcolor=djs_icolor('red'), errthick=thick, hatlength=!d.x_vsize/80.

; plotsym, 0
; oploterror, a[i_indep[0:4]].rp/1E3, a[i_indep[0:4]].ew*(1./a[i_indep[0:4]].line_ratio+1.)*3./2., a[i_indep[0:4]].sdev_nofit_2_s4_mc*1.5, $
;     psym=8, symsize=1.3, color=djs_icolor('blue'), thick=thick, errcolor=djs_icolor('blue'), errthick=thick, hatlength=!d.x_vsize/50.
; plotsym, 0
; oploterror, msm7_rp/1E3, msm7_Sigma/mass_factor*fbaryon, msm7_Sigma_err/mass_factor*fbaryon, $
;     psym=8, symsize=2, color=djs_icolor('gray'), thick=thick, errcolor=djs_icolor('gray'), errthick=thick
; djs_xyouts, 1.25, 0.475, $

; djs_xyouts, 0.0125, 3.00E-4, $

  djs_xyouts, 0.12, 5.00E-1*1.5^2, $
             'Si IV around LBGs at z\sim2.2 (Steidel et al. 2010)', color='magenta', charsize=charsize-0.4, charthick=charthick-0.1
; djs_xyouts, 0.99, 1.00E-1*1.5^5, $
;            'Si IV around LBGs at z\sim2.2', color='magenta', charsize=charsize-0.4, charthick=charthick-0.1
; djs_xyouts, 2.85, 1.00E-1*1.5^4, $
;            '(Steidel et al. 2010)', color='magenta', charsize=charsize-0.4, charthick=charthick-0.1

; djs_xyouts, 0.0125, 1.90E-4, $

  djs_xyouts, 0.12, 5.00E-1*1.5, $
             'Mg II around red galaxies at z\sim0.7 (Bordoloi et al. 2011)', color='orange', charsize=charsize-0.4, charthick=charthick-0.1
; djs_xyouts, 0.99, 1.00E-1*1.5^3, $
;            'Mg II around red galaxies at z\sim0.7', color='orange', charsize=charsize-0.4, charthick=charthick-0.1
; djs_xyouts, 2.85, 1.00E-1*1.5^2, $
;            '(Bordoloi et al. 2011)', color='orange', charsize=charsize-0.4, charthick=charthick-0.1

; djs_xyouts, 1.25, 0.850, $
;            'Zhu et al. (2013) [Mg II]', color='blue', charsize=charsize-0.2, charthick=charthick+0.3
; djs_xyouts, 0.32, 0.450, $
;            'around LRGs at z\sim0.5 (Zhu et al. 2013)', color='blue', charsize=charsize-0.2, charthick=charthick+0.3
; djs_xyouts, 1.50, 0.200, $
;            '(Zhu et al. 2013)', color='blue', charsize=charsize-0.2, charthick=charthick+0.3
; djs_xyouts, 5.00, 0.900, $
;            '(Zhu et al. 2013)', color='blue', charsize=charsize-0.3, charthick=charthick+0.3
; djs_xyouts, 0.43, 1.400, $
;            'Mg II', color='blue', charsize=charsize-0.2, charthick=charthick+0.3
 ;djs_xyouts, 0.43, 0.80, $
 ;           'Mg I', color='blue', charsize=charsize-0.2, charthick=charthick+0.3
; djs_xyouts, 0.018, 4.00E-4, $
;            'Zhu et al. (2013) [Mg I]', color='blue', charsize=charsize-0.2, charthick=charthick+0.3
; djs_xyouts, 0.0125, 1.2E-4, $

  djs_xyouts, 0.12, 5.00E-1, $
             ANSI_value('Ca II around all galaxies at z\sim0.1 (Zhu & M!3'+string("351B)+'!Xnard 2013)'), color='light blue', charsize=charsize-0.4, charthick=charthick-0.1
;            ANSI_value('Zhu & M!3'+string("351B)+'!Xnard (2013,')+' Ca II)', color='light blue', charsize=charsize-0.4, charthick=charthick-0.1
; djs_xyouts, 0.99, 1.00E-1*1.5, $
;            'Ca II around all galaxies at z\sim0.1', color='light blue', charsize=charsize-0.4, charthick=charthick-0.1
; djs_xyouts, 2.85, 1.00E-1, $
;            ANSI_value('(Zhu & M!3'+string("351B)+'!Xnard 2013)'), color='light blue', charsize=charsize-0.4, charthick=charthick-0.1
;            ANSI_value('Zhu & M!3'+string("351B)+'!Xnard (2013,')+' Ca II)', color='light blue', charsize=charsize-0.4, charthick=charthick-0.1

; djs_xyouts, 0.04, 1.360, $
  djs_xyouts, 0.04, 2.30, $
             'Mg II around luminous red galaxies at z\sim0.5 (Zhu et al. 2013)', color='blue', charsize=charsize-0.3, charthick=charthick+0.6



; djs_xyouts, 300./1E3, 0.4, $
;            "Mandelbaum & SDSS (Dark Matter)", color='gray', charsize=charsize-0.3, charthick=charthick
; djs_xyouts, 300./1E3, 0.7, $
;            "Zhu, Menard & SDSS ('Dark' Baryons)", color='blue', charsize=charsize-0.3, charthick=charthick

; text_relation = '<W_0>='+string(out_intercept, format='(f5.3)')+'\times(r_p/500 kpc)^{'+string(slope, format='(f5.2)')+'}'
; djs_oplot, [15,25], [1,1]*4.5E-4, thick=thick, linestyle=1, color='gray'
; djs_xyouts, 30, 4.E-4, text_relation, charsize=charsize-0.00, charthick=charthick-0.0, color='black'

;  ytitle1 = '<\Sigma_{Mg II}> (M_\odot pc^{-2})'
   ;; see Eq. 9.15 in Bruce Drain
;  djs_axis, yaxis=1, yra=yra/3.*factor*Mmg24/Msolar*pccm^2, /ylog, ythick=thick, ytitle=ytitle1, $
;      charsiz=charsize, charthick=charthick

k_end_print

k_print, filename=repstr(psfile4, '.ps', '_3.ps'), axis_char_scale=1.3, xsize=10, ysize=6
  djs_plot, a[i_indep].rp/1E3, a[i_indep].ew_nofit_2, xtitle=xtitle, ytitle=ytitle, $
      xra=xra, yra=yra, thick=thick, xthick=xthick, ythick=ythick, $
      charsize=charsize, charthick=charthick, /nodata, /xlog, /ylog, pos=pos, $
      ytickformat='jhusdss_tick_exponent', ystyle=1, xstyle=5

  if (DoScale) then begin
      djs_axis, xaxis=0, xtitle='r_p/r_{vir}', charsize=charsize, charthick=charthick, xra=xra/rvir*1E3, xthick=xthick
      djs_axis, xaxis=1, charsize=charsize, charthick=charthick, xra=xra/rvir*1E3, xthick=xthick, xtickformat='(A1)'
  endif else begin
      xtickv = [0.01, 0.1, 1., 10.]
      xticknames = ['10 kpc', '100 kpc', '1 Mpc', '10 Mpc']
      djs_axis, xaxis=0, xtitle=xtitle, charsize=charsize, charthick=charthick, xra=xra, xthick=xthick, $
          xtickv=xtickv, xtickname=xticknames;, xticks=n_elements(xtickv)
      djs_axis, xaxis=1, xtitle='', charsize=charsize, charthick=charthick, xra=xra, xthick=xthick, xtickformat='(A1)'
  endelse
  djs_xyouts, 0.02, 5.7, 'The Galaxy-Gas (Metal) Correlation', charsize=charsize, charthick=charthick

; plotsym, 3
; oploterror, rp_bor/1E3, ew_bor, ew_bor_err, psym=8, symsize=1., color=djs_icolor('orange'), $
;    errcolor=djs_icolor('orange')
; djs_oplot, [0.011, 0.011], [2.2E-4, 2.2E-4], psym=8, symsiz=1., color=djs_icolor('orange')
; djs_oplot, [0.10, 0.10], [5.3E-1, 5.3E-1]*1.5, psym=8, symsiz=1., color=djs_icolor('orange')

; plotsym, 4
; oploterror, rp_cs/1E3, ew_cs_siiv1393, ew_cs_siiv1393_err, psym=8, symsize=1., color=djs_icolor('magenta'), $
;    errcolor=djs_icolor('magenta')
; djs_oplot, [0.011, 0.011], [3.5E-4, 3.5E-4], psym=8, symsiz=1., color=djs_icolor('magenta')
; djs_oplot, [0.10, 0.10], [5.3E-1, 5.3E-1]*1.5^2, psym=8, symsiz=1., color=djs_icolor('magenta')

; plotsym, 0
; oploterror, a_caii[a_caii_i_fit].rp/1E3, a_caii[a_caii_i_fit].ew_nofit_2, a_caii[a_caii_i_fit].sdev_nofit_2_mc, $
;    psym=8, symsize=1., color=djs_icolor('light blue'), errcolor=djs_icolor('light blue')
; djs_oplot, [0.011, 0.011], [1.35E-4, 1.35E-4], psym=8, symsiz=1., color=djs_icolor('light blue')
; djs_oplot, [0.10, 0.10], [5.3E-1, 5.3E-1], psym=8, symsiz=1., color=djs_icolor('light blue')

  plotsym, 0, /fill
; oploterror, a[i_indep].rp/1E3, a[i_indep].ew*(1./a[i_indep].line_ratio+1.), a[i_indep].sdev_nofit_lr_1_s5_mc*2.0*bootstrap_corr[i_indep], $
  oploterror, a[i_indep].rp/1E3, a[i_indep].ew*(1./a[i_indep].line_ratio+1.), a_boot[i_indep].err_ew, $
      psym=8, symsize=1.6, color=djs_icolor('blue'), thick=thick, errcolor=djs_icolor('blue'), errthick=thick, hatlength=!d.x_vsize/80.
  djs_oplot, [0.033, 0.033], [2.56, 2.56], psym=8, symsize=1.3, color=djs_icolor('blue'), thick=thick

  ;; Mg I
; plotsym, 8, /fill
; oploterror, a[i_indep[0]].rp/1E3, [0.0831], [0.0645], $
;     psym=8, symsize=1.2, color=djs_icolor('blue'), thick=2, errcolor=djs_icolor('blue'), errthick=2, hatlength=!d.x_vsize/80.
; djs_oplot, [0.25, 0.25], [1.6, 1.6], psym=8, symsiz=1., color=djs_icolor('blue')

; djs_oplot, [0.35, 0.35], [0.9, 0.9], psym=8, symsiz=1., color=djs_icolor('blue')
; djs_oplot, [0.015, 0.015], [4.5E-4, 4.5E-4], psym=8, symsiz=1.2, color=djs_icolor('blue'), symthick=2

; plotsym, 0
; oploterror, a[i_indep].rp/1E3, a[i_indep].ew_nofit_2_s4*1.5, a[i_indep].sdev_nofit_2_s4_mc*1.5, $
;     psym=8, symsize=1.4, color=djs_icolor('dark green'), thick=thick, errcolor=djs_icolor('dark green'), errthick=thick, hatlength=!d.x_vsize/80.
; plotsym, 0
; oploterror, a[i_indep].rp/1E3, a[i_indep].ew_nofit_lr_1_s5*2.0, a[i_indep].sdev_nofit_lr_1_s5_mc*2.0, $
;     psym=8, symsize=1.4, color=djs_icolor('red'), thick=thick, errcolor=djs_icolor('red'), errthick=thick, hatlength=!d.x_vsize/80.

; plotsym, 0
; oploterror, a[i_indep[0:4]].rp/1E3, a[i_indep[0:4]].ew*(1./a[i_indep[0:4]].line_ratio+1.)*3./2., a[i_indep[0:4]].sdev_nofit_2_s4_mc*1.5, $
;     psym=8, symsize=1.3, color=djs_icolor('blue'), thick=thick, errcolor=djs_icolor('blue'), errthick=thick, hatlength=!d.x_vsize/50.
; plotsym, 0
; oploterror, msm7_rp/1E3, msm7_Sigma/mass_factor*fbaryon, msm7_Sigma_err/mass_factor*fbaryon, $
;     psym=8, symsize=2, color=djs_icolor('gray'), thick=thick, errcolor=djs_icolor('gray'), errthick=thick
; djs_xyouts, 1.25, 0.475, $

; djs_xyouts, 0.0125, 3.00E-4, $
; djs_xyouts, 0.12, 5.00E-1*1.5^2, $
;            'Si IV around LBGs at z\sim2.2 (Steidel et al. 2010)', color='magenta', charsize=charsize-0.4, charthick=charthick-0.1
; djs_xyouts, 0.0125, 1.90E-4, $
; djs_xyouts, 0.12, 5.00E-1*1.5, $
;            'Mg II around red galaxies at z\sim0.7 (Bordoloi et al. 2011)', color='orange', charsize=charsize-0.4, charthick=charthick-0.1
; djs_xyouts, 1.25, 0.850, $
;            'Zhu et al. (2013) [Mg II]', color='blue', charsize=charsize-0.2, charthick=charthick+0.3
; djs_xyouts, 0.32, 0.450, $
;            'around LRGs at z\sim0.5 (Zhu et al. 2013)', color='blue', charsize=charsize-0.2, charthick=charthick+0.3
; djs_xyouts, 1.50, 0.200, $
;            '(Zhu et al. 2013)', color='blue', charsize=charsize-0.2, charthick=charthick+0.3
; djs_xyouts, 5.00, 0.900, $
;            '(Zhu et al. 2013)', color='blue', charsize=charsize-0.3, charthick=charthick+0.3
; djs_xyouts, 0.43, 1.400, $
;            'Mg II', color='blue', charsize=charsize-0.2, charthick=charthick+0.3
 ;djs_xyouts, 0.43, 0.80, $
 ;           'Mg I', color='blue', charsize=charsize-0.2, charthick=charthick+0.3
; djs_xyouts, 0.018, 4.00E-4, $
;            'Zhu et al. (2013) [Mg I]', color='blue', charsize=charsize-0.2, charthick=charthick+0.3
; djs_xyouts, 0.0125, 1.2E-4, $
; djs_xyouts, 0.12, 5.00E-1, $
;            ANSI_value('Ca II around all galaxies at z\sim0.1 (Zhu & M!3'+string("351B)+'!Xnard 2013)'), color='light blue', charsize=charsize-0.4, charthick=charthick-0.1
;            ANSI_value('Zhu & M!3'+string("351B)+'!Xnard (2013,')+' Ca II)', color='light blue', charsize=charsize-0.4, charthick=charthick-0.1
; djs_xyouts, 0.04, 1.360, $
  djs_xyouts, 0.04, 2.30, $
             'Mg II around luminous red galaxies at z\sim0.5 (Zhu et al. 2013)', color='blue', charsize=charsize-0.3, charthick=charthick+0.6

; djs_xyouts, 300./1E3, 0.4, $
;            "Mandelbaum & SDSS (Dark Matter)", color='gray', charsize=charsize-0.3, charthick=charthick
; djs_xyouts, 300./1E3, 0.7, $
;            "Zhu, Menard & SDSS ('Dark' Baryons)", color='blue', charsize=charsize-0.3, charthick=charthick

; text_relation = '<W_0>='+string(out_intercept, format='(f5.3)')+'\times(r_p/500 kpc)^{'+string(slope, format='(f5.2)')+'}'
; djs_oplot, [15,25], [1,1]*4.5E-4, thick=thick, linestyle=1, color='gray'
; djs_xyouts, 30, 4.E-4, text_relation, charsize=charsize-0.00, charthick=charthick-0.0, color='black'

;  ytitle1 = '<\Sigma_{Mg II}> (M_\odot pc^{-2})'
   ;; see Eq. 9.15 in Bruce Drain
;  djs_axis, yaxis=1, yra=yra/3.*factor*Mmg24/Msolar*pccm^2, /ylog, ythick=thick, ytitle=ytitle1, $
;      charsiz=charsize, charthick=charthick

k_end_print


if keyword_set(thatsgood) then begin
xra=[8,48000]/1E3
yra=[6.E-4, 5.E-0]
pos = [0.20, 0.20, 0.90, 0.90]
plotsym, 0, /fill
k_print, filename=repstr(psfile4, '.ps', '_3.ps'), axis_char_scale=1.3, xsize=10, ysize=6
  djs_plot, a[i_indep].rp/1E3, a[i_indep].ew_nofit_2, xtitle=xtitle, ytitle=ytitle, $
      xra=xra, yra=yra, thick=thick, xthick=xthick, ythick=ythick, $
      charsize=charsize, charthick=charthick, /nodata, /xlog, /ylog, pos=pos, $
      ytickformat='jhusdss_tick_exponent', ystyle=1, xstyle=5

  if (DoScale) then begin
      djs_axis, xaxis=0, xtitle='r_p/r_{vir}', charsize=charsize, charthick=charthick, xra=xra/rvir*1E3, xthick=xthick
      djs_axis, xaxis=1, charsize=charsize, charthick=charthick, xra=xra/rvir*1E3, xthick=xthick, xtickformat='(A1)'
  endif else begin
      xtickv = [0.01, 0.1, 1., 10.]
      xticknames = ['10 kpc', '100 kpc', '1 Mpc', '10 Mpc']
      djs_axis, xaxis=0, xtitle=xtitle, charsize=charsize, charthick=charthick, xra=xra, xthick=xthick, $
          xtickv=xtickv, xtickname=xticknames;, xticks=n_elements(xtickv)
      djs_axis, xaxis=1, xtitle='', charsize=charsize, charthick=charthick, xra=xra, xthick=xthick, xtickformat='(A1)'
  endelse
  djs_xyouts, 0.10, 3.5, 'Galaxy-Gas Correlation', charsize=charsize, charthick=charthick

; ixx = where(xx lt 250., comp=jxx)
; djs_oplot, xx[ixx], 10.^(alog10(xx[ixx]/rvir*2.)*slope+intercept), thick=thick, color='blue'
; djs_oplot, xx[jxx], 10.^(alog10(xx[jxx]/rvir*2.)*slope+intercept), thick=thick, color='blue', linestyle=1

; djs_oplot, xx, 10.^(alog10(xx/rvir*2.)*slope+intercept), thick=thick, linestyle=1, color='gray'
; djs_oplot, xx, 10.^(alog10(xx)*slope_caii+intercept_caii), thick=thick, linestyle=1
; plotsym, 3
; oploterror, rp_bor/1E3, ew_bor, ew_bor_err, psym=8, symsize=1., color=djs_icolor('orange'), $
;    errcolor=djs_icolor('orange')
; djs_oplot, [0.011, 0.011], [2.2E-4, 2.2E-4], psym=8, symsiz=1., color=djs_icolor('orange')
; plotsym, 4
; oploterror, rp_cs/1E3, ew_cs_siiv1393, ew_cs_siiv1393_err, psym=8, symsize=1., color=djs_icolor('magenta'), $
;    errcolor=djs_icolor('magenta')
; djs_oplot, [0.011, 0.011], [3.5E-4, 3.5E-4], psym=8, symsiz=1., color=djs_icolor('magenta')
; plotsym, 0
; oploterror, a_caii[a_caii_i_fit].rp/1E3, a_caii[a_caii_i_fit].ew_nofit_2, a_caii[a_caii_i_fit].sdev_nofit_2_mc, $
;    psym=8, symsize=1., color=djs_icolor('light blue'), errcolor=djs_icolor('light blue')
; djs_oplot, [0.011, 0.011], [1.35E-4, 1.35E-4], psym=8, symsiz=1., color=djs_icolor('light blue')
; plotsym, 4
; oploterror, rp_cs/1E3, ew_cs_civ1549, ew_cs_civ1549_err, psym=8, symsize=1., color=djs_icolor('light blue'), $
;    errcolor=djs_icolor('light blue')
; djs_oplot, [0.015, 0.015], [4.5E-4, 4.5E-4], psym=8, symsiz=1., color=djs_icolor('light blue')
; oploterror, rp_bor, ew_bor_blue, ew_bor_blue_err, psym=8, symsize=1.8, thick=thick-1, color=djs_icolor('gray'), $
;    errcolor=djs_icolor('gray'), errthick=thick
  plotsym, 0, /fill
  oploterror, a[i_indep].rp/1E3, a[i_indep].ew*(1./a[i_indep].line_ratio+1.), a[i_indep].sdev_nofit_lr_1_s5_mc*2.0*bootstrap_corr[i_indep], $
      psym=8, symsize=1.6, color=djs_icolor('blue'), thick=thick, errcolor=djs_icolor('blue'), errthick=thick, hatlength=!d.x_vsize/80.
; djs_oplot, [0.35, 0.35], [1.6, 1.6], psym=8, symsize=1.3, color=djs_icolor('blue'), thick=thick
; djs_oplot, [0.033, 0.033], [1.56, 1.56], psym=8, symsize=1.3, color=djs_icolor('blue'), thick=thick

  ;; Mg I
; plotsym, 8, /fill
; oploterror, a[i_indep[0]].rp/1E3, [0.0831], [0.0645], $
;     psym=8, symsize=1.2, color=djs_icolor('blue'), thick=2, errcolor=djs_icolor('blue'), errthick=2, hatlength=!d.x_vsize/80.
; djs_oplot, [0.25, 0.25], [1.6, 1.6], psym=8, symsiz=1., color=djs_icolor('blue')

; djs_oplot, [0.35, 0.35], [0.9, 0.9], psym=8, symsiz=1., color=djs_icolor('blue')
; djs_oplot, [0.015, 0.015], [4.5E-4, 4.5E-4], psym=8, symsiz=1.2, color=djs_icolor('blue'), symthick=2

; plotsym, 0
; oploterror, a[i_indep].rp/1E3, a[i_indep].ew_nofit_2_s4*1.5, a[i_indep].sdev_nofit_2_s4_mc*1.5, $
;     psym=8, symsize=1.4, color=djs_icolor('dark green'), thick=thick, errcolor=djs_icolor('dark green'), errthick=thick, hatlength=!d.x_vsize/80.
; plotsym, 0
; oploterror, a[i_indep].rp/1E3, a[i_indep].ew_nofit_lr_1_s5*2.0, a[i_indep].sdev_nofit_lr_1_s5_mc*2.0, $
;     psym=8, symsize=1.4, color=djs_icolor('red'), thick=thick, errcolor=djs_icolor('red'), errthick=thick, hatlength=!d.x_vsize/80.

; plotsym, 0
; oploterror, a[i_indep[0:4]].rp/1E3, a[i_indep[0:4]].ew*(1./a[i_indep[0:4]].line_ratio+1.)*3./2., a[i_indep[0:4]].sdev_nofit_2_s4_mc*1.5, $
;     psym=8, symsize=1.3, color=djs_icolor('blue'), thick=thick, errcolor=djs_icolor('blue'), errthick=thick, hatlength=!d.x_vsize/50.
; plotsym, 0
; oploterror, msm7_rp/1E3, msm7_Sigma/mass_factor*fbaryon, msm7_Sigma_err/mass_factor*fbaryon, $
;     psym=8, symsize=2, color=djs_icolor('gray'), thick=thick, errcolor=djs_icolor('gray'), errthick=thick
; djs_xyouts, 1.25, 0.475, $

; djs_xyouts, 0.0125, 3.00E-4, $
;            'Si IV around LBGs at z\sim2.2 (Steidel et al. 2010)', color='magenta', charsize=charsize-0.3, charthick=charthick
; djs_xyouts, 0.0125, 1.90E-4, $
;            'Mg II around red galaxies at z\sim0.7 (Bordoloi et al. 2011)', color='orange', charsize=charsize-0.3, charthick=charthick
; djs_xyouts, 1.25, 0.850, $
;            'Zhu et al. (2013) [Mg II]', color='blue', charsize=charsize-0.2, charthick=charthick+0.3
; djs_xyouts, 0.32, 0.450, $
;            'around LRGs at z\sim0.5 (Zhu et al. 2013)', color='blue', charsize=charsize-0.2, charthick=charthick+0.3
; djs_xyouts, 1.50, 0.200, $
;            '(Zhu et al. 2013)', color='blue', charsize=charsize-0.2, charthick=charthick+0.3
; djs_xyouts, 0.04, 1.360, $
;            'Mg II around luminous red galaxies at z\sim0.5 (Zhu et al. 2013)', color='blue', charsize=charsize-0.3, charthick=charthick+0.3
; djs_xyouts, 5.00, 0.900, $
;            '(Zhu et al. 2013)', color='blue', charsize=charsize-0.3, charthick=charthick+0.3
; djs_xyouts, 0.43, 1.400, $
;            'Mg II', color='blue', charsize=charsize-0.2, charthick=charthick+0.3
 ;djs_xyouts, 0.43, 0.80, $
 ;           'Mg I', color='blue', charsize=charsize-0.2, charthick=charthick+0.3
; djs_xyouts, 0.018, 4.00E-4, $
;            'Zhu et al. (2013) [Mg I]', color='blue', charsize=charsize-0.2, charthick=charthick+0.3
; djs_xyouts, 0.0125, 1.2E-4, $
;            ANSI_value('Ca II around all galaxies at z\sim0.1 (Zhu & M!3'+string("351B)+'!Xnard 2013)'), color='light blue', charsize=charsize-0.3, charthick=charthick
;            ANSI_value('Zhu & M!3'+string("351B)+'!Xnard (2013,')+' Ca II)', color='light blue', charsize=charsize-0.3, charthick=charthick



; djs_xyouts, 300./1E3, 0.4, $
;            "Mandelbaum & SDSS (Dark Matter)", color='gray', charsize=charsize-0.3, charthick=charthick
; djs_xyouts, 300./1E3, 0.7, $
;            "Zhu, Menard & SDSS ('Dark' Baryons)", color='blue', charsize=charsize-0.3, charthick=charthick

; text_relation = '<W_0>='+string(out_intercept, format='(f5.3)')+'\times(r_p/500 kpc)^{'+string(slope, format='(f5.2)')+'}'
; djs_oplot, [15,25], [1,1]*4.5E-4, thick=thick, linestyle=1, color='gray'
; djs_xyouts, 30, 4.E-4, text_relation, charsize=charsize-0.00, charthick=charthick-0.0, color='black'

;  ytitle1 = '<\Sigma_{Mg II}> (M_\odot pc^{-2})'
   ;; see Eq. 9.15 in Bruce Drain
;  djs_axis, yaxis=1, yra=yra/3.*factor*Mmg24/Msolar*pccm^2, /ylog, ythick=thick, ytitle=ytitle1, $
;      charsiz=charsize, charthick=charthick

k_end_print
endif

xra=[20,30000]/1E3
xtitle='r_p' 
;ytitle='W_0^{\lambda 2796} (Mg II) [\AA]'
ytitle='\sigma (km/s)'
title='Double Gaussian Line Profile Measurement'

if (keyword_set(freedom)) then begin
yra=[60., 2.E3]
k_print, filename=psfile5, axis_char_scale=1.3, xsize=8, ysize=6
  djs_plot, a[i_indep].rp/1E3, (a[i_indep].sigma1+a[i_indep].sigma2)/2., xtitle=xtitle, ytitle=ytitle, $
      xra=xra, yra=yra, thick=thick, xthick=xthick, ythick=ythick, $
      charsize=charsize, charthick=charthick, /nodata, /xlog, ylog=1, pos=pos, $
      ystyle=5, xstyle=5

  if (DoScale) then begin
      djs_axis, xaxis=0, xtitle='r_p/r_{vir}', charsize=charsize, charthick=charthick, xra=xra/rvir*1E3, xthick=xthick
      djs_axis, xaxis=1, charsize=charsize, charthick=charthick, xra=xra/rvir*1E3, xthick=xthick, xtickformat='(A1)'
  endif else begin
;     xtickv = [0.01, 0.1, 1., 10.]
;     xticknames = ['10 kpc', '100 kpc', '1 Mpc', '10 Mpc']
      xtickv = [0.1, 1., 10.]
      xticknames = ['100 kpc', '1 Mpc', '10 Mpc']
      djs_axis, xaxis=0, xtitle=xtitle, charsize=charsize, charthick=charthick, xra=xra, xthick=xthick, $
          xtickv=xtickv, xtickname=xticknames;, xticks=n_elements(xtickv)
      djs_axis, xaxis=1, xtitle='', charsize=charsize, charthick=charthick, xra=xra, xthick=xthick, xtickformat='(A1)'
      ytickv = [100., 300., 600.]
      yticknames = ['100', '300', '600']
      print, ytickv, yticknames
      djs_axis, yaxis=0, ytitle=ytitle, charsize=charsize, charthick=charthick, yra=yra, ythick=ythick;, ytickformat='(A1)';$
;         ytickv=ytickv, ytickname=yticknames;, ytickformat='(A)'
      djs_axis, yaxis=1, ytitle='', charsize=charsize, charthick=charthick, yra=yra, ythick=ythick, ytickformat='(A1)'
      djs_xyouts, 0.010, 300*0.95, '300', charsize=charsize+0.3, charthick=charthick
      djs_xyouts, 0.010, 600*0.95, '600', charsize=charsize+0.3, charthick=charthick
  endelse
  djs_xyouts, 0.15, 807, 'Velocity Dispersion', charsize=charsize, charthick=charthick

  sigma_factor = 1./alog(10.)/2800.*1.D+4*69.03
; sigma = 69.03*sqrt((sigma^2-1.)>1.E-6)

  err_sigma1 = a[i_indep].err_sigma1
  err_sigma2 = a[i_indep].err_sigma2
  for ierr=0L, n_elements(err_sigma1)-1L do begin
      if (err_sigma1[ierr] le 0.) then err_sigma1[ierr] = 1.;max(err_sigma1)
      if (err_sigma2[ierr] le 0.) then err_sigma2[ierr] = 1.;max(err_sigma2)
  endfor

  plotsym, 0, /fill
  oploterror, a[i_indep].rp/1E3, sqrt(((a[i_indep].sigma1+a[i_indep].sigma2)/2.*sigma_factor)^2-69.03^2), sqrt(err_sigma1^2+err_sigma2^2)/2.*sigma_factor, $
; oploterror, a[i_indep].rp/1E3, (a[i_indep].sigma1+a[i_indep].sigma2)/2.*sigma_factor, sqrt(err_sigma1^2+err_sigma2^2)/2.*sigma_factor, $
      psym=8, symsize=1.4, color=djs_icolor('blue'), thick=thick, errcolor=djs_icolor('blue'), errthick=thick, hatlength=!d.x_vsize/80.
  djs_oplot, [1., 1.], [0.925, 0.925], psym=8, symsiz=1.2, color=djs_icolor('blue'), thick=thick
k_end_print
endif else begin
yra=[50., 2.E3]
k_print, filename=psfile5, axis_char_scale=1.3, xsize=8, ysize=6
  djs_plot, a[i_indep].rp/1E3, a[i_indep].sigma, xtitle=xtitle, ytitle=ytitle, $
      xra=xra, yra=yra, thick=thick, xthick=xthick, ythick=ythick, $
      charsize=charsize, charthick=charthick, /nodata, /xlog, ylog=1, pos=pos, $
      ystyle=5, xstyle=5

  if (DoScale) then begin
      djs_axis, xaxis=0, xtitle='r_p/r_{vir}', charsize=charsize, charthick=charthick, xra=xra/rvir*1E3, xthick=xthick
      djs_axis, xaxis=1, charsize=charsize, charthick=charthick, xra=xra/rvir*1E3, xthick=xthick, xtickformat='(A1)'
  endif else begin
;     xtickv = [0.01, 0.1, 1., 10.]
;     xticknames = ['10 kpc', '100 kpc', '1 Mpc', '10 Mpc']
      xtickv = [0.1, 1., 10.]
      xticknames = ['100 kpc', '1 Mpc', '10 Mpc']
      djs_axis, xaxis=0, xtitle=xtitle, charsize=charsize, charthick=charthick, xra=xra, xthick=xthick, $
          xtickv=xtickv, xtickname=xticknames;, xticks=n_elements(xtickv)
      djs_axis, xaxis=1, xtitle='', charsize=charsize, charthick=charthick, xra=xra, xthick=xthick, xtickformat='(A1)'
      ytickv = [100., 300., 600.]
      yticknames = ['100', '300', '600']
      print, ytickv, yticknames
      djs_axis, yaxis=0, ytitle=ytitle, charsize=charsize, charthick=charthick, yra=yra, ythick=ythick;, ytickformat='(A1)';$
;         ytickv=ytickv, ytickname=yticknames;, ytickformat='(A)'
      djs_axis, yaxis=1, ytitle='', charsize=charsize, charthick=charthick, yra=yra, ythick=ythick, ytickformat='(A1)'
      djs_xyouts, 0.010, 300*0.95, '300', charsize=charsize+0.3, charthick=charthick
      djs_xyouts, 0.010, 600*0.95, '600', charsize=charsize+0.3, charthick=charthick
  endelse
  djs_xyouts, 0.20, 735, 'Velocity Dispersion', charsize=charsize, charthick=charthick

  sigma_factor = 1./alog(10.)/2800.*1.D+4*69.03
; sigma = 69.03*sqrt((sigma^2-1.)>1.E-6)

  err_sigma = a_boot[i_indep].err_sigma
; for ierr=0L, n_elements(err_sigma)-1L do begin
;     if (err_sigma[ierr] le 0.) then err_sigma[ierr] = 1.;max(err_sigma1)
; endfor

  plotsym, 0, /fill
  oploterror, a[i_indep].rp/1E3, sqrt((a[i_indep].sigma*sigma_factor)^2-69.03^2), $
      err_sigma*sigma_factor*(a[i_indep].sigma*sigma_factor)/sqrt((a[i_indep].sigma*sigma_factor)^2-69.03^2), $
      psym=8, symsize=1.4, color=djs_icolor('blue'), thick=thick, errcolor=djs_icolor('blue'), errthick=thick, hatlength=!d.x_vsize/80.
; oploterror, a[i_indep].rp/1E3, a[i_indep].sigma*sigma_factor, $
;     err_sigma*sigma_factor, $
;     psym=8, symsize=1.4, color=djs_icolor('blue'), thick=thick, errcolor=djs_icolor('blue'), errthick=thick, hatlength=!d.x_vsize/80.
; djs_oplot, [1., 1.], [0.925, 0.925], psym=8, symsiz=1.2, color=djs_icolor('blue'), thick=thick

rr = findgen(100)*4./100.+1.
conc_ref = 10.
rvir_ref = 0.8E3
bbeta = 0.3
sigma_vir_ref = 200.*sqrt((1.-bbeta)/(3.-2.*bbeta))
yy = (alog(1.+conc_ref*10.^rr/rvir_ref)-conc_ref*10.^rr/rvir_ref/(1.+conc_ref*10.^rr/rvir_ref))/10.^rr/(alog(1.+conc_ref)-conc_ref/(1.+conc_ref))*rvir_ref*sigma_vir_ref
ss = 10.^(rr-4.0)
aa = ss/(ss+1.)
bb = 1./(ss+1.)
sigma_2h = 400.
djs_oplot, 10.^rr/1E3, sqrt(yy*yy*bb+sigma_2h*sigma_2h*aa), thick=thick, color='blue'
djs_oplot, 10.^rr/1E3, yy*sqrt(bb), linestyle=2, color='dark green', thick=thick
djs_oplot, 10.^rr/1E3, sigma_2h*sqrt(aa), linestyle=1, color='dark green', thick=thick

k_end_print

endelse

DR_individual = mrdfits('/home/menard/DATA/SDSS/Absorber/MgII/mgii_lineratio.fits',1)
xra=[1E-3,7.0E0]
xtitle='<W_0 (Mg II)> [\AA]' 
title = 'Doublet Ratio'
yra = [0.2, 10]
;ytitle='Doublet Ratio (W_0^{\lambda2796}/W_0^{\lambda2803})'
ytitle='<W_0^{\lambda2796}>/<W_0^{\lambda2803}>'
pos = [0.20, 0.20, 0.90, 0.90]
k_print, filename=psfile6, axis_char_scale=1.3, xsize=8, ysize=6
  djs_plot, a[i_indep].ew*(1.+1./a[i_indep].line_ratio), a[i_indep].line_ratio, xtitle=xtitle, ytitle=ytitle, $
      xra=xra, yra=yra, xst=1, yst=1, thick=thick, xthick=xthick, ythick=ythick, title=title, $
      charsize=charsize, charthick=charthick, /nodata, /xlog, ylog=1, xtickformat='jhusdss_tick_exponent', ytickformat='jhusdss_tick_exponent', pos=pos
  ;djs_oplot, [1E-4, 1.0E-1, 1.0E-1, 2E0], [2,2,1,1], color='orange', thick=12
  ;djs_oplot, [1E-4, 1.0E-1, 1.0E-1, 2E0], [2,2,1,1], color='orange', thick=12
  djs_oplot, xra, [2,2], color='blue', thick=6, linestyle=1
  djs_oplot, xra, [1,1], color='blue', thick=6, linestyle=1
  lflux1 = a.ew
  lflux2 = a.ew/a.line_ratio
  err_line_ratio = sqrt((1./lflux2)^2+(lflux1/lflux2^2)^2)*a.sdev_nofit_lr_1_s5_mc*sqrt(2.)
  lower_limit_line_ratio = (lflux1-a.sdev_nofit_lr_1_s5_mc*sqrt(2))/(lflux2+a.sdev_nofit_lr_1_s5_mc*sqrt(2))
; low_err_line_ratio = a.line_ratio - (lflux1-a.sdev_nofit_lr_1_s5_mc*sqrt(2))/(lflux2+a.sdev_nofit_lr_1_s5_mc*sqrt(2))
; high_err_line_ratio = (lflux1+a.sdev_nofit_lr_1_s5_mc)/lflux2
; low_err_line_ratio = (lflux1-a.sdev_nofit_lr_1_s5_mc)/lflux2
  ;err_line_ratio = sqrt((a.line_ratio/a.ew)^2+(a.line_ratio^2/a.ew)^2)*a.sdev_nofit_lr_1_s5_mc
  oploterror, a[i_indep[0:n_indep-3]].ew*(1.+1./a[i_indep[0:n_indep-3]].line_ratio), a[i_indep[0:n_indep-3]].line_ratio, err_line_ratio[i_indep[0:n_indep-3]], $
      psym=4, symsize=1.2, color=djs_icolor('blue'), thick=thick-1, errcolor=djs_icolor('blue'), errthick=thick-2
  plotsym, 0, 1.2, /fill
  djs_oplot, DR_individual.w_mean, DR_individual.dr_mean, psym=8, color='orange', symsize=1.8
  djs_xyouts, 0.8E-2, 130, 'Doublet Ratio', charsize=charsize, charthick=charthick
  dr_fid = 1.75
  djs_oplot, [xra[0], 0.15], [dr_fid, dr_fid], thick=thick+3, color='dark green'
  xtmp = findgen(50)*0.05+alog10(0.15)
  ytmp = -0.15*(xtmp-alog10(0.15))+alog10(dr_fid)
  xtmp = 10.^xtmp
; ytmp = -0.3*(xtmp-0.15)+1.7
  ytmp = 10.^ytmp
  djs_oplot, xtmp, ytmp, thick=thick+3, color='dark green'
; djs_oplot, dr_individual.w_mean, 10.^(-0.3*(alog10(dr_individual.w_mean)-alog10(0.15))+alog10(1.7)), thick=thick, color='brown'
  colors=[djs_icolor('orange'), djs_icolor('blue'), djs_icolor('dark green')]
  legend, ['Individual Systems', 'Stacking Measurements', 'Adopted Model'], psym=[8,4,0], linestyle=[0,0,0], color=colors, /right, box=0, thick=thick, $
     charsize=charsize, charthick=charthick, textcolors=colors, number=2, pspacing=0.8

; djs_oplot, a[i_indep[n_indep-2:n_indep-1]].ew*(1.+1./a[i_indep[n_indep-2:n_indep-1]].line_ratio), lower_limit_line_ratio[i_indep[n_indep-2:n_indep-1]]/2, $
;     psym=4, symsize=1.2, color=djs_icolor('blue'), thick=thick
  plotsym, 2, 2.25, thick=thick
  djs_oplot, a[i_indep[n_indep-2:n_indep-1]].ew*(1.+1./a[i_indep[n_indep-2:n_indep-1]].line_ratio), lower_limit_line_ratio[i_indep[n_indep-2:n_indep-1]], $
      psym=8, symsize=1.2, color=djs_icolor('blue'), thick=thick
  for itmp=0,1 do begin
      djs_oplot, [0.87,1/0.87]*a[i_indep[n_indep-2+itmp]].ew*(1.+1./a[i_indep[n_indep-2+itmp]].line_ratio), [0.5,0.5]*lower_limit_line_ratio[i_indep[n_indep-2+itmp]], $
          color='blue', thick=thick
      djs_oplot, [1,1]*a[i_indep[n_indep-2+itmp]].ew*(1.+1./a[i_indep[n_indep-2+itmp]].line_ratio), [1,0.5]*lower_limit_line_ratio[i_indep[n_indep-2+itmp]], $
          color='blue', thick=thick
  endfor
k_end_print

xra=[20,20000]/1E3
;xra=[20,30000]/1E3
xtitle='r_p' 
;yra = [0., 5]
yra = [0.2, 10]
;yra = [0.1, 22]
;ytitle='Doublet Ratio (W_0^{\lambda2796}/W_0^{\lambda2803})'
ytitle='<W_0^{\lambda2796}>/<W_0^{\lambda2803}>'
pos = [0.20, 0.20, 0.90, 0.90]
k_print, filename=psfile7, axis_char_scale=1.3, xsize=8, ysize=6
  djs_plot, a[i_indep].ew*(1.+1./a[i_indep].line_ratio), a[i_indep].line_ratio, xtitle=xtitle, ytitle=ytitle, $
      xra=xra, yra=yra, xst=5, yst=1, thick=thick, xthick=xthick, ythick=ythick, $
;     charsize=charsize, charthick=charthick, /nodata, /xlog, ylog=0, xtickformat='jhusdss_tick_exponent', ytickformat='jhusdss_tick_exponent', pos=pos
      charsize=charsize, charthick=charthick, /nodata, /xlog, /ylog, xtickformat='jhusdss_tick_exponent', ytickformat='jhusdss_tick_exponent', pos=pos


; xpoly = [xra[0], xra[1], xra[1], xra[0], xra[0]]
; ypoly = [yra[0], yra[0], 1, 1, yra[0]]
; polyfill, xpoly, ypoly, color=djs_icolor('light gray')
; xpoly = [xra[0], xra[1], xra[1], xra[0], xra[0]]
; ypoly = [yra[1], yra[1], 2, 2, yra[1]]
; polyfill, xpoly, ypoly, color=djs_icolor('light gray')
; djs_axis, yaxis=0, ythick=ythick, ytickformat='(A1)'
; djs_axis, yaxis=1, ythick=ythick, ytickformat='(A1)'

  djs_oplot, xra, [2,2], color='blue', thick=6, linestyle=1
  djs_oplot, xra, [0.98,0.98], color='blue', thick=6, linestyle=1
  lflux1 = a.ew
  lflux2 = a.ew/a.line_ratio
  err_line_ratio = sqrt((1./lflux2)^2+(lflux1/lflux2^2)^2)*a.sdev_nofit_lr_1_s5_mc*sqrt(2.)
  ;err_line_ratio = sqrt((a.line_ratio/a.ew)^2+(a.line_ratio^2/a.ew)^2)*a.sdev_nofit_lr_1_s5_mc
; tt = (lflux1-a.line_ratio*lflux2)^2/(1.+a.line_ratio^2)/a.sdev_nofit_lr_1_s5_mc^2
; qq = 1.-tt*a.sdev_nofit_lr_1_s5_mc^2*2./lflux2^2
; err_line_ratio2 = sqrt(1.+a.line_ratio^2-tt/lflux2^2*a.sdev_nofit_lr_1_s5_mc^2*2.)/lflux2/qq
; high_err_line_ratio = (lflux1+a.sdev_nofit_lr_1_s5_mc*sqrt(2))/lflux2 - a.line_ratio
  low_err_line_ratio = a.line_ratio - (lflux1-a.sdev_nofit_lr_1_s5_mc*sqrt(2))/(lflux2+a.sdev_nofit_lr_1_s5_mc*sqrt(2))
; oploterror, a[i_indep[0:n_indep-3].rp/1E3, a[i_indep[0:n_indep-3].line_ratio, err_line_ratio[i_indep[0:n_indep-3], /hibar, $
  oploterror, a[i_indep[0:n_indep-3]].rp/1E3, a[i_indep[0:n_indep-3]].line_ratio, err_line_ratio[i_indep[0:n_indep-3]], $
      psym=4, symsize=1.2, color=djs_icolor('blue'), thick=thick-1, errcolor=djs_icolor('blue'), errthick=thick-2

; err_line_ratio[i_indep[n_elements(i_indep)-2:n_elements(i_indep)-1]] = low_err_line_ratio[i_indep[n_elements(i_indep)-2:n_elements(i_indep)-1]]
; oploterror, a[i_indep].rp/1E3, a[i_indep].line_ratio, err_line_ratio[i_indep], /lobar, $
;     psym=4, symsize=1.2, color=djs_icolor('blue'), thick=thick-1, errcolor=djs_icolor('blue'), errthick=thick-2

; oploterror, a[i_indep].rp/1E3, a[i_indep].line_ratio, err_line_ratio[i_indep], /hibar, $
;     psym=4, symsize=1.2, color=djs_icolor('blue'), thick=thick-1, errcolor=djs_icolor('blue'), errthick=thick-2
; oploterror, a[i_indep].rp/1E3, a[i_indep].line_ratio, low_err_line_ratio[i_indep], /lobar, $
;     psym=4, symsize=1.2, color=djs_icolor('blue'), thick=thick-1, errcolor=djs_icolor('blue'), errthick=thick-2
  djs_oplot, [xra[0], 1.5E-1, 1.5E-1, xra[1]], [1,1,1.5,1.5], color='orange', thick=12
  djs_xyouts, 3E-1, 24, 'Doublet Ratio', charsize=charsize, charthick=charthick
      xtickv = [0.1, 1., 10.]
      xticknames = ['100 kpc', '1 Mpc', '10 Mpc']
      djs_axis, xaxis=0, xtitle=xtitle, charsize=charsize, charthick=charthick, xra=xra, xthick=xthick, $
          xtickv=xtickv, xtickname=xticknames;, xticks=n_elements(xtickv)
      djs_axis, xaxis=1, xtitle='', charsize=charsize, charthick=charthick, xra=xra, xthick=xthick, xtickformat='(A1)'

  rdisk = 1.5E-1
  djs_oplot, rdisk*[1,1], [5E0,yra[1]], linestyle=2, color='orange', thick=thick+1
; djs_xyouts, rdisk*0.79, 4.0E-3, '<r_{90}>', charsize=charsize-0.3, charthick=charthick, color='gray'
  djs_xyouts, rdisk*1.25, 0.73E1, '<W(Mg II)>', charsize=charsize, charthick=charthick+1, color='orange'
  djs_xyouts, rdisk*1.68, 0.58E1, '<0.15 \AA', charsize=charsize, charthick=charthick+1, color='orange'
  djs_xyouts, rdisk*0.18, 0.73E1, '<W(Mg II)>', charsize=charsize, charthick=charthick+1, color='orange'
  djs_xyouts, rdisk*0.24, 0.58E1, '>0.15 \AA', charsize=charsize, charthick=charthick+1, color='orange'

; djs_plot, alog10(xra), alog10(yra), /nodata, /noerase, xtickformat='(A1)', ytickformat='(A1)', charsize=charsize, charthick=charthick, pos=pos, xst=5, yst=5
; one_arrow, alog10(rdisk*1.2), alog10(8.5E0), 0., ' ', arrowsize=[50,15,40], $
;     charsize=charsize, thick=thick, color=djs_icolor('gray'), /data
; one_arrow, alog10(rdisk/1.2), alog10(8.5E0), 180., ' ', arrowsize=[50,15,40], $
;     charsize=charsize, thick=thick, color=djs_icolor('gray'), /data

k_end_print



end
