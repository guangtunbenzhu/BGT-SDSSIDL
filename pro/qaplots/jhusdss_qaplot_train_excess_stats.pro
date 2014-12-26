pro jhusdss_qaplot_train_excess_stats, nmfver

path=jhusdss_get_path(/nmfqso)+'/'+string(nmfver, format='(I3.3)')
jhu_nomfile = path+'/Absorbers/'+'JHU_nomatched_'+string(nmfver, format='(I3.3)')+'.fits'

jhu_nomatch0 = mrdfits(jhu_nomfile, 1)


;; trim
cutmin = 1550.
nabsmax = 10
sdevcut = 0.07

;; no match in window
ijhu_nomatch = jhusdss_absorber_trim(jhu_nomatch0)
jhu_nomatch = jhu_nomatch0[ijhu_nomatch]

;; only Mg II
ijhu_nomatch_mgii = where(jhu_nomatch.criterion_mgii eq 1b)
jhu_nomatch_mgii = jhu_nomatch[ijhu_nomatch_mgii]

;ijhu_nomatch_mgii = where(jhu_nomatch0.spec_snr_median gt 3. $
;           and jhu_nomatch0.med_sdeviation_red gt 0.00 $
;           and jhu_nomatch0.med_sdeviation_red le sdevcut $
;           and jhu_nomatch0.zqso ge 0.4 $
;           and jhu_nomatch0.nabs le nabsmax $
;          and jhu_nomatch0.zabs gt cutmin*(1.+jhu_nomatch0.zqso)/2803.53+0.02-1 $
;          and jhu_nomatch0.zabs lt jhu_nomatch0.zqso-0.04 $
;          and jhu_nomatch0.zabs gt 4000./2796.35-1. $
;          and jhu_nomatch0.criterion_mgii eq 1b)
;jhu_nomatch_mgii = jhu_nomatch0[ijhu_nomatch_mgii]

absorber = jhu_nomatch_mgii

;; init
thick=6
charsize=1.3
charthick=3
;; path
path=jhusdss_get_path(/nmfqso)+'/'+string(nmfver, format='(I3.3)')+'/Absorbers'
qapath = path+'/QAplots'
if (jhusdss_direxist(qapath) eq 0) then spawn, 'mkdir -p '+qapath

;; doublet ratio
psfile = qapath+'/'+'Excess_doublet_ratio_'+string(nmfver, format='(I3.3)')+'.ps'
k_print, filename=psfile, axis_char_scale=1.3, xsize=10, ysize=8

  xra = [0,6]
  yra = [0.0, 4]
  xtitle = textoidl('W_0^{2796} (\AA)')
  ytitle = textoidl('Doublet Ratio')
  x = absorber.rew_mgii_2796
  y = absorber.rew_mgii_2796/absorber.rew_mgii_2803

  flevels = [0.40, 0.70, 0.80]
  hogg_scatterplot, x, y, xra=xra, yra=yra, $
       xnpix=30L, ynpix=30L, exponent=0.80, $
       thick=thick, xthick=thick, ythick=thick, $
       levels=flevels, ytitle=ytitle, xst=1, yst=1,$
       charsize=charsize, charthick=charthick, xtitle=xtitle, $
       /outliers, outsymsize=0.2
   djs_oplot, !x.crange, [1.0, 1.0], linestyle=2, thick=3, color='blue'
   djs_oplot, !x.crange, [2.0, 2.0], linestyle=2, thick=3, color='blue'
k_end_print

;; vdisp
psfile = qapath+'/'+'Excess_vdisp_'+string(nmfver, format='(I3.3)')+'.ps'
k_print, filename=psfile, axis_char_scale=1.3, xsize=10, ysize=8

  xra = [0,6]
  yra = [0, 400]
  xtitle = textoidl('W_0^{2796} (\AA)')
  ytitle = textoidl('Gaussian Velocity Dispersion (km/s)')
  x = absorber.rew_mgii_2796
  y = absorber.vdisp_mgii_2796

  ii = where(y gt 0.1)
  flevels = [0.40, 0.70, 0.80]
; djs_plot, x, y, xra=xra, yra=yra, psym=3, $
;      thick=thick, xthick=thick, ythick=thick, $
;      xtitle=xtitle, ytitle=ytitle, xst=1, yst=1
  hogg_scatterplot, x[ii], y[ii], xra=xra, yra=yra, $
       xnpix=30L, ynpix=30L, exponent=0.90, $
       thick=thick, xthick=thick, ythick=thick, $
       levels=flevels, ytitle=ytitle, xst=1, yst=1,$
       charsize=charsize, charthick=charthick, xtitle=xtitle, $
       /outliers, outsymsize=0.2
;  djs_oplot, !x.crange, [1.0, 1.0], linestyle=2, thick=3, color='blue'
;  djs_oplot, !x.crange, [2.0, 2.0], linestyle=2, thick=3, color='blue'
k_end_print

stop
end
