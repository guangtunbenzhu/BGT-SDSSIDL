pro jhusdss_qaplot_stats, nmfver, boss=boss

qso = jhusdss_qso_readin(boss=boss)

absorber = jhusdss_absorber_readin(nmfver, boss=boss)
absorber_feii = jhusdss_absorber_readin(nmfver, boss=boss, /feii)
absorber_asso = jhusdss_absorber_readin(nmfver, boss=boss, /keepasso)

;; init
thick=8
charsize=1.4
charthick=3
;; path
path=jhusdss_get_path(/nmfqso)+'/'+string(nmfver, format='(I3.3)')
qapath = path+'/QAplots'
if (jhusdss_direxist(qapath) eq 0) then spawn, 'mkdir -p '+qapath

;; redshift distribution
psfile = qapath+'/'+'zabs_vs_zqso_'+string(nmfver, format='(I3.3)')+'.ps'
if keyword_set(boss) then $
psfile = qapath+'/BOSS_'+'zabs_vs_zqso_'+string(nmfver, format='(I3.3)')+'.ps'

k_print, filename=psfile, axis_char_scale=1.5, xsize=10, ysize=8

  zqso = absorber_asso.zqso
  zabs = absorber_asso.zabs
  xra = [0.2, 5.0]
  yra = [0.1, 2.5]
  xtitle = textoidl('z_{QSO}')
  ytitle = textoidl('z_{ABS}')
  djs_plot, zqso, zabs, xra=xra, yra=yra, psym=4, symsize=0.001, $
      thick=thick, xthick=thick, ythick=thick, $
      xtitle=xtitle, ytitle=ytitle, $
      charsize=charsize, charthick=charthick, $
      xstyle=1, ystyle=1, linestyle=0

; djs_oplot, !x.crange, [1., 1.]*(4000./2796.35-1.), thick=thick, linestyle=1
  djs_oplot, !x.crange, !x.crange, thick=thick, linestyle=2
; djs_oplot, !x.crange, 1909.*(1.+!x.crange)/2803.53-1., thick=thick, linestyle=2
  djs_oplot, !x.crange, 1550.*(1.+!x.crange)/2803.53-1., thick=thick, linestyle=2

  djs_xyouts, 1.9, 2.35, 'Mg II', $
      charsize=charsize+0.1, charthick=charthick
; djs_xyouts, 3.5, 2.35, 'C III', $
;     charsize=charsize+0.1, charthick=charthick
  djs_xyouts, 4.45, 1.9, 'C IV', $
      charsize=charsize+0.1, charthick=charthick
; djs_xyouts, 3., 0.45, '4000 \AA cut-off', $
;     charsize=charsize, charthick=charthick
k_end_print

;; redshift histogram
psfile = qapath+'/'+'zhist_'+string(nmfver, format='(I3.3)')+'.ps'
if keyword_set(boss) then $
psfile = qapath+'/BOSS_'+'zhist_'+string(nmfver, format='(I3.3)')+'.ps'

k_print, filename=psfile, axis_char_scale=1.5, xsize=10, ysize=8

  zqso = qso.z
  ii = where(zqso gt 0.) 
  zqso = zqso[ii]
  zabs = absorber_asso.zabs
  xra = [0,5]
  yra = [0.01, 6500]
  xtitle = textoidl('z')
  ytitle = textoidl('\DeltaN (\Deltaz=0.1)')
  plothist, zqso, xra=xra, yra=yra, bin=0.1, $
      thick=thick, xthick=thick, ythick=thick, $
      xtitle=xtitle, ytitle=ytitle, $
      charsize=charsize, charthick=charthick, $
      xstyle=1, ystyle=1, linestyle=0, color=djs_icolor('blue')
  plothist, zabs, bin=0.1, /fill, fcolor=djs_icolor('orange'), /overplot
; plothist, zqso, bin=0.1, thick=thick, /overplot, linestyle=0, color=djs_icolor('blue')

  items = ['Quasars', 'Mg II Absorbers']
  colors = [djs_icolor('blue'), djs_icolor('orange')]
  linestyle = [0, 0]
  legend, items, linestyle=linestyle, textcolors=colors, colors=colors, $
     charsize=1.8, charthick=2.5, thick=thick, /top, /right, box=0, pspacing=1.6
k_end_print

;; doublet ratio
psfile = qapath+'/'+'doublet_ratio_'+string(nmfver, format='(I3.3)')+'.ps'
if keyword_set(boss) then $
psfile = qapath+'/BOSS_'+'doublet_ratio_'+string(nmfver, format='(I3.3)')+'.ps'

k_print, filename=psfile, axis_char_scale=1.5, xsize=9, ysize=8

  ict = 3

  xra = [0,6]
  yra = [0.0, 4]
  xtitle = textoidl('W^{\lambda2796}_0 (\AA)')
  ytitle = textoidl('Doublet Ratio (W^{\lambda2796}_0/W^{\lambda2803}_0)')
  x = absorber.rew_mgii_2796
  y = absorber.rew_mgii_2796/absorber.rew_mgii_2803

  flevels = [0.50, 0.80, 0.95]
  loadct, ict
  hogg_scatterplot, x, y, xra=xra, yra=yra, $
       xnpix=45L, ynpix=45L, exponent=0.50, $
       thick=thick, xthick=thick, ythick=thick, $
       levels=flevels, ytitle=ytitle, xst=1, yst=1,$
       charsize=charsize, charthick=charthick, xtitle=xtitle, $
       /outliers, outsymsize=0.1, outcolor=djs_icolor('gray')
   loadct, 0
   djs_oplot, !x.crange, [1.0, 1.0], linestyle=2, thick=3, color='blue'
   djs_oplot, !x.crange, [2.0, 2.0], linestyle=2, thick=3, color='blue'
k_end_print

;; vdisp
psfile = qapath+'/'+'vdisp_'+string(nmfver, format='(I3.3)')+'.ps'
if keyword_set(boss) then $
psfile = qapath+'/BOSS_'+'vdisp_'+string(nmfver, format='(I3.3)')+'.ps'
k_print, filename=psfile, axis_char_scale=1.5, xsize=9, ysize=8

  xra = [0,6]
  yra = [0, 300]
  xtitle = textoidl('W^{\lambda2796}_0 (\AA)')
  ytitle = textoidl('Velocity Dispersion (\sigma km/s)')
  x = absorber.rew_mgii_2796
  y = absorber.vdisp_mgii_2796

  ii = where(y gt 0.1)
  flevels = [0.50, 0.80, 0.95]
; djs_plot, x, y, xra=xra, yra=yra, psym=3, $
;      thick=thick, xthick=thick, ythick=thick, $
;      xtitle=xtitle, ytitle=ytitle, xst=1, yst=1
  loadct, ict 
  hogg_scatterplot, x[ii], y[ii], xra=xra, yra=yra, $
       xnpix=45L, ynpix=45L, exponent=0.50, $
       thick=thick, xthick=thick, ythick=thick, $
       levels=flevels, ytitle=ytitle, xst=1, yst=1,$
       charsize=charsize, charthick=charthick, xtitle=xtitle, $
       /outliers, outsymsize=0.1, outcolor=djs_icolor('gray')
;  djs_oplot, !x.crange, [1.0, 1.0], linestyle=2, thick=3, color='blue'
;  djs_oplot, !x.crange, [2.0, 2.0], linestyle=2, thick=3, color='blue'
  loadct, 0
k_end_print

;; w(2796) hist
psfile = qapath+'/'+'w2796_hist_'+string(nmfver, format='(I3.3)')+'.ps'
if keyword_set(boss) then $
psfile = qapath+'/BOSS_'+'w2796_hist_'+string(nmfver, format='(I3.3)')+'.ps'

k_print, filename=psfile, axis_char_scale=1.5, xsize=10, ysize=8

  xra = [0.2, 7]
  yra = [0.8, 5000]
; xtitle = textoidl('log_{10}(W_0^{2796}) (\AA)')
  xtitle = textoidl('W^{\lambda2796}_0 (\AA)')
  ytitle = textoidl('Number')
  x = absorber.rew_mgii_2796

  plothist, x, xra=xra, yra=yra, bin=0.10, $
      thick=thick, xthick=thick, ythick=thick, $
      xtitle=xtitle, ytitle=ytitle, $
      charsize=charsize, charthick=charthick, $
      xstyle=1, ystyle=1, /ylog

k_end_print


stop
end
