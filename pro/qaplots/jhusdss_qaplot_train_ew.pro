pro jhusdss_qaplot_train_ew, nmfver

path=jhusdss_get_path(/nmfqso)+'/'+string(nmfver, format='(I3.3)')
mfile = path+'/Absorbers/'+'JHU_Pitts_matched_'+string(nmfver, format='(I3.3)')+'.fits'
nomfile = path+'/Absorbers/'+'Pitts_nomatched_'+string(nmfver, format='(I3.3)')+'.fits'
jhu_nomfile = path+'/Absorbers/'+'JHU_nomatched_'+string(nmfver, format='(I3.3)')+'.fits'

pitts0 = mrdfits(mfile, 1)
jhu0 = mrdfits(mfile, 2)
pitts_nomatch0 = mrdfits(nomfile, 1)
jhu_nomatch0 = mrdfits(jhu_nomfile, 1)

;; trim
inwindow_all = bytarr(n_elements(pitts0))
pitts0.criterion_mgii = 1b
pitts0.criterion_mgii_feii = 1b
pitts0.criterion_feii = 1b
ipitts_all = jhusdss_absorber_trim(pitts0)
inwindow_all[ipitts_all] = 1b

;; real matches
ireal_matches = bytarr(n_elements(pitts0))
ijhu_all = jhusdss_absorber_trim(jhu0)
pitts_all = pitts0[ijhu_all]
jhu_all = jhu0[ijhu_all]
ireal_matches[ijhu_all] = 1b

;; not real matches
ipitts_all_nomatch = where(inwindow_all and ~ireal_matches)

;; only Mg II
ipitts_mgii = where(jhu_all.criterion_mgii eq 1b)
pitts_mgii = pitts_all[ipitts_mgii]
jhu_mgii = jhu_all[ipitts_mgii]

;; no match in window
pitts_nomatch0.criterion_mgii = 1b
pitts_nomatch0.criterion_mgii_feii = 1b
pitts_nomatch0.criterion_feii = 1b
ipitts_nomatch = jhusdss_absorber_trim(pitts_nomatch0)
pitts_nomatch = [pitts_nomatch0[ipitts_nomatch], pitts0[ipitts_all_nomatch]]

;; no match in window
ijhu_nomatch = jhusdss_absorber_trim(jhu_nomatch0)
jhu_nomatch = jhu_nomatch0[ijhu_nomatch]

;; only Mg II
ijhu_nomatch_mgii = where(jhu_nomatch.criterion_mgii eq 1b)
jhu_nomatch_mgii = jhu_nomatch[ijhu_nomatch_mgii]


pitts = pitts_mgii
jhu = jhu_mgii

;; init
thick=6
charsize=1.3
charthick=3

qapath = path+'/QAplots'
if (jhusdss_direxist(qapath) eq 0) then spawn, 'mkdir -p '+qapath
psfile = qapath+'/'+'w2796_w2803_matched_'+string(nmfver, format='(I3.3)')+'.ps'
k_print, filename=psfile, axis_char_scale=1.4, xsize=12, ysize=8
; !p.multi = [0, 2, 1]
; !y.margin = 0
; !x.margin = 0

  ict = 3
  ;; w2796 vs. w2796
  x = pitts.rew_mgii_2796
  y = jhu.rew_mgii_2796
  xra = [0., 5.49]
  yra = [0., 5.49]
  ytitle = textoidl('W^{\lambda2796}_{0, JHU} (\AA)')
  xtitle = textoidl('W^{\lambda2796}_{0, Pitt} (\AA)')
  pos = [0.10, 0.35, 0.50, 0.90]
  title = "Comparison between JHU and Pittsburgh Catalog"

; flevels = [0.7, 0.85, 0.95]
  flevels = [0.7, 0.85, 0.95]
  loadct, ict 
  hogg_scatterplot, x, y, xra=xra, yra=yra, $
       xnpix=50L, ynpix=50L, exponent=0.35, $
       thick=thick, xthick=thick, ythick=thick, $
       levels=flevels, ytitle=ytitle, xst=1, yst=1,$
       charsize=charsize, charthick=charthick, xtickformat='(A1)', $
       /outliers, outsymsize=0.1, outcolor=djs_icolor('gray'), $
       position=pos, ytickformat='(i1.1)'
  loadct, 0
; djs_axis, xaxis=1, xtitle=xtitle, charsize=charsize, charthick=charthick, xtickformat='(i1.1)'
  djs_oplot, !x.crange, !x.crange, thick=thick

  x = pitts.rew_mgii_2796
  y = jhu.rew_mgii_2796 - pitts.rew_mgii_2796
  ymean = (moment(y, sdev=ysdev))[0]
  xra = [0., 5.49]
  yra = [-1.49, 1.49]
; ytitle = textoidl('W(2796)_{JHU}-W(2796)_{Pitt}')
  ytitle = textoidl('\Delta W^{\lambda2796}_0 (\AA)')
  xtitle = textoidl('W^{\lambda2796}_{0, Pitt} (\AA)')
  pos = [0.10, 0.10, 0.50, 0.35]

  legend = textoidl('<\Delta>='+string(ymean,format='(f6.3)')+' \AA')
  djs_xyouts, !x.crange[0]+0.6*(!x.crange[1]-!x.crange[0]), $
              !y.crange[0]+0.16*(!y.crange[1]-!y.crange[0]), $
              legend, charsize=1.5, charthick=charthick
  legend = textoidl('  \sigma_{\Delta}='+string(ysdev,format='(f4.2)')+' \AA')
  djs_xyouts, !x.crange[0]+0.6*(!x.crange[1]-!x.crange[0]), $
              !y.crange[0]+0.1*(!y.crange[1]-!y.crange[0]), $
              legend, charsize=1.5, charthick=charthick

  djs_xyouts, !x.crange[0]+0.30*(!x.crange[1]-!x.crange[0]), $
              !y.crange[0]+1.02*(!y.crange[1]-!y.crange[0]), $
              title, charsize=1.8, charthick=charthick

  loadct, ict 
  hogg_scatterplot, x, y, xra=xra, yra=yra, $
       xnpix=50L, ynpix=50L, exponent=0.35, $
       thick=thick, xthick=thick, ythick=thick, $
       levels=flevels, ytitle=ytitle, xst=1, yst=1,$
       charsize=charsize, charthick=charthick, xtitle=xtitle, $
       /outliers, outsymsize=0.1, outcolor=djs_icolor('gray'), $
       position=pos, /noerase, xtickformat='(i1.1)'
  loadct, 0
  djs_oplot, !x.crange, [0., 0.], thick=thick

  print, ymean

  ;; w2803 vs. w2803
  x = pitts.rew_mgii_2803
  y = jhu.rew_mgii_2803
  xra = [0., 5.49]
  yra = [0., 5.49]
  ytitle = textoidl('W^{\lambda2803}_{0, JHU} (\AA)')
  xtitle = textoidl('W^{\lambda2803}_{0, Pitt} (\AA)')
  pos = [0.50, 0.35, 0.90, 0.90]

  loadct, ict
  hogg_scatterplot, x, y, xra=xra, yra=yra, $
       xnpix=50L, ynpix=50L, exponent=0.35, $
       thick=thick, xthick=thick, ythick=thick, $
       levels=flevels, ytickformat='(A1)', xst=1, yst=1,$
       charsize=charsize, charthick=charthick, xtickformat='(A1)', $
       /outliers, outsymsize=0.1, outcolor=djs_icolor('gray'), $
       position=pos, /noerase
  loadct, 0
; djs_axis, xaxis=1, xtitle=xtitle, charsize=charsize, charthick=charthick, xtickformat='(i1.1)'
  djs_axis, yaxis=1, ytitle=ytitle, charsize=charsize, charthick=charthick, ytickformat='(i1.1)'
  djs_oplot, !x.crange, !x.crange, thick=thick

  x = pitts.rew_mgii_2803
  y = jhu.rew_mgii_2803 - pitts.rew_mgii_2803
  ymean = (moment(y, sdev=ysdev))[0]
  xra = [0., 5.49]
  yra = [-1.49, 1.49]
  ytitle = textoidl('\Delta W^{\lambda2803}_0 (\AA)')
  xtitle = textoidl('W^{\lambda2803}_{0, Pitt} (\AA)')
  pos = [0.50, 0.10, 0.90, 0.35]

  legend = textoidl('<\Delta>='+string(ymean,format='(f6.3)')+' \AA')
  djs_xyouts, !x.crange[0]+0.6*(!x.crange[1]-!x.crange[0]), $
              !y.crange[0]+0.16*(!y.crange[1]-!y.crange[0]), $
              legend, charsize=1.5, charthick=charthick
  legend = textoidl('  \sigma_{\Delta}='+string(ysdev,format='(f4.2)')+' \AA')
  djs_xyouts, !x.crange[0]+0.6*(!x.crange[1]-!x.crange[0]), $
              !y.crange[0]+0.1*(!y.crange[1]-!y.crange[0]), $
              legend, charsize=1.5, charthick=charthick

  loadct, ict
  hogg_scatterplot, x, y, xra=xra, yra=yra, $
       xnpix=50L, ynpix=50L, exponent=0.35, $
       thick=thick, xthick=thick, ythick=thick, $
       levels=flevels, ytickformat='(A1)', xst=1, yst=1,$
       charsize=charsize, charthick=charthick, xtitle=xtitle, $
       /outliers, outsymsize=0.1, outcolor=djs_icolor('gray'), $
       position=pos, /noerase, xtickformat='(i1.1)'
  loadct, 0
  djs_axis, yaxis=1, ytitle=ytitle, charsize=charsize, charthick=charthick
  djs_oplot, !x.crange, [0., 0.], thick=thick

k_end_print

end
