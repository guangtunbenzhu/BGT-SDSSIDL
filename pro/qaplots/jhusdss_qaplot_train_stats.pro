pro jhusdss_qaplot_train_stats, nmfver, nowlog=nowlog

path=jhusdss_get_path(/nmfqso)+'/'+string(nmfver, format='(I3.3)')

mfile = path+'/Absorbers/'+'JHU_Pitts_matched_'+string(nmfver, format='(I3.3)')+'.fits'
nomfile = path+'/Absorbers/'+'Pitts_nomatched_'+string(nmfver, format='(I3.3)')+'.fits'
jhu_nomfile = path+'/Absorbers/'+'JHU_nomatched_'+string(nmfver, format='(I3.3)')+'.fits'

pitts0 = mrdfits(mfile, 1)
jhu0 = mrdfits(mfile, 2)
pitts_nomatch0 = mrdfits(nomfile, 1)
jhu_nomatch0 = mrdfits(jhu_nomfile, 1)

;;##########
;; trim
;; This is not accurate -- Guangtun 04-15-2013
;; Let's use ijhu_all instead
inwindow_all = bytarr(n_elements(pitts0))
pitts0.criterion_mgii = 1b
pitts0.criterion_mgii_feii = 1b
pitts0.criterion_feii = 1b
pitts0.snr_mgii_2796 = jhu0.snr_mgii_2796
ipitts_all = jhusdss_absorber_trim(pitts0)
inwindow_all[ipitts_all] = 1b

;; real matches
ireal_matches = bytarr(n_elements(pitts0))
ijhu_all = jhusdss_absorber_trim(jhu0)
pitts_all = pitts0[ijhu_all]
jhu_all = jhu0[ijhu_all]
ireal_matches[ijhu_all] = 1b

;; in window but not real matches
;; Not accurate -- Guangtun 04-15-2013
;; ipitts_all_nomatch = where(inwindow_all and ~ireal_matches)

;; only Mg II
;; pitts_all includes Fe II, pitts_mgii only Mg II
;; Add comp=ipitts_all_nomatch -- Guangtun 04-15-2013
ipitts_mgii = where(jhu_all.criterion_mgii eq 1b, comp=ipitts_all_nomatch)
pitts_mgii = pitts_all[ipitts_mgii]
jhu_mgii = jhu_all[ipitts_mgii]

;; no match in window
pitts_nomatch0.criterion_mgii = 1b
pitts_nomatch0.criterion_mgii_feii = 1b
pitts_nomatch0.criterion_feii = 1b
pitts_nomatch0.snr_mgii_2796 = 0. ;; Low-z shouldn't be a reference, Guangtun 04-15-2013
ipitts_nomatch = jhusdss_absorber_trim(pitts_nomatch0)
;; Change to pitts_all[ipitts_all_nomatch] -- Guangtun 04-15-2013
;; pitts_nomatch = [pitts_nomatch0[ipitts_nomatch], pitts0[ipitts_all_nomatch]]
pitts_nomatch = [pitts_nomatch0[ipitts_nomatch], pitts_all[ipitts_all_nomatch]]

;; no match in window
ijhu_nomatch = jhusdss_absorber_trim(jhu_nomatch0)
jhu_nomatch = jhu_nomatch0[ijhu_nomatch]

;; only Mg II
ijhu_nomatch_mgii = where(jhu_nomatch.criterion_mgii eq 1b)
jhu_nomatch_mgii = jhu_nomatch[ijhu_nomatch_mgii]
;;##########

itmp = where(pitts_all.rew_mgii_2796 gt 4, ntmp)
print, ntmp
itmp = where(pitts_nomatch.rew_mgii_2796 gt 4, ntmp)
print, ntmp
itmp = where(jhu_nomatch_mgii.rew_mgii_2796 gt 4, ntmp)
print, ntmp

;; init
thick=8
charsize=1.3
charthick=2

qapath = path+'/QAplots'
if (jhusdss_direxist(qapath) eq 0) then spawn, 'mkdir -p '+qapath
psfile = qapath+'/'+'Frecover_Fexcess_'+string(nmfver, format='(I3.3)')+'.ps'

k_print, filename=psfile, axis_char_scale=1.4, xsize=12, ysize=8

  if (~keyword_set(nowlog)) then begin
      binmin=alog10(0.3)
      binmax=alog10(10.)
      binsize=alog10(1.5)/2.

      xra = alog10([0.25, 8.])

      ;; include matches with criterion_mgii_feii
      pitts_hist_mgfe = histogram(alog10(pitts_all.rew_mgii_2796), binsize=binsize, min=binmin, max=binmax)
      jhu_hist_mgfe = histogram(alog10(jhu_all.rew_mgii_2796), binsize=binsize, min=binmin, max=binmax)
      ;; only matches with criterion_mgii
      pitts_hist_mg = histogram(alog10(pitts_mgii.rew_mgii_2796), binsize=binsize, min=binmin, max=binmax)
      jhu_hist_mg = histogram(alog10(jhu_mgii.rew_mgii_2796), binsize=binsize, min=binmin, max=binmax)

      ;; nomatch if matches include criterion_mgii_feii
      pitts_nomatch_hist_mgfe = histogram(alog10(pitts_nomatch.rew_mgii_2796), binsize=binsize, min=binmin, max=binmax)
      ;; nomatch for criterion_mgii and criterion_mgii_feii
      jhu_nomatch_hist_mgfe = histogram(alog10(jhu_nomatch.rew_mgii_2796), binsize=binsize, min=binmin, max=binmax)
      ;; nomatch for criterion_mgii
      jhu_nomatch_hist_mg = histogram(alog10(jhu_nomatch_mgii.rew_mgii_2796), binsize=binsize, min=binmin, max=binmax)
  endif else begin
      binmin=0.0
      binmax=10.
      binsize=0.5

      xra = [0.00, 6.0]

      ;; include matches with criterion_mgii_feii
      pitts_hist_mgfe = histogram(pitts_all.rew_mgii_2796, binsize=binsize, min=binmin, max=binmax)
      jhu_hist_mgfe = histogram(jhu_all.rew_mgii_2796, binsize=binsize, min=binmin, max=binmax)
      ;; only matches with criterion_mgii
      pitts_hist_mg = histogram(pitts_mgii.rew_mgii_2796, binsize=binsize, min=binmin, max=binmax)
      jhu_hist_mg = histogram(jhu_mgii.rew_mgii_2796, binsize=binsize, min=binmin, max=binmax)

      ;; nomatch if matches include criterion_mgii_feii
      pitts_nomatch_hist_mgfe = histogram(pitts_nomatch.rew_mgii_2796, binsize=binsize, min=binmin, max=binmax)
      ;; nomatch for criterion_mgii and criterion_mgii_feii
      jhu_nomatch_hist_mgfe = histogram(jhu_nomatch.rew_mgii_2796, binsize=binsize, min=binmin, max=binmax)
      ;; nomatch for criterion_mgii
      jhu_nomatch_hist_mg = histogram(jhu_nomatch_mgii.rew_mgii_2796, binsize=binsize, min=binmin, max=binmax)
  endelse


  ;; all pitts catalog
  pitts_all_hist = pitts_hist_mgfe+pitts_nomatch_hist_mgfe
  ;; all jhu catalog for criterion_mgii and criterion_mgii_feii
  jhu_all_hist = jhu_hist_mgfe+jhu_nomatch_hist_mgfe
  ;; all jhu catalog for criterion_mgii only
  jhu_all_hist_mg = jhu_hist_mg+jhu_nomatch_hist_mg

  bins = findgen(n_elements(pitts_hist_mg))*binsize+binmin+binsize/2.

  ;; w2796 vs. w2796
  ii = where(pitts_all_hist gt 0)
  x = [binmin-binsize/2., bins[ii]]
  y = [0., float(pitts_hist_mg[ii])/float(pitts_all_hist[ii])]
  y1 = [0., float(pitts_hist_mgfe[ii])/float(pitts_all_hist[ii])]
  y_err = [0., sqrt(float(pitts_hist_mg[ii])-float(pitts_hist_mg[ii]*pitts_hist_mg[ii])/float(pitts_all_hist[ii]))/float(pitts_all_hist[ii])]


  yra = [0., 1.05]
  ytitle = textoidl('N(JHU in Pitt)/N(Pitt)')
  if (~keyword_set(nowlog)) then xtitle = textoidl('log_{10} (W^{\lambda2796}_{0, Pitt}/\AA)') else xtitle=textoidl('W^{\lambda2796}_{0, Pitt} (\AA)')
  xtitle1 = textoidl('W^{\lambda2796}_{0, Pitt} (\AA)')
  pos = [0.10, 0.10, 0.48, 0.90]

  djs_plot, x, y, psym=10, xra=xra, yra=yra, $
      ytitle=ytitle, $
      thick=thick, xthick=thick, ythick=thick, $
      charsize=charsize, charthick=charthick, $
      xstyle=5, ystyle=1, position=pos, xtickformat='(A1)'
  djs_oplot, x, y1, linestyle=2, thick=5, color='red', psym=10
  oploterror, x, y, y_err, psym=3, thick=thick, color=djs_icolor('black')
  djs_oplot, !x.crange, [1., 1.], linestyle=1, thick=5, color='gray'
  djs_oplot, !x.crange, [0.9, 0.9], linestyle=1, thick=5, color='gray'

  if (~keyword_set(nowlog)) then begin
     xtickv = [0.3, 0.4, 0.5, 0.6,0.7, 0.8, 0.9, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0]
     xtickname = [' ', ' ', '0.5', ' ',' ', ' ', ' ', '1', '2', '3', '4', '5', '6']
     axis, xaxis=0, xtitle=xtitle1, xra=10.^(!x.crange), /xlog, $
         xthick=thick, charsize=charsize, charthick=charthick, $
         xticks=n_elements(xtickv), xtickv=xtickv, xtickname=xtickname, xticklen=!p.ticklen/2.
     axis, xaxis=1, xtitle='', xra=10.^(!x.crange), /xlog, $
         xthick=thick, charsize=charsize, charthick=charthick, xtickformat='(A1)'
   endif else begin
     axis, xaxis=1, xtitle=xtitle, $
         xthick=thick, charsize=charsize, charthick=charthick
   endelse

  items = ['Only Mg II detections', 'Include Fe II detections']
  colors = [djs_icolor('black'), djs_icolor('red')]
  legend, items, colors=colors, textcolors=colors, box=0, /right, /bottom, $
     linestyle=[0, 2], thick=thick, charsize=1.3, charthick=1.5, pspacing=1.6

  djs_xyouts, -0.5, 0.30, 'Completeness of the JHU catalog', $
      charsize=1.3, charthick=1.6
  djs_xyouts, -0.5, 0.25, 'compared to the Pittsburgh catalog', $
      charsize=1.3, charthick=1.6

  ii = where(jhu_all_hist_mg gt 0)
  x = [binmin-binsize/2., bins[ii]]
  y = [0., float(jhu_hist_mg[ii])/float(jhu_all_hist_mg[ii])]
  y1 = [0., float(jhu_hist_mgfe[ii])/float(jhu_all_hist[ii])]
  y_err = [0., sqrt(float(jhu_hist_mg[ii])-float(jhu_hist_mg[ii]*jhu_hist_mg[ii])/float(jhu_all_hist_mg[ii]))/float(jhu_all_hist_mg[ii])]
  ytitle = textoidl('N(Pitt in JHU)/N(JHU)')
  if (~keyword_set(nowlog)) then xtitle = textoidl('log_{10} (W^{\lambda2796}_{0, JHU}/\AA)') else xtitle=textoidl('W^{\lambda2796}_{0, JHU} (\AA)')
  xtitle1 = textoidl('W^{\lambda2796}_{0, JHU} (\AA)')
  pos = [0.58, 0.10, 0.96, 0.90]

  djs_plot, x, y, psym=10, xra=xra, yra=yra, $
      ytitle=ytitle, $
      thick=thick, xthick=thick, ythick=thick, $
      charsize=charsize, charthick=charthick, $
      xstyle=5, ystyle=1, position=pos, /noerase
  djs_oplot, x, y1, linestyle=2, thick=5, color='red', psym=10
  oploterror, x, y, y_err, psym=3, thick=thick, color=djs_icolor('black')
  djs_oplot, !x.crange, [1., 1.], linestyle=1, thick=5, color='gray'
  djs_oplot, !x.crange, [0.9, 0.9], linestyle=1, thick=5, color='gray'

  djs_xyouts, -0.5, 0.30, 'Completeness of the Pittsburgh catalog', $
      charsize=1.3, charthick=1.6
  djs_xyouts, -0.5, 0.25, 'compared to the JHU catalog', $
      charsize=1.3, charthick=1.6

  if (~keyword_set(nowlog)) then begin
     axis, xaxis=0, xtitle=xtitle1, xra=10.^(!x.crange), /xlog, $
         xthick=thick, charsize=charsize, charthick=charthick, $
         xticks=n_elements(xtickv), xtickv=xtickv, xtickname=xtickname, xticklen=!p.ticklen/2.
     axis, xaxis=1, xtitle='', xra=10.^(!x.crange), /xlog, $
         xthick=thick, charsize=charsize, charthick=charthick, xtickformat='(A1)'
   endif else begin
     axis, xaxis=1, xtitle=xtitle, $
         xthick=thick, charsize=charsize, charthick=charthick
   endelse

k_end_print

stop
end
