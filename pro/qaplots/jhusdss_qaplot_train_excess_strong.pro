pro jhusdss_qaplot_train_excess_strong, nmfver

lines = jhusdss_train_lines()
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

jhu_nomatch = jhu_nomatch_mgii
ijhu_high = where(jhu_nomatch.rew_mgii_2796 le 8.0 and jhu_nomatch.rew_mgii_2796 ge 4.0, njhu_high)
print, njhu_high
ijhu_high = reverse(ijhu_high[bsort(jhu_nomatch[ijhu_high].rew_mgii_2796)])

ipitts_high = where(pitts_nomatch.rew_mgii_2796 le 8.0 and pitts_nomatch.rew_mgii_2796 ge 3.0, npitts_high)
print, npitts_high
ipitts_high = reverse(ipitts_high[bsort(pitts_nomatch[ipitts_high].rew_mgii_2796)])

;; init
thick=6
charsize=1.1
charthick=2.5

qapath = path+'/QAplots'
if (jhusdss_direxist(qapath) eq 0) then spawn, 'mkdir -p '+qapath
psfile = qapath+'/'+'JHU_nomatched_strong_'+string(nmfver, format='(I3.3)')+'.ps'
k_print, filename=psfile, axis_char_scale=1.2, xsize=12, ysize=8
   for i=0L, njhu_high-1L do begin
       spec = jhusdss_decompose_loadspec(jhu_nomatch[ijhu_high[i]].plate, jhu_nomatch[ijhu_high[i]].fiber, nmfver)
       zabs = jhu_nomatch[ijhu_high[i]].zabs
       nabs = jhu_nomatch[ijhu_high[i]].nabs
       ew2796 = jhu_nomatch[ijhu_high[i]].rew_mgii_2796
       sdev = jhu_nomatch[ijhu_high[i]].med_sdeviation_red
       x = spec.wave*(1.+spec.z)
       yflux = smooth(spec.flux, 5)
       yresi = smooth(spec.residual, 5)
       ynmf = smooth(spec.nmf_continuum*spec.med_continuum, 5)
       ymed = smooth(spec.med_continuum, 5)

       xtitle = textoidl('Observer-Frame \lambda (\AA)')
       ytitle = textoidl('Flux (arbitrary unit)')
       xra = [3700, 9500]
       yra = [0, 6]
       pos = [0.10, 0.35, 0.50, 0.90]
       djs_plot, x, yflux, xra=xra, xstyle=1, $
           yra=yra, ystyle=1, position=pos, xtickformat='(A1)', ytitle=ytitle, $
           thick=thick-3, xthick=thick, ythick=thick, $
           charsize=charsize, charthick=charthick
       djs_oplot, x, ynmf, color='red', thick=thick

       for iline=0L, n_elements(lines)-1L do begin
           djs_oplot, replicate(lines[iline].wave*(1.+zabs),2), $
              !y.crange[0]+[0.6, 1.0]*(!y.crange[1]-!y.crange[0]), $
              color='dark green', thick=thick, linestyle=2
       endfor

       djs_xyouts, 6800, !y.crange[0]+0.51*(!y.crange[1]-!y.crange[0]), 'Z(QSO)='+$
           string(spec.z, format='(f5.3)'), color='red', charsize=charsize, charthick=charthick
       djs_xyouts, 6800, !y.crange[0]+0.46*(!y.crange[1]-!y.crange[0]), 'CIV(QSO)='+$
           string(1550.*(1.+spec.z), format='(f6.1)'), color='red', charsize=charsize, charthick=charthick
       djs_xyouts, 6800, !y.crange[0]+0.41*(!y.crange[1]-!y.crange[0]), 'Z(ABS)='+$
           string(zabs, format='(f5.3)'), color='blue', charsize=charsize, charthick=charthick
       djs_xyouts, 6800, !y.crange[0]+0.36*(!y.crange[1]-!y.crange[0]), 'MgII(ABS)='+$
           string(2800.*(1.+zabs), format='(f6.1)'), color='blue', charsize=charsize, charthick=charthick
        djs_xyouts, 6800, !y.crange[0]+0.30*(!y.crange[1]-!y.crange[0]), 'Plate='+$
            string(spec.plate, format='(i4.4)')+', Fiber='+string(spec.fiber, format='(i4.4)'), $
            color='dark green', charsize=1.0, charthick=2.0
        djs_xyouts, 6800, !y.crange[0]+0.24*(!y.crange[1]-!y.crange[0]), 'N(ABS)='+$
            string(nabs, format='(i1.1)'), color='dark green', charsize=1.0, charthick=2.0

;      djs_xyouts, 7800, !y.crange[0]+0.48*(!y.crange[1]-!y.crange[0]), 'Ly\alpha (1216)='+$
;          string(1216.*(1.+spec.z), format='(f6.1)'), color='red', charsize=charsize, charthick=charthick
;      djs_xyouts, 7800, !y.crange[0]+0.36*(!y.crange[1]-!y.crange[0]), 'MgII (2800)='+$
;          string(2800.*(1.+spec.z), format='(f7.1)'), color='red', charsize=charsize, charthick=charthick

       ytitle = textoidl('Residual')
       pos = [0.10, 0.10, 0.50, 0.35]
       yra = [-0.1, 1.8]
       djs_plot, x, yresi, xra=xra, xstyle=1, $
           yra=yra, ystyle=1, position=pos, xtitle=xtitle, ytitle=ytitle, $
           thick=thick-3, xthick=thick, ythick=thick, $
           charsize=charsize, charthick=charthick, /noerase
       djs_oplot, !x.crange, [1,1], color='blue', thick=thick-3

       xtitle = textoidl('Observer-Frame \lambda (\AA)')
       ytitle = textoidl('Flux (arbitrary unit)')
       pos = [0.50, 0.35, 0.90, 0.90]
       xra = [2000, 3000]*(1.+zabs)
       yra = [0., 6.]
       djs_plot, x, yflux, xra=xra, xstyle=1, $
           yra=yra, ystyle=1, position=pos, xtickformat='(A1)', ytickformat='(A1)', $
           thick=thick-3, xthick=thick, ythick=thick, $
           charsize=charsize, charthick=charthick, /noerase, ytickinterval=400.
       djs_oplot, x, ynmf, color='red', thick=thick

        djs_xyouts, 2050.*(1.+zabs), !y.crange[0]+0.90*(!y.crange[1]-!y.crange[0]), 'W(2796)='+$
            string(ew2796, format='(f4.1)'), color='blue', charsize=1.2, charthick=2.0
        djs_xyouts, 2050.*(1.+zabs), !y.crange[0]+0.85*(!y.crange[1]-!y.crange[0]), 'SDEV='+$
            string(sdev, format='(f6.4)'), color='blue', charsize=1.2, charthick=2.0


      for iline=0L, n_elements(lines)-1L do begin
           djs_oplot, replicate(lines[iline].wave*(1.+zabs),2), $
              !y.crange[0]+[0.0, 1.0]*(!y.crange[1]-!y.crange[0]), $
              color='dark green', thick=thick, linestyle=2
       endfor

       ytitle = textoidl('Residual')
       pos = [0.50, 0.10, 0.90, 0.35]
       yra = [-0.1, 1.8]
       djs_plot, x, yresi, xra=xra, xstyle=1, $
           yra=yra, ystyle=1, position=pos, xtitle=xtitle, ytickformat='(A1)', $
           thick=thick-3, xthick=thick, ythick=thick, $
           charsize=charsize, charthick=charthick, /noerase, ytickinterval=400.
       djs_oplot, !x.crange, [1,1], color='blue', thick=thick-3

       for iline=0L, n_elements(lines)-1L do begin
           djs_oplot, replicate(lines[iline].wave*(1.+zabs),2), $
              !y.crange[0]+[0.0, 1.0]*(!y.crange[1]-!y.crange[0]), $
              color='dark green', thick=thick, linestyle=2
       endfor

   endfor

k_end_print

psfile = qapath+'/'+'Pitts_nomatched_strong_'+string(nmfver, format='(I3.3)')+'.ps'
k_print, filename=psfile, axis_char_scale=1.2, xsize=12, ysize=8
   for i=0L, npitts_high-1L do begin
       spec = jhusdss_decompose_loadspec(pitts_nomatch[ipitts_high[i]].plate, pitts_nomatch[ipitts_high[i]].fiber, nmfver)
       zabs = pitts_nomatch[ipitts_high[i]].zabs
       nabs = pitts_nomatch[ipitts_high[i]].nabs
       ew2796 = pitts_nomatch[ipitts_high[i]].rew_mgii_2796
       sdev = pitts_nomatch[ipitts_high[i]].med_sdeviation_red
       x = spec.wave*(1.+spec.z)
       yflux = smooth(spec.flux, 5)
       yresi = smooth(spec.residual, 5)
       ynmf = smooth(spec.nmf_continuum*spec.med_continuum, 5)
       ymed = smooth(spec.med_continuum, 5)

       xtitle = textoidl('Observer-Frame \lambda (\AA)')
       ytitle = textoidl('Flux (arbitrary unit)')
       xra = [3700, 9500]
       yra = [0, 6]
       pos = [0.10, 0.35, 0.50, 0.90]
       djs_plot, x, yflux, xra=xra, xstyle=1, $
           yra=yra, ystyle=1, position=pos, xtickformat='(A1)', ytitle=ytitle, $
           thick=thick-3, xthick=thick, ythick=thick, $
           charsize=charsize, charthick=charthick
       djs_oplot, x, ynmf, color='red', thick=thick

       for iline=0L, n_elements(lines)-1L do begin
           djs_oplot, replicate(lines[iline].wave*(1.+zabs),2), $
              !y.crange[0]+[0.6, 1.0]*(!y.crange[1]-!y.crange[0]), $
              color='dark green', thick=thick, linestyle=2
       endfor

       djs_xyouts, 6800, !y.crange[0]+0.51*(!y.crange[1]-!y.crange[0]), 'Z(QSO)='+$
           string(spec.z, format='(f5.3)'), color='red', charsize=charsize, charthick=charthick
       djs_xyouts, 6800, !y.crange[0]+0.46*(!y.crange[1]-!y.crange[0]), 'CIV(QSO)='+$
           string(1550.*(1.+spec.z), format='(f6.1)'), color='red', charsize=charsize, charthick=charthick
       djs_xyouts, 6800, !y.crange[0]+0.41*(!y.crange[1]-!y.crange[0]), 'Z(ABS)='+$
           string(zabs, format='(f5.3)'), color='blue', charsize=charsize, charthick=charthick
       djs_xyouts, 6800, !y.crange[0]+0.36*(!y.crange[1]-!y.crange[0]), 'MgII(ABS)='+$
           string(2800.*(1.+zabs), format='(f6.1)'), color='blue', charsize=charsize, charthick=charthick
        djs_xyouts, 6800, !y.crange[0]+0.30*(!y.crange[1]-!y.crange[0]), 'Plate='+$
            string(spec.plate, format='(i4.4)')+', Fiber='+string(spec.fiber, format='(i4.4)'), $
            color='dark green', charsize=1.0, charthick=2.0
        djs_xyouts, 6800, !y.crange[0]+0.24*(!y.crange[1]-!y.crange[0]), 'N(ABS)='+$
            string(nabs, format='(i1.1)'), color='dark green', charsize=1.0, charthick=2.0

;      djs_xyouts, 7800, !y.crange[0]+0.48*(!y.crange[1]-!y.crange[0]), 'Ly\alpha (1216)='+$
;          string(1216.*(1.+spec.z), format='(f6.1)'), color='red', charsize=charsize, charthick=charthick
;      djs_xyouts, 7800, !y.crange[0]+0.36*(!y.crange[1]-!y.crange[0]), 'MgII (2800)='+$
;          string(2800.*(1.+spec.z), format='(f7.1)'), color='red', charsize=charsize, charthick=charthick

       ytitle = textoidl('Residual')
       pos = [0.10, 0.10, 0.50, 0.35]
       yra = [-0.1, 1.8]
       djs_plot, x, yresi, xra=xra, xstyle=1, $
           yra=yra, ystyle=1, position=pos, xtitle=xtitle, ytitle=ytitle, $
           thick=thick-3, xthick=thick, ythick=thick, $
           charsize=charsize, charthick=charthick, /noerase
       djs_oplot, !x.crange, [1,1], color='blue', thick=thick-3

       xtitle = textoidl('Observer-Frame \lambda (\AA)')
       ytitle = textoidl('Flux (arbitrary unit)')
       pos = [0.50, 0.35, 0.90, 0.90]
       xra = [2000, 3000]*(1.+zabs)
       yra = [0., 6.]
       djs_plot, x, yflux, xra=xra, xstyle=1, $
           yra=yra, ystyle=1, position=pos, xtickformat='(A1)', ytickformat='(A1)', $
           thick=thick-3, xthick=thick, ythick=thick, $
           charsize=charsize, charthick=charthick, /noerase, ytickinterval=400.
       djs_oplot, x, ynmf, color='red', thick=thick

        djs_xyouts, 2050.*(1.+zabs), !y.crange[0]+0.90*(!y.crange[1]-!y.crange[0]), 'W(2796)='+$
            string(ew2796, format='(f4.1)'), color='blue', charsize=1.2, charthick=2.0
        djs_xyouts, 2050.*(1.+zabs), !y.crange[0]+0.85*(!y.crange[1]-!y.crange[0]), 'SDEV='+$
            string(sdev, format='(f6.4)'), color='blue', charsize=1.2, charthick=2.0


      for iline=0L, n_elements(lines)-1L do begin
           djs_oplot, replicate(lines[iline].wave*(1.+zabs),2), $
              !y.crange[0]+[0.0, 1.0]*(!y.crange[1]-!y.crange[0]), $
              color='dark green', thick=thick, linestyle=2
       endfor

       ytitle = textoidl('Residual')
       pos = [0.50, 0.10, 0.90, 0.35]
       yra = [-0.1, 1.8]
       djs_plot, x, yresi, xra=xra, xstyle=1, $
           yra=yra, ystyle=1, position=pos, xtitle=xtitle, ytickformat='(A1)', $
           thick=thick-3, xthick=thick, ythick=thick, $
           charsize=charsize, charthick=charthick, /noerase, ytickinterval=400.
       djs_oplot, !x.crange, [1,1], color='blue', thick=thick-3

       for iline=0L, n_elements(lines)-1L do begin
           djs_oplot, replicate(lines[iline].wave*(1.+zabs),2), $
              !y.crange[0]+[0.0, 1.0]*(!y.crange[1]-!y.crange[0]), $
              color='dark green', thick=thick, linestyle=2
       endfor

   endfor

k_end_print


stop
end
