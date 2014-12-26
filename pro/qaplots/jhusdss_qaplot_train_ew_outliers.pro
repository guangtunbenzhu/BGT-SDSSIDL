pro jhusdss_qaplot_train_ew_outliers, nmfver

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

pitts = pitts_mgii
jhu = jhu_mgii

delta_ew = abs(jhu.rew_mgii_2796 - pitts.rew_mgii_2796)
ioutlier = where(delta_ew gt 0.8, noutlier)

ioutlier = reverse(ioutlier[bsort(jhu[ioutlier].rew_mgii_2796)])

print, noutlier, n_elements(jhu), n_elements(pitts)

stop

;; init
thick=6
charsize=1.2
charthick=2.5

qapath = path+'/QAplots/JHU_Pitts_Outliers'
if (jhusdss_direxist(qapath) eq 0) then spawn, 'mkdir -p '+qapath

nexam = 10L
npage = noutlier/nexam
for ipage=0, npage do begin

psfile = qapath+'/JHU_matched_outlier_'+string(nmfver, format='(I3.3)')+'_page'+string(ipage, format='(I2.2)')+'.ps'
k_print, filename=psfile, axis_char_scale=2.0, xsize=10, ysize=12
   !p.multi = [0, 2, nexam]
   !x.margin = 0
   !y.margin = 0
   for i=0L, nexam-1L do begin
       j = ipage*nexam+i
       if (j gt noutlier-1) then continue
       spec = jhusdss_decompose_loadspec(jhu[ioutlier[j]].plate, jhu[ioutlier[j]].fiber, nmfver)
;      name = 'SDSS'+hogg_iau_name(spec.ra, spec.dec, 'SDSS')
       name = hogg_iau_name(spec.ra, spec.dec, '')
       zabs = jhu[ioutlier[j]].zabs
       nabs = jhu[ioutlier[j]].nabs
       ew2796 = jhu[ioutlier[j]].rew_mgii_2796
       sigma2796 = sqrt(((jhu[ioutlier[j]].vdisp/69.03)^2+1.))*1.D-4*alog(10.)*2800.
       ew2803 = jhu[ioutlier[j]].rew_mgii_2803
       err_ew2796 = jhu[ioutlier[j]].err_rew_mgii_2796
       ew2796_pitts = pitts[ioutlier[j]].rew_mgii_2796
       ew2803_pitts = pitts[ioutlier[j]].rew_mgii_2803
       sdev = jhu[ioutlier[j]].med_sdeviation_red
       x = spec.wave*(1.+spec.z)
       yflux = spec.flux
       yresi = spec.residual
       ynmf = spec.nmf_continuum*spec.med_continuum
       ymed = spec.med_continuum
;      yflux = smooth(spec.flux, 5)
;      yresi = smooth(spec.residual, 5)
;      ynmf = smooth(spec.nmf_continuum*spec.med_continuum, 5)
;      ymed = smooth(spec.med_continuum, 5)

       ;;line profile fit
       par_jhu   =[2796.35, 2803.53, ew2796, ew2803, sigma2796]
       par_pitts =[2796.35, 2803.53, ew2796_pitts, ew2803_pitts, sigma2796]
       out_jhu = doublegaussian_line_func(x/(1.+zabs), par_jhu)
       out_pitts = doublegaussian_line_func(x/(1.+zabs), par_pitts)

;      xtitle = textoidl('Observer-frame \lambda (\AA)')
       xtitle = textoidl("Absorber Rest-frame \lambda (\AA)")
;      ytitle = textoidl('Flux (arbitrary unit)')
;      ytitle = textoidl('Quasar Flux')
       ytitle = textoidl('JHU')
;      xra = [3700, 9500]
;      yra = [0, 4]
       xra = [2650, 2899]
       yra = [0.0, 1.8]
       pos = [0.10, 0.35, 0.50, 0.90]
;      djs_plot, x, yflux, xra=xra, xstyle=1, $
       djs_plot, x/(1.+zabs), yresi, xra=xra, xstyle=1, $
           yra=yra, ystyle=1, xtickformat='(A1)', $
           thick=thick-3, xthick=thick, ythick=thick, $
           charsize=charsize, charthick=charthick, ytickformat='(A1)', ytickinterval=10, /nodata
;      djs_oplot, x, ynmf, color='red', thick=thick
       djs_oplot, x/(1.+zabs), 1.-out_jhu, color='red', thick=thick, linestyle=0

       ii = where(spec.ivar gt 0., kk)
       dwave = jhusdss_dwave(spec.wave[ii]*(1.+spec.z))
       jj = where(dwave gt 10., ll)
;      if (i eq 12) then print, x[ii[jj]]
       if (ll gt 0) then begin
;         print, ll
          jj = [[0], jj, [kk-1]]
          for j=0L, ll, 2L do begin
;             djs_oplot, x[ii[jj[j]:jj[j+1]]], yflux[ii[jj[j]:jj[j+1]]], thick=thick-3
              djs_oplot, x[ii[jj[j]:jj[j+1]]]/(1.+zabs), yresi[ii[jj[j]:jj[j+1]]], thick=thick-2
          endfor
       endif else begin
;         djs_oplot, x[ii], yflux[ii], thick=thick-3
          djs_oplot, x[ii]/(1.+zabs), yresi[ii], thick=thick-2
       endelse
       for iline=0L, n_elements(lines)-2L do begin
;          djs_oplot, replicate(lines[iline].wave*(1.+zabs),2), $
              djs_oplot, replicate(lines[iline].wave,2), $
              !y.crange[0]+[0.8, 1.0]*(!y.crange[1]-!y.crange[0]), $
              color='dark green', thick=thick, linestyle=2
       endfor

        djs_xyouts, 2660, !y.crange[0]+0.04*(!y.crange[1]-!y.crange[0]), name, $
             color='black', charsize=1.2, charthick=2.5
        djs_xyouts, 2660, !y.crange[0]+0.72*(!y.crange[1]-!y.crange[0]), $
             string(ew2796, format='(F3.1)')+'\pm'+string(err_ew2796, format='(F3.1)')+ $
             ' \AA', $
             color='red', charsize=1.2, charthick=2.5

       djs_oplot, !x.crange, [1,1], linestyle=2, color='gray'


        if (i eq 0) then djs_axis, xaxis=1, xtitle=ytitle, thick=thick, ythick=thick, $
             charsize=charsize, charthick=charthick, xtickformat='(A1)'

       if (i eq nexam-1) then djs_axis, axis=0, xtitle=xtitle, thick=thick, xthick=thick, $
           charsize=charsize, charthick=charthick

       xtitle = textoidl("Absorber Rest-frame \lambda (\AA)")
;      ytitle = textoidl('Residual')
       ytitle = textoidl('Pittsburgh')
;      xra = [2500, 3000]
       xra = [2650, 2899]
       yra = [0.0, 1.8]

       djs_plot, x/(1.+zabs), yresi, xra=xra, xstyle=1, $
           xtickformat='(A1)', yra=yra, ystyle=1, ytickformat='(A1)', $
           thick=thick-3, xthick=thick, ythick=thick, $
           charsize=charsize, charthick=charthick, ytickinterval=2., /nodata
       djs_oplot, x/(1.+zabs), 1.-out_pitts, color='red', thick=thick, linestyle=0

       ii = where(spec.ivar gt 0., kk)
       dwave = jhusdss_dwave(spec.wave[ii]*(1.+spec.z))
       jj = where(dwave gt 10., ll)
       if (ll gt 0) then begin
          jj = [[0], jj, [kk-1]]
          for j=0L, ll, 2 do begin
              djs_oplot, x[ii[jj[j]:jj[j+1]]]/(1.+zabs), yresi[ii[jj[j]:jj[j+1]]], thick=thick-2
          endfor
       endif else begin
          djs_oplot, x[ii]/(1.+zabs), yresi[ii], thick=thick-2
       endelse

       if i eq 0 then begin
          for iline=0L, n_elements(lines)-2L do begin
              djs_oplot, replicate(lines[iline].wave,2), $
                 !y.crange[0]+[0.8, 1.0]*(!y.crange[1]-!y.crange[0]), $
                 color='dark green', thick=thick, linestyle=2
          endfor
       endif
        djs_xyouts, 2660, !y.crange[0]+0.72*(!y.crange[1]-!y.crange[0]), $
             string(ew2796_pitts, format='(F3.1)')+'\pm'+string(err_ew2796, format='(F3.1)')+ $
             ' \AA', $
             color='red', charsize=1.2, charthick=2.5

       djs_oplot, !x.crange, [1,1], linestyle=2, color='gray'

        if (i eq 0) then djs_axis, xaxis=1, xtitle=ytitle, thick=thick, ythick=thick, $
             charsize=charsize, charthick=charthick, xtickformat='(A1)'

   endfor
   djs_axis, axis=0, xtitle=xtitle, thick=thick, xthick=thick, $
      charsize=charsize, charthick=charthick
k_end_print

endfor
stop

end
