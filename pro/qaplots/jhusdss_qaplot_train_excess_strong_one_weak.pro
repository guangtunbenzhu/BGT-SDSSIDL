pro jhusdss_qaplot_train_excess_strong_one_weak, nmfver

lines = jhusdss_train_lines()
path=jhusdss_get_path(/nmfqso)+'/'+string(nmfver, format='(I3.3)')

mfile = path+'/Absorbers/'+'JHU_Pitts_matched_'+string(nmfver, format='(I3.3)')+'.fits'
nomfile = path+'/Absorbers/'+'Pitts_nomatched_'+string(nmfver, format='(I3.3)')+'.fits'
jhu_nomfile = path+'/Absorbers/'+'JHU_nomatched_'+string(nmfver, format='(I3.3)')+'.fits'

;;;; below is for 106
;mfile = path+'/Absorbers_bk/'+'JHU_Pitts_matched_'+string(nmfver, format='(I3.3)')+'.fits'
;nomfile = path+'/Absorbers_bk/'+'Pitts_nomatched_'+string(nmfver, format='(I3.3)')+'.fits'
;jhu_nomfile = path+'/Absorbers_bk/'+'JHU_nomatched_'+string(nmfver, format='(I3.3)')+'.fits'
;;;;

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

jhu_nomatch = jhu_nomatch_mgii
ijhu_high = where(jhu_nomatch.rew_mgii_2796 le 4.0 and jhu_nomatch.rew_mgii_2796 ge 2.0, njhu_high)
print, njhu_high
ijhu_high = reverse(ijhu_high[bsort(jhu_nomatch[ijhu_high].rew_mgii_2796)])

;; init
thick=6
charsize=1.2
charthick=2.5

qapath = path+'/QAplots'
if (jhusdss_direxist(qapath) eq 0) then spawn, 'mkdir -p '+qapath

nexam = 20L
npage = njhu_high/nexam
for ipage=0, npage do begin

psfile = qapath+'/JHU_nomatched_strong_one_weak_'+string(nmfver, format='(I3.3)')+'_page'+string(ipage, format='(I2.2)')+'.ps'
k_print, filename=psfile, axis_char_scale=2.0, xsize=10, ysize=12
;  nexam=15
;  nexam=njhu_high
   !p.multi = [0, 2, nexam]
   !x.margin = 0
   !y.margin = 0
;  for i=0L, njhu_high-1L do begin
   for i=0L, nexam-1L do begin
       j = ipage*nexam+i
       if (j gt njhu_high-1) then continue
       spec = jhusdss_decompose_loadspec(jhu_nomatch[ijhu_high[j]].plate, jhu_nomatch[ijhu_high[j]].fiber, nmfver)
;      name = 'SDSS'+hogg_iau_name(spec.ra, spec.dec, 'SDSS')
       name = hogg_iau_name(spec.ra, spec.dec, '')
       zabs = jhu_nomatch[ijhu_high[j]].zabs
       nabs = jhu_nomatch[ijhu_high[j]].nabs
       ew2796 = jhu_nomatch[ijhu_high[j]].rew_mgii_2796
       sdev = jhu_nomatch[ijhu_high[j]].med_sdeviation_red
       x = spec.wave*(1.+spec.z)
       yflux = smooth(spec.flux, 5)
       yresi = smooth(spec.residual, 5)
       ynmf = smooth(spec.nmf_continuum*spec.med_continuum, 5)
       ymed = smooth(spec.med_continuum, 5)

       xtitle = textoidl('Observer-frame \lambda (\AA)')
;      ytitle = textoidl('Flux (arbitrary unit)')
       ytitle = textoidl('Quasar Flux')
       xra = [3700, 9500]
       yra = [0, 4]
       pos = [0.10, 0.35, 0.50, 0.90]
       djs_plot, x, yflux, xra=xra, xstyle=1, $
           yra=yra, ystyle=1, xtickformat='(A1)', $
           thick=thick-3, xthick=thick, ythick=thick, $
           charsize=charsize, charthick=charthick, ytickformat='(A1)', ytickinterval=10, /nodata
;      djs_oplot, x, ynmf, color='red', thick=thick

       ii = where(spec.ivar gt 0., kk)
       dwave = jhusdss_dwave(spec.wave[ii]*(1.+spec.z))
       jj = where(dwave gt 10., ll)
;      if (i eq 12) then print, x[ii[jj]]
       if (ll gt 0) then begin
;         print, ll
          jj = [[0], jj, [kk-1]]
          for j=0L, ll, 2L do begin
              djs_oplot, x[ii[jj[j]:jj[j+1]]], yflux[ii[jj[j]:jj[j+1]]], thick=thick-3
          endfor
       endif else begin
          djs_oplot, x[ii], yflux[ii], thick=thick-3
       endelse


       for iline=0L, n_elements(lines)-2L do begin
           djs_oplot, replicate(lines[iline].wave*(1.+zabs),2), $
              !y.crange[0]+[0.8, 1.0]*(!y.crange[1]-!y.crange[0]), $
              color='dark green', thick=thick, linestyle=2
       endfor

        djs_xyouts, 6400, !y.crange[0]+0.64*(!y.crange[1]-!y.crange[0]), name, $
             color='black', charsize=1.2, charthick=2.5

        if (i eq 0) then djs_axis, xaxis=1, xtitle=ytitle, thick=thick, ythick=thick, $
             charsize=charsize, charthick=charthick, xtickformat='(A1)'

       if (i eq nexam-1) then djs_axis, axis=0, xtitle=xtitle, thick=thick, xthick=thick, $
           charsize=charsize, charthick=charthick

       xtitle = textoidl("Absorber Rest-frame \lambda (\AA)")
       ytitle = textoidl('Residual')
       xra = [2100, 3000]
       yra = [0.0, 1.8]

       djs_plot, x/(1.+zabs), yresi, xra=xra, xstyle=1, $
           xtickformat='(A1)', yra=yra, ystyle=1, ytickformat='(A1)', $
           thick=thick-3, xthick=thick, ythick=thick, $
           charsize=charsize, charthick=charthick, ytickinterval=2., /nodata

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

        djs_xyouts, 2850, !y.crange[0]+0.72*(!y.crange[1]-!y.crange[0]), $
;            'W^{\lambda2796}_0='+string(ew2796, format='(F3.1)')+' \AA', $
             string(ew2796, format='(F3.1)')+' \AA', $
             color='black', charsize=1.2, charthick=2.5

        djs_oplot, !x.crange, [1,1], linestyle=2, color='gray'

        if (i eq 0) then djs_axis, xaxis=1, xtitle=ytitle, thick=thick, ythick=thick, $
             charsize=charsize, charthick=charthick, xtickformat='(A1)'

   endfor
   djs_axis, axis=0, xtitle=xtitle, thick=thick, xthick=thick, $
      charsize=charsize, charthick=charthick
k_end_print

endfor

end
