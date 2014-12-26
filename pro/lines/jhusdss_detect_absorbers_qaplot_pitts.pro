;+
; I just have to record this: I wrote this piece of code on 12/31/2011.
;  -- Guangtun Benjamin Zhu, postdoc at the physics and astronomy department, Johns Hopkins University
;-

pro jhusdss_detect_absorbers_qaplot_pitts, spec, absorbers0, mgii, lines, writeout=writeout, path=path, boss=boss

if (keyword_set(writeout)) then begin
   thick=6
   charthick=1.7
   charsize=1.2
   if (n_elements(path) eq 0) then path='.'
   if (~keyword_set(boss)) then begin
      filename = 'Pitts_Absorber_'+string(spec.plate, format='(i4.4)')+'_'+$
                          string(spec.fiber, format='(i4.4)')+'.ps'
   endif else begin
      filename = 'Pitts_Absorber_BOSS_'+string(spec.plate, format='(i4.4)')+'_'+$
                          string(spec.fiber, format='(i4.4)')+'.ps'
   endelse
   print, filename
   filename = path+'/'+filename
   k_print, filename=filename, axis_char_scale=1.3, xsize=10, ysize=8
endif else begin
   thick=2
   charthick=1.8
   charsize=1.5
endelse


ii = where(spec.ivar gt 0., mm)

for i=0, 1 do begin
case i of
     0: tmp = absorbers0
     1: tmp = mgii
endcase

window, i
xtitle1 = 'Observer-Frame Wavelength (\AA)'
xtitle2 = 'QSO Rest-Frame Wavelength (\AA)'
ytitle1 = 'Flux (Arbitrary Unit)'
ytitle2 = 'Residual'

djs_plot, spec.wave[ii]*(1.+spec.z), spec.flux[ii], xra=[3700., 9500.], xstyle=1, $
    yra=[0,6], ystyle=1, position=[0.10, 0.55, 0.95, 0.93], xtitle=xtitle1, ytitle=ytitle1
djs_oplot, spec.wave[ii]*(1.+spec.z), spec.nmf_continuum[ii], color='red', thick=thick

for iline=0L, n_elements(lines)-1L do djs_oplot, $
    replicate(lines[iline].wave*(1.+tmp.zabs[0]),2), $
    !y.crange[0]+[0.6, 1.0]*(!y.crange[1]-!y.crange[0]), $
    color='red', thick=thick, linestyle=0
if (tmp.zabs[1] gt 0.) then begin
   for iline=0L, n_elements(lines)-1L do djs_oplot, $
       replicate(lines[iline].wave*(1.+tmp.zabs[1]),2), $
       !y.crange[0]+[0.6, 1.0]*(!y.crange[1]-!y.crange[0]), $
       color='dark green', thick=thick, linestyle=2
endif
if (tmp.zabs[2] gt 0.) then begin
   for iline=0L, n_elements(lines)-1L do djs_oplot, $
       replicate(lines[iline].wave*(1.+tmp.zabs[2]),2), $
       !y.crange[0]+[0.6, 1.0]*(!y.crange[1]-!y.crange[0]), $
       color='blue', thick=thick, linestyle=1
endif
if (tmp.zabs[3] gt 0.) then begin
   for iline=0L, n_elements(lines)-1L do djs_oplot, $
       replicate(lines[iline].wave*(1.+tmp.zabs[3]),2), $
       !y.crange[0]+[0.6, 1.0]*(!y.crange[1]-!y.crange[0]), $
       color='cyan', thick=thick, linestyle=1
endif

items = ['Absorber_1 z='+string(tmp.zabs[0], format='(f5.3)'), $
         'Absorber_2 z='+string(tmp.zabs[1], format='(f5.3)'), $
         'Absorber_3 z='+string(tmp.zabs[2], format='(f5.3)'), $
         'Absorber_4 z='+string(tmp.zabs[3], format='(f5.3)')]
colors = [djs_icolor('red'), djs_icolor('dark green'), djs_icolor('blue'), djs_icolor('cyan')]
legend, items, textcolors=colors, box=0, /right, charsize=charsize, charthick=charthick

djs_xyouts, 7800, !y.crange[0]+0.54*(!y.crange[1]-!y.crange[0]), 'Z(QSO)='+$
    string(spec.z, format='(f5.3)'), color='red', charsize=charsize, charthick=charthick
djs_xyouts, 7800, !y.crange[0]+0.48*(!y.crange[1]-!y.crange[0]), 'Ly\alpha (1216)='+$
    string(1216.*(1.+spec.z), format='(f6.1)'), color='red', charsize=charsize, charthick=charthick
djs_xyouts, 7800, !y.crange[0]+0.42*(!y.crange[1]-!y.crange[0]), 'CIV (1550)='+$
    string(1550.*(1.+spec.z), format='(f6.1)'), color='red', charsize=charsize, charthick=charthick
djs_xyouts, 7800, !y.crange[0]+0.36*(!y.crange[1]-!y.crange[0]), 'MgII (2800)='+$
    string(2800.*(1.+spec.z), format='(f7.1)'), color='red', charsize=charsize, charthick=charthick


djs_plot, spec.wave[ii], smooth(spec.residual[ii], 5), xra=[3700., 9500.]/(1.+spec.z), $
          yra=[0.5, 1.3], xstyle=1, ystyle=1, position=[0.10, 0.10, 0.95, 0.48], /noerase, $
          xtitle=xtitle2, ytitle=ytitle2
for iline=0L, n_elements(lines)-1L do djs_oplot, $
    replicate(lines[iline].wave*(1.+tmp.zabs[0])/(1.+spec.z),2), $
    !y.crange[0]+[0.6, 1.0]*(!y.crange[1]-!y.crange[0]), $
    color='red', thick=thick, linestyle=0
if (tmp.zabs[1] gt 0.) then $
   for iline=0L, n_elements(lines)-1L do djs_oplot, $
       replicate(lines[iline].wave*(1.+tmp.zabs[1])/(1.+spec.z),2), $
       !y.crange[0]+[0.6, 1.0]*(!y.crange[1]-!y.crange[0]), $
       color='dark green', thick=thick, linestyle=2
if (tmp.zabs[2] gt 0.) then $
   for iline=0L, n_elements(lines)-1L do djs_oplot, $
       replicate(lines[iline].wave*(1.+tmp.zabs[2])/(1.+spec.z),2), $
       !y.crange[0]+[0.6, 1.0]*(!y.crange[1]-!y.crange[0]), $
       color='blue', thick=thick, linestyle=1
if (tmp.zabs[3] gt 0.) then $
   for iline=0L, n_elements(lines)-1L do djs_oplot, $
       replicate(lines[iline].wave*(1.+tmp.zabs[3])/(1.+spec.z),2), $
       !y.crange[0]+[0.6, 1.0]*(!y.crange[1]-!y.crange[0]), $
       color='cyan', thick=thick, linestyle=1

djs_xyouts, lines[1].wave*(1.+tmp.zabs[0])/(1.+spec.z)-30., $
    !y.crange[0]+0.05*(!y.crange[1]-!y.crange[0]), 'Mg II', $
    color='red', charsize=charsize, charthick=charthick
djs_xyouts, lines[3].wave*(1.+tmp.zabs[0])/(1.+spec.z)-30., $
    !y.crange[0]+0.05*(!y.crange[1]-!y.crange[0]), 'Fe II', $
    color='red', charsize=charsize, charthick=charthick
djs_xyouts, lines[5].wave*(1.+tmp.zabs[0])/(1.+spec.z)-30., $
    !y.crange[0]+0.05*(!y.crange[1]-!y.crange[0]), 'Fe II', $
    color='red', charsize=charsize, charthick=charthick

;window, 1
;djs_plot, spec.wave[ii], smooth((1.-spec.residual[ii])*(sqrt(spec.ivar[ii])*spec.nmf_continuum[ii]*spec.med_continuum[ii]), 5), xra=[3700., 9500.]/(1.+spec.z), $
;          yra=[-1.5, 4.5], xstyle=1, ystyle=1, /noerase, $
;          xtitle=xtitle2, ytitle=ytitle2
;for iline=0L, n_elements(lines)-1L do djs_oplot, $
;    replicate(lines[iline].wave*(1.+absorbers0.zabs[0])/(1.+spec.z),2), $
;    !y.crange[0]+[0.6, 1.0]*(!y.crange[1]-!y.crange[0]), $
;    color='red', thick=thick, linestyle=0
;if (absorbers0.zabs[1] gt 0.) then $
;   for iline=0L, n_elements(lines)-1L do djs_oplot, $
;       replicate(lines[iline].wave*(1.+absorbers0.zabs[1])/(1.+spec.z),2), $
;       !y.crange[0]+[0.6, 1.0]*(!y.crange[1]-!y.crange[0]), $
;       color='dark green', thick=thick, linestyle=2
;if (absorbers0.zabs[2] gt 0.) then $
;   for iline=0L, n_elements(lines)-1L do djs_oplot, $
;       replicate(lines[iline].wave*(1.+absorbers0.zabs[2])/(1.+spec.z),2), $
;       !y.crange[0]+[0.6, 1.0]*(!y.crange[1]-!y.crange[0]), $
;       color='blue', thick=thick, linestyle=1
;if (absorbers0.zabs[3] gt 0.) then $
;   for iline=0L, n_elements(lines)-1L do djs_oplot, $
;       replicate(lines[iline].wave*(1.+absorbers0.zabs[3])/(1.+spec.z),2), $
;       !y.crange[0]+[0.6, 1.0]*(!y.crange[1]-!y.crange[0]), $
;       color='cyan', thick=thick, linestyle=1
;

endfor

if (keyword_set(writeout)) then k_end_print

end
