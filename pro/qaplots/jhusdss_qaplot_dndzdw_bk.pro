FUNCTION tick_exponent, axis, index, number

     ; A special case.
     IF number EQ 0 THEN RETURN, '0'

     ; Assuming multiples of 10 with format.
     ex = String(number, Format='(e8.0)')
     pt = StrPos(ex, '.')

     first = StrMid(ex, 0, pt)
     sign = StrMid(ex, pt+2, 1)
     thisExponent = StrMid(ex, pt+3)

     ; Shave off leading zero in exponent
     WHILE StrMid(thisExponent, 0, 1) EQ '0' DO thisExponent = StrMid(thisExponent, 1)

     ; Fix for sign and missing zero problem.
     IF (Long(thisExponent) EQ 0) THEN BEGIN
        sign = ''
        thisExponent = '0'
     ENDIF

     ; Make the exponent a superscript.
     IF sign EQ '-' THEN BEGIN
;       RETURN, first + 'x10!U' + sign + thisExponent + '!N'
        RETURN, '10!U' + sign + thisExponent + '!N'
     ENDIF ELSE BEGIN
;       RETURN, first + 'x10!U' + thisExponent + '!N'
        RETURN, '10!U' + thisExponent + '!N'
     ENDELSE

END

pro jhusdss_qaplot_dndzdw, nmfver

if (n_elements(nmfver) eq 0) then begin
   splog, 'nmfver required.'
   return
endif

dndz = jhusdss_montecarlo_dndzdw_readin(nmfver)

path=jhusdss_get_path(/nmfqso)+'/'+string(nmfver, format='(I3.3)')
qapath = path+'/QAplots'
if (jhusdss_direxist(qapath) eq 0) then spawn, 'mkdir -p '+qapath

;; init
thick=6
charsize=1.5
charthick=3

psfile = qapath+'/'+'dndzdw_z_'+string(nmfver, format='(I3.3)')+'.ps'
k_print, filename=psfile, axis_char_scale=2., xsize=12, ysize=8
  !p.multi=[0, 3, 2]
  !x.margin=0
  !y.margin=0
  z_nbin = n_elements(dndz.z)
  xra=[0., 5.9]
  yra=[2.E-5, 9.]
  xmodel=findgen(100)/10.
  xtitle=textoidl('W_0^{MgII 2796} (\AA)')
  ytitle=textoidl('!9D!6^2N/!9D!6W!9D!6z')

  for i=0, z_nbin-1 do begin
      ii = where(dndz.phi[*,i] gt 0.)
      x = 10.^dndz.median_log10w[ii,i]
      y = dndz.phi[ii,i]
      yerr = dndz.phi_poisson_err[ii,i]
      if (i eq 0 or i eq 3 or i eq 6 or i eq 9) then begin
          djs_plot, x, y, psym=4, xra=xra, yra=yra, $
             thick=thick, xthick=thick, ythick=thick, $
             xtickformat='(A1)', ytickformat='tick_exponent', ytitle=ytitle, /ylog, $
;            xtickformat='(A1)', ytitle=ytitle, /ylog, $
             charsize=charsize, charthick=charthick
      endif else begin
          djs_plot, x, y, psym=4, xra=xra, yra=yra, $
             thick=thick, xthick=thick, ythick=thick, $
             xtickformat='(A1)', ytickformat='(A1)', /ylog, $
             charsize=charsize, charthick=charthick
      endelse
      oploterror, x, y, yerr, thick=thick, psym=4
      z_min=dndz.z_min[i]
      z_max=dndz.z_max[i]
      if (i eq 0) then z_min=0.43
      legends = string(z_min, format='(f4.2)')+'<z<'+ string(z_max, format='(f4.2)')
      djs_xyouts, !x.crange[0]+0.6*(!x.crange[1]-!x.crange[0]), $
                  10.^(!y.crange[0]+0.85*(!y.crange[1]-!y.crange[0])), $
                  legends, charsize=1.3, charthick=charthick

      p = [dndz.n_star_strong[i]/dndz.w_star_strong[i], dndz.w_star_strong[i], $
           dndz.n_star_weak[i]/dndz.w_star_weak[i], dndz.w_star_weak[i]]
      ymodel = jhusdss_dndzdw_noz_func2(xmodel, p)
      djs_oplot, xmodel, ymodel, color='red', thick=4
      if (i eq 0) then  ymodel1 = ymodel $
       else djs_oplot, xmodel, ymodel1, color='blue', linestyle=2, thick=4

;     if (i eq 0 or i eq 3) then $
;        djs_axis, yaxis=0, ytitle=ytitle, thick=thick, $
;          charsize=charsize, charthick=charthick

      if ((i ge 3 and i lt 6) or (i ge 9 and i lt 12)) then $
         djs_axis, xaxis=0, xtitle=xtitle, thick=thick, $
           charsize=charsize, charthick=charthick
   ;  if (i eq 1) then $
   ;     djs_axis, xaxis=1, xtitle=ytitle, thick=thick, $
   ;       charsize=charsize, charthick=charthick, xtickformat='(A1)'
  endfor
k_end_print

psfile = qapath+'/'+'Cum_dndzdw_z_'+string(nmfver, format='(I3.3)')+'.ps'
k_print, filename=psfile, axis_char_scale=2., xsize=12, ysize=8
  !p.multi=[0, 3, 2]
  !x.margin=0
  !y.margin=0
  z_nbin = n_elements(dndz.z)
  xra=[0.4, 2.4]
  xmodel=findgen(100)/3.
  xtitle=textoidl('Redshift z')
  ytitle=textoidl('dN/dz (>W_{min})')

  w_min = [0.3, 0.7, 1.0, 2.0, 3.0, 4.0]
  for i=0, n_elements(w_min)-1 do begin
 
      if (i ge 3) then yra=[2.E-4, 0.3] else yra=[0.05, 2]
      tmp = min(abs(10.^dndz.log10w_min-w_min[i]), iw) 
      w_min_tmp = 10.^dndz.log10w_min[iw]
      print, w_min[i], w_min_tmp
      ii = where(dndz.cum_phi[iw,*] gt 0.)
      x = dndz.cum_median_z[iw,ii]
      y = dndz.cum_phi[iw,ii]
      yerr = dndz.cum_phi_poisson_err[iw,ii]
      if (i eq 0 or i eq 3) then begin
          djs_plot, x, y, psym=4, xra=xra, yra=yra, $
             thick=thick, xthick=thick, ythick=thick, $
             xtickformat='(A1)', ytickformat='', ytitle='', /ylog, $
;            xtickformat='(A1)', ytitle=ytitle, /ylog, $
             charsize=charsize, charthick=charthick
      endif else begin
          djs_plot, x, y, psym=4, xra=xra, yra=yra, $
             thick=thick, xthick=thick, ythick=thick, $
             xtickformat='(A1)', ytickformat='(A1)', /ylog, $
             charsize=charsize, charthick=charthick
      endelse
      oploterror, x, y, yerr, thick=thick, psym=4

;     z_min=dndz.z_min[i]
;     z_max=dndz.z_max[i]
;     if (i eq 0) then z_min=0.43
      legends = 'W>'+string(w_min_tmp, format='(f3.1)')+' \AA'
      if (i ge 3) then ypos=10.^(!y.crange[0]+0.2*(!y.crange[1]-!y.crange[0])) $
      else ypos=10.^(!y.crange[0]+0.2*(!y.crange[1]-!y.crange[0]))
      djs_xyouts, !x.crange[0]+0.6*(!x.crange[1]-!x.crange[0]), ypos, $
                  legends, charsize=1.3, charthick=charthick

;     p = [dndz.n_star_strong[i]/dndz.w_star_strong[i], dndz.w_star_strong[i], $
;          dndz.n_star_weak[i]/dndz.w_star_weak[i], dndz.w_star_weak[i]]
      tmp_n_strong = dndz.n_star_strong_all*(1.+xmodel)^dndz.alpha_strong_all
      tmp_w_strong = dndz.w_star_strong_all*(1.+xmodel)^dndz.bbeta_strong_all
      tmp_n_weak = dndz.n_star_weak_all*(1.+xmodel)^dndz.alpha_weak_all
      tmp_w_weak = dndz.w_star_weak_all*(1.+xmodel)^dndz.bbeta_weak_all
      ymodel = tmp_n_strong*exp(-w_min_tmp/tmp_w_strong)+tmp_n_weak*exp(-w_min_tmp/tmp_w_weak)
      djs_oplot, xmodel, ymodel, color='red', thick=4
;     if (i eq 0) then  ymodel1 = ymodel $
;      else djs_oplot, xmodel, ymodel1, color='blue', linestyle=2, thick=4

;     if (i eq 0 or i eq 3) then $
;        djs_axis, yaxis=0, ytitle=ytitle, thick=thick, $
;          charsize=charsize, charthick=charthick

      if (i ge 3) then $
         djs_axis, xaxis=0, xtitle=xtitle, thick=thick, $
           charsize=charsize, charthick=charthick
      if (i eq 1) then $
         djs_axis, xaxis=1, xtitle=ytitle, thick=thick, $
           charsize=charsize, charthick=charthick, xtickformat='(A1)'
  endfor
k_end_print

psfile = qapath+'/'+'Cum_dndzdw_z_w1.0_'+string(nmfver, format='(I3.3)')+'.ps'
k_print, filename=psfile, axis_char_scale=1.3, xsize=8, ysize=8
  !p.multi=[0, 1, 1]
  !x.margin=0
  !y.margin=0
  z_nbin = n_elements(dndz.z)
  xra=[0., 5.9]
  xmodel=findgen(100)/3.
  xtitle=textoidl('Redshift z')
  ytitle=textoidl('dN/dz (W>1.0 \AA)')

  w_min = 1.0
  yra = [0.03, 1.5]
  for i=0, n_elements(w_min)-1 do begin
 
      tmp = min(abs(10.^dndz.log10w_min-w_min[i]), iw) 
      w_min_tmp = 10.^dndz.log10w_min[iw]
      print, w_min[i], w_min_tmp
      ii = where(dndz.cum_phi[iw,*] gt 0.)
      x = dndz.cum_median_z[iw,ii]
      y = dndz.cum_phi[iw,ii]
      yerr = dndz.cum_phi_poisson_err[iw,ii]
      djs_plot, x, y, psym=5, xra=xra, yra=yra, $
          thick=thick, xthick=thick, ythick=thick, $
          xtitle=xtitle, ytitle=ytitle, /ylog, $
          charsize=charsize, charthick=charthick
      oploterror, x, y, yerr, thick=thick, psym=5

;     legends = 'W>'+string(w_min_tmp, format='(f3.1)')+' \AA'
;     ypos=10.^(!y.crange[0]+0.2*(!y.crange[1]-!y.crange[0]))
;     djs_xyouts, !x.crange[0]+0.6*(!x.crange[1]-!x.crange[0]), ypos, $
;                 legends, charsize=1.3, charthick=charthick

      tmp_n_strong = dndz.n_star_strong_all*(1.+xmodel)^dndz.alpha_strong_all
      tmp_w_strong = dndz.w_star_strong_all*(1.+xmodel)^dndz.bbeta_strong_all
      tmp_n_weak = dndz.n_star_weak_all*(1.+xmodel)^dndz.alpha_weak_all
      tmp_w_weak = dndz.w_star_weak_all*(1.+xmodel)^dndz.bbeta_weak_all
      ymodel = tmp_n_strong*exp(-w_min_tmp/tmp_w_strong)+tmp_n_weak*exp(-w_min_tmp/tmp_w_weak)
      djs_oplot, xmodel, ymodel, color='red', thick=4

      zt = [2.231, 2.722, 3.459, 4.776]
      ztlow = [1.947, 2.461, 3.150, 4.345]
      ztup = [2.461, 2.975, 3.805, 5.350]
      yt = [0.509, 0.578, 0.461, 0.325]
      yterr = [0.185, 0.170, 0.135, 0.164]
 
      oploterror, zt, yt, zt-ztlow, yterr, /lobar, thick=thick, symsize=1.5, $
        psym=4, color=djs_icolor('green'), errcolor=djs_icolor('green')
      oploterror, zt, yt, ztup-zt, yterr, /hibar, thick=thick, symsize=1.5, $
        psym=4, color=djs_icolor('green'), errcolor=djs_icolor('green')

      zp = [0.41, 0.50, 0.60, 0.70, 0.79, 0.89, 0.98, 1.08, 1.18, 1.27, 1.37, 1.46, 1.56, 1.66, 1.75, 1.85, 1.94, 2.04, 2.14, 2.23]
      yp0 = [0.056, 0.059, 0.066, 0.081, 0.088, 0.089, 0.080, 0.091, 0.092, 0.089, 0.102, 0.088, 0.097, 0.100, 0.084, 0.089, 0.094, 0.076, 0.070, 0.065]
      yperr0 = [0.006, 0.004, 0.004, 0.004, 0.004, 0.004, 0.004, 0.005, 0.005, 0.004, 0.005, 0.005, 0.005, 0.006, 0.007, 0.009, 0.010, 0.013, 0.020, 0.024]
      dxdz = (1+zp)^2/sqrt(0.3*(1.+zp)^3+0.7)
      yp = yp0*dxdz
      yperr = yperr0*dxdz
      
      oploterror, zp, yp, yperr, /lobar, thick=thick, symsize=1., $
        psym=6, color=djs_icolor('cyan'), errcolor=djs_icolor('cyan')
      oploterror, zp, yp, yperr, /hibar, thick=thick, symsize=1., $
        psym=6, color=djs_icolor('cyan'), errcolor=djs_icolor('cyan')

  endfor
k_end_print

psfile = qapath+'/'+'Wstar_dndzdw_z_'+string(nmfver, format='(I3.3)')+'.ps'
k_print, filename=psfile, axis_char_scale=1.4, xsize=8, ysize=8
   xra=[0., 6.]
   yra=[0.2, 1.2]
   x = dndz.zbin_median
   y = dndz.w_star_strong
   xerr = dndz.zbin_sdev
   yerr = dndz.w_star_strong_err
   xmodel=findgen(100)/3.
   xtitle='Redshift z' 
   ytitle=textoidl('W^*_{strong} \AA')

   djs_plot, x, y, psym=5, xra=xra, yra=yra, symsize=1.5, $
        thick=thick, xthick=thick, ythick=thick, $
        xtitle=xtitle, ytitle=ytitle, $
        charsize=charsize, charthick=charthick
   oploterror, x, y, xerr, yerr, thick=thick, psym=4
       
   ymodel = dndz.w_star_strong_all*(1.+xmodel)^dndz.bbeta_strong_all
   djs_oplot, xmodel, ymodel, color='red', thick=4

   zt = [2.51, 3.46, 4.78]
   ztlow = [1.947, 3.150, 4.345]
   ztup = [2.975, 3.805, 5.350]
   wt = [0.935, 0.766, 0.700]
   wterr = [0.150, 0.152, 0.180]

   oploterror, zt, wt, zt-ztlow, wterr, /lobar, thick=thick, symsize=1.5, $
      psym=4, color=djs_icolor('green'), errcolor=djs_icolor('green')
   oploterror, zt, wt, ztup-zt, wterr, /hibar, thick=thick, symsize=1.5, $
      psym=4, color=djs_icolor('green'), errcolor=djs_icolor('green')

   items = ['Zhu & Menard this work', 'Matejek & Simcoe 2012']
   colors = [djs_icolor('black'), djs_icolor('green')]
   legend, items, psym=[5,4], colors=colors, number=1, thick=thick, $
     box=0, /right, /bottom, charsize=1.3, charthick=charthick
k_end_print

stop
end
