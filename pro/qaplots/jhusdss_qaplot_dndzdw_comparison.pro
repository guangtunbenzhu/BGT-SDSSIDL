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

pro jhusdss_qaplot_dndzdw_comparison, nmfver, log=log

if (n_elements(nmfver) eq 0) then begin
   splog, 'nmfver required.'
   return
endif

readcol, '/home/gz323/Code/jhu-sdss/pro/literature/Nestor2005/w0.6.txt',  $
	nestor06_z, nestor06_zlow, nestor06_zhigh, nestor06_logdndz, nestor06_logdndz_high
readcol, '/home/gz323/Code/jhu-sdss/pro/literature/Nestor2005/w1.0.txt',  $
	nestor10_z, nestor10_zlow, nestor10_zhigh, nestor10_logdndz, nestor10_logdndz_high
readcol, '/home/gz323/Code/jhu-sdss/pro/literature/Nestor2005/w1.5.txt',  $
	nestor15_z, nestor15_zlow, nestor15_zhigh, nestor15_logdndz, nestor15_logdndz_high
readcol, '/home/gz323/Code/jhu-sdss/pro/literature/Nestor2005/w2.0.txt',  $
	nestor20_z, nestor20_zlow, nestor20_zhigh, nestor20_logdndz, nestor20_logdndz_high
;readcol, '/home/gz323/Code/jhu-sdss/pro/literature/Nestor2005/w3.0.txt',  $
;	nestor30_z, nestor30_zlow, nestor30_zhigh, nestor30_logdndz, nestor30_logdndz_high
readcol, '/home/gz323/Code/jhu-sdss/pro/literature/Nestor2005/w2.5.txt',  $
	nestor30_z, nestor30_zlow, nestor30_zhigh, nestor30_logdndz, nestor30_logdndz_high
readcol, '/home/gz323/Code/jhu-sdss/pro/literature/Lundgren2009/table2.txt',  $
	lundgren_z, lundgren08_dndz, lundgren08_dndz_err, lundgren10_dndz, lundgren10_dndz_err, lundgren15_dndz, lundgren15_dndz_err, lundgren20_dndz, lundgren20_dndz_err
readcol, '/home/gz323/Code/jhu-sdss/pro/literature/Prochter2006/w1.0.txt',  $
	prochter_z, prochter10_dndx, prochter10_dndx_err, prochter10_n, prochter14_dndx, prochter14_dndx_err, prochter14_n, prochter18_dndx, prochter18_dndx_err, prochter18_n
prochter_dxdz = jhusdss_dxdz(prochter_z)/0.7

dndz = jhusdss_montecarlo_dndzdw_readin(nmfver)

path=jhusdss_get_path(/nmfqso)+'/'+string(nmfver, format='(I3.3)')
qapath = path+'/QAplots'
if (jhusdss_direxist(qapath) eq 0) then spawn, 'mkdir -p '+qapath

;; init
thick=8
charsize=1.5
charthick=3

psfile = qapath+'/'+'Comparison_Cum_dndzdw_z_one_'+string(nmfver, format='(I3.3)')+'.ps'
k_print, filename=psfile, axis_char_scale=1.4, xsize=18, ysize=8

  z_nbin = n_elements(dndz.z)
  xra=[0.3, 3.0]
  xmodel=(findgen(100)+1.)/100.*2.4
  xtitle=textoidl('z')
  ytitle=textoidl('dN/dz (W_0^{\lambda2796}>W_{min})')
  yra=[8.E-3, 1E0]

  w_min = [0.6, 1.0, 1.5, 2.0, 2.5]

  psyms = [8, 8, 8, 8, 8]
  colors = [djs_icolor('grey'), djs_icolor('red'), djs_icolor('magenta'), djs_icolor('brown'), djs_icolor('brown')]

; pos = [0.1, 0.1, 0.48, 0.9]
  pos = [0.1, 0.12, 0.36, 0.9]
  loadct, 5
  for i=0L, n_elements(w_min)-1L do colors[i] = (175.-(i+1)/float(n_elements(w_min))*150.)

  for i=0L, n_elements(w_min)-1L do begin
 
      tmp = min(abs(dndz.w_cum_min-w_min[i]), iw) 
      w_min_tmp = dndz.w_cum_min[iw]
      print, w_min[i], w_min_tmp
      ii = where(dndz.cum_phi[iw,*] gt 0.)
      x = dndz.cum_median_z[iw,ii]
      y = dndz.cum_phi[iw,ii]
      yerr = dndz.cum_phi_poisson_err[iw,ii]
      if (i eq 0) then $
          djs_plot, x, y, psym=4, xra=xra, yra=yra, $
             thick=thick, xthick=thick, ythick=thick, $
             xtitle=xtitle, ytitle=ytitle, /ylog, ytickformat='tick_exponent', $
             charsize=charsize, charthick=charthick, /nodata, pos=pos

      case i of
         0: begin
            z_nestor = nestor06_z
            zlow_nestor = nestor06_zlow
            zhigh_nestor = nestor06_zhigh
            logdndz_nestor = nestor06_logdndz - 1.
            logdndz_high_nestor = nestor06_logdndz_high - 1.
	    psymbol = 8
            end
         1: begin
            z_nestor = nestor10_z
            zlow_nestor = nestor10_zlow
            zhigh_nestor = nestor10_zhigh
            logdndz_nestor = nestor10_logdndz - 1.
            logdndz_high_nestor = nestor10_logdndz_high - 1.
	    psymbol = 4
            end
         2: begin
            z_nestor = nestor15_z
            zlow_nestor = nestor15_zlow
            zhigh_nestor = nestor15_zhigh
            logdndz_nestor = nestor15_logdndz - 1.
            logdndz_high_nestor = nestor15_logdndz_high - 1.
	    psymbol = 0
            end
         3: begin
            z_nestor = nestor20_z
            zlow_nestor = nestor20_zlow
            zhigh_nestor = nestor20_zhigh
            logdndz_nestor = nestor20_logdndz - 1.
            logdndz_high_nestor = nestor20_logdndz_high - 1.
	    psymbol = 5
            end
         4: begin
            z_nestor = nestor30_z
            zlow_nestor = nestor30_zlow
            zhigh_nestor = nestor30_zhigh
            logdndz_nestor = nestor30_logdndz
            logdndz_high_nestor = nestor30_logdndz_high
	    psymbol = 3
            end

      endcase

;     plotsym, psymbol, 0.8, /fill
;     oploterror, x, y, yerr, thick=thick, psym=psyms[i], $
;         color=djs_icolor('black'), errcolor=djs_icolor('black'), symsize=2.5;, /nohat

      plotsym, psymbol, 0.8
      oploterror, z_nestor, 10.^logdndz_nestor, z_nestor-zlow_nestor, 10.^logdndz_high_nestor-10.^logdndz_nestor, $
;     oploterror, z_nestor, 10.^logdndz_nestor, 10.^logdndz_high_nestor-10.^logdndz_nestor, $
          thick=thick-2, psym=psyms[i], color=colors[i], errcolor=colors[i], symsize=2, /lobar
      oploterror, z_nestor, 10.^logdndz_nestor, zhigh_nestor-z_nestor, 10.^logdndz_high_nestor-10.^logdndz_nestor, $
;     oploterror, z_nestor, 10.^logdndz_nestor, 10.^logdndz_high_nestor-10.^logdndz_nestor, $
          thick=thick-2, psym=psyms[i], color=colors[i], errcolor=colors[i], symsize=2, /hibar

      plotsym, psymbol, 0.8, /fill
      oploterror, x, y, yerr, thick=thick, psym=psyms[i], $
          color=colors[i], errcolor=colors[i], symsize=2.8;, /nohat

      ymodel = dndz.f0*(1.+xmodel)^(dndz.f0_alpha+dndz.w0_alpha)/(1.+(xmodel/dndz.f0_bbeta)^dndz.f0_ggamma)*dndz.w0/(1.+(xmodel/dndz.w0_bbeta)^dndz.w0_ggamma)*exp(-w_min_tmp/dndz.w0*(1.+(xmodel/dndz.w0_bbeta)^dndz.w0_ggamma)/(1.+xmodel)^dndz.w0_alpha)
;     if i eq 0 then lstyle=2 else lstyle=0
;     if i ge 0 then djs_oplot, xmodel, ymodel, thick=thick, color=colors[i], linestyle=1

      legends = 'W_0>'+string(w_min_tmp, format='(f3.1)')+' \AA'

      if (i eq 0) then djs_xyouts, !x.crange[0]+0.78*(!x.crange[1]-!x.crange[0]), ymodel[n_elements(xmodel)-1], $
                  legends, charsize=1.4, charthick=charthick, color=colors[i]
      if (i gt 0) then djs_xyouts, !x.crange[0]+0.78*(!x.crange[1]-!x.crange[0]), ymodel[n_elements(xmodel)-1], $
                  legends, charsize=1.4, charthick=charthick, color=colors[i]
  endfor

; djs_xyouts, !x.crange[0]+0.30*(!x.crange[1]-!x.crange[0]), 0.015, $
;     "Solid Symbols: Zhu13", charsize=1.4, charthick=charthick, color='black'
  djs_xyouts, !x.crange[0]+0.30*(!x.crange[1]-!x.crange[0]), 0.011, $
      "Open Symbols: Nestor05", charsize=1.4, charthick=charthick, color='black'

  w_min = [1.0, 1.8]
  pos = [0.41, 0.12, 0.67, 0.90]
  xra=[0.3, 3.0]
; yra=[2.E-2, 5E-1]
; yra=[2.E-2, 1E0]
  loadct, 5
  for i=0L, n_elements(w_min)-1L do colors[i] = (200.-(i+1)/float(n_elements(w_min))*150.)

  for i=0L, n_elements(w_min)-1L do begin
 
      tmp = min(abs(dndz.w_cum_min-w_min[i]), iw) 
      w_min_tmp = dndz.w_cum_min[iw]
      print, w_min[i], w_min_tmp
      ii = where(dndz.cum_phi[iw,*] gt 0.)
      x = dndz.cum_median_z[iw,ii]
      y = dndz.cum_phi[iw,ii]
      yerr = dndz.cum_phi_poisson_err[iw,ii]
      if (i eq 0) then $
          djs_plot, x, y, psym=4, xra=xra, yra=yra, $
             thick=thick, xthick=thick, ythick=thick, $
             xtitle=xtitle, ytitle='', /ylog, ytickformat='tick_exponent', $
             charsize=charsize, charthick=charthick, /nodata, pos=pos, /noerase

      case i of
         0: begin
            dndz_prochter = prochter10_dndx*prochter_dxdz;/0.7
            dndz_err_prochter = prochter10_dndx_err*prochter_dxdz;/0.7
	    psymbol = 4
            end
         1: begin
            dndz_prochter = prochter18_dndx*prochter_dxdz;/0.7
            dndz_err_prochter = prochter18_dndx_err*prochter_dxdz;/0.7
	    psymbol = 5
            end
      endcase

;     plotsym, psymbol, 0.8, /fill
;     oploterror, x, y, yerr, thick=thick, psym=psyms[i], $
;         color=djs_icolor('black'), errcolor=djs_icolor('black'), symsize=2.5;, /nohat

;     print, lundgren_z, dndz_lundgren, dndz_err_lundgren
      plotsym, psymbol, 0.8
      oploterror, prochter_z, dndz_prochter, dndz_err_prochter, $
;     oploterror, prochter_z, dndz_prochter, 0.05*(fltarr(n_elements(prochter_z))+1.), dndz_err_prochter, $
;         thick=thick, psym=psyms[i], color=djs_icolor('grey'), errcolor=djs_icolor('grey'), symsize=2, /lobar
          thick=thick-2, psym=psyms[i], color=colors[i], errcolor=colors[i], symsize=2, /lobar
;     oploterror, prochter_z, dndz_prochter, 0.05*(fltarr(n_elements(prochter_z))+1.), dndz_err_prochter, $
      oploterror, prochter_z, dndz_prochter, dndz_err_prochter, $
;         thick=thick, psym=psyms[i], color=djs_icolor('grey'), errcolor=djs_icolor('grey'), symsize=2, /hibar
          thick=thick-2, psym=psyms[i], color=colors[i], errcolor=colors[i], symsize=2, /hibar

      plotsym, psymbol, 0.8, /fill
      oploterror, x, y, yerr, thick=thick, psym=psyms[i], $
          color=colors[i], errcolor=colors[i], symsize=2.8;, /nohat

      ymodel = dndz.f0*(1.+xmodel)^(dndz.f0_alpha+dndz.w0_alpha)/(1.+(xmodel/dndz.f0_bbeta)^dndz.f0_ggamma)*dndz.w0/(1.+(xmodel/dndz.w0_bbeta)^dndz.w0_ggamma)*exp(-w_min_tmp/dndz.w0*(1.+(xmodel/dndz.w0_bbeta)^dndz.w0_ggamma)/(1.+xmodel)^dndz.w0_alpha)
;     if i eq 0 then lstyle=2 else lstyle=0
;     if i ge 0 then djs_oplot, xmodel, ymodel, thick=thick, color=colors[i], linestyle=1

      legends = 'W_0>'+string(w_min_tmp, format='(f3.1)')+' \AA'

      if (i eq 0) then djs_xyouts, !x.crange[0]+0.78*(!x.crange[1]-!x.crange[0]), ymodel[n_elements(xmodel)-1], $
                  legends, charsize=1.4, charthick=charthick, color=colors[i]
      if (i gt 0) then djs_xyouts, !x.crange[0]+0.78*(!x.crange[1]-!x.crange[0]), ymodel[n_elements(xmodel)-1], $
                  legends, charsize=1.4, charthick=charthick, color=colors[i]


  endfor
  djs_xyouts, !x.crange[0]+0.30*(!x.crange[1]-!x.crange[0]), 0.011, $
      "Open Symbols: Prochter06", charsize=1.4, charthick=charthick, color='black'

  w_min = [0.8, 1.0, 1.5, 2.0]
  pos = [0.72, 0.12, 0.98, 0.90]
  xra=[0.3, 1.00]
  yra=[1.5E-2, 5E-1]
  loadct, 5
  for i=0L, n_elements(w_min)-1L do colors[i] = (200.-(i+1)/float(n_elements(w_min))*150.)

  for i=0L, n_elements(w_min)-1L do begin
 
      tmp = min(abs(dndz.w_cum_min-w_min[i]), iw) 
      w_min_tmp = dndz.w_cum_min[iw]
      print, w_min[i], w_min_tmp
      ii = where(dndz.cum_phi[iw,*] gt 0.)
      x = dndz.cum_median_z[iw,ii]
      y = dndz.cum_phi[iw,ii]
      yerr = dndz.cum_phi_poisson_err[iw,ii]
      if (i eq 0) then $
          djs_plot, x, y, psym=4, xra=xra, yra=yra, $
             thick=thick, xthick=thick, ythick=thick, $
             xtitle=xtitle, ytitle='', /ylog, ytickformat='tick_exponent', $
             charsize=charsize, charthick=charthick, /nodata, pos=pos, /noerase

      case i of
         0: begin
            dndz_lundgren = lundgren08_dndz
            dndz_err_lundgren = lundgren08_dndz_err
	    psymbol = 8
            end
         1: begin
            dndz_lundgren = lundgren10_dndz
            dndz_err_lundgren = lundgren10_dndz_err
	    psymbol = 4
            end
         2: begin
            dndz_lundgren = lundgren15_dndz
            dndz_err_lundgren = lundgren15_dndz_err
	    psymbol = 0
            end
         3: begin
            dndz_lundgren = lundgren20_dndz
            dndz_err_lundgren = lundgren20_dndz_err
	    psymbol = 5
            end
      endcase

;     plotsym, psymbol, 0.8, /fill
;     oploterror, x[0:2], y[0:2], yerr[0:2], thick=thick, psym=psyms[i], $
;         color=djs_icolor('black'), errcolor=djs_icolor('black'), symsize=2.5;, /nohat

;     print, lundgren_z, dndz_lundgren, dndz_err_lundgren
      plotsym, psymbol, 0.8
      oploterror, lundgren_z, dndz_lundgren, 0.025*(fltarr(n_elements(lundgren_z))+1.), dndz_err_lundgren, $
;         thick=thick, psym=psyms[i], color=djs_icolor('grey'), errcolor=djs_icolor('grey'), symsize=2, /lobar
          thick=thick-2, psym=psyms[i], color=colors[i], errcolor=colors[i], symsize=2, /lobar
      oploterror, lundgren_z, dndz_lundgren, 0.025*(fltarr(n_elements(lundgren_z))+1.), dndz_err_lundgren, $
;         thick=thick, psym=psyms[i], color=djs_icolor('grey'), errcolor=djs_icolor('grey'), symsize=2, /hibar
          thick=thick-2, psym=psyms[i], color=colors[i], errcolor=colors[i], symsize=2, /hibar

      plotsym, psymbol, 0.8, /fill
      oploterror, x[0:2], y[0:2], yerr[0:2], thick=thick, psym=psyms[i], $
          color=colors[i], errcolor=colors[i], symsize=2.8;, /nohat

      ymodel = dndz.f0*(1.+xmodel)^(dndz.f0_alpha+dndz.w0_alpha)/(1.+(xmodel/dndz.f0_bbeta)^dndz.f0_ggamma)*dndz.w0/(1.+(xmodel/dndz.w0_bbeta)^dndz.w0_ggamma)*exp(-w_min_tmp/dndz.w0*(1.+(xmodel/dndz.w0_bbeta)^dndz.w0_ggamma)/(1.+xmodel)^dndz.w0_alpha)
;     if i eq 0 then lstyle=2 else lstyle=0
;     if i ge 0 then djs_oplot, xmodel, ymodel, thick=thick, color=colors[i], linestyle=1

      legends = 'W_0>'+string(w_min_tmp, format='(f3.1)')+' \AA'

      if (i eq 0) then djs_xyouts, !x.crange[0]+0.78*(!x.crange[1]-!x.crange[0]), y[2], $
                  legends, charsize=1.4, charthick=charthick, color=colors[i]
      if (i gt 0) then djs_xyouts, !x.crange[0]+0.78*(!x.crange[1]-!x.crange[0]), y[2], $
                  legends, charsize=1.4, charthick=charthick, color=colors[i]



  endfor
  djs_xyouts, !x.crange[0]+0.30*(!x.crange[1]-!x.crange[0]), 0.018, $
      "Open Symbols: Lundgren09", charsize=1.4, charthick=charthick, color='black'



k_end_print

end
