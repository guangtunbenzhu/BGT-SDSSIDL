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

pro jhusdss_qaplot_dndzdw_nestor, nmfver, log=log

if (n_elements(nmfver) eq 0) then begin
   splog, 'nmfver required.'
   return
endif

dndz = jhusdss_montecarlo_dndzdw_readin(nmfver)

path=jhusdss_get_path(/nmfqso)+'/'+string(nmfver, format='(I3.3)')
qapath = path+'/QAplots'
if (jhusdss_direxist(qapath) eq 0) then spawn, 'mkdir -p '+qapath

;; init
thick=8
charsize=1.5
charthick=3

thick=8
psfile = qapath+'/'+'Cum_dndzdw_z_one_Nestor_'+string(nmfver, format='(I3.3)')+'.ps'
k_print, filename=psfile, axis_char_scale=1.4, xsize=16, ysize=8

  pos = [0.1, 0.1, 0.48, 0.9]
; w_min = [0.3, 0.6, 1.0, 1.5, 2.0, 3.0, 4.0]
  z_nbin = n_elements(dndz.z)
  xra=[0.3, 3.2]
  xmodel=(findgen(100)+1.)/100.*2.4
  xtitle=textoidl('Redshift z')
  ytitle=textoidl('dN/dz (W_0^{\lambda2796}>W_{min})')
  yra=[2.E-4, 3.]

  wvector = make_vector(0.6, 4.2, 6)
  w_min = [0.2, wvector.bound_min]
  w_max = [0.6, wvector.bound_max]
; w_min = [wvector.bound_min]
; w_max = [wvector.bound_max]

  psyms = [8, 8, 8, 8, 8, 8, 8]
  colors = [djs_icolor('grey'), djs_icolor('red'), djs_icolor('magenta'), djs_icolor('brown'), $
            djs_icolor('dark green'), djs_icolor('blue'), djs_icolor('black')]
  ypos = [1.0, 0.6, 0.3, 0.15, 0.09, 0.02, 0.006]

  loadct, 5
; colors[0] = getcolor('gray', 0)
  for i=0L, n_elements(w_min)-1L do colors[i] = (200.-(i+1)/float(n_elements(w_min))*150.)
  for i=0L, n_elements(w_min)-1L do begin
 
;     tmp = min(abs(10.^dndz.log10w_min-w_min[i]), iw) 
;     w_min_tmp = 10.^dndz.log10w_min[iw]
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

      plotsym, 8, 0.8, /fill
;     if (i gt 0) then plotsym, 8, 0.8, /fill else plotsym, 8, 0.8
;     if (i gt 0) then plotsym, 8, 0.8, /fill else plotsym, 8, 0.6
;     oploterror, x[1:n_elements(x)-2], y[1:n_elements(x)-2], yerr[1:n_elements(x)-2], thick=thick, psym=psyms[i], $
;     if (i gt 0) then oploterror, x, y, yerr, thick=thick, psym=psyms[i], $
      oploterror, x, y, yerr, thick=thick, psym=psyms[i], $
          color=colors[i], errcolor=colors[i], symsize=2;, /nohat
;         color=djs_icolor(colors[i]), errcolor=djs_icolor(colors[i]), symsize=2
;     plotsym, 8, 0.8
;     oploterror, x[0], y[0], yerr[0], thick=3, psym=psyms[i], $
;         color=djs_icolor(colors[i]), errcolor=djs_icolor(colors[i]), symsize=2
;     oploterror, x[n_elements(x)-1], y[n_elements(x)-1], yerr[n_elements(x)-1], thick=3, psym=psyms[i], $
;         color=djs_icolor(colors[i]), errcolor=djs_icolor(colors[i]), symsize=2

      ymodel = dndz.f0*(1.+xmodel)^(dndz.f0_alpha+dndz.w0_alpha)/(1.+(xmodel/dndz.f0_bbeta)^dndz.f0_ggamma)*dndz.w0/(1.+(xmodel/dndz.w0_bbeta)^dndz.w0_ggamma)*exp(-w_min_tmp/dndz.w0*(1.+(xmodel/dndz.w0_bbeta)^dndz.w0_ggamma)/(1.+xmodel)^dndz.w0_alpha)
      if i eq 0 then lstyle=2 else lstyle=0
;     if i gt 0 then djs_oplot, xmodel, ymodel, thick=thick, color=djs_icolor(colors[i]), linestyle=lstyle
      if i gt 0 then djs_oplot, xmodel, ymodel, thick=thick, color=colors[i], linestyle=lstyle

      ymodel_nestor = 1.001*(1.+xmodel)^0.226*exp(-w_min_tmp/0.443/(1.+xmodel)^0.634)
      djs_oplot, xmodel, ymodel_nestor, thick=thick, color='black', linestyle=1

;     z_min=dndz.z_min[i]
;     z_max=dndz.z_max[i]
;     if (i eq 0) then z_min=0.36
      legends = 'W_0>'+string(w_min_tmp, format='(f3.1)')+' \AA'

      if (i eq 0) then djs_xyouts, !x.crange[0]+0.75*(!x.crange[1]-!x.crange[0]), ymodel[n_elements(xmodel)-1], $;ypos[i], $
                  legends, charsize=1.4, charthick=charthick, color=colors[i]
      if (i gt 0) then djs_xyouts, !x.crange[0]+0.75*(!x.crange[1]-!x.crange[0]), ymodel[n_elements(xmodel)-1], $;ypos[i], $
                  legends, charsize=1.4, charthick=charthick, color=colors[i]
  endfor

  pos = [0.59, 0.1, 0.97, 0.90]
  z_nbin = n_elements(dndz.z)
  xra=[0.3, 3.2]
  xmodel=(findgen(100)+1.)/100.*2.4
  xtitle=textoidl('Redshift z')
  ytitle=textoidl('dN/dz (W_{min}<W_0^{\lambda2796}<W_{max})')
  yra=[2.E-4, 0.9]

;w_min = [0.3, 0.6, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0]
;w_max = [0.6, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 10.0]

  wvector = make_vector(0.6, 4.2, 6)
  w_min = [0.3, wvector.bound_min]
  w_max = [0.6, wvector.bound_max]
; w_max[n_elements(w_max)-1] = 10.0
  
; w_bin = jhusdss_make_bins(alog10(0.3), alog10(5.87), alog10(1.5))
; w_min = 10.^w_bin.min
; w_max = 10.^w_bin.max

; w_min = [0.2, 0.6, 1.0, 1.5, 2.2, 3.0, 4.0]
; w_max = [0.6, 1.0, 1.4, 2.2, 3.0, 4.0, 10.0]

; psyms = [1, 2, 4, 5, 6, 7, 8]
  psyms = [1, 8, 8, 8, 8, 8, 8]
; colors = ['grey', 'magenta', 'brown', 'dark yellow', 'dark green', 'blue', 'black']
; colors = ['grey', 'red', 'magenta', 'brown', 'dark green', 'blue', 'black']
; for i=0L, n_elements(w_min)-1L do colors[i] = (i+1)/float(n_elements(w_min))*200.-1.
  ypos = [0.3, 0.2, 0.09, 0.06, 0.03, 0.01,0.006]
  for i=0, n_elements(w_min)-1 do begin
 
      plotsym, 5
;     tmp = min(abs(10.^dndz.log10w_min-w_min[i]), iw) 
;     w_min_tmp = 10.^dndz.log10w_min[iw]
      tmp = min(abs(dndz.w_cum_min-w_min[i]), iw) 
      w_min_tmp = dndz.w_cum_min[iw]
;     maxtmp = min(abs(10.^dndz.log10w_min-w_max[i]), iwmax) 
;     w_max_tmp = 10.^dndz.log10w_min[iwmax]
      maxtmp = min(abs(dndz.w_cum_min-w_max[i]), iwmax) 
      w_max_tmp = dndz.w_cum_min[iwmax]
      print, w_min[i], w_min_tmp

      ii = where(dndz.cum_phi[iw,*] gt 0.)
      x = dndz.cum_median_z[iw,ii]
      y = dndz.cum_phi[iw,ii]-dndz.cum_phi[iwmax,ii]
      yerr = sqrt(dndz.cum_phi_poisson_err[iw,ii]^2+dndz.cum_phi_poisson_err[iwmax,ii]^2)

      if (i eq 0) then $
          djs_plot, x, y, psym=4, xra=xra, yra=yra, $
             thick=thick, xthick=thick, ythick=thick, $
             xtitle=xtitle, ytitle=ytitle, /ylog, ytickformat='tick_exponent', $
             charsize=charsize, charthick=charthick, /nodata, pos=pos, /noerase

      if (i gt 0) then begin
         plotsym, 8, 0.8, /fill 
;        oploterror, x[1:n_elements(x)-2], y[1:n_elements(x)-2], yerr[1:n_elements(x)-2], thick=thick, psym=psyms[i], $
         oploterror, x, y, yerr, thick=thick, psym=psyms[i], $
            color=djs_icolor(colors[i]), errcolor=djs_icolor(colors[i]), symsize=2
;        plotsym, 8, 0.8
;        oploterror, x[0], y[0], yerr[0], thick=3, psym=psyms[i], $
;            color=djs_icolor(colors[i]), errcolor=djs_icolor(colors[i]), symsize=2
;        oploterror, x[n_elements(x)-1], y[n_elements(x)-1], yerr[n_elements(x)-1], thick=3, psym=psyms[i], $
;            color=djs_icolor(colors[i]), errcolor=djs_icolor(colors[i]), symsize=2
      endif

;     ymodel = dndz.f0*(1.+xmodel)^4/(1.+(xmodel/dndz.f0_bbeta)^dndz.f0_ggamma)*dndz.w0/(1.+(xmodel/dndz.w0_bbeta)^dndz.w0_ggamma)*exp(-w_min_tmp/dndz.w0*(1.+(xmodel/dndz.w0_bbeta)^dndz.w0_ggamma)/(1.+xmodel))
;     ymodelmax = dndz.f0*(1.+xmodel)^4/(1.+(xmodel/dndz.f0_bbeta)^dndz.f0_ggamma)*dndz.w0/(1.+(xmodel/dndz.w0_bbeta)^dndz.w0_ggamma)*exp(-w_max_tmp/dndz.w0*(1.+(xmodel/dndz.w0_bbeta)^dndz.w0_ggamma)/(1.+xmodel))
      ymodel = dndz.f0*(1.+xmodel)^(dndz.f0_alpha+dndz.w0_alpha)/(1.+(xmodel/dndz.f0_bbeta)^dndz.f0_ggamma)*dndz.w0/(1.+(xmodel/dndz.w0_bbeta)^dndz.w0_ggamma)*exp(-w_min_tmp/dndz.w0*(1.+(xmodel/dndz.w0_bbeta)^dndz.w0_ggamma)/(1.+xmodel)^dndz.w0_alpha)
      ymodelmax = dndz.f0*(1.+xmodel)^(dndz.f0_alpha+dndz.w0_alpha)/(1.+(xmodel/dndz.f0_bbeta)^dndz.f0_ggamma)*dndz.w0/(1.+(xmodel/dndz.w0_bbeta)^dndz.w0_ggamma)*exp(-w_max_tmp/dndz.w0*(1.+(xmodel/dndz.w0_bbeta)^dndz.w0_ggamma)/(1.+xmodel)^dndz.w0_alpha)

      if i eq 0 then lstyle=2 else lstyle=0
      if i gt 0 then djs_oplot, xmodel, ymodel-ymodelmax, thick=thick, color=djs_icolor(colors[i]), linestyle=lstyle
;     djs_oplot, xmodel, ymodel-ymodelmax, thick=thick, color=djs_icolor(colors[i]), linestyle=lstyle

;     z_min=dndz.z_min[i]
;     z_max=dndz.z_max[i]
;     if (i eq 0) then z_min=0.36
;     if i eq n_elements(w_min)-1 then $
;        legends = 'W>'+string(w_min_tmp, format='(f3.1)')+' \AA' $
;     else $
         legends = string(w_min_tmp, format='(f3.1)')+'<W_0<'+string(w_max_tmp,format='(f3.1)')+' \AA'
;     if (i eq 0) then djs_xyouts, !x.crange[0]+0.75*(!x.crange[1]-!x.crange[0]), (ymodel-ymodelmax)[n_elements(xmodel)-1]+0.1, $;ypos[i], $
;                 legends, charsize=1.4, charthick=charthick, color=colors[i]
      if (i gt 0) then djs_xyouts, !x.crange[0]+0.74*(!x.crange[1]-!x.crange[0]), (ymodel-ymodelmax)[n_elements(xmodel)-1], $;ypos[i], $
                  legends, charsize=1.4, charthick=charthick, color=colors[i]
  endfor

k_end_print
loadct, 0

stop
end
