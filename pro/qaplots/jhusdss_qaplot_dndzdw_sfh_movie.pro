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

pro jhusdss_qaplot_dndzdw_sfh_movie, nmfver, log=log

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
charthick=4.0

  !p.multi=[0, 1, 1]
  !x.margin=0
  !y.margin=0
  z_nbin = n_elements(dndz.z)
  xra=[-0.1, 5.99]
; xra=[1.0, 10.]
  xmodel=(findgen(100)+1.)/100.*6.5
  xtitle=textoidl('Redshift z')
  ytitle=textoidl('dN/dz')

; w_min = [0.3, 0.6, 1.0]
  w_min = [1.0]
  w_max = [100.0]
  yra = [0.06, 1.99]

  ;; Simcoe 2012
;     zt = [[2.231, 2.727, 3.457, 4.786], [2.232, 2.723, 3.458, 4.779], [2.231, 2.722, 3.459, 4.776]]
;     ztlow = [[1.947, 2.461, 3.150, 4.345], [1.947, 2.461, 3.150, 4.345], [1.947, 2.461, 3.150, 4.345]]
;     ztup = [[2.461, 2.975, 3.805, 5.350], [2.461, 2.975, 3.805, 5.350], [2.461, 2.975, 3.805, 5.350]]
;     yt = [[0.362, 0.486, 0.453, 0.387], [0.343, 0.149, 0.432, 0.415], [0.509, 0.578, 0.461, 0.325]]
;     yterr = [[0.193, 0.176, 0.146, 0.195], [0.156, 0.087, 0.132, 0.187], [0.185, 0.170, 0.135, 0.164]]
      zt = [[2.232, 2.723, 3.458, 4.779], [2.231, 2.722, 3.459, 4.776]]
      ztlow = [[1.947, 2.461, 3.150, 4.345], [1.947, 2.461, 3.150, 4.345]]
      ztup = [[2.461, 2.975, 3.805, 5.350], [2.461, 2.975, 3.805, 5.350]]
      yt = [[0.343, 0.149, 0.432, 0.415], [0.509, 0.578, 0.461, 0.325]]
      yterr = [[0.156, 0.087, 0.132, 0.187], [0.185, 0.170, 0.135, 0.164]]

  for i=0, 2 do begin
 
psfile = qapath+'/'+'Cum_dndzdw_z_w1.0_movie_'+string(nmfver, format='(I3.3)')+'_pg'+string(i+1, format='(i1.1)')+'.ps'
;k_print, filename=psfile, axis_char_scale=1.2, xsize=10, ysize=7
k_print, filename=psfile, axis_char_scale=1.2, xsize=9, ysize=6
;     tmp = min(abs(10.^dndz.log10w_min-w_min[i]), iw_min) 
;     w_min_tmp = 10.^dndz.log10w_min[iw_min]
      tmp = min(abs(dndz.w_cum_min-w_min[0]), iw_min) 
      w_min_tmp = dndz.w_cum_min[iw_min]
;     tmp = min(abs(10.^dndz.log10w_min-w_max[i]), iw_max) 
;     w_max_tmp = 10.^dndz.log10w_min[iw_max]
      tmp = min(abs(dndz.w_cum_min-w_max[0]), iw_max) 
      w_max_tmp = dndz.w_cum_min[iw_max]

      print, w_min, w_min_tmp, w_max, w_max_tmp
      ii = where(dndz.cum_phi[iw_min,*] gt 0.)
      x = dndz.cum_median_z[iw_min,ii]
      y = dndz.cum_phi[iw_min,ii] - dndz.cum_phi[iw_max,ii]
      ;; use the dominant errors
      if (i ne 1) then yerr = sqrt(dndz.cum_phi_poisson_err[iw_min,ii]^2+dndz.cum_phi_poisson_err[iw_max, ii]^2) $
                  else yerr = dndz.cum_phi_poisson_err[iw_min,ii]

      djs_plot, x, y, psym=5, symsize=2.0, xra=xra, yra=yra, $
          thick=thick, xthick=thick, ythick=thick, $
          xtickformat='(A1)', ytitle=ytitle, /ylog,  $
          charsize=charsize, charthick=charthick, color='blue', $
          xstyle=1, ystyle=1
      oploterror, x, y, yerr, thick=thick, psym=5, symsize=1.5, color=djs_icolor('blue'), $
          errcolor=djs_icolor('blue')
;     djs_axis, xaxis=0, xtitle=xtitle, charsize=charsize, charthick=charthick
;     djs_axis, yaxis=0, ytitle=ytitle, charsize=charsize, charthick=charthick

;     ymodel = dndz.f0*(1.+xmodel)^4/(1.+(xmodel/dndz.f0_bbeta)^dndz.f0_ggamma)*dndz.w0/(1.+(xmodel/dndz.w0_bbeta)^dndz.w0_ggamma)*exp(-w_min_tmp/dndz.w0*(1.+(xmodel/dndz.w0_bbeta)^dndz.w0_ggamma)/(1.+xmodel)) $
;            - dndz.f0*(1.+xmodel)^4/(1.+(xmodel/dndz.f0_bbeta)^dndz.f0_ggamma)*dndz.w0/(1.+(xmodel/dndz.w0_bbeta)^dndz.w0_ggamma)*exp(-w_max_tmp/dndz.w0*(1.+(xmodel/dndz.w0_bbeta)^dndz.w0_ggamma)/(1.+xmodel)) 
      ymodel = dndz.f0*(1.+xmodel)^(dndz.f0_alpha+dndz.w0_alpha)/(1.+(xmodel/dndz.f0_bbeta)^dndz.f0_ggamma)*dndz.w0/(1.+(xmodel/dndz.w0_bbeta)^dndz.w0_ggamma)*exp(-w_min_tmp/dndz.w0*(1.+(xmodel/dndz.w0_bbeta)^dndz.w0_ggamma)/(1.+xmodel)^dndz.w0_alpha)
      ymodelmax = dndz.f0*(1.+xmodel)^(dndz.f0_alpha+dndz.w0_alpha)/(1.+(xmodel/dndz.f0_bbeta)^dndz.f0_ggamma)*dndz.w0/(1.+(xmodel/dndz.w0_bbeta)^dndz.w0_ggamma)*exp(-w_max_tmp/dndz.w0*(1.+(xmodel/dndz.w0_bbeta)^dndz.w0_ggamma)/(1.+xmodel)^dndz.w0_alpha)
;     djs_oplot, xmodel, ymodel-ymodelmax, thick=thick, color='blue',linestyle=0

      if (i eq 3) then legends = string(w_min,format='(f3.1)')+'<W_0<'+string(w_max, format='(f3.1)')+' \AA' $
                  else legends = 'W_0>'+string(w_min, format='(f3.1)')+' \AA' 
      djs_xyouts, !x.crange[0]+0.70*(!x.crange[1]-!x.crange[0]), 0.8, $
                  legends, charsize=1.6, charthick=charthick, color='blue'
 
     
      if (i ne 0) then begin
      oploterror, zt[*,1], yt[*,1], zt[*,1]-ztlow[*,1], yterr[*,1], /lobar, thick=thick, symsize=2.0, $
        psym=4, color=djs_icolor('light blue'), errcolor=djs_icolor('light blue')
;       psym=4, color=djs_icolor('blue'), errcolor=djs_icolor('blue')
      oploterror, zt[*,1], yt[*,1], ztup[*,1]-zt[*,1], yterr[*,1], /hibar, thick=thick, symsize=2.0, $
        psym=4, color=djs_icolor('light blue'), errcolor=djs_icolor('light blue')
;       psym=4, color=djs_icolor('blue'), errcolor=djs_icolor('blue')
      endif

      if (i eq 2) then begin
      asfh = 0.017
      bsfh = 0.13
      csfh = 3.3
      dsfh = 5.3
      ymodel_sfh = (asfh+bsfh*xmodel)/(1.+(xmodel/csfh)^dsfh)
      fudge = 1.65
      djs_oplot, xmodel, fudge*ymodel_sfh, color='red', linestyle=2, thick=thick

      apiece = -1.82
      bpiece = 3.28
      cpiece = -0.724
      dpiece = -0.26
      epiece = 4.99
      fpiece = -8.0
      zz1 = 1.04
      zz2 = 4.48
      ymodel_sfh = xmodel
      ii = where(xmodel lt zz1)
      ymodel_sfh[ii] = 10.^(bpiece*alog10(1.+xmodel[ii])+apiece)
      ii = where(xmodel ge zz1 and xmodel lt zz2)
      ymodel_sfh[ii] = 10.^(dpiece*alog10(1.+xmodel[ii])+cpiece)
      ii = where(xmodel gt zz2)
      ymodel_sfh[ii] = 10.^(fpiece*alog10(1.+xmodel[ii])+epiece)
      fudge = 2.0
;     djs_oplot, xmodel, fudge*ymodel_sfh, color='dark grey', linestyle=1, thick=thick
      endif

;     zp = [0.41, 0.50, 0.60, 0.70, 0.79, 0.89, 0.98, 1.08, 1.18, 1.27, 1.37, 1.46, 1.56, 1.66, 1.75, 1.85, 1.94, 2.04, 2.14, 2.23]
;     yp0 = [0.056, 0.059, 0.066, 0.081, 0.088, 0.089, 0.080, 0.091, 0.092, 0.089, 0.102, 0.088, 0.097, 0.100, 0.084, 0.089, 0.094, 0.076, 0.070, 0.065]
;     yperr0 = [0.006, 0.004, 0.004, 0.004, 0.004, 0.004, 0.004, 0.005, 0.005, 0.004, 0.005, 0.005, 0.005, 0.006, 0.007, 0.009, 0.010, 0.013, 0.020, 0.024]
;     dxdz = (1+zp)^2/sqrt(0.3*(1.+zp)^3+0.7)
;     yp = yp0*dxdz
;     yperr = yperr0*dxdz
      
;     oploterror, zp, yp, yperr, /lobar, thick=thick, symsize=1., $
;       psym=6, color=djs_icolor('cyan'), errcolor=djs_icolor('cyan')
;     oploterror, zp, yp, yperr, /hibar, thick=thick, symsize=1., $
;       psym=6, color=djs_icolor('cyan'), errcolor=djs_icolor('cyan')
;     if (i eq 0) then begin
;     !p.font=-1
;        items = [ANSI_Value('Zhu & M!31801!Xnard (2012)'), 'Matejek & Simcoe 2012']
      if (i ne 0) then begin
         items = ['Zhu & M!3'+string("351B)+'!Xnard (2012)', 'Matejek & Simcoe (2012)']
         psym1=[5,4]
      endif else begin
         items = ['Zhu & M!3'+string("351B)+'!Xnard (2012)', '  ']
         psym1=[5,3]
      endelse
         colors = [djs_icolor('blue'), djs_icolor('light blue')]
         legend, items, psym=psym1, colors=colors, number=1, thick=thick, $
           box=0, /left, /top, charsize=1.6, charthick=charthick, textcolors=colors, $
           pos=[0.2, 1.50]

;     endif else begin
      if (i eq 2) then begin
         items = ['Star Formation History (Hopkins & Beacom 2006)']
         colors = [djs_icolor('red')]
         legend, items, colors=colors, number=1, thick=thick, $
           box=0, /right, /top, charsize=1.3, charthick=charthick, textcolors=colors, $
           linestyle=[2], pspacing=1.5, $
           pos=[5.7, 0.06]
      endif
;     endif else begin
;        items = [ANSI_value('dN/dz Fit Zhu & M!3'+string("351B)+'!Xnard 2012')]
;        colors = [djs_icolor('blue')]
;        legend, items, colors=colors, number=1, thick=thick, $
;          box=0, /right, /top, charsize=1.3, charthick=charthick, textcolors=colors, $
;          linestyle=[0], pspacing=1.5, $
;          pos=[5.7, 0.06]
;     endelse

k_end_print
 spawn, 'convert -density 400 -quality 100% -transparent white '+psfile+' '+repstr(psfile, '.ps', '.png')
  endfor

;     items = ['Zhu & Menard this work', 'Zhu & Menard Fit (0.6<W_0^{\lambda2796}<5.0 \AA)', 'Matejek & Simcoe 2012']

;spawn, 'psselect -p1 '+psfile+' '+repstr(psfile, '.ps', '_pg1.ps')
;spawn, 'psselect -p2 '+psfile+' '+repstr(psfile, '.ps', '_pg2.ps')

end
