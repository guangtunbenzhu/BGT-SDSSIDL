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

pro jhusdss_qaplot_dndzdw_sfh, nmfver, log=log

if (n_elements(nmfver) eq 0) then begin
   splog, 'nmfver required.'
   return
endif

dndz = jhusdss_montecarlo_dndzdw_readin(nmfver)

sfrd_file='/home/gz323/DATA/SFRD/cosmic_SFRD.dat'
readcol, sfrd_file, indicator, sfrd_z, sfrd_z_higherr, sfrd_z_lowerr, sfrd, sfrd_higherr, sfrd_lowerr, sfrdcor, sfrcor_higherr, sfrdcor_lowerr, reference, $
         format='(a,f,f,f,f,f,f,f,f,f,a)'
path=jhusdss_get_path(/nmfqso)+'/'+string(nmfver, format='(I3.3)')
qapath = path+'/QAplots'
if (jhusdss_direxist(qapath) eq 0) then spawn, 'mkdir -p '+qapath

;; init
thick=8
charsize=1.5
charthick=3

  !p.multi=[0, 1, 1]
  !x.margin=0
  !y.margin=0
  z_nbin = n_elements(dndz.z)
  xra1=[-0.1, 5.99]
  xra2=[-0.1, 7.49]
; xra=[1.0, 10.]
  xmodel=(findgen(100)+1.)/100.*7.5
  xtitle=textoidl('z')
  ytitle=textoidl('dN/dz (Absorption)')

; w_min = [0.3, 0.6, 1.0]
  w_min = [0.6, 1.0, 2.0]
  w_max = [1.0, 100.0, 100.0]

  ;; Simcoe 2013
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

  for i=0, n_elements(w_min)-1 do begin
      if (i eq 0) then xra=xra1
      if (i eq 1) then xra=xra2
 
psfile = qapath+'/'+'Cum_dndzdw_z_w1.0_twofigures_'+string(nmfver, format='(I3.3)')+'_pg'+string(i+1, format='(i1.1)')+'.ps'
print, psfile
k_print, filename=psfile, axis_char_scale=1.2, xsize=10, ysize=7
      if (i le 1) then yra = [0.03, 3.0] else yra = [0.01, 0.5]
;     tmp = min(abs(10.^dndz.log10w_min-w_min[i]), iw_min) 
;     w_min_tmp = 10.^dndz.log10w_min[iw_min]
      tmp = min(abs(dndz.w_cum_min-w_min[i]), iw_min) 
      w_min_tmp = dndz.w_cum_min[iw_min]
;     tmp = min(abs(10.^dndz.log10w_min-w_max[i]), iw_max) 
;     w_max_tmp = 10.^dndz.log10w_min[iw_max]
      tmp = min(abs(dndz.w_cum_min-w_max[i]), iw_max) 
      w_max_tmp = dndz.w_cum_min[iw_max]

      print, w_min[i], w_min_tmp, w_max[i], w_max_tmp
      ii = where(dndz.cum_phi[iw_min,*] gt 0.)
      x = dndz.cum_median_z[iw_min,ii]
      y = dndz.cum_phi[iw_min,ii] - dndz.cum_phi[iw_max,ii]
      ;; use the dominant errors
      if (i ne 1) then yerr = sqrt(dndz.cum_phi_poisson_err[iw_min,ii]^2+dndz.cum_phi_poisson_err[iw_max, ii]^2) $
                  else yerr = dndz.cum_phi_poisson_err[iw_min,ii]

      djs_plot, x, y, psym=5, symsize=1.5, xra=xra, yra=yra, $
          thick=thick, xthick=thick, ythick=thick, $
          xtickformat='(A1)', ytitle='', /ylog,  yst=9, $
          charsize=charsize, charthick=charthick, color='blue', /nodata

      if (i eq 0) then begin
         dxdz = (1+xmodel)^2/sqrt(0.3*(1.+xmodel)^3+0.7)
         djs_oplot, xmodel, dxdz*0.1, thick=thick-0, color='gray',linestyle=2
         djs_axis, yaxis=1, ythick=thick, thick=thick, $
              yra=yra, ytitle='', charsize=charsize, charthick=charthick, ytickformat='(A1)'
      endif

      djs_oplot, x, y, psym=5, symsize=1.5, thick=thick, $
          charsize=charsize, charthick=charthick, color='blue'

      if (i ge 1) then begin
      if (i eq 1) then fudge = 1.00 else fudge = 0.3
      new_sfrd = 10.^sfrd*fudge*jhusdss_dxdz(sfrd_z)
      new_sfrd_higherr = (10.^(sfrd+sfrd_higherr)-10.^(sfrd))*fudge*jhusdss_dxdz(sfrd_z)
      new_sfrd_lowerr = (10.^(sfrd)-10.^(sfrd-sfrd_lowerr))*fudge*jhusdss_dxdz(sfrd_z)

      if (i eq 1) then begin
          sfh_ytitle='\rho_{SFR} dX/dz (M_\odot yr^{-1} Mpc^{-2}) (Emission)'
;         sfh_ytitle='\rho_{SFR} (M_\odot yr^{-1} Mpc^{-2})'
;         sfh_yra = yra/(fudge*jhusdss_dxdz(xra)*3E5/70.)
          sfh_yra = yra*3E5/70./fudge
          djs_xyouts, 8.3, 0.035, sfh_ytitle, orientation=90, color='gray', charsize=charsize, charthick=charthick
          djs_axis, yaxis=1, ythick=thick, thick=thick, $
              yra=sfh_yra, ytitle='', charsize=charsize, charthick=charthick, ytickformat='tick_exponent'
;         djs_axis, yaxis=1, ythick=thick, thick=thick, $
;             yra=sfh_yra, ytitle=sfh_ytitle, charsize=charsize, charthick=charthick, ytickformat='tick_exponent'
      endif

      plotsym, 0
      oploterror, sfrd_z, new_sfrd, sfrd_z_higherr, new_sfrd_higherr, /hibar, psym=8, color=djs_icolor('gray'), errcolor=djs_icolor('gray'), thick=2
      oploterror, sfrd_z, new_sfrd, sfrd_z_lowerr, new_sfrd_lowerr, /lobar, psym=8, color=djs_icolor('gray'), errcolor=djs_icolor('gray'), thick=2
      endif


      oploterror, x, y, yerr, thick=thick, psym=5, symsize=1.5, color=djs_icolor('blue'), $
          errcolor=djs_icolor('blue')
      djs_axis, xaxis=0, xtitle=xtitle, charsize=charsize, charthick=charthick
      djs_axis, yaxis=0, ytitle='', charsize=charsize, charthick=charthick
      djs_xyouts, -0.75, 0.09, ytitle, orientation=90, charsize=charsize+0.2, charthick=charthick, color='blue'

;     ymodel = dndz.f0*(1.+xmodel)^4/(1.+(xmodel/dndz.f0_bbeta)^dndz.f0_ggamma)*dndz.w0/(1.+(xmodel/dndz.w0_bbeta)^dndz.w0_ggamma)*exp(-w_min_tmp/dndz.w0*(1.+(xmodel/dndz.w0_bbeta)^dndz.w0_ggamma)/(1.+xmodel)) $
;            - dndz.f0*(1.+xmodel)^4/(1.+(xmodel/dndz.f0_bbeta)^dndz.f0_ggamma)*dndz.w0/(1.+(xmodel/dndz.w0_bbeta)^dndz.w0_ggamma)*exp(-w_max_tmp/dndz.w0*(1.+(xmodel/dndz.w0_bbeta)^dndz.w0_ggamma)/(1.+xmodel)) 
      ymodel = dndz.f0*(1.+xmodel)^(dndz.f0_alpha+dndz.w0_alpha)/(1.+(xmodel/dndz.f0_bbeta)^dndz.f0_ggamma)*dndz.w0/(1.+(xmodel/dndz.w0_bbeta)^dndz.w0_ggamma)*exp(-w_min_tmp/dndz.w0*(1.+(xmodel/dndz.w0_bbeta)^dndz.w0_ggamma)/(1.+xmodel)^dndz.w0_alpha)
      ymodelmax = dndz.f0*(1.+xmodel)^(dndz.f0_alpha+dndz.w0_alpha)/(1.+(xmodel/dndz.f0_bbeta)^dndz.f0_ggamma)*dndz.w0/(1.+(xmodel/dndz.w0_bbeta)^dndz.w0_ggamma)*exp(-w_max_tmp/dndz.w0*(1.+(xmodel/dndz.w0_bbeta)^dndz.w0_ggamma)/(1.+xmodel)^dndz.w0_alpha)
;     djs_oplot, xmodel, ymodel-ymodelmax, thick=thick-2, color='cyan blue',linestyle=0

      if (i lt 1) then legends = string(w_min[i],format='(f3.1)')+'<W_0<'+string(w_max[i], format='(f3.1)')+' \AA' $
                  else legends = 'W_0>'+string(w_min[i], format='(f3.1)')+' \AA' 
      djs_xyouts, !x.crange[0]+0.78*(!x.crange[1]-!x.crange[0]), 2.0, $
                  legends, charsize=1.4, charthick=charthick, color='blue'
 
      if (i le 1) then begin
      oploterror, zt[*,i], yt[*,i], zt[*,i]-ztlow[*,i], yterr[*,i], /lobar, thick=thick, symsize=1.5, $
;       psym=4, color=djs_icolor('dark green'), errcolor=djs_icolor('dark green')
        psym=5, color=djs_icolor('light blue'), errcolor=djs_icolor('light blue')
      oploterror, zt[*,i], yt[*,i], ztup[*,i]-zt[*,i], yterr[*,i], /hibar, thick=thick, symsize=1.5, $
;       psym=4, color=djs_icolor('dark green'), errcolor=djs_icolor('dark green')
        psym=5, color=djs_icolor('light blue'), errcolor=djs_icolor('light blue')
      endif

      if (i eq 1) then begin
      asfh = 0.017
      bsfh = 0.13
      csfh = 3.3
      dsfh = 5.3
      ymodel_sfh = (asfh+bsfh*xmodel)/(1.+(xmodel/csfh)^dsfh)
      fudge = 1.65
;     djs_oplot, xmodel, fudge*ymodel_sfh, color='black', linestyle=2, thick=thick

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
;     endif
      
;     oploterror, zp, yp, yperr, /lobar, thick=thick, symsize=1., $
;       psym=6, color=djs_icolor('cyan'), errcolor=djs_icolor('cyan')
;     oploterror, zp, yp, yperr, /hibar, thick=thick, symsize=1., $
;       psym=6, color=djs_icolor('cyan'), errcolor=djs_icolor('cyan')
;     if (i eq 0) then begin
;     !p.font=-1
;        items = [ANSI_Value('Zhu & M!31801!Xnard (2013)'), 'Matejek & Simcoe 2013']
         items = ['Zhu & M!3'+string("351B)+'!Xnard (2013)', 'Matejek & Simcoe (2012)']
;        colors = [djs_icolor('blue'), djs_icolor('dark green')]
         colors = [djs_icolor('blue'), djs_icolor('light blue')]
         legend, items, psym=[5,4], colors=colors, number=1, thick=thick, $
           box=0, /left, /top, charsize=1.3, charthick=charthick, textcolors=colors, $
           pos=[0.10, 2.20]
;     endif else begin
;if (keyword_set(ah)) then begin
      if (i ge 1) then begin
;        items = [ANSI_value('dN/dz Fit by Zhu & M!3'+string("351B)+'!Xnard (2013)'), '(Scaled) Cosmic SFH Fit by Hopkins & Beacom (2006)']
         items = [ANSI_value('dN/dz Fit by Zhu & M!3'+string("351B)+'!Xnard (2013)')]
;        items = [ANSI_value('dN/dz Fit by Zhu & M!3'+string("351B)+'!Xnard (2013)')]
;        colors = [djs_icolor('blue'), djs_icolor('gray')]
         colors = [djs_icolor('cyan blue')]

;         legend, items, colors=colors, number=1, thick=thick, $
;           box=0, /right, /top, charsize=1.3, charthick=charthick, textcolors=colors, $
;           linestyle=[0], pspacing=1.5, $
;;          linestyle=[0, 2], pspacing=1.5, $
;           pos=[6.1, 0.06]

         items = [ANSI_value('(Scaled) Star Formation Rate Density')]
;        colors = [djs_icolor('gray')]
         colors = [djs_icolor('gray')]
         plotsym, 0, 1, thick=5 
          legend, items, colors=colors, number=1, thick=thick, $
            box=0, /right, /top, charsize=1.3, charthick=charthick, textcolors=colors, psym=[8], $
            pos=[5.9, 0.05]
      endif else begin
;        items = [ANSI_value('dN/dz Fit Zhu & M!3'+string("351B)+'!Xnard 2013')]
         items = [ANSI_value('dN/dz Fit by Zhu & M!3'+string("351B)+'!Xnard (2013)'), 'Constant Comoving Density']
;        items = [ANSI_value('dN/dz Fit by Zhu & M!3'+string("351B)+'!Xnard (2013)'), 'Constant Comoving Density']
;        items = [ANSI_value('dN/dz Fit')]
         colors = [djs_icolor('cyan blue'), djs_icolor('gray')]
;        legend, items, colors=colors, number=1, thick=thick, $
;          box=0, /right, /top, charsize=1.3, charthick=charthick, textcolors=colors, $
;          linestyle=[0,2], pspacing=1.5, $
;          pos=[5.7, 0.06]

         items = ['Constant Comoving Density']
         colors = [djs_icolor('gray')]
         plotsym, 0, 1, thick=5 
         legend, items, colors=colors, number=1, thick=thick, $
           box=0, /right, /top, charsize=1.3, charthick=charthick, textcolors=colors, $
           linestyle=[2], pspacing=1.5, $
           pos=[5.7, 0.06]
      endelse
;endif

k_end_print
  endfor

;     items = ['Zhu & Menard this work', 'Zhu & Menard Fit (0.6<W_0^{\lambda2796}<5.0 \AA)', 'Matejek & Simcoe 2012']

;spawn, 'psselect -p1 '+psfile+' '+repstr(psfile, '.ps', '_pg1.ps')
;spawn, 'psselect -p2 '+psfile+' '+repstr(psfile, '.ps', '_pg2.ps')

end
