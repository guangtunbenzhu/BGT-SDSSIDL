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

pro jhusdss_qaplot_dndzdw, nmfver, log=log

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

psfile = qapath+'/'+'dndzdw_z_onepanel_'+string(nmfver, format='(I3.3)')+'.ps'
k_print, filename=psfile, axis_char_scale=1.4, xsize=8, ysize=12
  !p.multi=[0, 1, 1]
  !x.margin=0
  !y.margin=0
  z_nbin = n_elements(dndz.z)
  yra=[1.E-10, 8.]
  xmodel=findgen(100)/10.
  xtitle=textoidl('W_0^{\lambda2796} (\AA)')
  ytitle=textoidl('!9D!6^2N/!9D!6z!9D!6W^{\lambda2796}_0')

; zcolor = [djs_icolor('gray'), djs_icolor('black'), djs_icolor('pink'), djs_icolor('navy'), djs_icolor('orange'), $
  zcolor = [djs_icolor('black'), djs_icolor('pink'), djs_icolor('navy'), $
            djs_icolor('orange'), djs_icolor('blue'), djs_icolor('magenta red'), $
            djs_icolor('dark green'), djs_icolor('magenta'), djs_icolor('green'),  $
            djs_icolor('brown'), djs_icolor('dark cyan'), djs_icolor('red')]

  nmax = (z_nbin < 12)
  loadct, 4
  for i=0L, nmax-1L do zcolor[i] = (i+1)/float(nmax)*200.-1.
  for i=0L, nmax-1L do begin
      ii = where(dndz.phi[*,i] gt 0.)
;     x = 10.^dndz.median_log10w[ii,i]
      x = dndz.median_w[ii,i]
      y = dndz.phi[ii,i]
      yerr = dndz.phi_poisson_err[ii,i]
      plotsym, 0, 1.5, /fill 
      if (i mod 12 eq 0) then begin
          if (i eq 0) then begin  
             xra=[0.00, 5.00]
             ytickformat_tmp='tick_exponent'
             ytitle_tmp = ytitle
          endif else begin
             xra=[0.01, 5.00]
             ytickformat_tmp='(A1)'
             ytitle_tmp = ''
          endelse
          djs_plot, x, y, psym=8, xra=xra, yra=yra, $
             thick=thick, xthick=thick, ythick=thick, $
             xtitle=xtitle, ytickformat=ytickformat_tmp, ytitle=ytitle_tmp, /ylog, $
;            xtickformat='(A1)', ytitle=ytitle, /ylog, $
             charsize=charsize, charthick=charthick, color=zcolor[i], /nodata
      endif 
      dexp = 0.5
;     djs_oplot, x[2:*], y[2:*]*(1E-1)^((i mod 6)*dexp), thick=thick, psym=8, color=zcolor[i]
      plotsym, 0, 1.5, /fill 
      oploterror, x[2:*], y[2:*]*(1E-1)^((12- (i mod 12))*dexp), yerr[2:*]*(1E-1)^((12-(i mod 12))*dexp), $
          thick=thick, psym=8, /nohat, errthick=2, $
          color=zcolor[i], errcolor=zcolor[i]
      plotsym, 0, 1.5
      oploterror, x[0:1], y[0:1]*(1E-1)^((12-(i mod 12))*dexp), $
          yerr[0:1]*(1E-1)^((12-(i mod 12))*dexp), $
          thick=thick+4, psym=8, /nohat, errthick=2, $
          color=zcolor[i], errcolor=zcolor[i]

      z_min=dndz.z_min[i]
      z_max=dndz.z_max[i]

;     if (i eq 0) then z_min=0.37
      legends = string(z_min, format='(f4.2)')+'<z<'+ string(z_max, format='(f4.2)')
;     if (i mod 3 eq 0) then begin 
;        textyp = 0.85
;        textcolor = 'black'
;     endif
;     if (i mod 3 eq 1) then begin 
;        textyp = 0.8
;        textcolor = 'magenta'
;     endif
;     if (i mod 3 eq 2) then begin 
;        textyp = 0.75
;        textcolor = 'blue'
;     endif
      djs_xyouts, !x.crange[0]+0.77*(!x.crange[1]-!x.crange[0]), $
          10.^(!y.crange[0]+(0.98-(12-(i mod 12))*0.021)*(!y.crange[1]-!y.crange[0])), $
          legends, charsize=1.1, charthick=charthick, color=zcolor[i]

      p = [dndz.n_star[i]/dndz.w_star[i], dndz.w_star[i]]
      ymodel = jhusdss_dndzdw_noz_func(xmodel, p)
      if (i mod 12 eq 0) then begin
         ymodel1=ymodel
         mcolor=zcolor[i]
;        legend, 'Exponential Fit at '+legends, color=zcolor[i], $
;           charsize=1.2, charthick=charthick, $
;           linestyle=2, thick=thick, box=0, /left, /top, pspacing=4.0, textcolor=zcolor[i]
      endif
      djs_oplot, xmodel, ymodel*(1.E-1)^((12- (i mod 12))*dexp), color=zcolor[i], thick=10, linestyle=2
;     djs_oplot, xmodel, ymodel1*(1.E-1)^((i mod 6)*dexp), color=mcolor, thick=8, linestyle=2
  endfor

   ;; inset
   xra=[0.2, 2.5]
   yra=[0.41, 0.85]
   x = dndz.zbin_median
   y = dndz.w_star
   xerr = dndz.zbin_sdev
   yerr = dndz.w_star_err
   xmodel=(findgen(100)+1.)/100.*6.5
   xtitle='z'
   title=''
;  ytitle=textoidl('(W_0^{\lambda2796})^* \AA')
   ytitle=textoidl('W^* (\AA)')

   pos = [0.26, 0.16, 0.44, 0.28]
   djs_plot, x, y, psym=5, xra=xra, yra=yra, symsize=1.0, $
        thick=thick, xthick=thick, ythick=thick, $
        xtitle=xtitle, ytitle=ytitle, $
        charsize=0.8, charthick=1.3, /nodata, /noerase, $
        pos=pos
   for i=0L, nmax-1L do begin
        oploterror, x[i], y[i], yerr[i], thick=thick, psym=5, $
        color=zcolor[i], errcolor=zcolor[i]
   endfor

   ymodel = dndz.w0/(1.+(xmodel/dndz.w0_bbeta)^dndz.w0_ggamma)*(1.+xmodel)^dndz.w0_alpha
   djs_oplot, xmodel, ymodel, color='red', thick=6, linestyle=0

;  zt = [2.51, 3.46, 4.78]
;  ztlow = [1.947, 3.150, 4.345]
;  ztup = [2.975, 3.805, 5.350]
;  wt = [0.935, 0.766, 0.700]
;  wterr = [0.150, 0.152, 0.180]

;  oploterror, zt, wt, zt-ztlow, wterr, /lobar, thick=thick, symsize=1.5, $
;     psym=4, color=djs_icolor('dark green'), errcolor=djs_icolor('dark green')
;  oploterror, zt, wt, ztup-zt, wterr, /hibar, thick=thick, symsize=1.5, $
;     psym=4, color=djs_icolor('dark green'), errcolor=djs_icolor('dark green')

;  items = ['Zhu & Menard this work', 'Matejek & Simcoe 2012']
;  colors = [djs_icolor('black'), djs_icolor('dark green')]
;  legend, items, psym=[5,4], colors=colors, number=1, thick=thick, $
;    box=0, /right, /bottom, charsize=1.3, charthick=charthick, textcolor=colors

k_end_print
loadct, 0


psfile = qapath+'/'+'dndzdw_z_one_'+string(nmfver, format='(I3.3)')+'.ps'
k_print, filename=psfile, axis_char_scale=1.4, xsize=12, ysize=10
  !p.multi=[0, 2, 1]
  !x.margin=0
  !y.margin=0
  z_nbin = n_elements(dndz.z)
  yra=[1.E-8, 8.]
  xmodel=findgen(100)/10.
  xtitle=textoidl('W_0^{\lambda2796} (\AA)')
  ytitle=textoidl('!9D!6^2N/!9D!6z!9D!6W^{\lambda2796}_0')

; zcolor = [djs_icolor('gray'), djs_icolor('black'), djs_icolor('pink'), djs_icolor('navy'), djs_icolor('orange'), $
  zcolor = [djs_icolor('black'), djs_icolor('pink'), djs_icolor('navy'), djs_icolor('orange'), $
            djs_icolor('blue'), djs_icolor('magenta red'), djs_icolor('dark green'), djs_icolor('magenta'), $
            djs_icolor('green'),  djs_icolor('brown'), djs_icolor('dark cyan'), djs_icolor('red')]

  nmax = (z_nbin < 12)
  for i=0, nmax-1L do begin
      ii = where(dndz.phi[*,i] gt 0.)
;     x = 10.^dndz.median_log10w[ii,i]
      x = dndz.median_w[ii,i]
      y = dndz.phi[ii,i]
      yerr = dndz.phi_poisson_err[ii,i]
      plotsym, 0, 1.5, /fill 
      if (i mod 6 eq 0) then begin
          if (i eq 0) then begin  
             xra=[0.00, 5.00]
             ytickformat_tmp='tick_exponent'
             ytitle_tmp = ytitle
          endif else begin
             xra=[0.01, 5.00]
             ytickformat_tmp='(A1)'
             ytitle_tmp = ''
          endelse
          djs_plot, x, y, psym=8, xra=xra, yra=yra, $
             thick=thick, xthick=thick, ythick=thick, $
             xtitle=xtitle, ytickformat=ytickformat_tmp, ytitle=ytitle_tmp, /ylog, $
;            xtickformat='(A1)', ytitle=ytitle, /ylog, $
             charsize=charsize, charthick=charthick, color=zcolor[i], /nodata
      endif 
      dexp = 1.0
;     djs_oplot, x[2:*], y[2:*]*(1E-1)^((i mod 6)*dexp), thick=thick, psym=8, color=zcolor[i]
      plotsym, 0, 1.5, /fill 
      oploterror, x[2:*], y[2:*]*(1E-1)^((i mod 6)*dexp), yerr[2:*]*(1E-1)^((i mod 6)*dexp), $
          thick=thick, psym=8, /nohat, errthick=2, $
          color=zcolor[i], errcolor=zcolor[i]
      plotsym, 0, 1.5
      oploterror, x[0:1], y[0:1]*(1E-1)^((i mod 6)*dexp), yerr[0:1]*(1E-1)^((i mod 6)*dexp), $
          thick=thick+4, psym=8, /nohat, errthick=2, $
          color=zcolor[i], errcolor=zcolor[i]

      z_min=dndz.z_min[i]
      z_max=dndz.z_max[i]

;     if (i eq 0) then z_min=0.37
      legends = string(z_min, format='(f4.2)')+'<z<'+ string(z_max, format='(f4.2)')
;     if (i mod 3 eq 0) then begin 
;        textyp = 0.85
;        textcolor = 'black'
;     endif
;     if (i mod 3 eq 1) then begin 
;        textyp = 0.8
;        textcolor = 'magenta'
;     endif
;     if (i mod 3 eq 2) then begin 
;        textyp = 0.75
;        textcolor = 'blue'
;     endif
      djs_xyouts, !x.crange[0]+0.7*(!x.crange[1]-!x.crange[0]), $
          10.^(!y.crange[0]+(0.90-(i mod 6)*0.03)*(!y.crange[1]-!y.crange[0])), $
          legends, charsize=1.3, charthick=charthick, color=zcolor[i]

      p = [dndz.n_star[i]/dndz.w_star[i], dndz.w_star[i]]
      ymodel = jhusdss_dndzdw_noz_func(xmodel, p)
      if (i mod 6 eq 0) then begin
         ymodel1=ymodel
         mcolor=zcolor[i]
;        legend, 'Exponential Fit at '+legends, color=zcolor[i], $
;           charsize=1.2, charthick=charthick, $
;           linestyle=2, thick=thick, box=0, /left, /top, pspacing=4.0, textcolor=zcolor[i]
      endif
      djs_oplot, xmodel, ymodel*(1.E-1)^((i mod 6)*dexp), color=zcolor[i], thick=10, linestyle=2
;     djs_oplot, xmodel, ymodel1*(1.E-1)^((i mod 6)*dexp), color=mcolor, thick=8, linestyle=2
  endfor
k_end_print


psfile = qapath+'/'+'dndzdw_z_new_'+string(nmfver, format='(I3.3)')+'.ps'
k_print, filename=psfile, axis_char_scale=1.4, xsize=10, ysize=10
  !p.multi=[0, 2, 2]
  !x.margin=0
  !y.margin=0
  z_nbin = n_elements(dndz.z)
  xra=[0., 6.9]
  yra=[1.E-6, 5.]
  xmodel=findgen(100)/10.
  xtitle=textoidl('W_0^{\lambda2796} (\AA)')
  ytitle=textoidl('!9D!6^2N/!9D!6z!9D!6W^{\lambda2796}_0')

  nmax = (z_nbin < 12)
  for i=0, nmax-1 do begin
      ii = where(dndz.phi[*,i] gt 0.)
;     x = 10.^dndz.median_log10w[ii,i]
      x = dndz.median_w[ii,i]
      y = dndz.phi[ii,i]
      yerr = dndz.phi_poisson_err[ii,i]
      plotsym, 0, 1.5, /fill 
      if (i eq 0) then begin
          djs_plot, x, y, psym=8, xra=xra, yra=yra, $
             thick=thick, xthick=thick, ythick=thick, $
             xtickformat='(A1)', ytickformat='tick_exponent', ytitle=ytitle, /ylog, $
;            xtickformat='(A1)', ytitle=ytitle, /ylog, $
             charsize=charsize, charthick=charthick
      endif 
      if (i eq 6) then begin
          djs_plot, x, y, psym=8, xra=xra, yra=yra, $
             thick=thick, xthick=thick, ythick=thick, $
             xtitle=xtitle, ytickformat='tick_exponent', ytitle=ytitle, /ylog, $
;            xtickformat='(A1)', ytitle=ytitle, /ylog, $
             charsize=charsize, charthick=charthick
      endif 

      if (i eq 3 or i eq 9) then begin
          djs_plot, x, y, psym=8, xra=xra, yra=yra, $
             thick=thick, xthick=thick, ythick=thick, $
             xtickformat='(A1)', ytickformat='(A1)', /ylog, $
             charsize=charsize, charthick=charthick
      endif 

;     if (i eq 1 or i eq 4 or i eq 7 or i eq 10) then begin
;         djs_oplot, x, y*1.E-1, psym=4, xra=xra, yra=yra, $
;             thick=thick, xthick=thick, ythick=thick, $
;             xtickformat='(A1)', ytickformat='(A1)', /ylog, $
;             charsize=charsize, charthick=charthick, color='magenta'
;     endif
;     if (i eq 2 or i eq 5 or i eq 8 or i eq 11) then begin
;         djs_oplot, x, y*1.E-2, psym=4, xra=xra, yra=yra, $
;             thick=thick, xthick=thick, ythick=thick, $
;             xtickformat='(A1)', ytickformat='(A1)', /ylog, $
;             charsize=charsize, charthick=charthick, color='blue'
;     endif

      if (i eq 0 or i eq 3 or i eq 6 or i eq 9) then $
         oploterror, x, y, yerr, thick=thick, psym=8, /nohat, errthick=2
      if (i eq 1 or i eq 4 or i eq 7 or i eq 10) then $
         oploterror, x, y*1.E-1, yerr*1.E-1, thick=thick, psym=8, /nohat, errthick=3, $
            color=djs_icolor('magenta'),  errcolor=djs_icolor('magenta')
      if (i eq 2 or i eq 5 or i eq 8 or i eq 11) then $
         oploterror, x, y*1.E-2, yerr*1.E-2, thick=thick, psym=8, /nohat, errthick=3, $
            color=djs_icolor('blue'), errcolor=djs_icolor('blue')

      z_min=dndz.z_min[i]
      z_max=dndz.z_max[i]
;     if (i eq 0) then z_min=0.36
      legends = string(z_min, format='(f4.2)')+'<z<'+ string(z_max, format='(f4.2)')
      if (i mod 3 eq 0) then begin 
         textyp = 0.85
         textcolor = 'black'
      endif
      if (i mod 3 eq 1) then begin 
         textyp = 0.8
         textcolor = 'magenta'
      endif
      if (i mod 3 eq 2) then begin 
         textyp = 0.75
         textcolor = 'blue'
      endif
      djs_xyouts, !x.crange[0]+0.6*(!x.crange[1]-!x.crange[0]), $
          10.^(!y.crange[0]+textyp*(!y.crange[1]-!y.crange[0])), $
          legends, charsize=1.3, charthick=charthick, color=textcolor

   ;  p = [dndz.n_star_strong[i]/dndz.w_star_strong[i], dndz.w_star_strong[i], $
   ;       dndz.n_star_weak[i]/dndz.w_star_weak[i], dndz.w_star_weak[i]]
   ;  ymodel = jhusdss_dndzdw_noz_func2(xmodel, p)
      p = [dndz.n_star[i]/dndz.w_star[i], dndz.w_star[i]]
      ymodel = jhusdss_dndzdw_noz_func(xmodel, p)
      if (i mod 3 eq 0) then djs_oplot, xmodel, ymodel, color='black', thick=8, linestyle=2
      if (i mod 3 eq 1) then djs_oplot, xmodel, ymodel*1.E-1, color='magenta', thick=8
      if (i mod 3 eq 2) then djs_oplot, xmodel, ymodel*1.E-2, color='blue', thick=8
      if (i mod 3 eq 0) then ymodel1=ymodel
;     if (i eq 0) then  ymodel1 = ymodel $
;      else djs_oplot, xmodel, ymodel1, color='blue', linestyle=2, thick=4
      if (i mod 3 eq 1) then djs_oplot, xmodel, ymodel1*1.E-1, color='black', thick=8, linestyle=2
      if (i mod 3 eq 2) then djs_oplot, xmodel, ymodel1*1.E-2, color='black', thick=8, linestyle=2

;     if (i eq 0 or i eq 3) then $
;        djs_axis, yaxis=0, ytitle=ytitle, thick=thick, $
;          charsize=charsize, charthick=charthick

      if ((i eq 6) or (i eq 9)) then $
         djs_axis, xaxis=0, xtitle=xtitle, thick=thick, $
           charsize=charsize, charthick=charthick
   ;  if (i eq 1) then $
   ;     djs_axis, xaxis=1, xtitle=ytitle, thick=thick, $
   ;       charsize=charsize, charthick=charthick, xtickformat='(A1)'
  endfor
k_end_print


psfile = qapath+'/'+'dndzdw_z_'+string(nmfver, format='(I3.3)')+'.ps'
k_print, filename=psfile, axis_char_scale=2., xsize=12, ysize=8
  !p.multi=[0, 3, 2]
  !x.margin=0
  !y.margin=0
  z_nbin = n_elements(dndz.z)
  xra=[0., 5.9]
  yra=[2.E-5, 9.]
  xmodel=findgen(100)/10.
  xtitle=textoidl('W_0^{\lambda2796} (\AA)')
  ytitle=textoidl('!9D!6^2N/!9D!6z!9D!6W^{\lambda2796}_0')

  nmax = (z_nbin < 12)
  for i=0, nmax-1 do begin
      ii = where(dndz.phi[*,i] gt 0.)
;     x = 10.^dndz.median_log10w[ii,i]
      x = dndz.median_w[ii,i]
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
      oploterror, x, y, yerr, thick=thick, psym=4, symsize=1.5
      z_min=dndz.z_min[i]
      z_max=dndz.z_max[i]
;     if (i eq 0) then z_min=0.36
      legends = string(z_min, format='(f4.2)')+'<z<'+ string(z_max, format='(f4.2)')
      djs_xyouts, !x.crange[0]+0.6*(!x.crange[1]-!x.crange[0]), $
                  10.^(!y.crange[0]+0.85*(!y.crange[1]-!y.crange[0])), $
                  legends, charsize=1.3, charthick=charthick

   ;  p = [dndz.n_star_strong[i]/dndz.w_star_strong[i], dndz.w_star_strong[i], $
   ;       dndz.n_star_weak[i]/dndz.w_star_weak[i], dndz.w_star_weak[i]]
   ;  ymodel = jhusdss_dndzdw_noz_func2(xmodel, p)
      p = [dndz.n_star[i]/dndz.w_star[i], dndz.w_star[i]]
      ymodel = jhusdss_dndzdw_noz_func(xmodel, p)
      djs_oplot, xmodel, ymodel, color='red', thick=4
      if (i mod 3 eq 0) then  ymodel1 = ymodel $
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

thick=8
psfile = qapath+'/'+'Cum_dndzdw_z_one_'+string(nmfver, format='(I3.3)')+'.ps'
k_print, filename=psfile, axis_char_scale=1.4, xsize=16, ysize=8

  pos = [0.1, 0.1, 0.48, 0.9]
; w_min = [0.3, 0.6, 1.0, 1.5, 2.0, 3.0, 4.0]
  z_nbin = n_elements(dndz.z)
  xra=[0.3, 3.2]
  xmodel=(findgen(100)+1.)/100.*2.4
  xtitle=textoidl('z')
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
  xtitle=textoidl('z')
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


thick=8
psfile = qapath+'/'+'Cum_dndzdw_z_'+string(nmfver, format='(I3.3)')+'.ps'
k_print, filename=psfile, axis_char_scale=1.4, xsize=9, ysize=9
; !p.multi=[0, 3, 2]
; !x.margin=0
; !y.margin=0
  z_nbin = n_elements(dndz.z)
  xra=[0.3, 3.2]
  xmodel=(findgen(100)+1.)/100.*2.4
  xtitle=textoidl('z')
  ytitle=textoidl('dN/dz (W_0^{\lambda2796}>W_{min})')
  yra=[2.E-4, 2.]

  wvector = make_vector(0.6, 4.2, 6)
  w_min = [0.2, wvector.bound_min]
  w_max = [0.6, wvector.bound_max]
; w_min = [0.3, 0.6, 1.0, 1.5, 2.0, 3.0, 4.0]

  psyms = [1, 8, 8, 8, 8, 8, 8]
  colors = ['grey', 'red', 'magenta', 'brown', 'dark green', 'blue', 'black']
  ypos = [1.0, 0.6, 0.3, 0.15, 0.09, 0.02, 0.006]

  for i=0, n_elements(w_min)-1 do begin
 
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
             charsize=charsize, charthick=charthick, /nodata

      plotsym, 8, 0.8, /fill 
;     oploterror, x[1:n_elements(x)-2], y[1:n_elements(x)-2], yerr[1:n_elements(x)-2], thick=thick, psym=psyms[i], $
      oploterror, x, y, yerr, thick=thick, psym=psyms[i], $
          color=djs_icolor(colors[i]), errcolor=djs_icolor(colors[i]), symsize=2
;     plotsym, 8, 0.8
;     oploterror, x[0], y[0], yerr[0], thick=3, psym=psyms[i], $
;         color=djs_icolor(colors[i]), errcolor=djs_icolor(colors[i]), symsize=2
;     oploterror, x[n_elements(x)-1], y[n_elements(x)-1], yerr[n_elements(x)-1], thick=3, psym=psyms[i], $
;         color=djs_icolor(colors[i]), errcolor=djs_icolor(colors[i]), symsize=2

      ymodel = dndz.f0*(1.+xmodel)^(dndz.f0_alpha+dndz.w0_alpha)/(1.+(xmodel/dndz.f0_bbeta)^dndz.f0_ggamma)*dndz.w0/(1.+(xmodel/dndz.w0_bbeta)^dndz.w0_ggamma)*exp(-w_min_tmp/dndz.w0*(1.+(xmodel/dndz.w0_bbeta)^dndz.w0_ggamma)/(1.+xmodel)^dndz.w0_alpha)
      if i eq 0 then lstyle=2 else lstyle=0
      if i gt 0 then djs_oplot, xmodel, ymodel, thick=thick, color=djs_icolor(colors[i]), linestyle=lstyle

;     z_min=dndz.z_min[i]
;     z_max=dndz.z_max[i]
;     if (i eq 0) then z_min=0.36
      legends = 'W_0>'+string(w_min_tmp, format='(f3.1)')+' \AA'

      if (i eq 0) then djs_xyouts, !x.crange[0]+0.75*(!x.crange[1]-!x.crange[0]), ymodel[n_elements(xmodel)-1], $;ypos[i], $
                  legends, charsize=1.4, charthick=charthick, color=colors[i]
      if (i gt 0) then djs_xyouts, !x.crange[0]+0.75*(!x.crange[1]-!x.crange[0]), ymodel[n_elements(xmodel)-1], $;ypos[i], $
                  legends, charsize=1.4, charthick=charthick, color=colors[i]
  endfor
k_end_print

psfile = qapath+'/'+'Cum_dndzdw_z_diff_'+string(nmfver, format='(I3.3)')+'.ps'
k_print, filename=psfile, axis_char_scale=1.4, xsize=9, ysize=9
; !p.multi=[0, 3, 2]
; !x.margin=0
; !y.margin=0
  z_nbin = n_elements(dndz.z)
  xra=[0.3, 3.2]
  xmodel=(findgen(100)+1.)/100.*2.4
  xtitle=textoidl('z')
  ytitle=textoidl('dN/dz (W_{min}<W_0^{\lambda2796}<W_{max})')
  yra=[2.E-4, 0.8]

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
  colors = ['grey', 'red', 'magenta', 'brown', 'dark green', 'blue', 'black']
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
             charsize=charsize, charthick=charthick, /nodata

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
         legends = string(w_min_tmp, format='(f3.1)')+'<W<'+string(w_max_tmp,format='(f3.1)')+' \AA'
;     if (i eq 0) then djs_xyouts, !x.crange[0]+0.75*(!x.crange[1]-!x.crange[0]), (ymodel-ymodelmax)[n_elements(xmodel)-1]+0.1, $;ypos[i], $
;                 legends, charsize=1.4, charthick=charthick, color=colors[i]
      if (i gt 0) then djs_xyouts, !x.crange[0]+0.75*(!x.crange[1]-!x.crange[0]), (ymodel-ymodelmax)[n_elements(xmodel)-1], $;ypos[i], $
                  legends, charsize=1.4, charthick=charthick, color=colors[i]
  endfor
k_end_print


psfile = qapath+'/'+'Cum_dndzdw_z_w1.0_'+string(nmfver, format='(I3.3)')+'.ps'
k_print, filename=psfile, axis_char_scale=1.2, xsize=8, ysize=8
  !p.multi=[0, 1, 2]
  !x.margin=0
  !y.margin=0
  z_nbin = n_elements(dndz.z)
  xra=[0.1, 5.9]
; xra=[1.0, 10.]
  xmodel=(findgen(100)+1.)/100.*6.5
  xtitle=textoidl('z')
  ytitle=textoidl('dN/dz')

; w_min = [0.3, 0.6, 1.0]
  w_min = [0.6, 1.0]
  w_max = [1.0, 100.0]
  yra = [0.01, 1.9]

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

  for i=0, n_elements(w_min)-1 do begin
 
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

      djs_plot, x, y*1, psym=5, symsize=1.5, xra=xra, yra=yra, $
          thick=thick, xthick=thick, ythick=thick, $
          xtickformat='(A1)', ytitle='', /ylog,  $
          charsize=charsize, charthick=charthick, color='blue'
      oploterror, x, y*1, yerr, thick=thick, psym=5, symsize=1.5, color=djs_icolor('blue'), $
          errcolor=djs_icolor('blue')
      if (i eq 1) then djs_axis, xaxis=0, xtitle=xtitle, charsize=charsize, charthick=charthick
      djs_axis, yaxis=0, ytitle=ytitle, charsize=charsize, charthick=charthick

;     ymodel = dndz.f0*(1.+xmodel)^4/(1.+(xmodel/dndz.f0_bbeta)^dndz.f0_ggamma)*dndz.w0/(1.+(xmodel/dndz.w0_bbeta)^dndz.w0_ggamma)*exp(-w_min_tmp/dndz.w0*(1.+(xmodel/dndz.w0_bbeta)^dndz.w0_ggamma)/(1.+xmodel)) $
;            - dndz.f0*(1.+xmodel)^4/(1.+(xmodel/dndz.f0_bbeta)^dndz.f0_ggamma)*dndz.w0/(1.+(xmodel/dndz.w0_bbeta)^dndz.w0_ggamma)*exp(-w_max_tmp/dndz.w0*(1.+(xmodel/dndz.w0_bbeta)^dndz.w0_ggamma)/(1.+xmodel)) 
      ymodel = dndz.f0*(1.+xmodel)^(dndz.f0_alpha+dndz.w0_alpha)/(1.+(xmodel/dndz.f0_bbeta)^dndz.f0_ggamma)*dndz.w0/(1.+(xmodel/dndz.w0_bbeta)^dndz.w0_ggamma)*exp(-w_min_tmp/dndz.w0*(1.+(xmodel/dndz.w0_bbeta)^dndz.w0_ggamma)/(1.+xmodel)^dndz.w0_alpha)
      ymodelmax = dndz.f0*(1.+xmodel)^(dndz.f0_alpha+dndz.w0_alpha)/(1.+(xmodel/dndz.f0_bbeta)^dndz.f0_ggamma)*dndz.w0/(1.+(xmodel/dndz.w0_bbeta)^dndz.w0_ggamma)*exp(-w_max_tmp/dndz.w0*(1.+(xmodel/dndz.w0_bbeta)^dndz.w0_ggamma)/(1.+xmodel)^dndz.w0_alpha)
      djs_oplot, xmodel, (ymodel-ymodelmax)*1., thick=thick, color='blue',linestyle=0

      if (i ne 1) then legends = string(w_min[i],format='(f3.1)')+'<W_0(Mg II)<'+string(w_max[i], format='(f3.1)')+' \AA' $
                  else legends = 'W_0(Mg II)>'+string(w_min[i], format='(f3.1)')+' \AA' 
      djs_xyouts, !x.crange[0]+0.60*(!x.crange[1]-!x.crange[0]), 0.8, $
                  legends, charsize=1.4, charthick=charthick, color='black'
 
      oploterror, zt[*,i], yt[*,i], zt[*,i]-ztlow[*,i], yterr[*,i], /lobar, thick=thick, symsize=1.5, $
        psym=4, color=djs_icolor('dark green'), errcolor=djs_icolor('dark green')
;       psym=4, color=djs_icolor('blue'), errcolor=djs_icolor('blue')
      oploterror, zt[*,i], yt[*,i], ztup[*,i]-zt[*,i], yterr[*,i], /hibar, thick=thick, symsize=1.5, $
        psym=4, color=djs_icolor('dark green'), errcolor=djs_icolor('dark green')
;       psym=4, color=djs_icolor('blue'), errcolor=djs_icolor('blue')

      if (i eq 1) then begin
      asfh = 0.017
      bsfh = 0.13
      csfh = 3.3
      dsfh = 5.3
      ymodel_sfh = (asfh+bsfh*xmodel)/(1.+(xmodel/csfh)^dsfh)
;     fudge = 1.65
      fudge = 0.65
      djs_oplot, xmodel, fudge*ymodel_sfh*jhusdss_dxdz(xmodel), color='black', linestyle=2, thick=thick

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
      if (i eq 0) then begin
         items = ['Zhu & M!3'+string("351B)+'!Xnard (this work)', 'Matejek & Simcoe (2012)']
         colors = [djs_icolor('blue'), djs_icolor('dark green')]
;        colors = [djs_icolor('red'), djs_icolor('blue')]
         legend, items, psym=[5,4], colors=colors, number=1, thick=thick, $
           box=0, /right, /bottom, charsize=1.3, charthick=charthick, textcolors=colors
      endif else begin
         items = ['dN/dz Fit by Zhu & M!3'+string("351B)+'!Xnard (this work)', '(Scaled) SFH Fit by Hopkins & Beacom (2006)']
         colors = [djs_icolor('blue'), djs_icolor('black')]
         legend, items, colors=colors, number=1, thick=thick, $
           box=0, /right, /bottom, charsize=1.3, charthick=charthick, textcolors=colors, $
           linestyle=[0, 2], pspacing=1.5
      endelse

  endfor

;     items = ['Zhu & Menard this work', 'Zhu & Menard Fit (0.6<W_0^{\lambda2796}<5.0 \AA)', 'Matejek & Simcoe 2012']
k_end_print

psfile = qapath+'/'+'Wstar_dndzdw_z_'+string(nmfver, format='(I3.3)')+'.ps'
k_print, filename=psfile, axis_char_scale=1.4, xsize=8, ysize=8
   xra=[0., 6.]
   yra=[0.2, 1.2]
   x = dndz.zbin_median
   y = dndz.w_star
   xerr = dndz.zbin_sdev
   yerr = dndz.w_star_err
   xmodel=(findgen(100)+1.)/100.*6.5
   xtitle='z' 
   ytitle=textoidl('(W_0^{\lambda2796})^* \AA')

   djs_plot, x, y, psym=5, xra=xra, yra=yra, symsize=1.5, $
        thick=thick, xthick=thick, ythick=thick, $
        xtitle=xtitle, ytitle=ytitle, $
        charsize=charsize, charthick=charthick, /nodata
   oploterror, x, y, yerr, thick=thick, psym=5, color=djs_icolor('blue'), errcolor=djs_icolor('blue')
       

   ymodel = dndz.w0/(1.+(xmodel/dndz.w0_bbeta)^dndz.w0_ggamma)*(1.+xmodel)^dndz.w0_alpha
   djs_oplot, xmodel, ymodel, color='blue', thick=4, linestyle=2

   zt = [2.51, 3.46, 4.78]
   ztlow = [1.947, 3.150, 4.345]
   ztup = [2.975, 3.805, 5.350]
   wt = [0.935, 0.766, 0.700]
   wterr = [0.150, 0.152, 0.180]

   oploterror, zt, wt, zt-ztlow, wterr, /lobar, thick=thick, symsize=1.5, $
      psym=4, color=djs_icolor('dark green'), errcolor=djs_icolor('dark green')
   oploterror, zt, wt, ztup-zt, wterr, /hibar, thick=thick, symsize=1.5, $
      psym=4, color=djs_icolor('dark green'), errcolor=djs_icolor('dark green')

   items = ['Zhu & Menard 2012', 'Matejek & Simcoe 2012']
   colors = [djs_icolor('black'), djs_icolor('dark green')]
   legend, items, psym=[5,4], colors=colors, number=1, thick=thick, $
     box=0, /right, /bottom, charsize=1.3, charthick=charthick, textcolor=colors
k_end_print

for ifigure=0, 1 do begin

psfile = qapath+'/'+'Cum_dndzdw_z_w1.0_'+string(nmfver, format='(I3.3)')+'_'+string(ifigure, format='(i1.1)')+'.ps'

k_print, filename=psfile, axis_char_scale=1.2, xsize=10, ysize=7
  !p.multi=[0, 1, 1]
  !x.margin=0
  !y.margin=0
  z_nbin = n_elements(dndz.z)
  xra=[0.1, 5.9]
; xra=[1.0, 10.]
  xmodel=(findgen(100)+1.)/100.*7.5
  xtitle=textoidl('z')
  ytitle=textoidl('dN/dz')

  ms_color='light blue'
; w_min = [0.3, 0.6, 1.0]
; w_min = [0.6, 1.0]
; w_max = [1.0, 100.0]
  w_min = [0.6, 1.0]
  w_max = [1.0, 100.0]
  yra = [0.04, 1.9]

; dx/dz = (1+x)^2/sqrt(0.7+0.3*(1+x)^3)

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

  for i=1, n_elements(w_min)-1 do begin
 
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
;     yerr = dndz.cum_phi_poisson_err[iw_min,ii]

      djs_plot, x, y, psym=5, symsize=2.0, xra=xra, yra=yra, $
          thick=thick, xthick=thick, ythick=thick, $
          xtickformat='(A1)', ytitle='', /ylog,  $
          charsize=charsize, charthick=charthick, color='blue'
      oploterror, x, y, yerr, thick=thick, psym=5, symsize=2.0, color=djs_icolor('blue'), $
          errcolor=djs_icolor('blue')
      if (i eq 1) then djs_axis, xaxis=0, xtitle=xtitle, charsize=charsize, charthick=charthick
      djs_axis, yaxis=0, ytitle=ytitle, charsize=charsize, charthick=charthick

;     ymodel = dndz.f0*(1.+xmodel)^4/(1.+(xmodel/dndz.f0_bbeta)^dndz.f0_ggamma)*dndz.w0/(1.+(xmodel/dndz.w0_bbeta)^dndz.w0_ggamma)*exp(-w_min_tmp/dndz.w0*(1.+(xmodel/dndz.w0_bbeta)^dndz.w0_ggamma)/(1.+xmodel)) $
;            - dndz.f0*(1.+xmodel)^4/(1.+(xmodel/dndz.f0_bbeta)^dndz.f0_ggamma)*dndz.w0/(1.+(xmodel/dndz.w0_bbeta)^dndz.w0_ggamma)*exp(-w_max_tmp/dndz.w0*(1.+(xmodel/dndz.w0_bbeta)^dndz.w0_ggamma)/(1.+xmodel)) 
      ymodel = dndz.f0*(1.+xmodel)^(dndz.f0_alpha+dndz.w0_alpha)/(1.+(xmodel/dndz.f0_bbeta)^dndz.f0_ggamma)*dndz.w0/(1.+(xmodel/dndz.w0_bbeta)^dndz.w0_ggamma)*exp(-w_min_tmp/dndz.w0*(1.+(xmodel/dndz.w0_bbeta)^dndz.w0_ggamma)/(1.+xmodel)^dndz.w0_alpha)
      ymodelmax = dndz.f0*(1.+xmodel)^(dndz.f0_alpha+dndz.w0_alpha)/(1.+(xmodel/dndz.f0_bbeta)^dndz.f0_ggamma)*dndz.w0/(1.+(xmodel/dndz.w0_bbeta)^dndz.w0_ggamma)*exp(-w_max_tmp/dndz.w0*(1.+(xmodel/dndz.w0_bbeta)^dndz.w0_ggamma)/(1.+xmodel)^dndz.w0_alpha)
;     djs_oplot, xmodel, ymodel-ymodelmax, thick=thick, color='blue',linestyle=0

      if (i ne 1) then legends = string(w_min[i],format='(f3.1)')+'<W_0(Mg II)<'+string(w_max[i], format='(f3.1)')+' \AA' $
                  else legends = 'W_0(Mg II)>'+string(w_min[i], format='(f3.1)')+' \AA' 
      djs_xyouts, !x.crange[0]+0.80*(!x.crange[1]-!x.crange[0]), 1.2, $
                  legends, charsize=1.4, charthick=charthick, color='black'
 
      oploterror, zt[*,i], yt[*,i], zt[*,i]-ztlow[*,i], yterr[*,i], /lobar, thick=thick, symsize=2.0, $
        psym=4, color=djs_icolor(ms_color), errcolor=djs_icolor(ms_color)
;       psym=4, color=djs_icolor('blue'), errcolor=djs_icolor('blue')
      oploterror, zt[*,i], yt[*,i], ztup[*,i]-zt[*,i], yterr[*,i], /hibar, thick=thick, symsize=2.0, $
        psym=4, color=djs_icolor(ms_color), errcolor=djs_icolor(ms_color)
;       psym=4, color=djs_icolor('blue'), errcolor=djs_icolor('blue')

      if (ifigure eq 1) then begin
      asfh = 0.017
      bsfh = 0.13
      csfh = 3.3
      dsfh = 5.3
      ymodel_sfh = (asfh+bsfh*xmodel)/(1.+(xmodel/csfh)^dsfh)
      fudge = 1.65*1.05
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
      fudge = 2.0/3.0
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
;     if (i eq 1) then begin
         items = ['Zhu & Menard (this work)', 'Matejek & Simcoe (2012)']
         colors = [djs_icolor('blue'), djs_icolor(ms_color)]
         colors = [djs_icolor('blue'), djs_icolor(ms_color)]
         legend, items, psym=[5,4], colors=colors, number=1, thick=thick, $
           box=0, /left, /top, charsize=1.3, charthick=charthick, textcolors=colors
;     endif else begin
;        items = ['dN/dz Fit Zhu & Menard this work', '(Scaled) SFH Fit Hopkins & Beacom 2006']
;        colors = [djs_icolor('blue'), djs_icolor('black')]
;        legend, items, colors=colors, number=1, thick=thick, $
;          box=0, /right, /bottom, charsize=1.3, charthick=charthick, textcolors=colors, $
;          linestyle=[0, 2], pspacing=1.5
;     endelse

;     djs_xyouts, !x.crange[0]+0.2*(!x.crange[1]-!x.crange[0]), $
;                 !y.crange[0]+0.3*(!y.crange[1]-!y.crange[0]), $
;                 'Data Zhu & Menard this work', charsize=1.3, charthick=charthick,$
;                 color='blue'
;     djs_oplot, !x.crange[0]+[0.8, 0.9]*(!x.crange[1]-!x.crange[0]), $
;                 !y.crange[0]+[0.31, 0.31]*(!y.crange[1]-!y.crange[0]), $
;                 psym=3, thick=thick, color='blue'

;     djs_xyouts, !x.crange[0]+0.4*(!x.crange[1]-!x.crange[0]), $
;                 0.04,$
;                 'dN/dz Fit (Zhu & Menard 2012)', charsize=1.3, charthick=charthick,$
;                 color='blue'
;     djs_oplot, !x.crange[0]+[0.85, 0.944]*(!x.crange[1]-!x.crange[0]), $
;                 [0.04, 0.04]+0.001, $
;                 thick=thick, color='blue'

;     djs_xyouts, !x.crange[0]+0.2*(!x.crange[1]-!x.crange[0]), $
;                 !y.crange[0]+0.2*(!y.crange[1]-!y.crange[0]), $
;                 'Data Matejek & Simcoe 2012', charsize=1.3, charthick=charthick,$
;                 color='dark green'
;     djs_oplot, !x.crange[0]+[0.8, 0.9]*(!x.crange[1]-!x.crange[0]), $
;                 !y.crange[0]+[0.21, 0.21]*(!y.crange[1]-!y.crange[0]), $
;                 psym=4, thick=thick, color='dark green'

      if (ifigure eq 1) then begin
      djs_xyouts, !x.crange[0]+0.15*(!x.crange[1]-!x.crange[0]), $
                  0.052, $
;                 '(Scaled) SFH Fit Hopkins & Beacom 2006', charsize=1.3, charthick=charthick,$
                  'Star Formation History (Hopkins & Beacom 2006)', charsize=1.3, charthick=charthick,$
                  color='red'
      djs_oplot, !x.crange[0]+[0.85, 0.95]*(!x.crange[1]-!x.crange[0]), $
                 [0.052, 0.052]+0.001, $
                 linestyle=2, thick=thick, color='red'
      endif
 
  endfor

;     items = ['Zhu & Menard this work', 'Zhu & Menard Fit (0.6<W_0^{\lambda2796}<5.0 \AA)', 'Matejek & Simcoe 2012']
k_end_print

endfor

stop
end
