pro jhusdss_qaplot_lowz_caii, nmfver

if (n_elements(nmfver) eq 0) then message, 'nmfver required'

if (keyword_set(boss)) then begin
   lowzpath=jhusdss_get_path(/nmfqso)+'/'+string(nmfver, format='(I3.3)')+'/Lowz_BOSS'
endif else begin
   lowzpath=jhusdss_get_path(/nmfqso)+'/'+string(nmfver, format='(I3.3)')+'/Lowz'
endelse

infile = lowzpath+'/'+jhusdss_lowz_composite_filename(nmfver)

lowz = mrdfits(infile, 1)

;; init
thick=8
charsize=1.4
charthick=3
;; path
path=jhusdss_get_path(/nmfqso)+'/'+string(nmfver, format='(I3.3)')
qapath = path+'/QAplots'
if (jhusdss_direxist(qapath) eq 0) then spawn, 'mkdir -p '+qapath

all_lflux = fltarr(4)
all_lflux_err = fltarr(4)

all_sigma = fltarr(4)
all_sigma_err = fltarr(4)

psfile = qapath+'/'+'lowz_caii_'+string(nmfver, format='(I3.3)')+'.ps'
k_print, filename=psfile, axis_char_scale=1.8, xsize=12, ysize=12
   !p.multi=[0,1,4]
   !y.margin=0
   xtitle = 'Absorber-Frame \lambda (\AA)'
   xra = [3850, 4050]
   ierror = where(lowz.wave gt 3950.-300. and lowz.wave lt 3950.+300., nerror)
   iuse = where(lowz.wave gt 3850. and lowz.wave lt 4050., nuse)

   yra = [1.-0.08, 1+0.08]
   y = lowz.fluxmean[*,0]
;  y = lowz.fluxmean[*,0]
   terror = sqrt((moment(y[ierror]))[1])
   print, 'typical error= ', terror

   in_slope = 0.
   in_intercept = median(1.-y[ierror])
   in_center = 3934.78
   in_separation = 3969.59-3934.78
   in_lflux = 0.2
   in_ratio = 0.5
   in_sigma = 2.0

   jhusdss_lowz_doublet_fit, lowz.wave[iuse], 1.-y[iuse], fltarr(nuse)+1./terror^2, $
       in_slope, in_intercept, in_center, in_separation, in_lflux, in_ratio, in_sigma, $
       slope=slope, intercept=intercept, center=center, separation=separation, $
       lflux=lflux, ratio=ratio, sigma=sigma, $
       err_slope=err_slope, err_intercept=err_intercept, err_center=err_center, $
       err_separation=err_separation, err_lflux=err_lflux, err_ratio=err_ratio, err_sigma=err_sigma, $
       maxwidth=maxwidth

   p = [slope, intercept, center, separation, lflux, ratio, sigma]
   perror = [err_slope, err_intercept, err_center, err_separation, err_lflux, err_ratio, err_sigma]

   yfit = jhusdss_lowz_doublet_func(lowz.wave, p)
   print, p
   print, perror

   all_lflux[0] = lflux
   all_lflux_err[0] = err_lflux
   all_sigma[0] = sigma
   all_sigma_err[0] = err_sigma

   clevel = slope*(lowz.wave-center)+intercept
   djs_plot, lowz.wave, y+clevel, $
       xra=xra, yra=yra, thick=thick, xthick=thick, ythick=thick, $
       xtickformat='(A1)', charsize=charsize, charthick=charthick
   djs_oplot, lowz.wave, 1.-yfit+clevel, $
       thick=thick, color='red'

   djs_xyouts, !x.crange[0]+0.75*(!x.crange[1]-!x.crange[0]), $
               !y.crange[0]+0.90*(!y.crange[1]-!y.crange[0]), $
               ' b<15 kpc ', charsize=charsize, charthick=charthick
   djs_xyouts, !x.crange[0]+0.75*(!x.crange[1]-!x.crange[0]), $
               !y.crange[0]+0.80*(!y.crange[1]-!y.crange[0]), $
               '<b>='+string(lowz.rp[0],format='(i3)')+' kpc', charsize=charsize, charthick=charthick
   djs_xyouts, !x.crange[0]+0.05*(!x.crange[1]-!x.crange[0]), $
               !y.crange[0]+0.85*(!y.crange[1]-!y.crange[0]), $
               'W^{\lambda3934}='+string(lflux,format='(f5.3)')+'\pm'+string(err_lflux,format='(f5.3)')+' \AA', $
               charsize=1.6, charthick=charthick

   yra = [1.-0.08*0.25, 1.+0.08*0.25]
   y = lowz.fluxmean[*,1]
;  y = lowz.fluxmean[*,1]
   terror = sqrt((moment(y[ierror]))[1])
   print, 'typical error= ', terror

   in_intercept = median(1.-y[ierror])
   in_lflux = 0.05

   jhusdss_lowz_doublet_fit, lowz.wave[iuse], 1.-y[iuse], fltarr(nuse)+1./terror^2, $
       in_slope, in_intercept, in_center, in_separation, in_lflux, in_ratio, in_sigma, $
       slope=slope, intercept=intercept, center=center, separation=separation, $
       lflux=lflux, ratio=ratio, sigma=sigma, $
       err_slope=err_slope, err_intercept=err_intercept, err_center=err_center, $
       err_separation=err_separation, err_lflux=err_lflux, err_ratio=err_ratio, err_sigma=err_sigma, $
       maxwidth=maxwidth

   p = [slope, intercept, center, separation, lflux, ratio, sigma]
   perror = [err_slope, err_intercept, err_center, err_separation, err_lflux, err_ratio, err_sigma]

   yfit = jhusdss_lowz_doublet_func(lowz.wave, p)
   print, p
   print, perror

   all_lflux[1] = lflux
   all_lflux_err[1] = err_lflux
   all_sigma[1] = sigma
   all_sigma_err[1] = err_sigma

   clevel = slope*(lowz.wave-center)+intercept
   djs_plot, lowz.wave, y+clevel, $
       xra=xra, yra=yra, thick=thick, xthick=thick, ythick=thick, $
       xtickformat='(A1)', yminor=5 , charsize=charsize, charthick=charthick;ytickinterval=0.01
   djs_oplot, lowz.wave, 1.-yfit+clevel, $
       thick=thick, color='red'

   djs_xyouts, !x.crange[0]+0.75*(!x.crange[1]-!x.crange[0]), $
               !y.crange[0]+0.90*(!y.crange[1]-!y.crange[0]), $
               '15<b<45 kpc ', charsize=charsize, charthick=charthick
   djs_xyouts, !x.crange[0]+0.75*(!x.crange[1]-!x.crange[0]), $
               !y.crange[0]+0.80*(!y.crange[1]-!y.crange[0]), $
               '<b>='+string(lowz.rp[1],format='(i3)')+' kpc', charsize=charsize, charthick=charthick
   djs_xyouts, !x.crange[0]+0.05*(!x.crange[1]-!x.crange[0]), $
               !y.crange[0]+0.85*(!y.crange[1]-!y.crange[0]), $
               'W^{\lambda3934}='+string(lflux,format='(f5.3)')+'\pm'+string(err_lflux,format='(f5.3)')+' \AA', $
               charsize=1.6, charthick=charthick

   yra = [1.-0.08*0.25^2, 1.+0.08*0.25^2]
   y = lowz.fluxmean[*,2]
;  y = lowz.fluxmean[*,2]
   terror = sqrt((moment(y[ierror]))[1])
   print, 'typical error= ', terror

   in_intercept = median(1.-y[ierror])
   in_lflux = 0.04

   jhusdss_lowz_doublet_fit, lowz.wave[iuse], 1.-y[iuse], fltarr(nuse)+1./terror^2, $
       in_slope, in_intercept, in_center, in_separation, in_lflux, in_ratio, in_sigma, $
       slope=slope, intercept=intercept, center=center, separation=separation, $
       lflux=lflux, ratio=ratio, sigma=sigma, $
       err_slope=err_slope, err_intercept=err_intercept, err_center=err_center, $
       err_separation=err_separation, err_lflux=err_lflux, err_ratio=err_ratio, err_sigma=err_sigma, $
       maxwidth=maxwidth

   p = [slope, intercept, center, separation, lflux, ratio, sigma]
   perror = [err_slope, err_intercept, err_center, err_separation, err_lflux, err_ratio, err_sigma]

   yfit = jhusdss_lowz_doublet_func(lowz.wave, p)
   print, p
   print, perror

   all_lflux[2] = lflux
   all_lflux_err[2] = err_lflux
   all_sigma[2] = sigma
   all_sigma_err[2] = err_sigma

   clevel = slope*(lowz.wave-center)+intercept
   djs_plot, lowz.wave, y+clevel, $
       xra=xra, yra=yra, thick=thick, xthick=thick, ythick=thick, $
       xtickformat='(A1)', charsize=charsize, charthick=charthick
   djs_oplot, lowz.wave, 1.-yfit+clevel, $
       thick=thick, color='red'

   djs_xyouts, !x.crange[0]+0.75*(!x.crange[1]-!x.crange[0]), $
               !y.crange[0]+0.90*(!y.crange[1]-!y.crange[0]), $
               '45<b<135 kpc ', charsize=charsize, charthick=charthck
   djs_xyouts, !x.crange[0]+0.75*(!x.crange[1]-!x.crange[0]), $
               !y.crange[0]+0.80*(!y.crange[1]-!y.crange[0]), $
               '<b>='+string(lowz.rp[2],format='(i3)')+' kpc', charsize=charsize, charthick=charthick
   djs_xyouts, !x.crange[0]+0.05*(!x.crange[1]-!x.crange[0]), $
               !y.crange[0]+0.85*(!y.crange[1]-!y.crange[0]), $
               'W^{\lambda3934}='+string(lflux,format='(f5.3)')+'\pm'+string(err_lflux,format='(f5.3)')+' \AA', $
               charsize=1.6, charthick=charthick

   yra = [1.-0.15*0.25^3, 1.+0.15*0.25^3]
   y = lowz.fluxmean[*,3]
;  y = lowz.fluxmean[*,3]
   terror = sqrt((moment(y[ierror]))[1])
   print, 'typical error= ', terror

   in_intercept = median(1.-y[ierror])
   in_lflux = 0.02

   jhusdss_lowz_doublet_fit, lowz.wave[iuse], 1.-y[iuse], fltarr(nuse)+1./terror^2, $
       in_slope, in_intercept, in_center, in_separation, in_lflux, in_ratio, in_sigma, $
       slope=slope, intercept=intercept, center=center, separation=separation, $
       lflux=lflux, ratio=ratio, sigma=sigma, $
       err_slope=err_slope, err_intercept=err_intercept, err_center=err_center, $
       err_separation=err_separation, err_lflux=err_lflux, err_ratio=err_ratio, err_sigma=err_sigma, $
       maxwidth=maxwidth

   p = [slope, intercept, center, separation, lflux, ratio, sigma]
   perror = [err_slope, err_intercept, err_center, err_separation, err_lflux, err_ratio, err_sigma]

   yfit = jhusdss_lowz_doublet_func(lowz.wave, p)
   print, p
   print, perror


   all_lflux[3] = lflux
   all_lflux_err[3] = err_lflux
   all_sigma[3] = sigma
   all_sigma_err[3] = err_sigma

   clevel = slope*(lowz.wave-center)+intercept
   djs_plot, lowz.wave, y+clevel, $
       xra=xra, yra=yra, thick=thick, xthick=thick, ythick=thick, $
       xtitle=xtitle, yminor=5, charsize=charsize, charthick=charthick
   djs_oplot, lowz.wave, 1.-yfit+clevel, $
       thick=thick, color='red'

   djs_xyouts, !x.crange[0]+0.75*(!x.crange[1]-!x.crange[0]), $
               !y.crange[0]+0.90*(!y.crange[1]-!y.crange[0]), $
               '135<b<405 kpc ', charsize=charsize, charthick=charthick
   djs_xyouts, !x.crange[0]+0.75*(!x.crange[1]-!x.crange[0]), $
               !y.crange[0]+0.80*(!y.crange[1]-!y.crange[0]), $
               '<b>='+string(lowz.rp[3],format='(i3)')+' kpc', charsize=charsize, charthick=charthick
   djs_xyouts, !x.crange[0]+0.05*(!x.crange[1]-!x.crange[0]), $
               !y.crange[0]+0.85*(!y.crange[1]-!y.crange[0]), $
               'W^{\lambda3934}='+string(lflux,format='(f5.3)')+'\pm'+string(err_lflux,format='(f5.3)')+' \AA', $
               charsize=1.6, charthick=charthick

k_end_print

charsize=1.6
psfile = qapath+'/'+'lowz_caii_ew_density_'+string(nmfver, format='(I3.3)')+'.ps'
k_print, filename=psfile, axis_char_scale=1.5, xsize=8, ysize=8
;load_dp, /b
   xra = [5,500]
   yra = [2E-3,0.5]
   xtitle = '<b> (kpc)'
   ytitle = 'W^{\lambda3934} (\AA)'
   ytitle1 = 'N (Ca II) (cm^{-2})'
;  factor = 0.511*1E6*3E18/!dpi/3934.78/0.6267
   ;; see Eq. 9.15 in Bruce Drain
   factor = 1.13E12*1E8/0.6267/(3934.78)^2
   djs_plot, lowz.rp, all_lflux, psym=4, xra=xra, yra=yra, /xlog, /ylog, $
       thick=thick, xthick=thick, ythick=thick, charsize=charsize, charthick=charthick, $
       xtitle=xtitle, ytitle=ytitle, color='blue', ystyle=9
   djs_oplot, lowz.rp, all_lflux, thick=thick, color='blue'
   plotsym, 0, 2, /fill
   oploterror, lowz.rp, all_lflux, all_lflux_err, psym=8, thick=thick, $
       color=djs_icolor('blue'), errcolor=djs_icolor('blue'), symsize=1
   djs_axis, yaxis=1, yra=yra*factor, /ylog, ythick=thick, ytitle=ytitle1, $
       charsiz=charsize, charthick=charthick
k_end_print

;djs_plot, lowz.rp, all_sigma, psym=4, xra=[4, 300], yra=[0.3, 4], /xlog, /ylog
;djs_oplot, lowz.rp, all_sigma_
;oploterror, lowz.rp, all_sigma, all_sigma_err, psym=4

stop

end
