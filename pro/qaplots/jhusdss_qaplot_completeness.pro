pro jhusdss_qaplot_completeness, nmfver

if (n_elements(nmfver) eq 0) then begin
   splog, 'nmfver required.'
   return
endif

mc = jhusdss_montecarlo_readin(nmfver)
log10w = mc.log10w
z = mc.z
f_log10w = float(total(mc.ndetected, 2))/float(total(mc.ncovered, 2)+(total(mc.ncovered, 2) eq 0.))
;f_z = float(total(mc.ndetected, 1))/float(total(mc.ncovered, 1)+(total(mc.ncovered, 1) eq 0.))
f_z = median(float(mc.ndetected)/float(mc.ncovered+(mc.ncovered eq 0)), dimen=1);float(total(mc.ndetected, 1))/float(total(mc.ncovered, 1)+(total(mc.ncovered, 1) eq 0.))

;; init
thick=6
charsize=1.3
charthick=3
;; path
path=jhusdss_get_path(/nmfqso)+'/'+string(nmfver, format='(I3.3)')
qapath = path+'/QAplots'

;; redshift distribution
psfile = qapath+'/'+'Completeness_'+string(nmfver, format='(I3.3)')+'.ps'
k_print, filename=psfile, axis_char_scale=1.4, xsize=8, ysize=9

   xra = [0.E-1, 6.]
   yra = [0, 1]
   xtitle = textoidl('W_0^{\lambda2796} (\AA)')
   ytitle = textoidl('Completeness f')
   title = textoidl('Average over All Redshifts')
   pos = [0.10, 0.58, 0.90, 0.93]
   ii  = where(f_log10w gt 0.)
   x = 10.^log10w[ii]
   y = f_log10w[ii]
   djs_plot, x, y, xra=xra, yra=yra, title=title, $
       thick=thick+2, xthick=thick, ythick=thick, $
       xtitle=xtitle, ytitle=ytitle, xst=1, yst=1,$
       charsize=charsize, charthick=charthick, position=pos,  color='blue'
;  djs_oplot, x, y, psym=4, thick=thick

   xra = [0.36, 2.3]
   yra = [0, 1]
   xtitle = textoidl('z')
   ytitle = textoidl('Completeness f')
   title = textoidl('Average over All Rest Equivalent Widths')
   pos = [0.10, 0.10, 0.90, 0.45]
   ii  = where(f_z gt 0.)
   x = z[ii]
   y = f_z[ii]
   djs_plot, x, y, xra=xra, yra=yra, title=title, $
       thick=thick, xthick=thick, ythick=thick, $
       xtitle=xtitle, ytitle=ytitle, xst=1, yst=1,$
       charsize=charsize, charthick=charthick, position=pos, color='blue', /noerase
;  djs_oplot, x, y, psym=4, thick=thick
  
k_end_print

end
