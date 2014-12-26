;+
; Documentation Needed!
;-
pro jhusdss_composite_qaplot, objs, stack, qapath=qapath, qafile=qafile

if (n_elements(qapath) eq 0) then qapath='.'
if (n_elements(qafile) eq 0) then qafile='QSO_composite_qaplot.ps'

thick=3
charthick=3
charsize=1.5

ztag = jhusdss_get_tags(/ztag)
zindex = tag_indx(objs[0], ztag)

magtag = jhusdss_get_tags(/magtag)
magindex = tag_indx(objs[0], magtag)


k_print, filename=qapath+'/'+qafile, axis_char_scale=1.3, xsize=10, ysize=6
   xtitle = textoidl('Wavelength (\AA)')

   djs_plot, stack.wave, stack.fluxmean, thick=thick, xthick=thick, ythick=thick, $
       charthick=charthick, charsize=charsize, xra=xra, yra=yra, $
       xtitle=xtitle, ytitle=ytitle
   djs_oplot, stack.wave, stack.fluxmedian, thick=thick, color='red'
   djs_oplot, stack.wave, stack.fluxgeomean, thick=thick, color='blue'

   xra = [0, 5]
   yra = [-22., -29.]
   djs_plot, objs.(zindex), objs.(magindex), thick=thick, xthick=thick, ythick=thick, $
       charthick=charthick, charsize=charsize, xra=xra, yra=yra, $
       xtitle=xtitle, ytitle=ytitle, psym=4, symsize=0.2

k_end_print

end
