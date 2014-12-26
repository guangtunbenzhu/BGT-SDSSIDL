pro jhusdss_qaplot_train_excess_fit, nmfver

path=jhusdss_get_path(/nmfqso)+'/'+string(nmfver, format='(I3.3)')
infile = path+'/Absorbers/'+'JHU_nomatched_composite_'+string(nmfver, format='(I3.3)')+'.fits'

comp_high = mrdfits(infile, 1)
spec_high = {wave:comp_high.wave, z:0., residual:comp_high.fluxmedian, ivar:fltarr(n_elements(comp_high.wave))+1.}
comp_medium = mrdfits(infile, 2)
spec_medium = {wave:comp_medium.wave, z:0., residual:comp_medium.fluxmedian, ivar:fltarr(n_elements(comp_medium.wave))+1.}
comp_low = mrdfits(infile, 3)
spec_low = {wave:comp_low.wave, z:0., residual:comp_low.fluxmedian, ivar:fltarr(n_elements(comp_low.wave))+1.}
comp_verylow = mrdfits(infile, 4)
spec_verylow = {wave:comp_verylow.wave, z:0., residual:comp_verylow.fluxmedian, ivar:fltarr(n_elements(comp_verylow.wave))+1.}
comp_1 = mrdfits(infile, 5)
spec_1 = {wave:comp_1.wave, z:0., residual:comp_1.fluxmedian, ivar:fltarr(n_elements(comp_1.wave))+1.}
comp_2 = mrdfits(infile, 6)
spec_2 = {wave:comp_2.wave, z:0., residual:comp_2.fluxmedian, ivar:fltarr(n_elements(comp_2.wave))+1.}
comp_3 = mrdfits(infile, 7)
spec_3 = {wave:comp_3.wave, z:0., residual:comp_3.fluxmedian, ivar:fltarr(n_elements(comp_3.wave))+1.}
comp_4 = mrdfits(infile, 8)
spec_4 = {wave:comp_4.wave, z:0., residual:comp_4.fluxmedian, ivar:fltarr(n_elements(comp_4.wave))+1.}

jhusdss_finalpass_fit, spec_high, 0., newzabs=zabs_spec_high, err_newzabs=err_zabs_spec_high, $
       allzabs=allzabs_high, ew=ew_high, err_ew=err_ew_high, sigma=sigma_high, err_sigma=err_sigma_high, $
       meansigma=meansigma_high, err_meansigma=err_meansigma_high

jhusdss_finalpass_fit, spec_medium, 0., newzabs=zabs_spec_medium, err_newzabs=err_zabs_spec_medium, $
       allzabs=allzabs_medium, ew=ew_medium, err_ew=err_ew_medium, sigma=sigma_medium, err_sigma=err_sigma_medium, $
       meansigma=meansigma_medium, err_meansigma=err_meansigma_medium

jhusdss_finalpass_fit, spec_low, 0., newzabs=zabs_spec_low, err_newzabs=err_zabs_spec_low, $
       allzabs=allzabs_low, ew=ew_low, err_ew=err_ew_low, sigma=sigma_low, err_sigma=err_sigma_low, $
       meansigma=meansigma_low, err_meansigma=err_meansigma_low

jhusdss_finalpass_fit, spec_verylow, 0., newzabs=zabs_spec_verylow, err_newzabs=err_zabs_spec_verylow, $
       allzabs=allzabs_verylow, ew=ew_verylow, err_ew=err_ew_verylow, sigma=sigma_verylow, err_sigma=err_sigma_verylow, $
       meansigma=meansigma_verylow, err_meansigma=err_meansigma_verylow

jhusdss_finalpass_fit, spec_1, 0., newzabs=zabs_spec_1, err_newzabs=err_zabs_spec_1, $
       allzabs=allzabs_1, ew=ew_1, err_ew=err_ew_1, sigma=sigma_1, err_sigma=err_sigma_1, $
       meansigma=meansigma_1, err_meansigma=err_meansigma_1

jhusdss_finalpass_fit, spec_2, 0., newzabs=zabs_spec_2, err_newzabs=err_zabs_spec_2, $
       allzabs=allzabs_2, ew=ew_2, err_ew=err_ew_2, sigma=sigma_2, err_sigma=err_sigma_2, $
       meansigma=meansigma_2, err_meansigma=err_meansigma_2

jhusdss_finalpass_fit, spec_3, 0., newzabs=zabs_spec_3, err_newzabs=err_zabs_spec_3, $
       allzabs=allzabs_3, ew=ew_3, err_ew=err_ew_3, sigma=sigma_3, err_sigma=err_sigma_3, $
       meansigma=meansigma_3, err_meansigma=err_meansigma_3

jhusdss_finalpass_fit, spec_4, 0., newzabs=zabs_spec_4, err_newzabs=err_zabs_spec_4, $
       allzabs=allzabs_4, ew=ew_4, err_ew=err_ew_4, sigma=sigma_4, err_sigma=err_sigma_4, $
       meansigma=meansigma_4, err_meansigma=err_meansigma_4

bin_high    = [1.5, 4.0]
bin_medium  = [1.1, 1.5]
bin_low     = [0.8, 1.1]
bin_verylow = [0.2, 0.8]
bin_1       = [0.65, 0.8]
bin_2       = [0.5, 0.65]
bin_3       = [0.4, 0.5]
bin_4       = [0.2, 0.4]

bin_min = [1.5, 1.1, 0.8, 0.2, 0.65, 0.5, 0.4, 0.2]
bin_max = [4.0, 1.5, 1.1, 0.8, 0.8, 0.65, 0.5, 0.4]

median_ew= [1.96057, 1.24123, 0.924520, 0.523950, 0.715320, 0.574426, 0.452532, 0.264549]

;; init
thick=5
charsize=1.5
charthick=2.5

psfile = path+'/QAplots/'+'JHU_nomatched_composite_one_ew_'+string(nmfver, format='(I3.3)')+'.ps'
k_print, filename=psfile, axis_char_scale=1.4, xsize=8, ysize=6

   x = median_ew
   y_2796 = [ew_high[5], ew_medium[5], ew_low[5], ew_verylow[5], $
             ew_1[5], ew_2[5], ew_3[5], ew_4[5]]
   y_2803 = [ew_high[4], ew_medium[4], ew_low[4], ew_verylow[4], $
             ew_1[4], ew_2[4], ew_3[4], ew_4[4]]
   y_2600 = [ew_high[6], ew_medium[6], ew_low[6], ew_verylow[6], $
             ew_1[6], ew_2[6], ew_3[6], ew_4[6]]
   y_2586 = [ew_high[7], ew_medium[7], ew_low[7], ew_verylow[7], $
             ew_1[7], ew_2[7], ew_3[7], ew_4[7]]
   y_2383 = [ew_high[8], ew_medium[8], ew_low[8], ew_verylow[8], $
             ew_1[8], ew_2[8], ew_3[8], ew_4[8]]

   xra = [0.1, 3.0]
   yra = [0.01, 3.0]
   title = textoidl('')
   xtitle = textoidl('Median W_0^{\lambda 2796} (\AA)')
   ytitle = textoidl('Composite W_0 (\AA)')

   djs_plot, x, y_2796, xra=xra, yra=yra, $
       xtitle=xtitle, ytitle=ytitle, title=title, $
       thick=thick, xthick=thick, ythick=thick, $
       charsize=charsize, charthick=charthick, $
       xstyle=1, ystyle=1, psym=4, /xlog, /ylog

   djs_oplot, x, y_2803, $
       thick=thick, xthick=thick, ythick=thick, $
       xstyle=1, ystyle=1, psym=5, color='red'

   djs_oplot, x, y_2600, $
       thick=thick, xthick=thick, ythick=thick, $
       xstyle=1, ystyle=1, psym=6, color='blue'

   djs_oplot, x, y_2586, $
       thick=thick, xthick=thick, ythick=thick, $
       xstyle=1, ystyle=1, psym=2, color='dark green'

   djs_oplot, x, y_2383, $
       thick=thick, xthick=thick, ythick=thick, $
       xstyle=1, ystyle=1, psym=7, color='brown'

   djs_oplot, xra, xra, thick=thick
   djs_oplot, xra, 0.70*xra, thick=thick, color='red'
   djs_oplot, xra, 0.25*xra, thick=thick, color='blue'
   djs_oplot, xra, 0.12*xra, thick=thick, color='dark green'
;  djs_oplot, !x.crange, !x.crange, thick=thick
;  djs_oplot, !x.crange, 0.70*!x.crange, thick=thick
;  djs_oplot, !x.crange, 0.25*!x.crange, thick=thick

   colors = [djs_icolor('black'), djs_icolor('red'), djs_icolor('blue'), djs_icolor('dark green'), djs_icolor('brown')] 
   items = [textoidl('\lambda2796'), textoidl('\lambda2803'), textoidl('\lambda2600'), textoidl('\lambda2586'), textoidl('\lambda2383')]
   legend, items, color=colors, psym=[4, 5, 6, 2, 7], charsize=charsize, charthick=charthick, $
       box=0, /right, /bottom

k_end_print

stop
end
