;pro jhusdss_garching_galqso_stack, nmfver, boss=boss, ivarweight=ivarweight

;read,'BOSS? [1=yes, 0=no]: ',BOSS

;; galaxy redshift range
zgalmin = 0.030
zgalmax = 0.400
;zgalmin = 0.030
;zgalmax = 0.060
massive = 0b
quiescent = 0b

nmfver = 106
ivarweigth = 1b
overwrite = 1b
savespec = 0b
savemean = 1b
snrcut = 2.
sdevcut = 0.10
sigma_cut = 2.

;; ca II
ca_line_wave = [3934.79, 3969.59]
;;ca_line_wave = [3935.00, 3969.00]
ca_xra = [3800., 4300.]
;; Na I
na_line_wave = [5890.00, 5896.00]
na_xra = [5700., 6100.]

ha_line_wave = [6563.00, 6584.00]
ha_xra = [6400., 6800.]

line_wave = ca_line_wave
xra = ca_xra

stackpath = jhusdss_get_path(/nmfqso)+'/'+string(nmfver, format='(I3.3)')+'/Stack/'
outfile = stackpath + jhusdss_garching_galqso_stack_filename(nmfver, boss=boss)
psfile = repstr(outfile, '.fits', '_masssfr_depend.ps')

if (file_test(outfile) and ~overwrite) then begin
   splog, 'File already exists, use /overwrite to overwrite' 
   return
endif else begin
   splog, 'Will write into this file: '
   print, psfile
endelse

quiescent_mass = [10.^10.35, 10.^10.98]
sf_mass = [10.^9.57, 10.^10.52]
;quiescent_caii = [2809., 4854.]
;err_quiescent_caii = [1434., 1185.]
;sf_caii = [2332., 5099.]
;err_sf_caii = [1026., 1019.]

;; the data
;quiescent_caii = [693., 6011.]
;err_quiescent_caii = [2420., 1974.]
;sf_caii = [4150., 6773.]
;err_sf_caii = [1602, 1728.]

quiescent_caii = [881., 2016.]/0.62
err_quiescent_caii = [894., 911.]/0.62
sf_caii = [3227., 5742.]/0.62
err_sf_caii = [916, 914.]/0.62

;scale_quiescent_caii = [894., 2448.]
;err_scale_quiescent_caii = [1085., 1106.]
;scale_sf_caii = [2620., 6114.]
;err_scale_sf_caii = [1112., 1110.]
scale_quiescent_caii = [864., 2365.]/0.62
err_scale_quiescent_caii = [1049., 1069.]/0.62
scale_sf_caii = [2532., 5909.]/0.62
err_scale_sf_caii = [1075., 1073.]/0.62

all_mass = [10.^10.30, 10.^10.30]/0.62
all_caii = [3076., 3076.]/0.62
err_all_caii = [618., 618.]/0.62
;all_caii = [3183., 3183.]
;err_all_caii = [639., 639.]

thick=8
xthick=8
ythick=8
charsize=1.3
charthick=3

xra=[10.^9.3,10.^11.2]
yra=[10.^2.3, 10.^4.2]
xtitle=textoidl('M_* (M_\odot)')
;ytitle=textoidl('<M(Ca II)> within virial radius [M_\odot]')
ytitle=textoidl('<M(Ca II, <r_{vir})> [M_\odot]')

pos = [0.17, 0.20, 0.95, 0.95]
k_print, filename=psfile, axis_char_scale=1.5, xsize=7, ysize=7

  djs_plot, quiescent_mass, quiescent_caii, /xlog, /ylog, xra=xra, yra=yra, $
       thick=thick, xthick=thick, ythick=thick, $
       xtitle=xtitle, ytitle='', xst=1, yst=1, $
       charsize=charsize, charthick=charthick, /nodata, ytickformat='(A1)'

   djs_xyouts, 10.^8.95, 270., '3\times10^2', charsize=charsize+0.3, charthick=charthick+0.5
   djs_xyouts, 10.^8.95, 900., '1\times10^3', charsize=charsize+0.3, charthick=charthick+0.5
   djs_xyouts, 10.^8.95, 4500., '5\times10^3', charsize=charsize+0.3, charthick=charthick+0.5
   djs_xyouts, 10.^8.95, 9000., '1\times10^4', charsize=charsize+0.3, charthick=charthick+0.5
   djs_xyouts, 10.^8.90, 0370., ytitle, charsize=charsize+0.5, charthick=charthick+0.5, orientation=90.

;  plotsym, 4, /fill
;  plotsym, 4, thick=thick
;  oploterror, quiescent_mass, quiescent_caii, err_quiescent_caii, color=djs_icolor('red'), psym=8, thick=thick, symsize=2, errcolor=djs_icolor('red')
;  djs_oplot, quiescent_mass, quiescent_caii, color='red', linestyle=0, thick=thick
;  plotsym, 0, /fill
;  plotsym, 0, thick=thick
;  oploterror, sf_mass, sf_caii, err_sf_caii, color=djs_icolor('navy'), errcolor=djs_icolor('navy'), psym=8, symsize=2, thick=thick
;  djs_oplot, sf_mass, sf_caii, color='navy', linestyle=0, thick=thick

   plotsym, 4, /fill
;  plotsym, 4, thick=thick
   oploterror, quiescent_mass, scale_quiescent_caii, err_scale_quiescent_caii, color=djs_icolor('red'), psym=8, thick=thick, symsize=2, errcolor=djs_icolor('red')
   djs_oplot, quiescent_mass, scale_quiescent_caii, color='red', linestyle=0, thick=thick
   plotsym, 0, /fill
;  plotsym, 0, thick=thick
   oploterror, sf_mass, scale_sf_caii, err_scale_sf_caii, color=djs_icolor('navy'), errcolor=djs_icolor('navy'), psym=8, symsize=2, thick=thick
   djs_oplot, sf_mass, scale_sf_caii, color='navy', linestyle=0, thick=thick

;  plotsym, 3, /fill
;  oploterror, all_mass, all_caii, err_all_caii, color=djs_icolor('gray'), errcolor=djs_icolor('gray'), psym=8, symsize=2, thick=thick

   plotsym, 0, /fill
   djs_xyouts, 10.^9.4, 10.^2.52, 'Star-forming', charsize=charsize, charthick=charthick+0.5, color='navy'
   djs_oplot, 10.^10.12*[1,1], 10.^2.54*[1,1], psym=8, symsize=1.5, thick=thick, color='navy'
   plotsym, 4, /fill
   djs_xyouts, 10.^9.4, 10.^2.42, '   Quiescent', charsize=charsize, charthick=charthick+0.5, color='red'
   djs_oplot, 10.^10.12*[1,1], 10.^2.44*[1,1], psym=8, symsize=1.5, thick=thick, color='red'

;  plotsym, 3, /fill
;  djs_xyouts, 10.^9.4, 10.^2.422, '         All', charsize=charsize, charthick=charthick+0.5, color='gray'
;  djs_oplot, 10.^10.12*[1,1], 10.^2.442*[1,1], psym=8, symsize=1.5, thick=thick, color='gray'

k_end_print

end
