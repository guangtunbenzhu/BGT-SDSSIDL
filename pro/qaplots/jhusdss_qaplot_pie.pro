;jhusdss_qaplot_pie, nmfver

qso = jhusdss_qso_readin()
absorber = jhusdss_absorber_readin(107)

choice_ps = 0
read,'choice_ps=',choice_ps
load_dp,/b

if choice_ps then begin
   SET_PLOT,'ps'
   DEVICE,filename='pie_plot_abs.eps',/color,bits_per_pixel=32,$
     yoffset=0,xoffset=0,xsize=8,ysize=8,/inches,/encapsulated,/cmyk
   load_dp,/ps
   !p.charsize = 3
   !p.thick = 10
   !x.thick = 12
   !x.thick = 12
endif

my_z_max = 2.9

erase
loadct,1 ;blue
my_colors = 0*indgen(n_elements(absorber))+2. / my_z_max * 255
PLOTSYM,0,.3,/fill
PIE_PLOT,qso.z,qso.ra,psym=8,rotate=-98,rrange=[0,my_z_max],color=my_colors,$
    thaxes=98,rlabel='redshift'

;; -- absorbers
loadct,3  ;red
my_colors = 0*indgen(n_elements(absorber))+2 / my_z_max * 255
PLOTSYM,0,.2,/fill
PIE_PLOT,absorber.zabs,absorber.ra,psym=8,rotate=-98,rrange=[0,my_z_max],color=my_colors,/over

end
