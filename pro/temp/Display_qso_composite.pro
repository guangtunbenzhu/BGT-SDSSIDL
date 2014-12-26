;; Use BM_Display_composites.pro instead

image_choice = 1

choice_load_data = 0
read,'load data? [1=yes]: ',choice_load_data
if choice_load_data eq 1 then begin
    nmfver = 106
    qso_image =  jhusdss_qso_composite_image_readin(nmfver)
    qso_resi_image =  jhusdss_qso_composite_image_readin(nmfver, /residual)
    qso =  jhusdss_qso_readin()
    
    abs_image =  jhusdss_absorber_composite_image_readin(nmfver)
    abs_path=jhusdss_get_path(/nmfqso)+'/'+string(nmfver, format='(I3.3)')+'/Absorbers'
    abs_filename = abs_path+'/OnlyAbsorbers_Window_'+jhusdss_absorbers_filename(nmfver, /mgii, boss=boss)

;   absorber = jhusdss_absorber_readin(106)
    absorber = [mrdfits(abs_filename,1), mrdfits(abs_filename,2), mrdfits(abs_filename,3), mrdfits(abs_filename,4)]
    
    loglam = jhusdss_get_loglam(minwave=3700., maxwave=9200D0)
    wave = 10.^loglam
    iwave_min = value_locate(wave, 3800.)
    wave = wave[iwave_min:*]
endif


if image_choice eq 0 then begin
;; -- QSO composite
    Z = qso_image.image_median[iwave_min:*,*]
    my_z_min = .1
    my_z_max = 4.
    Z = (Z<my_z_max)>my_z_min
;    Z = Z - min(Z)
;    Z = sqrt(Z)  ;; increase contrast
    
    my_z_list = [min(QSO_image.z),0.5,1,2,3,4,5] ;[0.5,1,1.5,2,2.5,3,3.5,4,5]
    redshift_list = QSO_image.z
    loadct,3 ;20
endif


if image_choice eq 1 then begin
;; -- QSO flat
    Z = qso_resi_image.image_median[iwave_min:*,*]
    my_z_min = 0.8
    my_z_max = 1.2
    Z = (Z<my_z_max)>my_z_min
    Z = Z - min(Z)

;    Z = sqrt(Z)  ;; increase contrast
    
    my_z_list = [min(QSO_resi_image.z),0.5,1,2,3,4,5] ;[0.5,1,1.5,2,2.5,3,3.5,4,5]
    redshift_list = QSO_resi_image.z
    loadct,3
endif


;; ---------------------------------------------------------------
choice_ps = 0
;read,'choice_ps=',choice_ps
load_dp,/b

if choice_ps then begin
   SET_PLOT,'ps'
   DEVICE,filename='composite_image.eps',/color,bits_per_pixel=32,$
     yoffset=0,xoffset=0,xsize=8,ysize=8,/inches,/cmyk ;/encapsulated,/cmyk
   load_dp,/ps
   !p.charthick = 6
   !x.thick = 6
   !y.thick = 6
   my_thick = 10
endif
;; ---------------------------------------------------------------



n_lambda = (size(Z))(1)
n_z = (size(Z))(2)
y_axis = indgen(n_z)

;n = 1000
;W = CONGRID(Z,n,n,/INTERP)
;contour,(W),nlevels=20,/fill,$
;  xtit=textoidl('\lambda_{obs} [Angstrom]'),$
;  ytit='redshift   ',ytickname=replicate(' ',30),xmargin=[5,2],$
;  charsize=1.7,xstyle=1,ystyle=1
;pause

contour,(Z),wave,y_axis,nlevels=20,/fill,$
  xtit=textoidl('\lambda_{obs} [Angstrom]'),$
  ytit='redshift   ',ytickname=replicate(' ',30),xmargin=[5,2],$
  charsize=1.7,xstyle=1,ystyle=1

loadct,0
contour,0*(Z),wave,y_axis,$
  xtit=textoidl('\lambda_{obs} [Angstrom]'),$
  ytit='redshift   ',ytickname=replicate(' ',30),xmargin=[5,2],$
  charsize=1.7,xstyle=1,ystyle=1,/nodata,/noerase


my_x = 3600.
for i=0,n_elements(my_z_list)-1 do begin
    my_y = 1.d0*(where(redshift_list ge my_z_list[i]))(0)
    xyouts,my_x,my_y,string(format="(f3.1)",my_z_list[i]),charsize=1.2
endfor

stop
;pause

ici:
n_abs = n_elements(absorber)
PLOTSYM,0,.01
for i_abs=0,n_abs-1 do begin
    my_z = absorber[i_abs].zabs
    x = 2796.*(1.+my_z)
    y = INTERPOL(y_axis,QSO_image.z,absorber[i_abs].zqso)
    plots,x,y,psym=8,color=getcolor('black',1)
endfor



if choice_ps then begin
   device,/close
   set_plot,'x'
   print,'device closed'
endif


end
