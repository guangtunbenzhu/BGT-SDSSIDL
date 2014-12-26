

;; ------------------------------------------------
;;             USER PARAMETERS
;; ------------------------------------------------
choice_image = 1
print,'Image choice: [0] QSO'
print,'              [1] QSO flat'
read,'              [2] Absorber:  ',choice_image

;!p.font=-1
!p.charthick=2.
!p.thick=2.
!x.thick=1.8
!y.thick=1.8
my_charsize = 1.5

add_labels = 0b
add_absorbers = 0b

n_congrid_x = 7600L   ;; resolution of the final image
n_congrid_y = 3800L   ;; resolution of the final image
my_nlevels = 15L
choice_name = 1
;; ------------------------------------------------


;; -- Display choice (PS/PNG/Screen)
choice_ps = 0
read,'choice_ps [1=ps, 2=screen, 3=png]: ',choice_ps

;; -- Load data from fits files
choice_load_data = 0
read,'load data? [1=yes, 2=no]: ',choice_load_data
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
    ii = where(absorber.zabs lt absorber.zqso-0.04 and absorber.zabs gt 1250.*(absorber.zqso+1.)/2800.-1.)
    absorber = absorber[ii]

    
    loglam = jhusdss_get_loglam(minwave=3700., maxwave=9200D0)
    wave = 10.^loglam
    iwave_min = value_locate(wave, 3800.)
    wave = wave[iwave_min:*]
endif



if choice_image eq 0 then begin

  top_lines_name = ['Ly\beta', 'Ly\alpha', 'C_{ }IV'];, 'C III]', 'C II]']
  top_factor = (1.+4.0)
  top_lines_wave = [1025.72*top_factor-30., 1215.67*top_factor-10., 1549.06*top_factor-20.];, 1908.73*top_factor-40., 2326.44*top_factor-10.]

;; -- QSO composite
    my_filename = 'QSO_composite.eps'
    my_ytit = 'QSO redshift'

    Z = qso_image.image_median[iwave_min:*,*]
    my_z_min = .1
    my_z_max = 4.
    Z = (Z<my_z_max)>my_z_min
    
    my_z_list = [min(QSO_image.z),0.5,1,1.5,2,3,max(QSO_image.z)]
    redshift_list = QSO_image.z
    loadct,3

  right_lines_name = ['C_{ }III]', 'Mg_{ }II', 'H\beta', 'H\alpha']
  right_lines_z = [9200./1908.73-1., 9200./2798.75-1., 9200./4862.68-1., 9200./6564.-1.]-0.02
  right_lines_zgrid = right_lines_z
  for iline=0L, n_elements(right_lines_name)-1L do $
      right_lines_zgrid[iline] = value_locate(redshift_list, right_lines_z[iline])

endif


if choice_image eq 1 then begin

  top_lines_name = ['Ly\beta', 'Ly\alpha', 'C_{ }IV'];, 'C III]', 'C II]']
  top_factor = (1.+4.0)
  top_lines_wave = [1025.72*top_factor-30., 1215.67*top_factor-10., 1549.06*top_factor-20.];, 1908.73*top_factor-40., 2326.44*top_factor-10.]

;; -- QSO flat
    my_filename = 'Flat_composite.eps'
    my_ytit = 'QSO redshift'
    my_y_title = 0.7
    Z = qso_resi_image.image_median[iwave_min:*,*]
    my_z_min = 0.80
    my_z_max = 1.10
    Z = (Z<my_z_max)>my_z_min
    
    my_z_list = [min(QSO_image.z),0.5,1,1.5,2,3,max(QSO_image.z)]
    redshift_list = QSO_resi_image.z
    loadct,3

  right_lines_name = ['C_{ }III]', 'Mg_{ }II', 'H\beta', 'H\alpha']
  right_lines_z = [9200./1908.73-1., 9200./2798.75-1., 9200./4862.68-1., 9200./6564.-1.]-0.02
  right_lines_zgrid = right_lines_z
  for iline=0L, n_elements(right_lines_name)-1L do $
      right_lines_zgrid[iline] = value_locate(redshift_list, right_lines_z[iline])


endif


if choice_image eq 2 then begin

  n_congrid_y=750
  top_lines_name = ['Fe_{ }II', 'O_{ }I', 'C_{ }II', 'Si_{ }IV', 'Si_{ }II', 'C_{ }IV', 'Fe_{ }II ', 'Al_{ }II', 'Al_{ }III', $
        'Fe_{ }II  ', 'Fe_{ }II', 'Mg_{ }II']
  top_factor = (1.+2.2)
  top_lines_wave = [1260.53*top_factor-40., 1304.86*top_factor-10., 1334.53*top_factor+20., $
                    1402.76*top_factor-10., 1526.71*top_factor-40., $
                    1550.78*top_factor+10., 1608.45*top_factor-20., $
                    1670.79*top_factor-25., 1862.79*top_factor-30., $
                    2382.77*top_factor-60., 2600.17*top_factor-35., $
                    2803.53*top_factor-30]


;; -- absorber composite
    my_filename = 'Abs_composite.eps'
    my_ytit = 'Absorber redshift'

    Z = abs_image.image_median[iwave_min:*,*]
    my_z_min = 0.8
    my_z_max = 1.1
    Z = (Z<my_z_max)>my_z_min

    redshift_list = abs_image.z
    my_z_list = [0.5,1,1.5,2]
    loadct,3

  right_lines_name = ['Mg_{ }I', 'Ca_{ }II']
  right_lines_z = [9200./2803.53-1., 9200./3950.00-1.]-0.02
  right_lines_zgrid = right_lines_z
  for iline=0L, n_elements(right_lines_name)-1L do $
      right_lines_zgrid[iline] = value_locate(redshift_list, right_lines_z[iline])


;; -- we can create a sky composite and subtract it to the data:
;sky = CONGRID(REFORM(MEDIAN(b,dim=2),n_lambda,1),n_lambda,n_z)
;Z = ((b-sky+1.)>0.8)<1.1

endif


;; ---------------------------------------------------------------
load_dp,/b

if choice_ps eq 3 then begin
   SET_PLOT,'z'
   DEVICE,set_resolution=[1000,1000]
   load_dp,/ps
;  !p.charthick = 6
;  !x.thick = 6
;  !y.thick = 6
;  my_thick = 10
endif


if choice_ps eq 1 then begin
   xsize=16
   ysize=8

;  if (choice_image eq 2) then begin
;     xsize=10
;     ysize=5
;  endif

   SET_PLOT,'ps'
   DEVICE,filename=my_filename,/color,bits_per_pixel=32,$
     yoffset=0,xoffset=0,xsize=xsize,ysize=ysize,/inches,/cmyk ;/encapsulated,/cmyk
   load_dp,/ps
;  !p.charthick = 6
;  !x.thick = 6
;  !y.thick = 6
;  my_thick = 10
endif
;; ---------------------------------------------------------------


n_lambda = (size(Z))(1)
n_z = (size(Z))(2)
y_axis = indgen(n_z)


Z = CONGRID(Z,n_congrid_x,n_congrid_y,/INTERP)
x_axis_congrid = (make_vector(min(wave),max(wave),n_congrid_x,/log)).bound_min
y_axis_congrid = (make_vector(min(y_axis),max(y_axis),n_congrid_y)).bound_min

n_redshift = n_elements(redshift_list)
ratio = float(n_redshift)/float(n_congrid_y)


;; -- artificially flatten the Lyman-alpha forest
;if keyword_set(nolya) then begin
if choice_image eq 1 then begin

    image_median = MEDIAN(Z)
    minmax_original = minmax(Z)
    Z2 = Z
    for i_qso=0,n_congrid_y-1 do begin
        my_z = redshift_list[i_Qso*ratio]
        id = where(x_axis_congrid/(1.+my_z) le 1218., count)
        if count gt 10 then begin
            my_x_tilt = 1.*(indgen(n_elements(id))+2)/n_elements(id)
            tmpmoment = moment(Z[id, i_qso], sdev=sdev)
            tmpnorm = median(Z[id, i_qso])+sdev*0.5
            Z2[id, i_qso] = (((Z[id, i_qso]/tmpnorm)>0.8)<1.1)*exp(-0.04*(1.-my_x_tilt/(count+1.)))
;           my_overall_shift = .7*image_median*(my_z-2.) / (my_x_tilt)^(0.5)
;           Z2[id,i_qso] = ((Z[id,i_qso] - my_overall_shift)>0.8)<1.1 ;minmax_original[1]
        endif
    endfor
    Z2 = (Z2>0.8)
    Z = Z2
; for i=35,60 do begin
;     print,i
;     plot,Z[*,n_congrid-i],yr=[0,.6]
;     oplot,Z2[*,n_congrid-i],color=getcolor('red',1),linestyle=2
;     pause
; endfor
; pause
endif
;endif

contour,(Z),x_axis_congrid,y_axis_congrid,nlevels=my_nlevels,/fill,$
  xmargin=[6,4], ymargin=[6,3], charsize=my_charsize,xstyle=5,ystyle=5


if (add_labels) then begin

   axis, xaxis=0, xtitle=textoidl('\lambda_{obs} (\AA)'), charsize=my_charsize, xticklen=0.01, xstyle=1
   axis, yaxis=0, ytitle='', charsize=my_charsize, ytickformat='(A1)', yticklen=0.001, ystyle=1
   axis, xaxis=1, xtitle='', xtickformat='(A1)', charsize=my_charsize, xticklen=0.01, xstyle=1
   axis, yaxis=1, ytitle='', charsize=my_charsize, ytickformat='(A1)', yticklen=0.001, ystyle=1

   ;; -- take care of the left/right axis
   my_title_x = 3400.
   my_x = 3450.
   if (choice_image eq 2) then my_title_x = 3450.
   if (choice_image eq 2) then my_x = 3500.

   my_title_y_position = n_elements(redshift_list)/2.
   xyouts,my_title_x,my_title_y_position,my_ytit,orientation=90,align=0.5,charsize=my_charsize

   for i=0,n_elements(my_z_list)-1 do begin
       my_String = string(format="(f3.1)",my_z_list[i])
       my_y = value_locate(redshift_list, my_z_list[i])
       xyouts,my_x,my_y-0.015*(!y.crange[1]-!y.crange[0]),my_string,charsize=my_charsize-0.4
       oplot,min(x_axis_congrid)*[1,1.01],my_y*[1,1]
       oplot,max(x_axis_congrid)*[1,1.01],my_y*[1,1]
   endfor

;   ;; -- add vertical lines on borders
;   oplot,(min(x_axis_congrid))*[1,1],$
;;      minmax(y_axis_congrid)
;  oplot,(x_axis_congrid[n_elements(x_axis_congrid)-2])*[1,1],$
;     minmax(y_axis_congrid)

   for iline=0L, n_elements(top_lines_name)-1L do begin 
       djs_xyouts, top_lines_wave[iline], !y.crange[0]+1.015*(!y.crange[1]-!y.crange[0]), $
           top_lines_name[iline], charsize=1.0, orientation=45, charthick=2.0
   endfor

   for iline=0L, n_elements(right_lines_name)-1L do begin 
       djs_xyouts, !x.crange[0]+1.015*(!x.crange[1]-!x.crange[0]), $
           right_lines_zgrid[iline], $
           right_lines_name[iline], charsize=1.0, orientation=45, charthick=2.0
   endfor

endif 

;; -- add absorbers
if add_absorbers eq 1 then begin
    n_abs = n_elements(absorber)
    PLOTSYM,0,.01
    for i_abs=0,n_abs-1,2 do begin
        my_z = absorber[i_abs].zabs
        x = 2796.*(1.+my_z)
        y = INTERPOL(y_axis,QSO_image.z,absorber[i_abs].zqso)
        plots,x,y,psym=8,color=getcolor('black',1)
    endfor
endif


;; -- add our names
;if choice_image eq 2 AND choice_name eq 1 then begin
;    djs_xyouts,9060,!y.crange[0]+0.03*(!y.crange[1]-!y.crange[0]), $
;      '!8'+'Guangtun Ben Zhu & Brice Menard (2012)!X',align=1,$
;       ANSI_Value('Zhu & Ménard (2012)'),align=1,$
;      charsize=1.2
;endif


if choice_ps eq 3 then begin
    write_png,'my.png',tvrd()
endif

if choice_ps eq 1 then begin
    device,/close
    set_plot,'x'
    print,'device closed'
endif

for i=0,2 do BEEP
end
