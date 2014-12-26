

file = '/home/menard/composite_image/FIT/my_image.fit'
;file = 'FIT/MgII_composite_image.fit'
load_dp,/b

nmfver = 005
absorber = jhusdss_absorber_readin(nmfver) ;, /byabs, /trim, boss=boss)
absorber = absorber(SORT(absorber.zabs))

choice_ps = 1
read,'choice_ps=',choice_ps
load_dp,/b

if choice_ps then begin
   SET_PLOT,'ps'
   DEVICE,filename='composite_image.eps',/color,bits_per_pixel=32,$
;    yoffset=0,xoffset=0,xsize=3.840*3,ysize=0.603*6,/cmyk ;/encapsulated,/cmyk
     yoffset=0,xoffset=0,xsize=20,ysize=20,/inches,/cmyk ;/encapsulated,/cmyk
   load_dp,/ps
   !p.charthick = 8
   !x.thick = 8
   !y.thick = 8
   my_thick = 12
endif



a = MRDFITS(file)
a(where(a eq 0))=1.

z_min = 0.8
z_max = 1.1
b = (a>z_min)<z_max


n_lambda = (size(b))(1)
n_z = (size(b))(2)
image_square = CONGRID(b,n_lambda,n_lambda)


loadct,3 ;10

print,minmax(a)
print,moment(a)


;; observer-frame wavelength
minwave = 3800.
maxwave = 9200.
loglam = jhusdss_get_loglam(minwave=minwave, maxwave=maxwave)
n_lambda = n_elements(loglam)

;; redshift vector
n_in_group = 50
n_abs = n_elements(absorber)
n_z_bin = n_abs / n_in_group
my_z_median = fltarr(n_z_bin)
for i=0L,n_z_bin-1 do begin
    ind_min = i*n_in_group
    ind_max = ind_min + n_in_group
    my_z_median[i] = MEDIAN(absorber[ind_min:ind_max].zabs)
endfor

n_z = (size(a))(2)
y_axis = indgen(n_z)

;sky = CONGRID(REFORM(MEDIAN(b,dim=2),n_lambda,1),n_lambda,n_z)

;Z = ((b-sky+1.)>0.8)<1.1
Z = ((b)>0.8)<1.1

contour,Z,10^loglam,y_axis,nlevels=20,/fill,$
  xtit=textoidl('\lambda_{obs} [Angstrom]'),$
  ytit='redshift   ',ytickname=replicate(' ',30), xmargin=[5,2],$
  charsize=3.0, xst=9, yst=5;, xra=[minwave, maxwave], yra=minmax(my_z_median), xst=1, yst=1
;contour,b,10^loglam,my_z_median,nlevels=20,/fill,$
;  xtit=textoidl('\lambda_{obs} [Angstrom]')

my_z_list = [0.5,.6,.7,.8,.9,1.,1.1,1.2,1.3,1.4,1.5,1.6,1.8,2.2]
my_y_list = fltarr(n_elements(my_z_list))
for i_z=0,n_elements(my_z_list)-1 do begin
;   my_y_list[i_z] = 1.d0*(where(my_z_median ge my_z_list[i_z]))(0)
    tmp = min(abs(my_z_median-my_z_list[i_z]), imin)
    my_y_list[i_z] = 1.d0*imin/687.*602.
;   my_y_list[i_z] = (size(Z))[2]*(my_z_list[i_z]-min(my_z_median))/(max(my_z_median)-min(my_z_median))
endfor

my_x = 3550.
for i=0,n_elements(my_z_list)-1 do begin
    my_y = my_y_list[i]
    xyouts,my_x,my_y,string(format="(f3.1)",my_z_list[i]),charsize=3.0
endfor

;plothist,a,bin=0.001,/ylog,yr=[1,2e5],$
;  position=[.7,.2,.88,.4],/noerase,charsize=1,xr=[.5,1.24],xstyle=1


;xticks=[1260.53, 1302.17*3.29, 1304.86*3.29, 1334.53*3.29, $
;       1393.76, 1402.77*3.29, 1526.71*3.29, 1548.20*3.29, $
;       1550.78, 1608.45*3.29, 1670.79*3.29, 1854.72*3.29, $
;       1862.79, 2344.21*3.29, 2374.46*3.29, 2382.77*3.29,  $
;       2586.55, 2600.17*3.29, 2796.35*3.29, 2803.53*3.29]*3.20
;djs_axis, xaxis=1, xtickv=ticks, xthick=8, xtickformat='(A1)'
names = ['Fe II   ', 'O I', 'C II', 'Si IV', 'Si II', 'C IV', 'Fe II ', 'Al II', 'Al III', $
        'Fe II  ', 'Fe II', 'Mg II']
lines = [1260.53*3.18-30., 1304.86*3.18-10., 1334.53*3.18-10., 1402.76*3.18-20., 1526.71*3.18-40., $
         1550.78*3.18+10., 1608.45*3.18-20., 1670.79*3.18-25., 1862.79*3.18-30., 2382.77*3.18-60., $
         2600.17*3.18-35., 2803.53*3.18-30]

for i=0L, n_elements(names)-1L do begin
    djs_xyouts, lines[i], !y.crange[0]+1.005*(!y.crange[1]-!y.crange[0]), $
        names[i], charsize=2.2, charthick=6.0, color='black', orientation=45
endfor
        

if choice_ps then begin
   device,/close
   set_plot,'x'
   print,'device closed'
endif


end
