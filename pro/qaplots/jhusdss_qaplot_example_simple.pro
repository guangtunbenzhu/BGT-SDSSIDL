pro jhusdss_qaplot_example_simple, nmfver, choose=choose

linename=['Mg I', 'Mg II', 'Fe II', 'Fe II', 'Al II', 'C IV', 'Si II', 'Si IV', 'C II', 'O I'] 
linewave=[2860.,  2790.,   2590.,   2360.,   1661.,   1560.,  1507.,   1390., 1330.,  1290.]
linewave_all = [2852.96, 2803.53, 2796.35, 2600.17, 2586.55, 2382.77, 2374.46, 2344.21, 1670.79, $
                1550.78, 1548.20, 1526.71, 1402.77, 1393.76, 1334.53, 1304.86, 1302.17]

if keyword_set(choose) then begin
   lines = (jhusdss_train_lines())[0:5]
   qso_trim = jhusdss_absorber_readin(nmfver, /trim)
   ;; Select 3 examples
   i1 = where(qso_trim.nabs eq 1 and qso_trim.ew[5, 0] gt 2. and qso_trim.zqso gt 1.9, n1)
   i2 = where(qso_trim.nabs eq 2 and qso_trim.ew[5, 1] gt 1., n2)
   i3 = where(qso_trim.nabs eq 3 and qso_trim.ew[5, 2] gt 1., n3)
endif

;iran1 = floor(randomu(seed)*n1)
;iran2 = floor(randomu(seed)*n2)
;iran3 = floor(randomu(seed)*n3)

;;
plate1=389 & fiber1=178 & zabs1=1.97213

;; backup
;; plate1=1490 & fiber1=368


;;
plate2=542 & fiber2=208 & zabs21=1.24324 & zabs22=1.10576

;; backup
;; plate2=1915 & fiber2=445
;; plate2=396 & fiber2=452
;; plate2=390 & fiber2=443

;;
plate3=2594 & fiber3=247 & zabs31=1.74355 & zabs32=1.69355 & zabs33=1.33032

;; backup
;; plate3=572 & fiber3=190
;; plate3=2583 & fiber3=189
;; plate3=1941 & fiber3=183
;; plate3=1274 & fiber3=425
;; plate3=566 & fiber3=76
;; plate3=1298 & fiber3=274
;; plate3=756 & fiber3=164
;; plate3=389 & fiber3=239

if keyword_set(choose) then begin
;; go through the list one by one and select a good one
ichoose1 = -1L
colors = ['dark green', 'yellow', 'magenta']
for j=0L, n1-1L do begin
    plate = qso_trim[i1[j]].plate
    fiber = qso_trim[i1[j]].fiber
    spec = jhusdss_decompose_loadspec(plate, fiber, nmfver)

    x = spec.wave*(1.+spec.z)
    yflux = smooth(spec.flux, 5)
    yresi = smooth(spec.residual, 5)
    ynmf = smooth(spec.nmf_continuum*spec.med_continuum, 5)
    ymed = smooth(spec.med_continuum, 5)

    ytitle = textoidl('Flux (Normalized)')
    xra = [3700, 9500]
    yra = [0, 6]

    pos = [0.10, 0.35, 0.90, 0.90]
    djs_plot, x, yflux, xra=xra, xstyle=1, $
        yra=yra, ystyle=1, position=pos, xtickformat='(A1)', ytitle=ytitle, $
        charsize=charsize, charthick=charthick
    djs_oplot, x, ynmf, color='red', thick=thick
    djs_oplot, x, ynmf, color='red', thick=thick

    ytitle = textoidl('Residual')
    pos = [0.10, 0.10, 0.90, 0.35]
    yra = [-0.1, 1.9]
    djs_plot, x, yresi, xra=xra, xstyle=1, $
        yra=yra, ystyle=1, position=pos, xtitle=xtitle, ytitle=ytitle, $
        charsize=charsize, charthick=charthick, /noerase
    djs_oplot, !x.crange, [1,1], color='blue'

    for jabs=0L, 0L do begin
        zabs = qso_trim[i1[j]].zabs[jabs]
        for iline=0L, n_elements(lines)-1L do begin
            djs_oplot, replicate(lines[iline].wave*(1.+zabs),2), $
                !y.crange[0]+[0.0, 1.0]*(!y.crange[1]-!y.crange[0]), $
                color=colors[jabs], thick=thick, linestyle=2
        endfor
    endfor

    a = 'a'
    read, a
    if a eq 'q' then begin
       print, plate, fiber
       stop
    endif
    if a eq 'c' then begin
       ichoose1 = i1[j]
       break
    endif
endfor
endif

;; init
thick=4
charsize=1.1
charthick=2.5

path=jhusdss_get_path(/nmfqso)+'/'+string(nmfver, format='(I3.3)')
qapath = path+'/QAplots'
if (jhusdss_direxist(qapath) eq 0) then spawn, 'mkdir -p '+qapath
psfile = qapath+'/'+'JHU_examples_simple_'+string(nmfver, format='(I3.3)')+'.ps'

k_print, filename=psfile, axis_char_scale=1.2, xsize=14, ysize=8
  
;; page 1: 1 absorber
    plate = plate1
    fiber = fiber1
    zabs = zabs1
    spec = jhusdss_decompose_loadspec(plate, fiber, nmfver)
    objname = hogg_iau_name(spec.ra, spec.dec, 'SDSS')

    ii = where(spec.ivar gt 0.)
    x = spec.wave*(1.+spec.z)
;   yflux = smooth(spec.flux, 5)
;   yresi = smooth(spec.residual, 5)
;   ynmf = smooth(spec.nmf_continuum, 5)
;   ymed = smooth(spec.med_continuum, 5)
    yflux = spec.flux
    yresi = spec.residual
    ynmf = spec.nmf_continuum
    ymed = spec.med_continuum

    xtitle = textoidl('Observer-frame \lambda (\AA)')
    ytitle = textoidl('Flux (Normalized)')
    xra = [3500, 9490]
    yra = [0, 5]

    pos = [0.10, 0.50, 0.90, 0.85]
    djs_plot, x[ii], yflux[ii], xra=xra, xstyle=1, $
        yra=yra, ystyle=1, position=pos, xtickformat='(A1)', ytitle=ytitle, $
        charsize=charsize, charthick=charthick, $
        xthick=thick, ythick=thick, thick=1, title=objname
    djs_oplot, x, ynmf*ymed, color='red', thick=thick

    ;; overplot windows
    xpoly = [1550.*(1.+spec.z+0.02), 1550.*(1.+spec.z+0.02), 2803.*(1.+spec.z+0.04), 2803.*(1.+spec.z+0.04)]
    ypoly = !y.crange[0]+[1., 0.94, 0.94, 1.0]*(!y.crange[1]-!y.crange[0])
    polyfill, xpoly, ypoly, color=djs_icolor('dark green')
;   xpoly = [1550.*(1.+spec.z+0.02), 1550.*(1.+spec.z+0.02), 1909.*(1.+spec.z-0.03), 1909.*(1.+spec.z-0.03)]
;   ypoly = !y.crange[0]+[1., 0.9, 0.9, 1.0]*(!y.crange[1]-!y.crange[0])
;   polyfill, xpoly, ypoly, color=djs_icolor('dark green')

;   xpoly = [1909.*(1.+spec.z+0.01), 1909.*(1.+spec.z+0.01), 2803.*(1.+spec.z-0.04), 2803.*(1.+spec.z-0.04)]
;   ypoly = !y.crange[0]+[1., 0.9, 0.9, 1.0]*(!y.crange[1]-!y.crange[0])
;   polyfill, xpoly, ypoly, color=djs_icolor('dark green')

;   xpoly = [2803.*(1.+spec.z+0.04), 2803.*(1.+spec.z+0.04), 2803.*(1.+spec.z-0.04), 2803.*(1.+spec.z-0.04)]
;   ypoly = !y.crange[0]+[1., 0.94, 0.94, 1.0]*(!y.crange[1]-!y.crange[0])
;   polyfill, xpoly, ypoly, color=djs_icolor('dark green'), /line_fill, $
;     orientation=45., spacing=0.1

    djs_oplot, [1., 1.]*1550.*(1.+spec.z+0.02), $
        !y.crange[0]+[0.88, 1]*(!y.crange[1]-!y.crange[0]), $
        thick=8, color='dark green', linestyle=0
;   djs_oplot, [1., 1.]*2803.53*(1.+spec.z-0.04), $
;       !y.crange[0]+[0.88, 1]*(!y.crange[1]-!y.crange[0]), $
;       thick=8, color='dark green', linestyle=0
    djs_oplot, [1., 1.]*2803.53*(1.+spec.z+0.04), $
        !y.crange[0]+[0.88, 1]*(!y.crange[1]-!y.crange[0]), $
        thick=8, color='dark green', linestyle=0
    djs_oplot, 2803.53*(1.+spec.z+[-0.04, 0.04]), $
        !y.crange[0]+[0.943, 0.943]*(!y.crange[1]-!y.crange[0]), $
        thick=2, color='dark green', linestyle=0

;   djs_oplot, [1., 1.]*1909.*(1.+spec.z-0.02), $
;       !y.crange[0]+[0.8, 1]*(!y.crange[1]-!y.crange[0]), $
;       thick=8, color='dark green', linestyle=1
;   djs_oplot, [1., 1.]*1909.*(1.+spec.z+0.01), $
;       !y.crange[0]+[0.8, 1]*(!y.crange[1]-!y.crange[0]), $


    djs_xyouts, 1550.*(1.+spec.z-0.10), !y.crange[0]+0.77*(!y.crange[1]-!y.crange[0]), $
        'C IV', charsize=charsize, charthick=charthick, color='dark green'
    djs_xyouts, 1909.*(1.+spec.z-0.06), !y.crange[0]+0.45*(!y.crange[1]-!y.crange[0]), $
        'C III', charsize=charsize, charthick=charthick, color='dark green'
    djs_xyouts, 2800.*(1.+spec.z-0.04), !y.crange[0]+0.35*(!y.crange[1]-!y.crange[0]), $
        'Mg II', charsize=charsize, charthick=charthick, color='dark green'

;   djs_xyouts, 1580.*(1.+spec.z), !y.crange[0]+0.82*(!y.crange[1]-!y.crange[0]), $
;       'Search Window', charsize=charsize, charthick=charthick, color='dark green'
;   djs_xyouts, 2180.*(1.+spec.z), !y.crange[0]+0.82*(!y.crange[1]-!y.crange[0]), $
;       'Search Window', charsize=charsize, charthick=charthick, color='dark green'
    djs_xyouts, 2100.*(1.+spec.z), !y.crange[0]+0.88*(!y.crange[1]-!y.crange[0]), $
        'Search Window', charsize=charsize, charthick=charthick, color='dark green'

;   djs_xyouts, 8100., !y.crange[0]+0.32*(!y.crange[1]-!y.crange[0]), $
;       objname, charsize=charsize, charthick=charthick, color='black'

;   items = ['NMF Continuum', 'Median Continuum']
;   colors = [djs_icolor('red'), djs_icolor('blue')]
;   legend, items, textcolor=colors, colors=colors, $
;       box=0, /right, /top, charsize=charsize, charthick=charthick, $
;       linestyle=[0,0], thick=thick, pspacing=1.6


;   ytitle = textoidl('NMF Residual')
;   pos = [0.10, 0.35, 0.90, 0.60]
;   yra = [0.0, 1.9]
;   djs_plot, x[ii], yflux[ii]/ynmf[ii], xra=xra, xstyle=1, $
;       yra=yra, ystyle=1, position=pos, xtickformat='(A1)', ytitle=ytitle, $
;       charsize=charsize, charthick=charthick, /noerase, $
;       xthick=thick, ythick=thick, thick=1
;   djs_oplot, x[ii], ymed[ii], color='blue', thick=thick+2
;   items = ['NMF Continuum', 'Median Continuum']
;   colors = [djs_icolor('red'), djs_icolor('blue')]
;   legend, items, textcolor=colors, colors=colors, $
;       box=0, /right, /top, charsize=charsize, charthick=charthick, $
;       linestyle=[0,0], thick=thick, pspacing=1.6
;   djs_oplot, !x.crange, [1,1], color='gray', thick=thick, linestyle=2

    ytitle = textoidl('Residual')
    pos = [0.10, 0.15, 0.90, 0.50]
    yra = [0.0, 1.8]
    djs_plot, x[ii], yresi[ii], xra=xra, xstyle=1, $
        yra=yra, ystyle=1, position=pos, xtitle=xtitle, ytitle=ytitle, $
        charsize=charsize, charthick=charthick, /noerase, $
        xthick=thick, ythick=thick, thick=1
    djs_oplot, !x.crange, [1,1], color='gray', thick=thick, linestyle=2

    for jabs=0L, n_elements(zabs)-1L do begin
        zabstmp = zabs[jabs]
;       jhusdss_qaplot_oplotlines, zabstmp, thiscolor='dark green'
        for iline=0L, n_elements(linewave_all)-1L do begin
            djs_oplot, replicate(linewave_all[iline]*(1.+zabstmp),2), $
                !y.crange[0]+[0.72, 0.80]*(!y.crange[1]-!y.crange[0]), $
                color='dark green', thick=thick, linestyle=0
        endfor
        for iline=0L, n_elements(linewave)-1L do begin
            djs_xyouts, linewave[iline]*(1.+zabstmp-0.03), $
                !y.crange[0]+0.85*(!y.crange[1]-!y.crange[0]), $
                linename[iline], color='dark green', charsize=0.8, charthick=1.8
        endfor
    endfor

k_end_print

end
