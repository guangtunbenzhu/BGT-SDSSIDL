
makeplot1 = 0b
Coarse = 0b
read, 'Coarse 1/0?: ', Coarse 

nmfver=106
stackpath = jhusdss_get_path(/nmfqso)+'/'+string(nmfver, format='(I3.3)')+'/Stack/'
infile = stackpath + jhusdss_garching_galqso_stack_filename(nmfver, boss=boss)
if (Coarse) then infile=repstr(infile, '.fits', '_coarse.fits')
comp = mrdfits(infile,1)

line_wave = [3934.79, 3969.59]

fix_separation = long((line_wave[1]-line_wave[0])/(line_wave[1]+line_wave[0])/alog(10.)*2E+4)

fit_infile = repstr(infile, '.fits', '_fit.fits')
fit = mrdfits(fit_infile, 1)

ew_nofit_mc = fit.ew_nofit_mc
ew_nofit_2_mc = fit.ew_nofit_2_mc
;  for i=0,11 do begin
;      shuffle_infile = repstr(fit_infile, '_fit.fits', '_shuffle_'+string(i+1, format='(i2.2)')+'_fit.fits')
;      ashuffle = mrdfits(shuffle_infile, 1)
;      ew_nofit_mc = [ew_nofit_mc, ashuffle.ew_nofit_mc]
;      ew_nofit_2_mc = [ew_nofit_2_mc, ashuffle.ew_nofit_2_mc]
;  endfor

sdev_nofit_mc = fltarr(n_elements(fit))
sdev_nofit_2_mc = fltarr(n_elements(fit))
;njump = n_elements(ew_nofit_mc[*,0])
;ijump = lindgen(njump/5)*5
for irp=0,n_elements(fit)-1 do begin
    tmp = moment(ew_nofit_mc[*,irp], sdev=sdev)
;   tmp = moment(ew_nofit_mc[ijump,irp], sdev=sdev)
    sdev_nofit_mc[irp] = sdev
    tmp = moment(ew_nofit_2_mc[*,irp], sdev=sdev)
;   tmp = moment(ew_nofit_2_mc[ijump,irp], sdev=sdev)
    sdev_nofit_2_mc[irp] = sdev
endfor
;print, sdev_nofit_mc[i_indep]

if (makeplot1) then psfile = repstr(infile, '.fits', '_spec_1.ps') else psfile = repstr(infile, '.fits', '_spec_2.ps')
if (Coarse) then psfile = repstr(psfile, '.ps', '_coarse.ps')

;; manually check, independent bins
i_indep = [[0, 1], lindgen(12)*2+2]
if (Coarse) then i_indep = [[0], lindgen(3)*2+1]

thick=8
xthick=8
ythick=8
charsize=1.4
charthick=3

xra=[3800,4100]
xtitle='\lambda (\AA)' 
;xtitle1='(1-R_{\lambda,d})/\sigma(R_{\lambda,d})' 
xtitle1='(1-<R_c>)/\sigma(<R_c>)' 
;xtitle1='(1-!3'+string("303B)+'!XF_\lambda)/\sigma_F' 
;ytitle='Normalized Flux F_\lambda'
ytitle='Normalized Flux <R>';_\lambda'
;title='Single Gaussian Line Profile Measurement'
xra1=[-4.9,6.1]
yra1=[0,1.1]

xx = 10.^(findgen(1000)*0.0024+1.)

plotsym, 0, /fill

;k_print, filename=psfile, axis_char_scale=1.3, xsize=12, ysize=12
k_print, filename=repstr(psfile, '.ps', '_simple.ps'), axis_char_scale=1.3, xsize=12, ysize=8

nrp=4

pos = [0.15, 0.9-0.8/nrp*(0+1), 0.67, 0.9-0.8/nrp*0]
pos1 = [0.70, 0.9-0.8/nrp*(0+1), 0.95, 0.9-0.8/nrp*0]

djs_plot, xra, [1., 1.], psym=10, xra=xra, yra=[1.-0.001, 1.+0.001], $
    position=pos, /nodata, xtickformat='(A1)', ytickformat='(A1)', $ 
    thick=2, xst=5, yst=5, xthick=xthick, ythick=ythick

for i=0L, nrp-1 do begin
    rp = comp[i_indep[i]].rp
    wave = comp[i_indep[i]].wave
    tmp = min(abs(comp[i_indep[i]].wave-line_wave[0]), icaii)
    y = comp[i_indep[i]].fluxgeomean
    med_continuum = median(y, 21, /even)
    y = (y/med_continuum)

    iuse = where(wave gt xra[0] and wave lt xra[1], nuse)
;   iuse_jump = lindgen(nuse/3)*3
;   iuse_hist = iuse[iuse_jump]
    iuse_nocaii = where((wave gt xra[0] and wave lt (line_wave[0]-10.) $
                     or (wave gt (line_wave[0]+10.) and wave lt (line_wave[1]-10.))) $
                     or (wave gt (line_wave[1]+10.) and wave lt (xra[1])), nuse_nocaii)

    ;; first fit

    in_quadra = 0.
    in_slope = 0.
    in_center = line_wave[0]
    in_separation = line_wave[1]-line_wave[0]
    in_lflux = 0.2
    in_ratio = 0.5
    in_sigma = 2.0

    in_intercept = median(1.-y[iuse_nocaii])
    terror = sqrt((moment(1.-y[iuse_nocaii]))[1])

    for iter=0,1 do begin
        jhusdss_lowz_doublet_fit2, wave[iuse], 1.-y[iuse], fltarr(nuse)+1./terror^2, $
           in_quadra, in_slope, in_intercept, in_center, in_separation, in_lflux, in_ratio, in_sigma, $
           quadra=quadra, slope=slope, intercept=intercept, center=center, separation=separation, $
           lflux=lflux, ratio=ratio, sigma=sigma, $
           err_quadra=err_quadra, err_slope=err_slope, err_intercept=err_intercept, err_center=err_center, $
           err_separation=err_separation, err_lflux=err_lflux, err_ratio=err_ratio, err_sigma=err_sigma, $
           maxwidth=maxwidth

        ;; get continuum residuals
        p = [quadra, slope, intercept, center, separation, lflux, ratio, sigma]
        yfit = jhusdss_lowz_doublet_func2(wave, p)
        clevel = quadra*(wave-center)^2+slope*(wave-center)+intercept

        ;; new typical error
        terror = sqrt((moment(1.-y[iuse_nocaii]-clevel[iuse_nocaii]))[1])
    endfor

    newy = y+clevel

    jhusdss_singlet_smooth, 1.-newy, fltarr(n_elements(y))+1./terror^2, $
       outflux=single_newy, outivar=single_newyivar, sigma=1.0, /normalize, factor_norm=factor_norm1
    jhusdss_doublet_smooth, 1.-newy, fltarr(n_elements(y))+1./terror^2, $
       outflux=double_newy, outivar=double_newyivar, sigma=2.0, /normalize, factor_norm=factor_norm2, separation=fix_separation

    pos = [0.15, 0.9-0.8/nrp*(i+1), 0.67, 0.9-0.8/nrp*i]
    pos1 = [0.70, 0.9-0.8/nrp*(i+1), 0.95, 0.9-0.8/nrp*i]
    if (makeplot1) then begin
       djs_plot, wave, 1.-single_newy, xra=xra, xst=1, yra=[1.-6.*terror, 1.+4.*terror], $
           position=pos, /noerase, xtickformat='(A1)', thick=5, xthick=xthick, ythick=ythick
       djs_oplot, line_wave[0]*[1,1.], !y.crange[0]+[0.60,0.95]*(!y.crange[1]-!y.crange[0]), $
           thick=5, color='blue'
       djs_oplot, line_wave[1]*[1,1.], !y.crange[0]+[0.60,0.95]*(!y.crange[1]-!y.crange[0]), $
           thick=5, color='blue'
       ii = where(wave gt (line_wave[0]-6.) and wave lt (line_wave[0]+6.))
       djs_oplot, wave[ii], 1.-single_newy[ii], color='blue', thick=5
       ii = where(wave gt (line_wave[1]-6.) and wave lt (line_wave[1]+6.))
       djs_oplot, wave[ii], 1.-single_newy[ii], color='blue', thick=5
       legend = '<r_p>='+string(rp, format='(i3)')+' kpc'
       djs_xyouts, !x.crange[0]+0.7*(!x.crange[1]-!x.crange[0]), $
                   !y.crange[0]+0.1*(!y.crange[1]-!y.crange[0]), $
                   legend, charsize=1.4, charthick=3
       
   endif else begin
       case i of
           0: begin
              ticks=5 
              tickv=[0.85, 0.90, 0.95, 1.00, 1.05] 
              interval = 0.05
              end
           1: begin
              ticks=4 
              tickv=[0.92, 0.96, 1.00, 1.04] 
              interval = 0.02
              end
           2: begin
              ticks=5 
              tickv=[0.94, 0.96, 0.98, 1.00, 1.02] ;interval = 0.02
              interval = 0.005
              end
           3: begin
              ticks=4 
              tickv=[0.96, 0.98, 1.00, 1.02] ;interval = 0.02
              interval = 0.002
              end
           4: begin
              ticks=4 
              tickv=[0.98, 0.99, 1.00, 1.01] ;interval = 0.01
              interval = 0.01
              end
           5: begin
              ticks=4 
              tickv=[0.990, 0.995, 1.000, 1.005] ;interval = 0.005
              interval = 0.005
              end
           6: begin
              ticks=4 
              tickv=[0.992, 0.996, 1.000, 1.004] ;interval = 0.004
              interval = 0.004
              end
           7: begin
              ticks=5 
              tickv=[0.994, 0.996, 0.998, 1.000, 1.002] ;interval = 0.002
              interval = 0.002
              end
       endcase
       djs_plot, wave, 1.-single_newy, xra=xra, xst=1, yra=[1.-5.*terror, 1.+3.*terror], $
           position=pos, /noerase, xtickformat='(A1)', thick=5, xthick=xthick, ythick=ythick, color='gray', yminor=2, ytickinterval=interval;, ytickv=tickv;, psym=10
       djs_oplot, line_wave[0]*[1,1.], !y.crange[0]+[0.65,0.95]*(!y.crange[1]-!y.crange[0]), $
           thick=5, color='gray'
       djs_oplot, line_wave[1]*[1,1.], !y.crange[0]+[0.65,0.95]*(!y.crange[1]-!y.crange[0]), $
           thick=5, color='gray'
       djs_oplot, wave, 1.+fltarr(n_elements(wave)), color='light gray', thick=4, linestyle=2
       ii = where(wave gt (line_wave[0]-6.) and wave lt (line_wave[0]+6.))
       djs_oplot, wave[ii], 1.-single_newy[ii], color='blue', thick=10
;      polyfill, wave[ii], 1.-single_newy[ii], color=djs_icolor('gray')
       if (i le 1) then begin
          ii = where(wave gt (line_wave[1]-8.) and wave lt (line_wave[1]+8.))
          djs_oplot, wave[ii], 1.-single_newy[ii], color='blue', thick=10
       endif
       legend = '<r_p>='+string(rp, format='(i3)')+' kpc'
       djs_xyouts, !x.crange[0]+0.7*(!x.crange[1]-!x.crange[0]), $
                   !y.crange[0]+0.1*(!y.crange[1]-!y.crange[0]), $
                   legend, charsize=1.4, charthick=3

;      djs_oplot, wave, 1.-double_newy-3.*terror, thick=5, color='light blue';, psym=10
;      ii = where(wave gt (line_wave[0]-6.) and wave lt (line_wave[0]+6.))
;      djs_oplot, wave[ii], 1.-double_newy[ii]-3.*terror, color='blue', thick=10
;      polyfill, wave[ii], 1.-double_newy[ii]-3.*terror, color=djs_icolor('blue')
;      djs_oplot, wave, 1.+fltarr(n_elements(wave))-3.*terror, color='light gray', thick=5, linestyle=2
;      if (i eq 4) then djs_xyouts, 3760, !y.crange[0]-0.3*(!y.crange[1]-!y.crange[0]), $
;         ytitle, charsize=charsize+0.7, charthick=charthick, orientation=90
       if (i eq 2) then djs_xyouts, 3760, !y.crange[0]-0.40*(!y.crange[1]-!y.crange[0]), $
          ytitle, charsize=charsize+0.7, charthick=charthick, orientation=90
;      if (i eq 2) then djs_xyouts, 3760, !y.crange[0]+1.90*(!y.crange[1]-!y.crange[0]), $
;         '(<R_c>)', charsize=charsize+0.7, charthick=charthick, orientation=90, color='blue'
       if (i eq nrp-1) then djs_axis, xaxis=0, xtitle=xtitle, charsize=charsize, charthick=charthick
;      plothist, ew_nofit_2_mc[*, i_indep[i]]/sdev_nofit_2_mc[i_indep[i]], bin=0.1, $

;      if (i eq 0) then djs_xyouts, 3800., 1.15, 'Flux residuals (<R>)', charsize=charsize+0.2, charthick=charthick, color='gray'
;      if (i eq 0) then djs_xyouts, 3800., 1.10, 'Flux residuals convolved with doublet profile (<R_c>)', charsize=charsize+0.2, charthick=charthick, color='blue'

       plothist, double_newy[iuse]/factor_norm2*2./3./sdev_nofit_2_mc[i_indep[i]], bin=0.2, $
           xra=xra1, xst=1, yra=yra1, /fill, color=djs_icolor('light blue'), fcolor=djs_icolor('light blue'), peak=1, $
           position=pos1, /nodata, /noerase, ytickformat='(A1)', xtickformat='(A1)', thick=5, xthick=xthick, ythick=ythick
       djs_oplot, [1,1]*double_newy[icaii]/sdev_nofit_2_mc[i_indep[i]]/factor_norm2*2./3., !y.crange, thick=30, $
           color='blue'
       if (i eq nrp-1) then djs_axis, xaxis=0, xtitle=xtitle1, charsize=charsize, charthick=charthick

   endelse

endfor

k_end_print

end
