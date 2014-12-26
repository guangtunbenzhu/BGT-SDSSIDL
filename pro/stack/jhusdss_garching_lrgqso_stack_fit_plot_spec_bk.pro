
makeplot1 = 1b

nmfver=106
stackpath = jhusdss_get_path(/nmfqso)+'/'+string(nmfver, format='(I3.3)')+'/Stack/'
infile = stackpath + jhusdss_garching_lrgqso_stack_filename(nmfver, boss=boss)
comp = mrdfits(infile,1)

line_wave = [2796.35, 2803.53, 2853.35]

fit_infile = repstr(infile, '.fits', '_fit.fits')
fit = mrdfits(fit_infile, 1)

if (makeplot1) then psfile = repstr(infile, '.fits', '_spec_1.ps') else psfile = repstr(infile, '.fits', '_spec_2.ps')

;; manually check, independent bins
nrp = 6
i_indep = [[1], lindgen(nrp)*4+2]

thick=8
xthick=8
ythick=8
charsize=1.4
charthick=3

xra = [2700., 2950.]

xtitle='\lambda (\AA)' 
ytitle='Normalized Flux'
;title='Single Gaussian Line Profile Measurement'

xx = 10.^(findgen(1000)*0.0024+1.)

plotsym, 0, /fill

k_print, filename=psfile, axis_char_scale=1.3, xsize=8, ysize=12

pos = [0.15, 0.9-0.8/nrp*(0+1), 0.95, 0.9-0.8/nrp*0]
djs_plot, xra, [1., 1.], psym=10, xra=xra, yra=[1.-0.001, 1.+0.001], $
    position=pos, /nodata, xtickformat='(A1)', ytickformat='(A1)', $ 
    thick=2, xst=5, yst=5, xthick=xthick, ythick=ythick

for i=0L, nrp-1 do begin
    rp = comp[i_indep[i]].rp
    wave = comp[i_indep[i]].wave
    tmp = min(abs(comp[i_indep[i]].wave-line_wave[0]), icaii)
    y = comp[i_indep[i]].fluxmean

    iuse = where(wave gt xra[0] and wave lt xra[1], nuse)
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
       outflux=double_newy, outivar=double_newyivar, sigma=2.0, /normalize, factor_norm=factor_norm2

    pos = [0.15, 0.9-0.8/nrp*(i+1), 0.95, 0.9-0.8/nrp*i]
    if (makeplot1) then begin
       djs_plot, wave, 1.-single_newy, xra=xra, xst=1, yra=[1.-6.*terror, 1.+4.*terror], $
           position=pos, /noerase, xtickformat='(A1)', thick=5, xthick=xthick, ythick=ythick
       djs_oplot, line_wave[0]*[1,1.], !y.crange[0]+[0.60,0.95]*(!y.crange[1]-!y.crange[0]), $
           thick=5, color='red'
       djs_oplot, line_wave[1]*[1,1.], !y.crange[0]+[0.60,0.95]*(!y.crange[1]-!y.crange[0]), $
           thick=5, color='red'
       djs_oplot, line_wave[2]*[1,1.], !y.crange[0]+[0.60,0.95]*(!y.crange[1]-!y.crange[0]), $
           thick=5, color='red'
       ii = where(wave gt (line_wave[0]-6.) and wave lt (line_wave[0]+6.))
       djs_oplot, wave[ii], 1.-single_newy[ii], color='red', thick=5
       ii = where(wave gt (line_wave[1]-6.) and wave lt (line_wave[1]+6.))
       djs_oplot, wave[ii], 1.-single_newy[ii], color='red', thick=5
       ii = where(wave gt (line_wave[2]-6.) and wave lt (line_wave[2]+6.))
       djs_oplot, wave[ii], 1.-single_newy[ii], color='red', thick=5
       legend = '<rp>='+string(rp, format='(i4)')+' Kpc'
       djs_xyouts, !x.crange[0]+0.7*(!x.crange[1]-!x.crange[0]), $
                   !y.crange[0]+0.1*(!y.crange[1]-!y.crange[0]), $
                   legend, charsize=1.4, charthick=3
   endif else begin
       djs_plot, wave, 1.-single_newy, xra=xra, xst=1, yra=[1.-8.*terror, 1.+4.*terror], $
           position=pos, /noerase, xtickformat='(A1)', thick=5, xthick=xthick, ythick=ythick
       djs_oplot, line_wave[0]*[1,1.], !y.crange[0]+[0.65,0.95]*(!y.crange[1]-!y.crange[0]), $
           thick=5, color='red'
       djs_oplot, line_wave[1]*[1,1.], !y.crange[0]+[0.65,0.95]*(!y.crange[1]-!y.crange[0]), $
           thick=5, color='red'
       djs_oplot, line_wave[2]*[1,1.], !y.crange[0]+[0.65,0.95]*(!y.crange[1]-!y.crange[0]), $
           thick=5, color='red'
       ii = where(wave gt (line_wave[0]-6.) and wave lt (line_wave[0]+6.))
       djs_oplot, wave[ii], 1.-single_newy[ii], color='red', thick=5
       ii = where(wave gt (line_wave[1]-6.) and wave lt (line_wave[1]+6.))
       djs_oplot, wave[ii], 1.-single_newy[ii], color='red', thick=5
       ii = where(wave gt (line_wave[2]-6.) and wave lt (line_wave[2]+6.))
       djs_oplot, wave[ii], 1.-single_newy[ii], color='red', thick=5
       legend = '<rp>='+string(rp, format='(i4)')+' Kpc'
       djs_xyouts, !x.crange[0]+0.7*(!x.crange[1]-!x.crange[0]), $
                   !y.crange[0]+0.1*(!y.crange[1]-!y.crange[0]), $
                   legend, charsize=1.4, charthick=3
       djs_oplot, wave, 1.-double_newy-3.5*terror, thick=5, color='blue'
       ii = where(wave gt (line_wave[0]-6.) and wave lt (line_wave[0]+6.))
       djs_oplot, wave[ii], 1.-double_newy[ii]-3.5*terror, color='red', thick=5
       djs_oplot, line_wave[0]*[1,1.], !y.crange[0]+[0.40,0.50]*(!y.crange[1]-!y.crange[0]), $
           thick=5, color='red'

   endelse

endfor
djs_axis, xaxis=0, xtitle=xtitle, charsize=charsize, charthick=charthick

k_end_print

end
