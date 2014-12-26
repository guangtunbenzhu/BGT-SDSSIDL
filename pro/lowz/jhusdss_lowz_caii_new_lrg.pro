;pro jhusdss_lowz_caii_new

lrgver = 101

infile = jhusdss_get_path(/fitlrg)+'/'+string(lrgver,format='(i3.3)') $
       +'/' + 'Lowz_composite_spec_rp_lrg.fits'

comp = mrdfits(infile,1)
nrp = 8

qapath=jhusdss_get_path(/fitlrg)+'/'+string(lrgver, format='(I3.3)')+'/QAplots/'

xra = [3850., 4070]
thick = 5
xthick = 7
ythick = 7 
charthick = 3.5
charsize = 1.1
;title =  'ABSORBER COMPOSITE SPECTRA'
xtitle =  '\lambda (\AA)'

   in_quadra = 0.
   in_slope = 0.
   in_center = 3934.78
   in_separation = 3969.59-3934.78
   in_lflux = 0.2
   in_ratio = 0.5
   in_sigma = 2.0


psfile = qapath+'/LRG_CaII_mean_doublet.ps'

k_print, filename=psfile, axis_char_scale=1.2, xsize=12, ysize=10
  for irp=0L, nrp-1L do begin
      if (irp eq 0) then noerase=0 else noerase=1
      pos1 = [0.1, 0.1+(nrp-1-irp)*0.8/nrp, 0.59, 0.1+(nrp-1-irp+1)*0.8/nrp]
      pos2 = [0.61, 0.1+(nrp-1-irp)*0.8/nrp, 0.9, 0.1+(nrp-1-irp+1)*0.8/nrp]
      rp = comp[irp].rp
      wave = comp[irp].wave
      y = comp[irp].fluxmean
      iwave = where(wave gt 3934.78-3. and wave lt 3934.78+3.)
      jwave = where(wave gt 3969.59-3. and wave lt 3969.59+3.)


      iuse = where(wave gt 3800. and wave lt 4200., nuse)
      iuse_nocaii = where((wave gt 3800. and wave lt 3925.) $
                       or (wave gt 3945. and wave lt 3960.) $
                       or (wave gt 3980. and wave lt 4200.), nuse_nocaii)

      ;; first fit

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

      print, lflux, err_lflux, sigma, err_sigma

      jhusdss_singlet_smooth, 1.-clevel-y, fltarr(n_elements(y))+1./terror^2, $
        outflux=single_newy, outivar=single_newyivar, sigma=1
      jhusdss_doublet_smooth, 1.-clevel-y, fltarr(n_elements(y))+1./terror^2, $
        outflux=newy, outivar=newyivar, sigma=1
      
      single_newy = 1.-single_newy
      newy = 1.-newy

      yra = [1.-3.5*terror, 1.+3.5*terror]
      djs_plot, wave, single_newy, xra=xra, yra=yra, $
          xtickformat='(A1)', pos=pos1, $
          psym=10, thick=thick, xthick=xthick, ythick=ythick, $
          noerase=noerase, color='gray', charsize=charsize, charthick=charthick
      djs_oplot, wave, newy, color='blue', $
          psym=10, thick=thick

      djs_oplot, !x.crange, [1, 1], thick=3, color='gray'
      djs_oplot, wave[iwave], newy[iwave], thick=thick+2, color='red'
;     djs_oplot, wave[jwave], newy[jwave], thick=thick, color='red'
      djs_oplot, [3934.78, 3934.78], !y.crange, linestyle=1, thick=thick
      djs_oplot, [3969.59, 3969.59], !y.crange, linestyle=1, thick=thick
      legend = '<rp>='+string(rp, format='(f6.1)')+' (kpc)'
      djs_xyouts, !x.crange[0]+0.6*(!x.crange[1]-!x.crange[0]), $
                  !y.crange[0]+0.85*(!y.crange[1]-!y.crange[0]), $
                  legend, charsize=charsize, charthick=charthick

      if (irp eq nrp-1) then $
         djs_axis, xaxis=0, xtitle=xtitle, charsize=charsize, $
                   charthick=charthick

;     histxra = [1.-3.*terror, 1.+2.5*terror]
      histxra = [-2.99, 2.5]
      plothist, (newy-1.)/terror, binval, hist, bin=0.1, peak=1, xra=histxra, $
          pos=pos2, /noerase, xtickformat='(A1)', ytickformat='(A1)', $
          thick=thick, xthick=xthick, ythick=ythick
      jj = where(newy gt 1.-2.*terror and newy lt 1.+2.*terror)
      plothist, (newy[jj]-1.)/terror, bin=0.1, peak=1, $
          /overplot, /fill, fcolor=djs_icolor('blue'), thick=thick
      plothist, (newy[iwave]-1.)/terror, bin=0.1, peak=1, /overplot, $
          color=djs_icolor('red'), fcolor=djs_icolor('red'), /fill, thick=thick

      if (irp eq nrp-1) then $
         djs_axis, xaxis=0, xtitle='\Delta (\sigma)', charsize=charsize, $
                   charthick=charthick
  endfor
k_end_print

end
