pro jhusdss_lowz_lrgqso_caii, nmfver, lrgver

if (n_elements(nmfver) eq 0) then message, 'nmfver required'
if (n_elements(lrgver) eq 0) then message, 'lrgver required'
if (n_elements(nssfr) eq 0) then nssfr=2
if (n_elements(nmass) eq 0) then nmass=2

if (keyword_set(boss)) then begin
   lowzpath=jhusdss_get_path(/nmfqso)+'/'+string(nmfver, format='(I3.3)')+'/Lowz_BOSS'
endif else begin
   lowzpath=jhusdss_get_path(/nmfqso)+'/'+string(nmfver, format='(I3.3)')+'/Lowz'
endelse

infile = lowzpath+'/'+jhusdss_lowz_lrgqso_composite_filename(nmfver)
outfile = lowzpath+'/EW_'+jhusdss_lowz_lrgqso_composite_filename(nmfver)
lowz = mrdfits(infile, 1)
;lowz = reform(nssfr+1, nmass+1)

nrp = n_elements(lowz[0].rp)
str_tmp = {lflux:fltarr(nrp), err_lflux:fltarr(nrp), sigma:fltarr(nrp), err_sigma:fltarr(nrp)}
outstr = replicate(str_tmp, n_elements(lowz))

wave = lowz[0].wave
ierror = where(wave gt 3950.-300. and wave lt 3950.+300., nerror)
iuse = where(wave gt 3850. and wave lt 4050., nuse)
load_dp, /b
xra = [3850., 4050]
;for i=0L, n_elements(lowz)-1L do begin
for i=0L, 0L do begin

   for j=0L, nrp-1L do begin

;  y = lowz[i].fluxgeomean[j,*]
   y = lowz[i].fluxmean[j,*]
   terror = sqrt((moment(y[ierror]))[1])
   yivar = fltarr(nuse)+1./terror^2.

   print, 'typical error= ', terror
   print, 'radius= ', lowz[i].rp[j,*]

   in_slope = 0.
   in_intercept = median(1.-y[ierror])
   in_center = 3934.78
   in_separation = 3969.59-3934.78
   in_lflux = 0.2
   in_ratio = 0.5
   in_sigma = 2.0

   if (total(y[iuse]) eq 0.) then continue
   jhusdss_lowz_doublet_fit, wave[iuse], 1.-y[iuse], fltarr(nuse)+1./terror^2, $
       in_slope, in_intercept, in_center, in_separation, in_lflux, in_ratio, in_sigma, $
       slope=slope, intercept=intercept, center=center, separation=separation, $
       lflux=lflux, ratio=ratio, sigma=sigma, $
       err_slope=err_slope, err_intercept=err_intercept, err_center=err_center, $
       err_separation=err_separation, err_lflux=err_lflux, err_ratio=err_ratio, err_sigma=err_sigma, $
       maxwidth=maxwidth

   p = [slope, intercept, center, separation, lflux, ratio, sigma]
   perror = [err_slope, err_intercept, err_center, err_separation, err_lflux, err_ratio, err_sigma]

   yfit = jhusdss_lowz_doublet_func(wave, p)
   print, p
   print, perror

   outstr[i].lflux[j] = lflux
   outstr[i].err_lflux[j] = err_lflux
   outstr[i].sigma[j] = sigma
   outstr[i].err_sigma[j] = err_sigma

   yra = [1.-5.*terror, 1.+5.*terror]
   charsize = 1.6
   clevel = slope*(lowz.wave-center)+intercept

   jhusdss_doublet_smooth, 1.-clevel-reform(y), reform(yivar), outflux=newy, outivar=newyivar
   newy = 1.-newy-clevel

   djs_plot, wave, smooth(y+clevel, 3), $
       xra=xra, yra=yra, thick=thick, xthick=thick, ythick=thick, $
       xtickformat='(A1)', charsize=charsize, charthick=charthick
   djs_oplot, wave, newy+clevel, $
       thick=thick, color='green'
   djs_oplot, wave, 1.-yfit+clevel, $
       thick=thick, color='red'
   djs_oplot, [1.,1.]*3934.78, !y.crange[0]+[1.0,0.6]*(!y.crange[1]-!y.crange[0]), $
       thick=thick, color='cyan'
   djs_oplot, [1.,1.]*3969.59, !y.crange[0]+[1.0,0.6]*(!y.crange[1]-!y.crange[0]), $
       thick=thick, color='cyan'

   djs_xyouts, !x.crange[0]+0.05*(!x.crange[1]-!x.crange[0]), $
               !y.crange[0]+0.90*(!y.crange[1]-!y.crange[0]), $
               'mass='+string(lowz[i].mass[j], format='(f5.2)'), charsize=charsize, charthick=charthick
   djs_xyouts, !x.crange[0]+0.05*(!x.crange[1]-!x.crange[0]), $
               !y.crange[0]+0.80*(!y.crange[1]-!y.crange[0]), $
               'ssfr='+string(lowz[i].ssfr[j],format='(f6.2)'), charsize=charsize, charthick=charthick
   djs_xyouts, !x.crange[0]+0.05*(!x.crange[1]-!x.crange[0]), $
               !y.crange[0]+0.85*(!y.crange[1]-!y.crange[0]), $
               'W^{\lambda3934}='+string(lflux,format='(f5.3)')+'\pm'+string(err_lflux,format='(f5.3)')+' \AA', $
               charsize=1.6, charthick=charthick
   a = 'a'
   read, a
   if (a eq 'q') then stop
   endfor
endfor

mwrfits, outstr, outfile, /create
lowz=reform(lowz, 3, 3)
outstr=reform(outstr, 3, 3)
djs_plot, lowz[0,0].rp, outstr[0,0].lflux, psym=4, /xlog, /ylog, yra=[0.00005, 2.0], xra=[5,2000], xst=1, yst=1
oploterror, lowz[0,0].rp, outstr[0,0].lflux, outstr[0,0].err_lflux, psym=4, $
     color=djs_icolor('red'), errcolor=djs_icolor('red')
print, outstr[0,0].lflux/outstr[0,0].err_lflux
stop
i = 0
djs_plot, lowz.mass[i], outstr.lflux[i,*], psym=4, /ylog, yra=[7E-3, 8E-2], xra=[9.7, 11.2], xst=1, yst=1
oploterror, lowz[0,*].mass[i], outstr[0,*].lflux[i,*], outstr[0,*].err_lflux[i,*], psym=4, $
     color=djs_icolor('green'), errcolor=djs_icolor('green')
oploterror, lowz[1,*].mass[i], outstr[1,*].lflux[i,*], outstr[1,*].err_lflux[i,*], psym=4, $
     color=djs_icolor('red'), errcolor=djs_icolor('red')
oploterror, lowz[2,*].mass[i], outstr[2,*].lflux[i,*], outstr[2,*].err_lflux[i,*], psym=4, $
     color=djs_icolor('blue'), errcolor=djs_icolor('blue')
stop
end
