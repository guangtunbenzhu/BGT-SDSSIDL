nmfver = 106

choice_load_data = 0
read,'load data? [1=yes, 0=no]: ',choice_load_data

if choice_load_data eq 1 then begin
   qsoshen = mrdfits('/home/menard/DATA/SDSS/QSO/dr7_bh_Jun25_2010.fits.gz', 1)
   stat = jhusdss_qsostats_readin(nmfver, boss=boss)
   allspec = jhusdss_read_allqsospec(nmfver, /flux, boss=boss)
   allcont = jhusdss_read_allqsospec(nmfver, /continuum, boss=boss)
endif

ntry = 8
for i=0L, ntry-1 do begin
    zmin = 0.6+0.2*i
    zmax = 0.6+0.2*(i+1)
    tmp = execute('i'+string(i,format='(i2.2)')+'=where(stat.zqso gt '+string(zmin, format='(f4.2)')+ ' and stat.zqso lt '+string(zmax, format='(f4.2)')+' and stat.spec_snr_median gt 12.)')
endfor

qapath = jhusdss_get_path(/nmfqso)+'/'+string(nmfver, format='(I3.3)')+'/QAplots/'
psfile = qapath+'Every_Quasar_Is_Different_1.ps'
xra = [1400, 3200]
yra = [8E-1, 8]

charsize = 1.5
charthick = 3.0

xtitle = '\lambda (\AA)'
ytitle = 'Normalized Flux'

k_print, filename=psfile, axis_char_scale=1.3, xsize=10, ysize=6 
  djs_plot, allspec.wave, allcont.continuum[0,*], xtitle=xtitle, ytitle=ytitle, $
      xra=xra, yra=yra, charsize=charsize, charthick=charthick, $
      xst=1, yst=1, thick=6, xthick=6, ythick=6, /nodata, /ylog

  loadct, 5
  for i=0L, ntry-1L do begin
      icolor = (200.-(i+1)/float(ntry)*150.)
      tmp = execute('iran=i'+string(i,format='(i2.2)')+'[0]')
      y = allcont.continuum[iran, *]/allspec.norm_ratio[iran]
      ii = where(y gt 0. and allspec.wave gt 3850. and allspec.wave lt 9000.)
      y = smooth(y, 3);*allspec.zqso[iran]^2
      y = y[ii]
      x = allspec.wave/(1.+allspec.zqso[iran])
      x = x[ii]
      djs_oplot, x, y, thick=4, color=icolor
  endfor
  loadct, 0
k_end_print

end
