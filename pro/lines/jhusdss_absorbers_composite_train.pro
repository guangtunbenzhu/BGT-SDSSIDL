;+
; 01/01/2012
;-
pro jhusdss_absorbers_composite_train, nmfver=nmfver

if (n_elements(nmfver) eq 0) then nmfver = jhusdss_get_nmf_version()
path=jhusdss_get_path(/nmfqso)+'/'+string(nmfver, format='(I3.3)')+'/Absorbers'
filename=path+'/'+'Absorbers_NMF_'+string(nmfver, format='(i3.3)')+'_MgII.fits'
bfilename=path+'/'+'Absorbers_BOSS_NMF_'+string(nmfver, format='(i3.3)')+'_MgII.fits'
a = mrdfits(filename,  1)
b = mrdfits(bfilename,  1)

ii = where((a.zabs[0] le 1.0) and (a.zabs[1] eq 0.) and (a.snr[0,0] gt 4. and a.snr[1,0] gt 4. and (a.snr[2,0] gt 2 or a.snr[3,0] gt 2)), nn)
jj = ii[0:1999]
astack1 = jhusdss_absorbers_composite_engine(a[jj], nmfver=nmfver)
ii = where((a.zabs[0] gt 1.0) and (a.zabs[1] eq 0.) and (a.snr[0,0] gt 4. and a.snr[1,0] gt 4. and (a.snr[2,0] gt 2 or a.snr[3,0] gt 2)), nn)
jj = ii[0:1999]
astack2 = jhusdss_absorbers_composite_engine(a[jj], nmfver=nmfver)

ii = where((b.zabs[0] le 1.0) and (b.zabs[1] eq 0.) and (b.snr[0,0] gt 2. and b.snr[1,0] gt 2. and (b.snr[2,0] gt 1. or b.snr[3,0] gt 1.)), nn)
jj = ii
bstack1 = jhusdss_absorbers_composite_engine(b[jj], nmfver=nmfver, /boss)
ii = where((b.zabs[0] gt 1.0) and (b.zabs[1] eq 0.) and (b.snr[0,0] gt 2. and b.snr[1,0] gt 2. and (b.snr[2,0] gt 1. or b.snr[3,0] gt 1.)), nn)
jj = ii
bstack2 = jhusdss_absorbers_composite_engine(b[jj], nmfver=nmfver, /boss)


save, filename='composite.sav', astack1, astack2, bstack1, bstack2
help, /st, astack1
djs_plot, astack1.wave, astack1.fluxmedian
stop

thick=4
charsize=1.3
charthick=1.8

aflux1 = astack1.fluxmedian
ai1 = where(aflux1 lt 0.9, ni1)
aflux1[ai1] = !values.f_nan
fluxtmp = median(aflux1, 21)
aflux1 = astack1.fluxmedian/fluxtmp

aflux2 = astack2.fluxmedian
ai2 = where(aflux2 lt 0.9, ni2)
aflux2[ai2] = !values.f_nan
fluxtmp = median(aflux2, 21)
aflux2 = astack2.fluxmedian/fluxtmp

lines = jhusdss_abslines_all(/forplot)

k_print, filename='composite.ps', axis_char_scale=1.3, xsize=10, ysize=6
   djs_plot, astack1.wave, aflux1, position=[0.1, 0.60, 0.60, 0.95], $
       xra=[1500, 6000], yra=[0.4, 1.2], xstyle=1, ystyle=1, $
       thick=thick, xthick=thick, ythick=ythick
   djs_xyouts, !x.crange[0]-0.08*(!x.crange[1]-!x.crange[0]), $
               !y.crange[0]+1.07*(!y.crange[1]-!y.crange[0]), $
               'z_{abs}<1', $
               charsize=charsize, charthick=charthick, $
               color='red'
   djs_plot, astack1.wave, aflux1, position=[0.68, 0.60, 0.95, 0.95], $
       xra=[3650, 4000], yra=[0.96, 1.02], xstyle=1, ystyle=1, /noerase, $
       thick=thick, xthick=thick, ythick=ythick, xtickinterval=100.
   xyouts, 3727., $
           !y.crange[0]+0.80*(!y.crange[1]-!y.crange[0]), $
           'OII 3727', charsize=0.6, charthick=0.4, orient=90., $
           color=djs_icolor('red')
   ilines = where(lines.wave gt !x.crange[0] and lines.wave lt !x.crange[1], nlines)
   for i=0L, nlines-1L do $
       xyouts, lines[ilines[i]].wave-1., $
                   !y.crange[0]+0.05*(!y.crange[1]-!y.crange[0]), $
                   lines[ilines[i]].name, charsize=0.6, charthick=0.4, orient=90., $
                   color=djs_icolor('red')

   djs_plot, astack1.wave, aflux2, position=[0.45, 0.1, 0.95, 0.45], $
       xra=[1000, 5000], yra=[0.4, 1.2], xstyle=1, ysytle=1, /noerase, $
       thick=thick, xthick=thick, ythick=ythick
   djs_xyouts, !x.crange[0]-0.1*(!x.crange[1]-!x.crange[0]), $
               !y.crange[0]-0.2*(!y.crange[1]-!y.crange[0]), $
               'Rest-frame Wavelength (\AA)', $
               charsize=charsize, charthick=charthick
   djs_plot, astack1.wave, aflux2, position=[0.1, 0.1, 0.37, 0.45], $
       xra=[1200, 1600], yra=[0.4, 1.2], xstyle=1, ysytle=1, /noerase, $
       thick=thick, xthick=thick, ythick=ythick
   djs_xyouts, !x.crange[0]-0.15*(!x.crange[1]-!x.crange[0]), $
               !y.crange[0]+1.1*(!y.crange[1]-!y.crange[0]), $
               'z_{abs}>1', $
               charsize=charsize, charthick=charthick, $
               color=djs_icolor('red') 
   ilines = where(lines.wave gt !x.crange[0] and lines.wave lt !x.crange[1], nlines)
   for i=0L, nlines-1L do $
       xyouts, lines[ilines[i]].wave-1., $
                   !y.crange[0]+0.05*(!y.crange[1]-!y.crange[0]), $
                   lines[ilines[i]].name, charsize=0.6, charthick=0.4, orient=90., $
                   color=djs_icolor('red') 
k_end_print

end
