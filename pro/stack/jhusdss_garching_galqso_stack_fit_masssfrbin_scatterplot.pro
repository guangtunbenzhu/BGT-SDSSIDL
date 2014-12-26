;pro jhusdss_garching_galqso_stack, nmfver, boss=boss, ivarweight=ivarweight

;read,'BOSS? [1=yes, 0=no]: ',BOSS

;; galaxy redshift range
zgalmin = 0.030
zgalmax = 0.400
;zgalmin = 0.030
;zgalmax = 0.060
massive = 0b
quiescent = 0b

nmfver = 106
ivarweigth = 1b
overwrite = 1b
savespec = 0b
savemean = 1b
snrcut = 2.
sdevcut = 0.10
sigma_cut = 2.

;; ca II
ca_line_wave = [3934.79, 3969.59]
;;ca_line_wave = [3935.00, 3969.00]
ca_xra = [3800., 4300.]
;; Na I
na_line_wave = [5890.00, 5896.00]
na_xra = [5700., 6100.]

ha_line_wave = [6563.00, 6584.00]
ha_xra = [6400., 6800.]

line_wave = ca_line_wave
xra = ca_xra

;; radius bins
;minrad = [0.003, 0.010, 0.020, 0.040, 0.200];, 0.160, 0.320, 0.640];, 1.280, 2.560];, 0.270]
;maxrad = [0.010, 0.020, 0.040, 0.080, 0.300];, 0.320, 0.640, 1.280];, 2.560, 5.120];, 0.810]
;minrad = [0.010, 0.060, 0.120, 0.130, 0.140, 0.150, 0.160];, 0.170, 0.180];, 0.150]
;maxrad = [0.030, 0.120, 0.240, 0.260, 0.280, 0.300, 0.320];, 0.340, 0.360];, 0.300]
;minrad = [0.030, 0.060, 0.120, 0.135];, 0.150]
;maxrad = [0.060, 0.120, 0.240, 0.270];, 0.300]

;minrad = [0.010, 0.050, 0.100, 0.120, 0.150, 0.170];, 0.180];, 0.200];, 0.170, 0.180];, 0.150]
;maxrad = [0.050, 0.100, 0.200, 0.240, 0.300, 0.340];, 0.360];, 0.400];, 0.340, 0.360];, 0.300]

;minrad = [0.03, 0.10, 0.20, 0.25, 0.184, 0.30, 0.45, 0.0.010, 0.050, 0.100, 0.120, 0.150, 0.170];, 0.180];, 0.200];, 0.170, 0.180];, 0.150] ;maxrad = [0.10, 0.20, 0.30, 0.0.45, 0.675, 0.100, 0.200, 0.240, 0.300, 0.340];, 0.360];, 0.400];, 0.340, 0.360];, 0.300]

minrad_tmp = 10.^(alog10(0.020)+findgen(23)*(0.5*alog10(1.5)))
maxrad_tmp = minrad_tmp*1.5
;minrad = [[0.003, 0.010], minrad_tmp]
;maxrad = [[0.010, 0.020], maxrad_tmp]
minrad = [[0.003], minrad_tmp]
maxrad = [[0.020], maxrad_tmp]

;minrad_tmp = 10.^(alog10(0.015)+findgen(23)*(0.5*alog10(1.5)))
;maxrad_tmp = minrad_tmp*1.5
;minrad = [[0.003], minrad_tmp]
;maxrad = [[0.015], maxrad_tmp]

fix_separation = long((line_wave[1]-line_wave[0])/(line_wave[1]+line_wave[0])/alog(10.)*2E+4)

wave_exclude = 4360.

rpmean = fltarr(n_elements(minrad))
wmean = fltarr(n_elements(minrad))

; maxrad = minrad*3.
 print, 'radius (kpc): ', minrad*1E3
 print, 'radius (kpc): ', maxrad*1E3

choice_load_data = 0
read,'load data? [1=yes, 0=no]: ',choice_load_data

stackpath = jhusdss_get_path(/nmfqso)+'/'+string(nmfver, format='(I3.3)')+'/Stack/'
outfile = stackpath + jhusdss_garching_galqso_stack_filename(nmfver, boss=boss)
psfile = repstr(outfile, '.fits', '_masssfr_scatterplot.ps')
psfile1 = repstr(outfile, '.fits', '_masssfr_scatterplot_sfr.ps')
psfile2 = repstr(outfile, '.fits', '_masssfr_scatterplot_sfr1.ps')
psfile3 = repstr(outfile, '.fits', '_masssfr_scatterplot_sfr2.ps')
psfile4 = repstr(outfile, '.fits', '_masssfr_scatterplot_sfr3.ps')

if (file_test(outfile) and ~overwrite) then begin
   splog, 'File already exists, use /overwrite to overwrite' 
   return
endif else begin
   splog, 'Will write into this file: '
   print, psfile
endelse

if choice_load_data eq 1 then begin

   ;; foreground galaxies
   garching_path = jhusdss_get_path(/garching)
   garching_file = garching_path+'/'+'gal_info_dr7_v5_2.fit.gz'
   gal = mrdfits(garching_file, 1)
   uniq_file = garching_path+'/'+'gal_uniq_dr7_v5_2.fits'
   galuniq = mrdfits(uniq_file, 1)

   sfr_file = garching_path+'/'+'gal_totsfr_dr7_v5_2.fits.gz'
   sfr = mrdfits(sfr_file, 1)
   mass_file = garching_path+'/'+'totlgm_dr7_v5_2.fit.gz'
   mass = mrdfits(mass_file, 1)
   ssfr_file = garching_path+'/'+'gal_totspecsfr_dr7_v5_2.fits.gz'
   ssfr = mrdfits(ssfr_file, 1)

   ;; background quasars
;  print, 'Load SDSS DR7'
;  qso = jhusdss_qso_readin(boss=boss)
;  stat = jhusdss_qsostats_readin(nmfver, boss=boss)
;  allspec = jhusdss_read_allqsospec(nmfver, /normresi, boss=boss)

   ;; match
   match0 = jhusdss_galqso_match_readin(boss=boss)

    isub = where(match0.rp_mpc gt 0. $
             and match0.rp_mpc le 0.200, nmatch)
    if nmatch eq 0 then message, "Can't find any pair within the annulus"
    match = match0[isub]
    nmatch = n_elements(match)

;   sdev_red_tmp = stat[match.index_qso].med_sdeviation_red
;   sdev_blue_tmp = stat[match.index_qso].med_sdeviation_blue
;   snr_tmp = stat[match.index_qso].SPEC_SNR_MEDIAN
    ssfr_tmp = ssfr[match.index_gal].avg
    sfr_tmp = sfr[match.index_gal].avg
    mass_tmp = mass[match.index_gal].avg
    zgal = gal[match.index_gal].z
    zuniq = galuniq[match.index_gal].choose
;   zqso = qso[match.index_qso].z
    rp_tmp = match.rp_mpc*1E3

endif

thick=8
xthick=8
ythick=8
charsize=1.3
charthick=3


xx = 10.^(findgen(1500)*0.0024+0.1)
print, minmax(xx)


pos = [0.17, 0.20, 0.90, 0.90]
k_print, filename=psfile, axis_char_scale=1.3, xsize=7, ysize=5

   ;yra=[5.E-4, 8.E-1]
   xra=[7.1,13.0]
   yra=[-8.2, -13.5]
   xtitle=textoidl('log_{10} (M_*/M_\odot)')
   ytitle=textoidl('log_{10} (sSFR yr)')

  ii = where(zuniq and mass_tmp gt 7. and mass_tmp lt 13. and sfr_tmp gt -3. and sfr_tmp lt 3., nn)
; x = mass[ii].avg
; y = ssfr[ii].avg
; yy = sfr[ii].avg
  x = mass_tmp[ii]
  y = ssfr_tmp[ii]
  yy = sfr_tmp[ii]
  ihq = where(x gt 10.6 and yy lt -1.29+0.65*(x-10.)+0.5, nhq)
  ihsf = where(x gt 10.0 and yy gt -1.29+0.65*(x-10.)+0.5, nhsf)
  ilq = where(x lt 10.6 and yy lt -1.29+0.65*(x-10.)+0.5, nlq)
  ilsf = where(x lt 10.0 and yy gt -1.29+0.65*(x-10.)+0.5, nlsf)

  iq = where(yy lt -1.29+0.65*(x-10.)+0.5)
  isf = where(yy gt -1.29+0.65*(x-10.)+0.5)
  print, 'high-mass quiescent', median(x[ihq])
  print, 'high-mass star-forming', median(x[ihsf])
  print, 'low-mass quiescent', median(x[ilq])
  print, 'low-mass star-forming', median(x[ilsf])
  print, 'all', median(x)
  print, 'quiescent', median(x[iq])
  print, 'star-forming', median(x[isf])

  flevels = [0.30, 0.60, 0.85]
  hogg_scatterplot, x, y, xra=xra, yra=yra, $
       xnpix=45L, ynpix=45L, exponent=0.50, $
       thick=thick, xthick=thick, ythick=thick, $
       levels=flevels, ytitle=ytitle, xst=1, yst=1,$
       charsize=charsize, charthick=charthick, xtitle=xtitle, $
       outliers=0, outsymsize=0.1, outcolor=djs_icolor('gray')

   djs_oplot, !x.crange, -1.29+0.65*(!x.crange-10.)+0.5-!x.crange, color='dark green', thick=thick 
   djs_oplot, 10.7*[1,1], [-11.0, !y.crange[1]], color='red', thick=thick
   djs_oplot, 10.1*[1,1], [!y.crange[0],-10.8], color='blue', thick=thick

;  djs_xyouts, 11.2, -13.0, 'log_{10}<M_*>=11.3', charsize=charsize-0.1, charthick=charthick-1, color='red'
;  djs_xyouts, 11.2, -12.65, 'log_{10}<M_{CaII}>=3.4', charsize=charsize-0.1, charthick=charthick-1, color='red'
;  djs_xyouts, 11.2, -13.0, 'High-mass Quiescent (HQ)', charsize=charsize-0.3, charthick=charthick-1, color='red'
;  djs_xyouts, 11.2, -12.65, '<M_*>=10.9', charsize=charsize-0.3, charthick=charthick-1, color='red'
;  djs_xyouts, 11.2, -12.30, '<M_{CaII}>=3.68', charsize=charsize-0.3, charthick=charthick-1, color='red'
   djs_xyouts, 11.2, -13.0, 'High-mass', charsize=charsize-0.3, charthick=charthick-1, color='dark green'
   djs_xyouts, 11.2, -12.65, 'Quiescent', charsize=charsize-0.3, charthick=charthick-1, color='red'

;  djs_xyouts, 11.2, -9.5, 'log_{10}<M_*>=10.8', charsize=charsize-0.1, charthick=charthick-1, color='blue'
;  djs_xyouts, 11.2, -9.1, 'log_{10}<M_{CaII}>=3.6', charsize=charsize-0.1, charthick=charthick-1, color='blue'

   djs_xyouts, 11.0, -9.5, 'High-mass', charsize=charsize-0.3, charthick=charthick-1, color='dark green'
   djs_xyouts, 11.0, -9.1, 'Star-forming', charsize=charsize-0.3, charthick=charthick-1, color='blue'

;  djs_xyouts, 8.2, -12.7, 'log_{10}<M_*/M_\odot>=10.7', charsize=charsize-0.1, charthick=charthick-1, color='red'
;  djs_xyouts, 8.2, -12.35, 'log_{10}<M_{CaII}/M_\odot>=3.0', charsize=charsize-0.1, charthick=charthick-1, color='red'

   djs_xyouts, 8.2, -12.7, 'Low-mass', charsize=charsize-0.3, charthick=charthick-1, color='magenta'
   djs_xyouts, 8.2, -12.35, 'Quiescent', charsize=charsize-0.3, charthick=charthick-1, color='red'

;  djs_xyouts, 8.0, -8.9, 'log_{10}<M_*>=10.0', charsize=charsize-0.1, charthick=charthick-1, color='blue'
;  djs_xyouts, 8.0, -8.55, 'log_{10}<M_{CaII}>=3.2', charsize=charsize-0.1, charthick=charthick-1, color='blue'

   djs_xyouts, 7.6, -8.9, 'Low-mass', charsize=charsize-0.3, charthick=charthick-1, color='magenta'
   djs_xyouts, 7.6, -8.55, 'Star-forming', charsize=charsize-0.3, charthick=charthick-1, color='blue'

k_end_print

k_print, filename=psfile1, axis_char_scale=1.3, xsize=7, ysize=5

xra=[7.1,13.0]
yra=[-3.0, 2.0]
xtitle=textoidl('log_{10} (M_*/M_\odot)')
ytitle=textoidl('log_{10} (SFR/M_\odot yr^{-1})')

  flevels = [0.30, 0.60, 0.85]
  hogg_scatterplot, x, yy, xra=xra, yra=yra, $
       xnpix=45L, ynpix=45L, exponent=0.50, $
       thick=thick, xthick=thick, ythick=thick, $
       levels=flevels, ytitle=ytitle, xst=1, yst=1,$
       charsize=charsize, charthick=charthick, xtitle=xtitle, $
       outliers=0, outsymsize=0.1, outcolor=djs_icolor('gray')

   djs_oplot, !x.crange, -1.29+0.65*(!x.crange-10.)+0.5, color='dark green', thick=thick 
   djs_oplot, 10.1*[1,1], [-0.7, !y.crange[1]], color='blue', thick=thick
   djs_oplot, 10.7*[1,1], [!y.crange[0],-0.35], color='red', thick=thick
;  djs_oplot, !x.crange, [-0.5,-0.5], color='cyan blue', thick=thick, linestyle=2
;  djs_oplot, !x.crange, [-1.0,-1.0], color='cyan blue', thick=thick, linestyle=2

   djs_xyouts, 11.2, -2.0, 'High-mass', charsize=charsize-0.3, charthick=charthick-1, color='dark green'
   djs_xyouts, 11.2, -2.3, 'Quiescent', charsize=charsize-0.3, charthick=charthick-1, color='red'

   djs_xyouts, 10.9, 1.6, 'High-mass', charsize=charsize-0.3, charthick=charthick-1, color='dark green'
   djs_xyouts, 10.9, 1.3, 'Star-forming', charsize=charsize-0.3, charthick=charthick-1, color='blue'

   djs_xyouts, 8.2, -2.3, 'Low-mass', charsize=charsize-0.3, charthick=charthick-1, color='magenta'
   djs_xyouts, 8.2, -2.6, 'Quiescent', charsize=charsize-0.3, charthick=charthick-1, color='red'

   djs_xyouts, 7.6, 1.2, 'Low-mass', charsize=charsize-0.3, charthick=charthick-1, color='magenta'
   djs_xyouts, 7.6, 0.9, 'Star-forming', charsize=charsize-0.3, charthick=charthick-1, color='blue'

k_end_print

k_print, filename=psfile2, axis_char_scale=1.3, xsize=7, ysize=5

xra=[7.1,13.0]
yra=[-3.0, 2.0]
xtitle=textoidl('log_{10} (M_*/M_\odot)')
ytitle=textoidl('log_{10} (SFR/M_\odot yr^{-1})')

  flevels = [0.30, 0.60, 0.85]
  hogg_scatterplot, x, yy, xra=xra, yra=yra, $
       xnpix=45L, ynpix=45L, exponent=0.50, $
       thick=thick, xthick=thick, ythick=thick, $
       levels=flevels, ytitle=ytitle, xst=1, yst=1,$
       charsize=charsize, charthick=charthick, xtitle=xtitle, $
       outliers=0, outsymsize=0.1, outcolor=djs_icolor('gray')

   djs_oplot, !x.crange, [-1.0,-1.0], color='cyan blue', thick=thick, linestyle=2
;  djs_oplot, !x.crange, -1.29+0.65*(!x.crange-10.)+0.5, color='dark green', thick=thick 
;  djs_oplot, 10.1*[1,1], [-0.7, !y.crange[1]], color='blue', thick=thick
;  djs_oplot, 10.7*[1,1], [!y.crange[0],-0.35], color='red', thick=thick

;  djs_xyouts, 11.2, -2.0, 'High-mass', charsize=charsize-0.3, charthick=charthick-1, color='dark green'
;  djs_xyouts, 11.2, -2.3, 'Quiescent (HQ)', charsize=charsize-0.3, charthick=charthick-1, color='red'

;  djs_xyouts, 10.9, 1.6, 'High-mass', charsize=charsize-0.3, charthick=charthick-1, color='dark green'
;  djs_xyouts, 10.9, 1.3, 'Star-forming (HSF)', charsize=charsize-0.3, charthick=charthick-1, color='blue'

;  djs_xyouts, 8.2, -2.3, 'Low-mass', charsize=charsize-0.3, charthick=charthick-1, color='magenta'
;  djs_xyouts, 8.2, -2.6, 'Quiescent (LQ)', charsize=charsize-0.3, charthick=charthick-1, color='red'

;  djs_xyouts, 7.6, 1.2, 'Low-mass', charsize=charsize-0.3, charthick=charthick-1, color='magenta'
;  djs_xyouts, 7.6, 0.9, 'Star-forming (LSF)', charsize=charsize-0.3, charthick=charthick-1, color='blue'

k_end_print

k_print, filename=psfile3, axis_char_scale=1.3, xsize=7, ysize=5

xra=[7.1,13.0]
yra=[-3.0, 2.0]
xtitle=textoidl('log_{10} (M_*/M_\odot)')
ytitle=textoidl('log_{10} (SFR/M_\odot yr^{-1})')

  flevels = [0.30, 0.60, 0.85]
  hogg_scatterplot, x, yy, xra=xra, yra=yra, $
       xnpix=45L, ynpix=45L, exponent=0.50, $
       thick=thick, xthick=thick, ythick=thick, $
       levels=flevels, ytitle=ytitle, xst=1, yst=1,$
       charsize=charsize, charthick=charthick, xtitle=xtitle, $
       outliers=0, outsymsize=0.1, outcolor=djs_icolor('gray')

;  djs_oplot, !x.crange, -1.29+0.65*(!x.crange-10.)+0.5, color='dark green', thick=thick 
;  djs_oplot, 10.1*[1,1], [-0.7, !y.crange[1]], color='blue', thick=thick
;  djs_oplot, 10.7*[1,1], [!y.crange[0],-0.35], color='red', thick=thick
;  djs_oplot, !x.crange, [-0.5,-0.5], color='cyan blue', thick=thick, linestyle=2
;  djs_oplot, !x.crange, [-1.0,-1.0], color='cyan blue', thick=thick, linestyle=2

;  djs_xyouts, 11.2, -2.0, 'High-mass', charsize=charsize-0.3, charthick=charthick-1, color='dark green'
;  djs_xyouts, 11.2, -2.3, 'Quiescent', charsize=charsize-0.3, charthick=charthick-1, color='red'

;  djs_xyouts, 10.9, 1.6, 'High-mass', charsize=charsize-0.3, charthick=charthick-1, color='dark green'
;  djs_xyouts, 10.9, 1.3, 'Star-forming', charsize=charsize-0.3, charthick=charthick-1, color='blue'

;  djs_xyouts, 8.2, -2.3, 'Low-mass', charsize=charsize-0.3, charthick=charthick-1, color='magenta'
;  djs_xyouts, 8.2, -2.6, 'Quiescent', charsize=charsize-0.3, charthick=charthick-1, color='red'

;  djs_xyouts, 7.6, 1.2, 'Low-mass', charsize=charsize-0.3, charthick=charthick-1, color='magenta'
;  djs_xyouts, 7.6, 0.9, 'Star-forming', charsize=charsize-0.3, charthick=charthick-1, color='blue'

k_end_print

k_print, filename=psfile4, axis_char_scale=1.3, xsize=7, ysize=5

xra=[7.1,13.0]
yra=[2.0, -3.0]
xtitle=textoidl('log_{10} (M_*/M_\odot)')
ytitle=textoidl('log_{10} (SFR/M_\odot yr^{-1})')

  flevels = [0.30, 0.60, 0.85]
  hogg_scatterplot, x, yy, xra=xra, yra=yra, $
       xnpix=45L, ynpix=45L, exponent=0.50, $
       thick=thick, xthick=thick, ythick=thick, $
       levels=flevels, ytitle=ytitle, xst=1, yst=1,$
       charsize=charsize, charthick=charthick, xtitle=xtitle, $
       outliers=0, outsymsize=0.1, outcolor=djs_icolor('gray')

;  djs_oplot, !x.crange, -1.29+0.65*(!x.crange-10.)+0.5, color='dark green', thick=thick 
;  djs_oplot, 10.1*[1,1], [-0.7, !y.crange[1]], color='blue', thick=thick
;  djs_oplot, 10.7*[1,1], [!y.crange[0],-0.35], color='red', thick=thick
;  djs_oplot, !x.crange, [-0.5,-0.5], color='cyan blue', thick=thick, linestyle=2
;  djs_oplot, !x.crange, [-1.0,-1.0], color='cyan blue', thick=thick, linestyle=2

;  djs_xyouts, 11.2, -2.0, 'High-mass', charsize=charsize-0.3, charthick=charthick-1, color='dark green'
;  djs_xyouts, 11.2, -2.3, 'Quiescent', charsize=charsize-0.3, charthick=charthick-1, color='red'

;  djs_xyouts, 10.9, 1.6, 'High-mass', charsize=charsize-0.3, charthick=charthick-1, color='dark green'
;  djs_xyouts, 10.9, 1.3, 'Star-forming', charsize=charsize-0.3, charthick=charthick-1, color='blue'

;  djs_xyouts, 8.2, -2.3, 'Low-mass', charsize=charsize-0.3, charthick=charthick-1, color='magenta'
;  djs_xyouts, 8.2, -2.6, 'Quiescent', charsize=charsize-0.3, charthick=charthick-1, color='red'

;  djs_xyouts, 7.6, 1.2, 'Low-mass', charsize=charsize-0.3, charthick=charthick-1, color='magenta'
;  djs_xyouts, 7.6, 0.9, 'Star-forming', charsize=charsize-0.3, charthick=charthick-1, color='blue'

k_end_print



end
