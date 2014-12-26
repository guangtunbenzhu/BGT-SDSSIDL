;pro jhusdss_qaplot_qsolrg_residuals_stats

lrgver = 101
qsover = 106

charsize = 1.5
charthick = 2.
thick = 6

choice_load_data = 0
read,'load data? [1=yes]: ',choice_load_data
if choice_load_data eq 1 then begin
    qso = jhusdss_qso_readin()
    qso_spec = jhusdss_read_allqsospec(qsover, boss=boss)
    qso_stats = jhusdss_qsostats_readin(qsover, boss=boss)

    lrg = jhusdss_lrg_readin()
    lrg_spec = jhusdss_read_alllrgspec(lrgver, boss=boss)
    lrg_resi = jhusdss_read_alllrgspec(lrgver, boss=boss, /residual)

    ;; apply lrg correction
    cor_file = jhusdss_get_path(/fitlrg)+'/'+string(lrgver,format='(i3.3)') $
          +'/'+'lrg_correction_'+string(lrgver, format='(i3.3)')+'.fits'
    corr = mrdfits(cor_file, 1)

    ;; apply the correction
    inwave = lrg_spec.wave
    ninwave = n_elements(inwave)
    in_iwave = value_locate(inwave, 5000.*(lrg.z+1.))
    in_iwave_begin = value_locate(inwave, 3800.)
    in_iwave = in_iwave - in_iwave_begin

;   loglam = jhusdss_get_loglam(minwave=3700./1.6, maxwave=9200.)
    outwave = corr.wave
    noutwave = n_elements(outwave)
    out_iwave = value_locate(outwave, 5000.)

    lrg_residual_cor = lrg_resi.residual

    for i=0L, n_elements(lrg.z)-1L do begin
        if (lrg[i].z gt 0.6) then continue
        wave_begin = out_iwave-in_iwave[i]
        wave_end = out_iwave-in_iwave[i]+ninwave-1L -in_iwave_begin
        lrg_residual_cor[i, in_iwave_begin:*] = lrg_resi.residual[i, in_iwave_begin:*] $
                                         - corr.fluxmedian[wave_begin:wave_end] + 1.
    endfor
endif

iqso_trim = where(qso_stats.spec_snr_median gt 5. $
              and qso_stats.med_sdeviation_red gt 0.0 $
              and qso_stats.med_sdeviation_red le 0.1 $
              and qso.z ge 0.1 $
              and qso.z le 4.8, nqso_trim)

ilrg_trim = jhusdss_lrg_trim(lrg, sncut=5.)

yra=[0.0, 1.02E0]
xmin=0.3
xmax=1.7
xra=[xmin, xmax]
;yra=[8.E-0, 5.E7]

;qso
path=jhusdss_get_path(/nmfqso)+'/'+string(qsover, format='(I3.3)')
qapath = path+'/QAplots'
psfile = qapath+'/'+'QSO_residual_stats_'+string(qsover, format='(I3.3)')+'.ps'
k_print, filename=psfile, axis_char_scale=1.2, xsize=8, ysize=8
  qso_residual = qso_spec.residual[iqso_trim, *]
; qso_ivar = qso_spec.ivar[iqso_trim,*]
; qso_iuse = where(qso_ivar ne 0. and qso_residual gt 0.1 and qso_residual lt 2.)
  qso_iuse = where(qso_residual gt xmin and qso_residual lt xmax)

  xtitle='Residual'
  ytitle=textoidl('Number (\Delta=0.002)')
  title=textoidl('QSO Residual Distribution')
  plothist, qso_residual[qso_iuse], bin=0.002, $
     xra=xra, yra=yra, xst=1, yst=1, $
     title=title, xtitle=xtitle, ytitle=ytitle, $
     thick=thick, xthick=xthick, ythick=ythick, $
     charsize=charsize, charthick=charthick, peak=1
  tmp = moment(qso_residual[qso_iuse], sdev=sdev)
  legend = '\sigma='+string(sdev, format='(f5.3)')
  djs_xyouts, !x.crange[0] + 0.8*(!x.crange[1]-!x.crange[0]), $
              !y.crange[0] + 0.9*(!y.crange[1]-!y.crange[0]), $
              legend, charsize=charsize, charthick=charthick
k_end_print

;lrg
path=jhusdss_get_path(/fitlrg)+'/'+string(lrgver, format='(I3.3)')
qapath = path+'/QAplots'
psfile = qapath+'/'+'LRG_residual_stats_'+string(lrgver, format='(I3.3)')+'.ps'
k_print, filename=psfile, axis_char_scale=1.2, xsize=8, ysize=8
  lrg_residual = lrg_residual_cor[ilrg_trim, *]
; lrg_residual = lrg_resi.residual[ilrg_trim, *]
; lrg_ivar = lrg_spec.ivar[ilrg_trim,*]
; lrg_iuse = where(lrg_ivar ne 0. and lrg_residual gt 0.1 and lrg_residual lt 2.)
  lrg_iuse = where(lrg_residual gt xmin and lrg_residual lt xmax)

  xtitle='Residual'
  ytitle=textoidl('Number (\Delta=0.002)')
  title=textoidl('LRG Residual Distribution')
  plothist, lrg_residual[lrg_iuse], bin=0.002, $
     xra=xra, yra=yra, xst=1, yst=1, $
     title=title, xtitle=xtitle, ytitle=ytitle, $
     thick=thick, xthick=xthick, ythick=ythick, $
     charsize=charsize, charthick=charthick, peak=1
  tmp = moment(lrg_residual[lrg_iuse], sdev=sdev)
  legend = '\sigma='+string(sdev, format='(f5.3)')
  djs_xyouts, !x.crange[0] + 0.8*(!x.crange[1]-!x.crange[0]), $
              !y.crange[0] + 0.9*(!y.crange[1]-!y.crange[0]), $
              legend, charsize=charsize, charthick=charthick
k_end_print

end
