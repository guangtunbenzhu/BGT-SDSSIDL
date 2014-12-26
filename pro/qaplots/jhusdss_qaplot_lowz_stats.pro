;pro jhusdss_qaplot_lowz_stats

choice_load_data = 0
if (choice_load_data eq 1) then begin
    path = jhusdss_get_path(/garching)
    infile = path + '/gal_info_dr7_v5_2.fit.gz'

    info = mrdfits(infile, 1)
    ii = where(info.z gt 0.001 and info.z le 0.6)
endif

    psfile = 'lowz_hist.ps'
    k_print, filename=psfile, axis_char_scale=1.2, xsize=8, ysize=6
      plothist, info.z, bin=0.01, xra=[0., 0.6], thick=8, xthick=8, ythick=8, $
         xtitle='Redshift', charsize=1.5, charthick=2, peak=1, color=djs_icolor('blue'), $
         ytickformat='(A1)', yra=[0., 1.2], ytitle=''
    k_end_print
    
end
