FUNCTION tick_exponent, axis, index, number

     ; A special case.
     IF number EQ 0 THEN RETURN, '0'

     ; Assuming multiples of 10 with format.
     ex = String(number, Format='(e8.0)')
     pt = StrPos(ex, '.')

     first = StrMid(ex, 0, pt)
     sign = StrMid(ex, pt+2, 1)
     thisExponent = StrMid(ex, pt+3)

     ; Shave off leading zero in exponent
     WHILE StrMid(thisExponent, 0, 1) EQ '0' DO thisExponent = StrMid(thisExponent, 1)

     ; Fix for sign and missing zero problem.
     IF (Long(thisExponent) EQ 0) THEN BEGIN
        sign = ''
        thisExponent = '0'
     ENDIF

     ; Make the exponent a superscript.
     IF sign EQ '-' THEN BEGIN
;       RETURN, first + 'x10!U' + sign + thisExponent + '!N'
        RETURN, '10!U' + sign + thisExponent + '!N'
     ENDIF ELSE BEGIN
;       RETURN, first + 'x10!U' + thisExponent + '!N'
        RETURN, '10!U' + thisExponent + '!N'
     ENDELSE

END

nmfver = 106
;pro jhusdss_qaplot_completeness_gz, nmfver

;; path
if (~keyword_set(path)) then begin
   if (keyword_set(boss)) then begin
       path=jhusdss_get_path(/nmfqso)+'/'+$
          string(nmfver, format='(I3.3)')+'/Wmin_BOSS'
   endif else begin
       path=jhusdss_get_path(/nmfqso)+'/'+$
          string(nmfver, format='(I3.3)')+'/Wmin'
   endelse
endif

infile = path+'/'+jhusdss_montecarlo_wmin_gigantic_filename(nmfver)

choice_load_data = 0
read,'load data? [1=yes, 0=no]: ',choice_load_data
if (choice_load_data eq 1) then begin
    zgrid = mrdfits(infile, 1)
    outstr = mrdfits(infile, 2)
    absorbers = jhusdss_absorber_readin(107)
    qsos = jhusdss_absorber_readin(nmfver, /byqso)
    nz = n_elements(zgrid.zgrid)
    dzgrid = median(zgrid.zgrid[1:nz-2]-zgrid.zgrid[0:nz-1])

    print, "Mask out Ca II"
    ;; mask out Calcium II
;   ca_wave = [3934.7750, 3969.5901]
;   for ica=0L, n_elements(ca_wave)-1L do begin
;       ii = where(abs(zgrid.zgrid - (ca_wave[ica]/2796.35-1.)) le 0.005, nn)
;       outstr.isitcovered[ii] = 0b
;   endfor

    print, "CIV to Mg II"
    ;; CIV to MgII
    wave_limit = 1550.
    mgii_limit = 0.04
    for i=0L, n_elements(qsos)-1L do begin
        ii = where((zgrid.zgrid le wave_limit*(1.+qsos[i].zqso+0.02)/2796.35-1.) $
                or (zgrid.zgrid ge qsos[i].zqso-mgii_limit), nn)
        if (nn gt 0) then outstr[i].isitcovered[ii] = 0b
    endfor

    print, "Histogram..."

    icovered = where(outstr.isitcovered, ncovered)
;   ixcovered = icovered mod nz
;   histcovered = histogram(ixcovered, min=0, max=nz-1)
    total_Delta_z = n_elements(icovered)*dzgrid

    itmp = where(outstr.rewmin_mgii_2796 le 6.0 $
                and outstr.rewmin_mgii_2803 le 6.0/2. $
                and outstr.isitcovered)
    ixtmp = itmp mod nz
    hist60 = histogram(ixtmp, min=0, max=nz-1)

    itmp = where(outstr.rewmin_mgii_2796 le 3.0 $
                and outstr.rewmin_mgii_2803 le 3.0/2. $
                and outstr.isitcovered)
    ixtmp = itmp mod nz
    hist30 = histogram(ixtmp, min=0, max=nz-1)

    itmp = where(outstr.rewmin_mgii_2796 le 2.0 $
                and outstr.rewmin_mgii_2803 le 2.0/2. $
                and outstr.isitcovered)
    ixtmp = itmp mod nz
    hist20 = histogram(ixtmp, min=0, max=nz-1)

    itmp = where(outstr.rewmin_mgii_2796 le 1.0 $ 
                and outstr.rewmin_mgii_2803 le 1.0/2. $
                and outstr.isitcovered)
    ixtmp = itmp mod nz
    hist10 = histogram(ixtmp, min=0, max=nz-1)

    itmp = where(outstr.rewmin_mgii_2796 le 0.6 $
                and outstr.rewmin_mgii_2803 le 0.6/2. $
                and outstr.isitcovered)
    ixtmp = itmp mod nz
    hist06 = histogram(ixtmp, min=0, max=nz-1)

    itmp = where(outstr.rewmin_mgii_2796 le 0.3 $
                and outstr.rewmin_mgii_2803 le 0.3/2. $
                and outstr.isitcovered)
    ixtmp = itmp mod nz
    hist03 = histogram(ixtmp, min=0, max=nz-1)

    ;; ew(2796) bins
    all_ewmin = alog10(0.1)
    all_ewmax = alog10(10.0)
    ew_binsize = alog10(1.5)/4.
    ew_bin = jhusdss_make_bins(all_ewmin, all_ewmax, ew_binsize, nbin=ew_nbin)
    z_path = fltarr(ew_nbin)

    for i=0L, ew_nbin-1L do begin
        counter, i+1, ew_nbin
        ii = where(outstr.rewmin_mgii_2796 le 10.^ew_bin[i].min and $
                   outstr.rewmin_mgii_2803 le 10.^ew_bin[i].min/2. and $
                   outstr.isitcovered, nn)
        z_path[i] = float(nn)*dzgrid
;       ix = ii mod nz
;       histx = histogram(ix, min=0, max=nz-1)
    endfor
endif
stop

;; init
thick=6
charsize=1.3
charthick=3.0

qapath = jhusdss_get_path(/nmfqso)+'/'+string(nmfver, format='(I3.3)')+'/QAplots'
if (jhusdss_direxist(qapath) eq 0) then spawn, 'mkdir -p '+qapath

psfile = qapath+'/Wmin_gz_'+string(nmfver, format='(I3.3)')+'.ps'
k_print, filename=psfile, axis_char_scale=1.4, xsize=8, ysize=9
;  !p.multi = [0, 2, nexam]
;  !x.margin = 0
;  !y.margin = 0

   pos = [0.10, 0.10, 0.90, 0.45]
   xtitle = 'z'
   ytitle = 'Sightline Coverage'
   ytitle1 = 'Sightline Fraction'
   xra = [0.36, 2.3]
   yra = [0, 66000L]
;  yra = [0, 84533L]
   djs_plot, zgrid.zgrid, hist03, xra=xra, xstyle=1, $
           yra=yra, ystyle=1, $
           thick=thick, xthick=thick, ythick=thick, $
           charsize=charsize, charthick=charthick, ytickformat='(A1)', $
           xtitle=xtitle, ytitle='', pos=pos, /nodata

   colors = [djs_icolor('grey'), djs_icolor('red'), djs_icolor('magenta'), djs_icolor('brown')]

   loadct, 5
   for i=0L, 3L do colors[i] = (200.-(i+1)/float(4)*150.)

   djs_oplot, zgrid.zgrid, hist20, thick=thick, linestyle=0, color=colors[3]
   djs_xyouts, 0.63, 6.0E4, "W_0^{min}=2.0 \AA", charsize=charsize+0.0, charthick=charthick, color=colors[3]
   djs_oplot, zgrid.zgrid, hist10, thick=thick, linestyle=0, color=colors[2]
   djs_xyouts, 0.63, 5.1E4, "W_0^{min}=1.0 \AA", charsize=charsize+0.0, charthick=charthick, color=colors[2]
;  djs_oplot, zgrid.zgrid, hist20, thick=thick, linestyle=0
;  djs_xyouts, 0.55, 7.5E4, "W_0^{min}=2.0 \AA", charsize=charsize+0.3, charthick=charthick
;  djs_oplot, zgrid.zgrid, hist30, thick=thick, linestyle=0
   djs_oplot, zgrid.zgrid, hist06, thick=thick, linestyle=0, color=colors[1]
   djs_xyouts, 0.63, 3.2E4, "W_0^{min}=0.6 \AA", charsize=charsize+0.0, charthick=charthick, color=colors[1]
   djs_oplot, zgrid.zgrid, hist03, thick=thick, linestyle=0, color=colors[0]
   djs_xyouts, 0.63, 0.9E4, "W_0^{min}=0.3 \AA", charsize=charsize+0.0, charthick=charthick, color=colors[0]
   djs_axis, yaxis=0, charsize=charsize, charthick=charthick, ytitle=ytitle, yra=yra/1E4
   djs_xyouts, !x.crange[0], !y.crange[1]+0.02*(!y.crange[1]-!y.crange[0]), '\times1e4', charsize=charsize, charthick=charthick
   loadct, 0

;  djs_axis, yaxis=1, ytitle=ytitle1, yra=[0,1], $
;      charsize=charsize, charthick=charthick
;k_end_print

;psfile1 = qapath+'/Wmin_Deltaz_'+string(nmfver, format='(I3.3)')+'.ps'
;pos = [0.2, 0.15, 0.95, 0.95]
 thick = 8
;k_print, filename=psfile1, axis_char_scale=2.0, xsize=9, ysize=7.2
;  !p.multi = [0, 2, nexam]
;  !x.margin = 0
;  !y.margin = 0

   pos = [0.1, 0.58, 0.90, 0.93]
   xtitle = textoidl('W_0^{\lambda2796} (\AA)')
   ytitle = textoidl('\Deltaz')
   ytitle1 = 'Sightline Fraction'
   xra = [0.0, 6.0]
   yra = [0, max(z_path)/0.98]
;  yra = [0, 84533L]
   djs_plot, 10.^ew_bin.min, z_path, xra=xra, xstyle=1, $
           yra=yra, ystyle=1, $
           thick=thick, xthick=thick, ythick=thick, $
           charsize=charsize, charthick=charthick, ytickformat='(A1)', $
           xtitle=xtitle, ytitle='', pos=pos, color='blue', /noerase
   djs_axis, yaxis=0, charsize=charsize, charthick=charthick, ytitle=ytitle, yra=yra/1E4
   djs_xyouts, !x.crange[0], !y.crange[1]+0.02*(!y.crange[1]-!y.crange[0]), '\times1e4', charsize=charsize, charthick=charthick
k_end_print

psfile2 = qapath+'/Wmin_DeltaN_'+string(nmfver, format='(I3.3)')+'.ps'
pos = [0.2, 0.15, 0.95, 0.95]
thick = 8
k_print, filename=psfile2, axis_char_scale=2.0, xsize=9, ysize=7.2
;  !p.multi = [0, 2, nexam]
;  !x.margin = 0
;  !y.margin = 0

   xtitle = textoidl('W_0^{\lambda2796} (\AA)')
   ytitle = textoidl('\Deltaz')
   ytitle1 = 'Sightline Fraction'
   xra = [0.0, 7.0]
   yra = [1, 60000L]
;  yra = [0, 84533L]
   binsize=0.1
   hist2796 = histogram(absorbers.rew_mgii_2796, bin=binsize, min=0.2, max=7)
   hist2796_bin = 0.2+findgen(fix((7-0.2)/binsize)+1)*binsize
   zz_path = interpol(z_path, 10.^ew_bin.min, hist2796_bin)
   hist2796_corrected = hist2796*total_Delta_z/zz_path
   djs_plot, [0.0, hist2796_bin]+0.5*binsize, [yra[0], hist2796], xra=xra, xstyle=1, $
           yra=yra, ystyle=1, /ylog, psym=10, $
           thick=thick, xthick=thick, ythick=thick, $
           charsize=charsize, charthick=charthick, ytickformat='tick_exponent', $
           xtitle=xtitle, ytitle=ytitle, pos=pos, color='black'
   djs_oplot, [0.0, hist2796_bin]+0.5*binsize, [yra[0], hist2796_corrected], psym=10, $
           thick=thick, color='red'
k_end_print

psfile3 = qapath+'/Wmin_Completeness_'+string(nmfver, format='(I3.3)')+'.ps'
pos = [0.2, 0.15, 0.95, 0.95]
thick = 8
k_print, filename=psfile3, axis_char_scale=2.0, xsize=9, ysize=5.
;  !p.multi = [0, 2, nexam]
;  !x.margin = 0
;  !y.margin = 0

   xtitle = textoidl('W_0^{\lambda2796} (\AA)')
   ytitle = textoidl('Completeness f')
   ytitle1 = 'Sightline Fraction'
   xra = [0.0, 6.0]
   yra = [0, 1.]
;  yra = [0, 84533L]
   binsize=0.1
   hist2796 = histogram(absorbers.rew_mgii_2796, bin=binsize, min=0.2, max=7)
   hist2796_bin = 0.2+findgen(fix((7-0.2)/binsize)+1)*binsize
   zz_path = interpol(z_path, 10.^ew_bin.min, hist2796_bin)
   hist2796_corrected = hist2796*total_Delta_z/zz_path
   djs_plot, hist2796_bin+0.5*binsize, zz_path/total_Delta_z, xra=xra, xstyle=1, $
           yra=yra, ystyle=1, $
           thick=thick, xthick=thick, ythick=thick, $
           charsize=charsize, charthick=charthick, $
           xtitle=xtitle, ytitle=ytitle, pos=pos, color='blue'
k_end_print

stop

end
