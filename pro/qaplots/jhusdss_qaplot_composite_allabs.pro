pro jhusdss_qaplot_composite_allabs, nmfver
    compath = jhusdss_get_path(/nmfqso)+'/'+string(nmfver, format='(I3.3)')+'/Composite'
    filename = 'Absorbers_Composite_Allabs_0.8AA.fits'
    infile = compath + '/' + filename

    comp = mrdfits(infile, 1)
    
    thick=5
    xthick=8
    ythick=8
    charsize=1.3
    charthick=4.0
    xtitle='\lambda (\AA)'

    qapath = jhusdss_get_path(/nmfqso)+'/'+string(nmfver, format='(I3.3)')+'/QAplots'
    psfile = qapath+'/JHU_allabs_composite_'+string(nmfver, format='(I3.3)')+'.ps'

    k_print, filename=psfile, axis_char_scale=1.2, xsize=12, ysize=12
    for i=1, 4 do begin
        xra = [1200., 2200.]+1000.*(i-1)
;       if i eq 4 then xra=[4700., 6100.]
        iwave = where(comp.wave gt xra[0] and comp.wave lt xra[1], nwave)
        if (i eq 1) then iwave = where(comp.wave gt 1230 and comp.wave lt xra[1], nwave)
        case i of
           1: begin
                yra=[0.4, 1.4]
                pos=[0.1, 0.720, 0.9, 0.880]
                noerase=0
              end
           2: begin
                yra=[0.4, 1.4]
                pos=[0.1, 0.520, 0.9, 0.680]
                noerase=1
              end
           3: begin
                yra=[0.95, 1.04]
                pos=[0.1, 0.320, 0.9, 0.480]
                noerase=1
              end
           4: begin
                yra=[0.95, 1.04]
                pos=[0.1, 0.120, 0.9, 0.280]
                noerase=1
              end
        endcase

;       djs_plot, comp.wave[iwave], smooth(comp.fluxmedian[iwave],2), pos=pos, noerase=noerase, $
        djs_plot, comp.wave[iwave], comp.fluxmedian[iwave], pos=pos, noerase=noerase, $
            thick=thick, xthick=xthick, ythick=ythick, xra=xra, yra=yra, xst=1, yst=1, $
            charsize=charsize, charthick=charthick, yticklen=0.01, xticklen=0.04
    endfor
    djs_axis, xaxis=0, xtickformat='(A1)', xtitle=xtitle, $
        charsize=charsize+0.5, charthick=charthick+1
        
    k_end_print

    qapath = jhusdss_get_path(/nmfqso)+'/'+string(nmfver, format='(I3.3)')+'/QAplots'
    psfile = qapath+'/JHU_allabs_composite_onepanel_'+string(nmfver, format='(I3.3)')+'.ps'

    k_print, filename=psfile, axis_char_scale=1.2, xsize=12, ysize=6
      xra = [1200, 5200]
      yra = [0.4, 1.4]
      iwave = where(comp.wave gt 1230 and comp.wave lt xra[1])
      djs_plot, comp.wave[iwave], comp.fluxmedian[iwave], $
            thick=1.5, xthick=xthick, ythick=ythick, xra=xra, yra=yra, xst=1, yst=1, $
            charsize=charsize, charthick=charthick, yticklen=0.01, xtitle=xtitle
    k_end_print

end
