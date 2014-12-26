;; a naive matching routine
pro jhusdss_detect_absorbers_train_zmatch, pitts, jhu, index1, index2, nmatch

;; There could be multiple matches in jhu to pitts, but not vice versa because min(dz) = 0.01 in JHU
nmax = 7
index1 = lonarr(nmax)-1L
tmpindex2 = lonarr(nmax)-1L
imatch1 = bytarr(nmax)
imatch2 = bytarr(10)
nmatch = 0L
    
for j=0, pitts.nabs-1 do begin
    for k=0, jhu.nabs-1 do begin
        if (abs(jhu.zabs[k]-pitts.zabs[j]) le 0.005) then begin
            index1[nmatch] = j
            tmpindex2[nmatch] = k
            imatch1[j] = 1b
            imatch2[k] = 1b
            nmatch++
        endif
    endfor
;   if ((nmatch eq jhu.nabs) or (nmatch eq pitts.nabs)) then break
endfor

iuse1 = where(imatch1 eq 0b, nuse1)
if nuse1 gt 0 then index1[nmax-nuse1:nmax-1] = iuse1
iuse2 = where(imatch2 eq 0b, nuse2)
if nuse2 gt 0 then begin
   if (nmatch gt 0) then begin
       index2 = [tmpindex2[0:nmatch-1], iuse2]
   endif else begin
       index2 = iuse2
   endelse
endif

end

pro jhusdss_detect_absorbers_train_makeplot, pitts, jhu, nmfver, filename=filename
    
lines = jhusdss_abslines_all(/train)

thick=6
charthick=1.7
charsize=1.2
k_print, filename=filename, axis_char_scale=1.3, xsize=10, ysize=8

    xtitle1 = 'Observer-Frame Wavelength (\AA)'
    xtitle2 = 'QSO Rest-Frame Wavelength (\AA)'
    ytitle1 = 'Flux (Arbitrary Unit)'
    ytitle2 = 'Residual'
    xra = [3700., 9500.]
    colors = ['red', 'dark green', 'blue', 'yellow', 'green', 'cyan', 'cyan', 'cyan', 'cyan', 'cyan']

    for i=0L, n_elements(pitts)-1L do begin
        counter, i+1, n_elements(pitts)

        if (pitts[i].nabs gt 0) then index1 = indgen(pitts[i].nabs)
        if (jhu[i].nabs gt 0) then index2 = indgen(jhu[i].nabs)

        if (jhu[i].nabs gt 0) and (pitts[i].nabs gt 0) then $ 
            jhusdss_detect_absorbers_train_zmatch, pitts[i], jhu[i], index1, index2, nmatch

        ;; load the spectrum
        spec = jhusdss_decompose_loadspec(pitts[i].plate, pitts[i].fiber, $
                  nmfver, error=error)

        ii = where(spec.ivar gt 0., mm)
        djs_plot, spec.wave[ii]*(1.+spec.z), spec.flux[ii], xra=xra, xstyle=1, $
            yra=[0,6], ystyle=1, position=[0.10, 0.55, 0.95, 0.93], xtitle=xtitle1, ytitle=ytitle1
        djs_oplot, spec.wave[ii]*(1.+spec.z), spec.nmf_continuum[ii], color='red', thick=thick

        if (pitts[i].nabs gt 0) then begin
        for jabs=0, pitts[i].nabs-1 do begin
            ztmp = pitts[i].zabs[index1[jabs]]
            colortmp = colors[jabs]
            for iline=0L, n_elements(lines)-1L do djs_oplot, $
                replicate(lines[iline].wave*(1.+ztmp),2), $
                !y.crange[0]+[0.6, 1.0]*(!y.crange[1]-!y.crange[0]), $
                color=colortmp, thick=thick, linestyle=0
        endfor
        endif
        items = ['Absorber_1 z='+string(pitts[i].zabs[index1[0]], format='(f5.3)'), $
                 'Absorber_2 z='+string(pitts[i].zabs[index1[1]], format='(f5.3)'), $
                 'Absorber_3 z='+string(pitts[i].zabs[index1[2]], format='(f5.3)'), $
                 'Absorber_4 z='+string(pitts[i].zabs[index1[3]], format='(f5.3)')]
        textcolors = [djs_icolor('red'), djs_icolor('dark green'), djs_icolor('blue'), djs_icolor('cyan')]
        legend, items, textcolors=textcolors, box=0, /right, charsize=charsize, charthick=charthick

        djs_xyouts, 7800, !y.crange[0]+0.54*(!y.crange[1]-!y.crange[0]), 'Z(QSO)='+$
            string(spec.z, format='(f5.3)'), color='red', charsize=charsize, charthick=charthick
        djs_xyouts, 7800, !y.crange[0]+0.48*(!y.crange[1]-!y.crange[0]), 'Ly\alpha (1216)='+$
            string(1216.*(1.+spec.z), format='(f6.1)'), color='red', charsize=charsize, charthick=charthick
        djs_xyouts, 7800, !y.crange[0]+0.42*(!y.crange[1]-!y.crange[0]), 'CIV (1550)='+$
            string(1550.*(1.+spec.z), format='(f6.1)'), color='red', charsize=charsize, charthick=charthick
        djs_xyouts, 7800, !y.crange[0]+0.36*(!y.crange[1]-!y.crange[0]), 'MgII (2800)='+$
            string(2800.*(1.+spec.z), format='(f7.1)'), color='red', charsize=charsize, charthick=charthick
        djs_xyouts, 7800, !y.crange[0]+0.26*(!y.crange[1]-!y.crange[0]), 'N(abs)='+$
            string(pitts[i].nabs, format='(i1.1)'), color='Dark Green', charsize=charsize, charthick=charthick

        djs_plot, spec.wave[ii], smooth(spec.residual[ii], 5), xra=xra/(1.+spec.z), $
            yra=[0.5, 1.5], xstyle=1, ystyle=1, position=[0.10, 0.10, 0.95, 0.48], /noerase, $
            xtitle=xtitle2, ytitle=ytitle2

        if (jhu[i].nabs gt 0) then begin
        for jabs=0, jhu[i].nabs-1 do begin
            ztmp = jhu[i].zabs[index2[jabs]]
            colortmp = colors[jabs]
            for iline=0L, n_elements(lines)-1L do djs_oplot, $
                replicate(lines[iline].wave*(1.+ztmp)/(1.+spec.z),2), $
                !y.crange[0]+[0.6, 1.0]*(!y.crange[1]-!y.crange[0]), $
                color=colortmp, thick=thick, linestyle=0
        endfor
        endif
        items = ['Absorber_1 z='+string(jhu[i].zabs[index2[0]], format='(f5.3)'), $
                 'Absorber_2 z='+string(jhu[i].zabs[index2[1]], format='(f5.3)'), $
                 'Absorber_3 z='+string(jhu[i].zabs[index2[2]], format='(f5.3)'), $
                 'Absorber_4 z='+string(jhu[i].zabs[index2[3]], format='(f5.3)')]
        textcolors = [djs_icolor('red'), djs_icolor('dark green'), djs_icolor('blue'), djs_icolor('cyan')]
        legend, items, textcolors=textcolors, box=0, /right, charsize=charsize, charthick=charthick

        djs_xyouts, 7800/(1.+spec.z), !y.crange[0]+0.20*(!y.crange[1]-!y.crange[0]), 'Plate='+$
            string(spec.plate, format='(i4.4)')+' Fiber='+string(spec.fiber, format='(i4.4)'), $
            color='red', charsize=charsize, charthick=charthick
        djs_xyouts, 7800/(1.+spec.z), !y.crange[0]+0.10*(!y.crange[1]-!y.crange[0]), 'N(abs)='+$
            string(jhu[i].nabs, format='(i1.1)'), color='Dark Green', charsize=charsize, charthick=charthick
    endfor
    print, ''
k_end_print
end

pro jhusdss_detect_absorbers_train_plot_stats, pitts0, jhu0, nmfver, path=path

    itmp = where(pitts0.nabs le 6); and jhu0.nabs le 6)
    pitts = pitts0[itmp]
    jhu = jhu0[itmp]
    
    nall = n_elements(pitts)
    iabs = where(pitts.nabs ge 1, comp=inoabs, nabs)
    ;; find all matches

    pitts_w2796 = -1.
    pitts_w2803 = -1.
    jhu_w2796 = -1.
    jhu_w2803 = -1.

    no_pitts_w2796 = -1.
    no_pitts_w2803 = -1.
    no_jhu_z = -1.
    no_jhu_zqso = -1.
    no_jhu_w2796 = -1.
    no_jhu_w2803 = -1.

    no_jhu_plate = -1L
    no_jhu_fiber = -1L
    no_pitts_plate = -1L
    no_pitts_fiber = -1L

    all_pitts_plate = -1L
    all_pitts_fiber = -1L
    all_jhu_plate = -1L
    all_jhu_fiber = -1L

    all_pitts_w2796 = -1.
    all_pitts_w2803 = -1.
    all_jhu_w2796 = -1.
    all_jhu_w2803 = -1.

    for i=0L, nabs-1L do begin
        for j=0L, pitts[iabs[i]].nabs-1L do begin
;           if 1550.*(1.+pitts[iabs[i]].redshift)/2803. lt (1.+pitts[iabs[i]].zabs[j]) then begin
            all_pitts_w2796 = [all_pitts_w2796, pitts[iabs[i]].ewall[1,j]]
            all_pitts_w2803 = [all_pitts_w2803, pitts[iabs[i]].ewall[0,j]]
            all_pitts_plate = [all_pitts_plate, pitts[iabs[i]].plate]
            all_pitts_fiber = [all_pitts_fiber, pitts[iabs[i]].fiber]
;           endif
        endfor

        for j=0L, jhu[iabs[i]].nabs-1L do begin
            if (jhu[iabs[i]].zabs[j] gt 0.$
            and (jhu[iabs[i]].magic[0,j] or jhu[iabs[i]].magic[1,j])) then begin
;              if 1550.*(1.+jhu[iabs[i]].zqso)/2803. lt (1.+jhu[iabs[i]].zabs[j]) then begin
               all_jhu_w2796 = [all_jhu_w2796, jhu[iabs[i]].ew[1,j]]
               all_jhu_w2803 = [all_jhu_w2803, jhu[iabs[i]].ew[0,j]]
               all_jhu_plate = [all_jhu_plate, jhu[iabs[i]].plate]
               all_jhu_fiber = [all_jhu_fiber, jhu[iabs[i]].fiber]
;              endif
            endif
        endfor


        jhusdss_detect_absorbers_train_zmatch, pitts[iabs[i]], jhu[iabs[i]], index1, index2, nmatch
        if nmatch gt 0 then begin
           for j=0L, nmatch-1L do begin
;              if 1550.*(1.+pitts[iabs[i]].redshift)/2803. lt (1.+pitts[iabs[i]].zabs[index1[j]]) then begin
               pitts_w2796 = [pitts_w2796, pitts[iabs[i]].ewall[1,index1[j]]]
               pitts_w2803 = [pitts_w2803, pitts[iabs[i]].ewall[0,index1[j]]]
               jhu_w2796 = [jhu_w2796, jhu[iabs[i]].ew[1,index2[j]]]
               jhu_w2803 = [jhu_w2803, jhu[iabs[i]].ew[0,index2[j]]]
;              endif
           endfor
       endif
        if n_elements(index2) gt nmatch then begin
           for j=nmatch, n_elements(index2)-1 do begin
               if (jhu[iabs[i]].zabs[index2[j]] gt 0.$
               and (jhu[iabs[i]].magic[0,index2[j]] or jhu[iabs[i]].magic[1,index2[j]])) then begin
;                  if 1550.*(1.+jhu[iabs[i]].zqso)/2803. lt (1.+jhu[iabs[i]].zabs[index2[j]]) then begin
                   no_jhu_w2796 = [no_jhu_w2796, jhu[iabs[i]].ew[1,index2[j]]]
                   no_jhu_w2803 = [no_jhu_w2803, jhu[iabs[i]].ew[0,index2[j]]]
                   no_jhu_z = [no_jhu_z, jhu[iabs[i]].zabs[index2[j]]]
                   no_jhu_zqso = [no_jhu_zqso, jhu[iabs[i]].zqso]
                   no_jhu_plate = [no_jhu_plate, jhu[iabs[i]].plate]
                   no_jhu_fiber = [no_jhu_fiber, jhu[iabs[i]].fiber]
;                  endif
               endif
           endfor
        endif
        if n_elements(index1) gt nmatch then begin
           for j=nmatch, n_elements(index1)-1 do begin
               if (pitts[iabs[i]].zabs[index1[j]] gt 0.) then begin
;                  if 1550.*(1.+pitts[iabs[i]].redshift)/2803. lt (1.+pitts[iabs[i]].zabs[index1[j]]) then begin
                   no_pitts_w2796 = [no_pitts_w2796, pitts[iabs[i]].ewall[1,index1[j]]]
                   no_pitts_w2803 = [no_pitts_w2803, pitts[iabs[i]].ewall[0,index1[j]]]
                   no_pitts_plate = [no_pitts_plate, pitts[iabs[i]].plate]
                   no_pitts_fiber = [no_pitts_fiber, pitts[iabs[i]].fiber]
;                  endif
               endif
           endfor
        endif

    endfor

    nmatch_all = n_elements(pitts_w2796)-1L
    nabs_all0 = n_elements(all_pitts_w2796)-1L
    nabs_all = total(pitts[iabs].nabs)
    print, nmatch_all, nabs_all0, nabs_all, nmatch_all/nabs_all

    binmin=0.1
    binmax=8.
    binsize=0.2
    all_hist = histogram(all_pitts_w2796, binsize=binsize, min=binmin, max=binmax)
    hist = histogram(pitts_w2796, binsize=binsize, min=binmin, max=binmax)
    bins = findgen(n_elements(hist))*binsize+binmin+binsize/2.
    ii = where(all_hist gt 0)

    all_jhu_hist = histogram(all_jhu_w2796, binsize=binsize, min=binmin, max=binmax)
    jhu_hist = histogram(jhu_w2796, binsize=binsize, min=binmin, max=binmax)

    ;F(missing) stats
    print, total(hist), total(all_hist), float(total(hist))/total(all_hist)
    djs_plot, [binmin-binsize/2., bins[ii]], [0., float(hist[ii])/float(all_hist[ii])], psym=10, yra=[0.0,1.2], ystyle=1
    djs_oplot, !x.crange, [1.00, 1.00]

    ;; w(2796)>4 without match:
    ;; plate=402, fiber=246
    ;; plate=862, fiber=500

;   kk=where(all_pitts_w2796 ge 5.3 and all_pitts_w2796 le 5.5)
;   print, all_pitts_plate[kk]
;   print, all_pitts_fiber[kk]
;   kk=where(all_pitts_w2796 ge 5.7 and all_pitts_w2796 le 5.9)
;   print, all_pitts_plate[kk]
;   print, all_pitts_fiber[kk]
    stop
    ii = where(all_jhu_hist gt 0)
    print, total(jhu_hist), total(all_jhu_hist), float(total(jhu_hist))/total(all_jhu_hist)
    djs_plot, [binmin-binsize/2., bins[ii]], [0., float(jhu_hist[ii])/float(all_jhu_hist[ii])], psym=10, yra=[0.0,1.2], ystyle=1
    djs_oplot, !x.crange, [1.00, 1.00]
    stop

    djs_plot, no_jhu_w2796, no_jhu_w2796/no_jhu_w2803, psym=3, xra=[0.3, 6], yra=[-1, 6], xstyle=1, ystyle=1
;   kk = where(abs(no_jhu_z+1.-1550.*(1+no_jhu_zqso)/2803) gt 3)
;   kk = where(abs(no_jhu_z+1.-1550.*(1+no_jhu_zqso)/2803) le 0.02)
    kk = where(no_jhu_w2796 gt 3)
    print, no_jhu_plate[kk[0:20]]
    print, no_jhu_fiber[kk[0:20]]
    print, no_jhu_z[kk[0:20]]

    stop
    kk = where(no_jhu_z gt 1550.*(no_jhu_zqso+1.)/2803-1.)
    no_jhu_hist = histogram(no_jhu_w2796[kk], binsize=binsize, min=binmin, max=binmax)
    print, total(no_jhu_hist), total(all_jhu_hist), float(total(no_jhu_hist))/total(all_jhu_hist)
    ii = where(all_jhu_hist gt 0)
    djs_plot, [binmin-binsize/2., bins[ii]], [0., float(no_jhu_hist[ii])/float(all_jhu_hist[ii])], psym=10, yra=[0.0,1.0], ystyle=1
    stop
    djs_plot, [binmin-binsize/2., bins[ii]], [0., float(no_jhu_hist[ii])], psym=10, ystyle=1
    stop
    djs_plot, no_jhu_zqso, no_jhu_z, psym=3
;   djs_oplot, !x.crange, 1216.*(1.+!x.crange)/2796.-1
    stop
    x = (1.+no_jhu_zqso)/(1.+no_jhu_z)
    ix = where(x gt 0.1 and x lt 2)
    plothist, x[ix], bin=0.001
    stop
end

;; all pitts.nabs should be the same here
pro jhusdss_detect_absorbers_train_plot_nabs_jhu, nabs, pitts, jhu, nmfver, path=path

itmp = uniq(pitts.nabs)
if n_elements(itmp) ne 1 then message, 'All pitts.nabs should be the same here.'
if pitts[itmp].nabs ne nabs then message, 'nabs you give is not the same as in pitts.'

for i=0,1 do begin
    ewtext = ['Strong', 'Weak']
    case i of
         0: ii = where(pitts.ewall[1,0] ge 1.0, nn)
         1: ii = where(pitts.ewall[1,0] lt 1.0 and pitts.ewall[1,0] ge 0.3, nn)
    endcase

    tmppitts = pitts[ii]
    tmpjhu = jhu[ii]
    splog, 'There are '+string(nn, format='(i4.4)')+' '+string(nabs,format='(i1.1)')+'-absorber systems in this training set.'

    for j=0, nabs+1 do begin
        if (j le nabs) then begin
           icat = where(tmpjhu.nabs eq j, ncat) 
           kmax=j-1
        endif else begin
           icat = where(tmpjhu.nabs ge j, ncat)
           kmax=j-2
        endelse
        if ncat gt 0 then begin
           for k=0, kmax do begin
               subicat = -1L
               for l=0L, ncat-1L do begin
                   if (nabs eq 1 and tmppitts[icat[l]].zabs[0]-1550*(1.+tmppitts[icat[l]].redshift)/2803.+1. le 0.02) then continue
                   jhusdss_detect_absorbers_train_zmatch, tmppitts[icat[l]], tmpjhu[icat[l]], index1, index2, nmatch
                   if (nmatch eq k) then subicat = [subicat, icat[l]]
               endfor
               if (n_elements(subicat) gt 1) then begin
                  thiscat = subicat[1:*]
                  print, 'N(Pitts)='+string(nabs,format='(i1.1)') $
                        +' N(JHU)=' +string(j,format='(i1.1)') $
                        +' N(match)='+string(k,format='(i1.1)')
                  print, n_elements(thiscat), n_elements(icat), float(n_elements(thiscat))/n_elements(icat)
                  filename = path+'/'+ewtext[i]+'_Pitts_'+string(nabs,format='(i1.1)')+'_JHU_' $
                                              +string(j,format='(i1.1)')+'_Match_' $
                                              +string(k, format='(i1.1)')+'.ps'
                  jhusdss_detect_absorbers_train_makeplot, tmppitts[thiscat], tmpjhu[thiscat], nmfver, filename=filename
               endif
           endfor

        endif
    endfor
endfor
end


;; all pitts.nabs should be the same here
pro jhusdss_detect_absorbers_train_plot_nabs, nabs, pitts, jhu, nmfver, path=path

itmp = uniq(pitts.nabs)
if n_elements(itmp) ne 1 then message, 'All pitts.nabs should be the same here.'
if pitts[itmp].nabs ne nabs then message, 'nabs you give is not the same as in pitts.'

for i=0,1 do begin
    ewtext = ['Strong', 'Weak']
    case i of
         0: ii = where(pitts.ewall[1,0] ge 1.0, nn)
         1: ii = where(pitts.ewall[1,0] lt 1.0 and pitts.ewall[1,0] ge 0.3, nn)
    endcase

    tmppitts = pitts[ii]
    tmpjhu = jhu[ii]
    splog, 'There are '+string(nn, format='(i4.4)')+' '+string(nabs,format='(i1.1)')+'-absorber systems in this training set.'

    for j=0, nabs+1 do begin
        if (j le nabs) then begin
           icat = where(tmpjhu.nabs eq j, ncat) 
           kmax=j-1
        endif else begin
           icat = where(tmpjhu.nabs ge j, ncat)
           kmax=j-2
        endelse
        if ncat gt 0 then begin
           for k=0, kmax do begin
               subicat = -1L
               for l=0L, ncat-1L do begin
                   if (nabs eq 1 and tmppitts[icat[l]].zabs[0]-1550*(1.+tmppitts[icat[l]].redshift)/2803.+1. le 0.02) then continue
                   jhusdss_detect_absorbers_train_zmatch, tmppitts[icat[l]], tmpjhu[icat[l]], index1, index2, nmatch
                   if (nmatch eq k) then subicat = [subicat, icat[l]]
               endfor
               if (n_elements(subicat) gt 1) then begin
                  thiscat = subicat[1:*]
                  print, 'N(Pitts)='+string(nabs,format='(i1.1)') $
                        +' N(JHU)=' +string(j,format='(i1.1)') $
                        +' N(match)='+string(k,format='(i1.1)')
                  print, n_elements(thiscat), n_elements(icat), float(n_elements(thiscat))/n_elements(icat)
                  filename = path+'/'+ewtext[i]+'_Pitts_'+string(nabs,format='(i1.1)')+'_JHU_' $
                                              +string(j,format='(i1.1)')+'_Match_' $
                                              +string(k, format='(i1.1)')+'.ps'
                  jhusdss_detect_absorbers_train_makeplot, tmppitts[thiscat], tmpjhu[thiscat], nmfver, filename=filename
               endif
           endfor

        endif
    endfor
endfor
end

pro jhusdss_detect_absorbers_train_plot, nmfver=nmfver

if (n_elements(nmfver) eq 0) then nmfver = jhusdss_get_nmf_version()

;; the Pittsburgh master catalog
pittsfile = jhusdss_get_path(/absorber)+'/MgII/Master_Pitts_Catalog.fits'
pitts = mrdfits(pittsfile, 1)
nspec = n_elements(pitts)

;; the corresponding JHU catalog
path=jhusdss_get_path(/nmfqso)+'/'+string(nmfver, format='(I3.3)')+'/Absorbers'
jhufile = path+'/Pitts_'+jhusdss_absorbers_filename(nmfver, /mgii, boss=boss)
jhu = mrdfits(jhufile, 1)

iconv = where(jhu.isitconvolved eq 1b, nconv)
print, nconv, n_elements(jhu), float(nconv)/n_elements(jhu)

pitts = pitts[iconv]
jhu = jhu[iconv]

lines = jhusdss_abslines_all(/train)
outpath = path+'/Pitts/'
if (jhusdss_direxist(outpath) eq 0) then message, "Can't find the directory."

load_dp, /b
jhusdss_detect_absorbers_train_plot_stats, pitts, jhu, nmfver, path=outpath

i1 = where(pitts.nabs eq 1 and pitts.zabs[0] lt pitts.redshift-0.01, n1)
print, n1
jhusdss_detect_absorbers_train_plot_nabs, 1, pitts[i1], jhu[i1], nmfver, path=outpath

i2 = where(pitts.nabs eq 2 and pitts.zabs[0] lt pitts.redshift-0.01, n2)
print, n2
jhusdss_detect_absorbers_train_plot_nabs, 2, pitts[i2], jhu[i2], nmfver, path=outpath

i3 = where(pitts.nabs eq 3 and pitts.zabs[0] lt pitts.redshift-0.01, n3)
print, n3
jhusdss_detect_absorbers_train_plot_nabs, 3, pitts[i3], jhu[i3], nmfver, path=outpath

end


