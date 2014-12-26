pro jhusdss_qaplot_nmf_example_small, nmfver

linename=['Mg I', 'Mg II', 'Fe II', 'Fe II', 'Al II', 'C IV', 'Si II', 'C II', 'O I'] 
linewave=[2860.,  2800.,   2590.,   2360.,   1661.,   1560.,  1507.,   1340.,  1294.]
linewave_all = [2852.96, 2803.53, 2796.35, 2600.17, 2586.55, 2382.77, 2374.46, 2344.21, 1670.79, $
                1550.78, 1548.20, 1526.71, 1334.53, 1304.86, 1302.17]

lines = (jhusdss_train_lines())[0:5]
qso_trim = jhusdss_absorber_readin(nmfver, /trim)

zall_min = 1.0
zall_max = 2.5
dz=0.1
zbin = jhusdss_make_bins(zall_min, zall_max, dz, nbin=z_nbin)

;; init
thick=4
charsize=1.1
charthick=2.5

vector = [0, 8, 2, 1 3]

path=jhusdss_get_path(/nmfqso)+'/'+string(nmfver, format='(I3.3)')
qapath = path+'/QAplots'
if (jhusdss_direxist(qapath) eq 0) then spawn, 'mkdir -p '+qapath
psfile = qapath+'/'+'JHU_nmf_examples_'+string(nmfver, format='(I3.3)')+'.ps'

k_print, filename=psfile, axis_char_scale=1.2, xsize=12, ysize=8
  
;; page 1: 1 absorber
    for i=0L, z_nbin-1L do begin
        counter, i+1, z_nbin
        iz = where(qso_trim.nabs eq 1 and qso_trim.ew[5,0] gt 1.5 $
               and qso_trim.zqso gt zbin[i].min and qso_trim.zqso lt zbin[i].max $
               and qso_trim.spec_snr_median gt 10., nz)
        if nz gt 0 then j = floor(randomu(seed)*nz)
        index = iz[j]
        plate = qso_trim[index].plate
        fiber = qso_trim[index].fiber

        spec = jhusdss_decompose_loadspec(plate, fiber, nmfver)
        objname = hogg_iau_name(spec.ra, spec.dec, 'SDSS')

        x = spec.wave*(1.+spec.z)
        yflux = smooth(spec.flux, 5)
        yresi = smooth(spec.residual, 5)
        ynmf = smooth(spec.nmf_continuum, 5)
        ymed = smooth(spec.med_continuum, 5)

        xtitle = textoidl('Observer-frame \lambda (\AA)')
        ytitle = textoidl('Flux (Normalized)')
        xra = [3500, 9990]
        yra = [0, 6]

        pos = [0.10, 0.40, 0.90, 0.90]
        djs_plot, x, yflux, xra=xra, xstyle=1, $
            yra=yra, ystyle=1, position=pos, xtickformat='(A1)', ytitle=ytitle, $
            charsize=charsize, charthick=charthick, $
            xthick=thick, ythick=thick, thick=1, /nodata

        ii = where(spec.ivar gt 0., kk)
        dwave = jhusdss_dwave(spec.wave[ii]*(1.+spec.z))
        jj = where(dwave gt 10., ll)
        if (ll gt 0) then begin
           jj = [[0], jj, [kk-1]]
           for j=0L, ll do begin
               djs_oplot, x[ii[jj[j]:jj[j+1]]], yflux[ii[jj[j]:jj[j+1]]], thick=thick
           endfor
        endif else begin
           djs_oplot, x[ii], yflux[ii], thick=thick
        endelse
        
        djs_oplot, x, ynmf, color='red', thick=thick

        ytitle = textoidl('Residual')
        pos = [0.10, 0.10, 0.90, 0.4]
        yra = [-0.4, 1.6]
        djs_plot, x, yresi, xra=xra, xstyle=1, $
            yra=yra, ystyle=1, position=pos, xtitle=xtitle, ytitle=ytitle, $
            charsize=charsize, charthick=charthick, /noerase, $
            xthick=thick, ythick=thick, thick=1, /nodata
        if (ll gt 0) then begin
           for j=0L, ll do begin
               djs_oplot, x[ii[jj[j]:jj[j+1]]], yresi[ii[jj[j]:jj[j+1]]], thick=thick
           endfor
        endif else begin
           djs_oplot, x[ii], yresi[ii], thick=thick
        endelse
        
        djs_oplot, !x.crange, [1,1], color='gray', thick=thick, linestyle=2
    endfor

k_end_print

end
