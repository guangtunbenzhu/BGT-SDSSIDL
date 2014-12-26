;; See jhusdss_create_nmf_basis
;pro jhusdss_lrg_nmf_all, version, path=path

reject = 0
maxreject = 10L
lrgver = 101
outfile = jhusdss_get_path(/nmflrg)+'/'+string(lrgver, format='(i3.3)')$
          +'/'+jhusdss_lrg_nmf_basis_name(lrgver)


if (file_test(outfile)) then message, 'File already exists!'

;; -- Load data from fits files
choice_load_data = 0
read,'load data? [1=yes, others=no]: ',choice_load_data
if choice_load_data eq 1 then begin
    lrg0 = jhusdss_lrg_readin()
    allspec = jhusdss_read_alllrgspec(lrgver)
endif

    ;;see /export/scratch1/menard/gz323/SDSS/Garching/garching_lrg_selection.pro
    itrim1 = where(lrg0.sn_median gt 5. $
              and lrg0.mass gt 10.7 and lrg0.ssfr lt -0.40*(lrg0.mass-11.)-11.4 $
              and lrg0.mass lt 12.1 and lrg0.ssfr gt -13., ntrim1)

    itrim2 = where(lrg0.sn_median gt 5. and lrg0.z gt 0.32999 and lrg0.z lt 0.6, ntrim2)

    tmp_index = lindgen(ntrim1/10L)*10L
    iuse = [itrim1[tmp_index], itrim2]
    nuse = n_elements(iuse)

    lrg = lrg0[iuse]
    influx = allspec.flux[iuse,*]
    inivar = allspec.ivar[iuse,*]

    inwave = allspec.wave
    ninwave = n_elements(inwave)
    in_iwave = value_locate(inwave, 5000.*(lrg.z+1.))
    in_iwave_begin = value_locate(inwave, 3800.)
    in_iwave = in_iwave - in_iwave_begin

    loglam = jhusdss_get_loglam(minwave=3700./1.6, maxwave=9200.)
    outwave = 10.^loglam
    noutwave = n_elements(outwave)
    out_iwave = value_locate(outwave, 5000.)


    allflux = fltarr(nuse, noutwave)
    allivar = fltarr(nuse, noutwave)

    for i=0L, nuse-1L do begin
        wave_begin = out_iwave-in_iwave[i]
        wave_end = out_iwave-in_iwave[i]+ninwave-1L -in_iwave_begin
        allflux[i,wave_begin:wave_end] = influx[i,in_iwave_begin:*]
        allivar[i,wave_begin:wave_end] = inivar[i,in_iwave_begin:*]
    endfor

    stop
    ;; assign 0 weight to non-positive flux.
    izero = where(allflux le 0., nzero)
    if (nzero gt 0) then allivar[izero]=0.

    splog, "Sanity checking ..."
    jhusdss_nmf_sanity_check, lrg, loglam, flux=allflux, ivar=allivar, $
       spec_index=spec_index, pix_index=pix_index

    splog, "Performing NMF ..."
    jhusdss_nmf_engine, allflux, weight=allivar, n_dimension=n_dimension, $
    eigen_vectors=w, eigen_values=h, maxiters=maxiters

;; reject weird ones, sigma-clipping, 
;; see CREATE_NMF_EIGEN_VECTORS.pro by Brice
if (keyword_set(reject)) then begin
    n_sigma_cut = 5. ;; hardcoded
    iterreject=0L

    totaloutlier = 0L
    repeat begin
       isitoutlier = bytarr((size(allflux))[1])
       ;; w[n_dimension, n_pix], h[n_spec, n_dimension], wh[n_spec, n_pix]
       wh = w##h
       htmp = h
       ;; note we did not normalize spectra when loading, so we need to normalize here temporarily
       ;; Warning, this requires normalization before hand
;      for i=0L, (size(allflux))[1]-1 do htmp[i,*] = h[i,*]/median(wh[i,*])
       for i=0L, n_dimension-1L do begin
           average = median(htmp[*,i])
           sigma = stddev(htmp[*,i])
           ii = where(abs(htmp[*,i]-average) gt sigma*n_sigma_cut, nn)
           if (nn gt 0) then isitoutlier[ii] = 1b ;;Logical OR
       endfor
       ioutliers = where(isitoutlier, comp=index_use, noutliers)
       if (noutliers gt 0) then begin
          totaloutlier = totaloutlier+noutliers
          splog, "I found "+strtrim(string(noutliers),2)+" outliers!"

          ;; chuck outliers, remember to apply it to h as well
          lrg = lrg[index_use]
          allflux = allflux[index_use, *]
          allivar = allivar[index_use, *]
          if (n_elements(index_use) ne (size(h))[1]) then h = h[index_use, *]
    
          ;; sanity check
          jhusdss_nmf_sanity_check, lrg, loglam, flux=allflux, ivar=allivar, $
             spec_index=spec_index, pix_index=pix_index

          ;; reinitialization
          if (n_elements(spec_index) ne (size(h))[1]) then h = h[spec_index, *]
          if (n_elements(pix_index) ne (size(w))[2]) then w = w[*, pix_index]

          ;; ready for nmf
          jhusdss_nmf_engine, allflux, weight=allivar, n_dimension=n_dimension, $
              eigen_vectors=w, eigen_values=h, maxiters=maxiters
       endif
       iterreject++
    endrep until (noutliers eq 0 or iterreject gt maxreject)
    splog, "Total outliers found "+strtrim(string(totaloutlier),2)+" out of "+strtrim(string(noriginal),2)
    splog, "Total rejection iterations: "+strtrim(string(iterreject),2)
endif

mask = (allivar eq 0.)
nobjuse = total(~mask, 1)
basis = {eigen_vectors:w, wave:10.^loglam, nobjuse:nobjuse}

mwrfits, basis, outfile, /create

end
