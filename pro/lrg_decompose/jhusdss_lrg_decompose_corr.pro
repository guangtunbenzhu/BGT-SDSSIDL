;; See jhusdss_create_nmf_basis
;pro jhusdss_lrg_nmf_all, version, path=path

lrgver = 101
outfile = jhusdss_get_path(/fitlrg)+'/'+string(lrgver,format='(i3.3)') $
          +'/'+'lrg_correction_101.fits'

;; -- Load data from fits files
choice_load_data = 0
read,'load data? [1=yes, others=no]: ',choice_load_data
if choice_load_data eq 1 then begin
    lrg0 = jhusdss_lrg_readin()
    allspec = jhusdss_read_alllrgspec(lrgver)
    allresi = jhusdss_read_alllrgspec(lrgver, /residual)
endif

    ;;see /export/scratch1/menard/gz323/SDSS/Garching/garching_lrg_selection.pro
    itrim = where((lrg0.sn_median gt 5. $
              and lrg0.mass gt 10.7 and lrg0.ssfr lt -0.40*(lrg0.mass-11.)-11.4 $
              and lrg0.mass lt 12.1 and lrg0.ssfr gt -13. and lrg0.z lt 0.60) $
               or (lrg0.sn_median gt 5. and lrg0.z gt 0.32999 and lrg0.z lt 0.60), ntrim)

;   tmp_index = lindgen(ntrim/5L)*5L
;   iuse = itrim[tmp_index]
    iuse = itrim
    nuse = n_elements(iuse)

    lrg = lrg0[iuse]
    influx = allresi.residual[iuse,*]
    ;; just to use the mask
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
;   allivar = allflux
;   allivar[*,*] = 1.

    for i=0L, nuse-1L do begin
        wave_begin = out_iwave-in_iwave[i]
        wave_end = out_iwave-in_iwave[i]+ninwave-1L -in_iwave_begin
        allflux[i,wave_begin:wave_end] = influx[i,in_iwave_begin:*]
        allivar[i,wave_begin:wave_end] = inivar[i,in_iwave_begin:*]
    endfor

    jhusdss_composite_engine, allflux, allivar, fmean=fluxmean, fmedian=fluxmedian, fgeomean=fluxgeomean

    outstr = {wave:outwave, fluxmean:fluxmean, fluxmedian:fluxmedian, fluxgeomean:fluxgeomean}
    mwrfits, outstr, outfile, /create
end
