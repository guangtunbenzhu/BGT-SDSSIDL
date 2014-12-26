;+
; Documentation needed!
; The wrapper
; No need to use so many bins eventually
; we can normalize the spectra using mean flux within a certain range but then re-normzalize them
; to another certain range using mean conversion factor obtained with quasars with common coverage.
;-
pro jhusdss_get_fluxratio, boss=boss, overwrite=overwrite

qsopath = jhusdss_get_path(/qso)

if (keyword_set(boss)) then begin
   filename = 'qso_norm_ratio.fits'
endif else begin
   filename = 'BOSS_qso_norm_ratio.fits'
endelse
outfile = qsopath+'/'+filename

if (file_test(outfile) and ~keyword_set(overwrite)) then begin
   splog, 'File exists. Use Overwrite to overwrite!'
   return
endif else begin
   splog, 'Will write into this file:'
   splog, outfile
endelse

qso = jhusdss_qso_readin(boss=boss)
zmin = [0.4, 0.8, 2.0]
zmax = [1.0, 1.8, 2.8]
outratio = replicate({ratio:fltarr(3)}, n_elements(qso))

minwave = 3800D0
maxwave = 9200D0
loglam = jhusdss_get_loglam(minwave=minwave, maxwave=maxwave, nwave=nwave)
wave = 10.^loglam

for i=1L, 2L do begin
    iz = where(qso.z gt zmin[i] and qso.z lt zmax[i], nz)
;   iz = iz[0:400]
;   nz = n_elements(iz)
    print, nz
    normwave1 = jhusdss_normwave_minmax(option=i+1)
    normwave2 = jhusdss_normwave_minmax(option=i+2)
    jhusdss_load_interp_spec, qso[iz], loglam=loglam, $
            zuse=fltarr(nz), boss=boss, allflux=allflux, allivar=allivar

    iwavemin1 = value_locate(wave, normwave1[0]*(1.+qso[iz].z))
    iwavemax1 = value_locate(wave, normwave1[1]*(1.+qso[iz].z))
    iwavemin2 = value_locate(wave, normwave2[0]*(1.+qso[iz].z))
    iwavemax2 = value_locate(wave, normwave2[1]*(1.+qso[iz].z))

    tmpflux1 = fltarr(nz)
    tmpflux2 = fltarr(nz)
    for j=0L, nz-1L do begin
        tmpflux1[j] = median(allflux[j, iwavemin1[j]:iwavemax1[j]])
        tmpflux2[j] = median(allflux[j, iwavemin2[j]:iwavemax2[j]])
    endfor
    outratio[iz].ratio[i] = tmpflux1/tmpflux2
;   plothist, outratio[iz].ratio[i], bin=0.01, xra=[-5,5]
    print, median(outratio[iz].ratio[i])
;   stop
endfor


mwrfits, outratio, outfile, /create
;; median(ratio) = [0.546970, 0.685592, 0.560925]

end
