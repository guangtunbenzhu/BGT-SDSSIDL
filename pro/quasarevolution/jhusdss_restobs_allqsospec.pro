;+
; load and interpolate *ALL* quasar spectra
; output at both restframe and observer framer, with common central wavelength grid for each 
;-
pro jhusdss_restobs_allqsospec, boss=boss, dr12=dr12, overwrite=overwrite, $
       doobserver=doobserver, dorest=dorest


;; qsopath
parent_path='~/SDATA/SDSS/AllInOne'
outfile_base = parent_path+'/AIO_QSO_SDSS_DR07'
wave_outfile = parent_path+'/AIO_CommonWave.fits'

;; no header files, use explicit file names
 observer_outfile = outfile_base+'_ObserverFrame_Wave03650_10400A.fits'
rest_outfile_base = outfile_base+'_HWzzRestFrame.fits'

qsofile = '~/SDATA/SDSS/QSO/HW_dr7qso_newz.fits'
qso = mrdfits(qsofile, 1)
nqso = n_elements(qso)

if (file_test(wave_outfile)) then begin
   splog, "Warning: wavelength grid already exists. Use the existing one."
   wave = (mrdfits(wave_outfile, 1)).wave
   loglam = alog10(wave)
endif else begin
   loglam = jhusdss_get_loglam(minwave=448., maxwave=10404.)
   wave = 10.^loglam
   mwrfits, {wave:wave}, wave_outfile, /create
endelse
nwave = n_elements(loglam)

delta_z = 0.1
nzbin = fix((max(qso.z)+0.101)/delta_z)
print, "Memory used: ", memory(/current)/1024./1024., ' MB'

;; #######################
;; observer frame
if keyword_set(doobserver) then begin
;; #######################
print, "I am now working on the observer frame."

if (file_test(observer_outfile) and ~keyword_set(overwrite)) then begin
   splog, 'File exists. Use Overwrite to overwrite!'
   return
endif else begin
   splog, 'Will write into these files:'
   print, observer_outfile 
endelse

observer_wave_range = [3650., 10400.]
observer_loc = value_locate(wave, observer_wave_range)
observer_loglam = loglam[observer_loc[0]:observer_loc[1]+1]
observer_nwave = n_elements(observer_loglam)
observer_outstr = {ra:qso.ra, dec:qso.dec, z:qso.z, $
                   mjd:qso.mjd, plate:qso.plate, fiber:qso.fiber, $
                   wave:wave[observer_loc[0]:observer_loc[1]+1], $
                   flux:fltarr(nqso, observer_nwave), $
                   ivar:fltarr(nqso, observer_nwave)}

observer_allflux_out = fltarr(nqso, observer_nwave)
observer_allivar_out = fltarr(nqso, observer_nwave)
print, "Memory used: ", memory(/current)/1024./1024., ' MB'

;; observer frame
for ibin=0L, nzbin-1L do begin
    zbegin = ibin*delta_z
    zend = (ibin+1)*delta_z
    print, 'redshift = ', zbegin, zend
    iz = where(qso.z gt zbegin and qso.z le zend, nz)
    if (nz gt 0) then begin
        ;; observer frame
        zuse = fltarr(nz)
        jhusdss_load_interp_spec, qso[iz], loglam=observer_loglam, $
                zuse=zuse, boss=boss, dr12=dr12, allflux=allflux, allivar=allivar
        observer_allflux_out[iz, *] = allflux
        observer_allivar_out[iz, *] = allivar
    endif
endfor

observer_outstr.flux = observer_allflux_out
observer_outstr.ivar = observer_allivar_out
mwrfits, observer_outstr, observer_outfile, /create
delvar, observer_allflux_out
delvar, observer_allvar_out
delvar, observer_outstr

;; #######################
endif ;; observer frame
;; #######################

;; #######################
;; rest frame
if keyword_set(dorest) then begin
;; #######################
print, "I am now working on the rest frame."

wave_range_aa = [450., 900.]
rest_loc_aa = value_locate(wave, wave_range_aa)
rest_outfile_aa = repstr(rest_outfile_base, '.fits', '_Wave00450_00900A.fits')
wave_range_bb = [900., 1800.]
rest_loc_bb = value_locate(wave, wave_range_bb)
rest_outfile_bb = repstr(rest_outfile_base, '.fits', '_Wave00900_01800A.fits')
wave_range_cc = [1800., 3600.]
rest_loc_cc = value_locate(wave, wave_range_cc)
rest_outfile_cc = repstr(rest_outfile_base, '.fits', '_Wave01800_03600A.fits')
wave_range_dd = [3600., 7200.]
rest_loc_dd = value_locate(wave, wave_range_dd)
rest_outfile_dd = repstr(rest_outfile_base, '.fits', '_Wave03600_07200A.fits')

if (file_test(rest_outfile_aa) and ~keyword_set(overwrite)) then begin
   splog, 'File exists. Use Overwrite to overwrite!'
   return
endif else begin
   splog, 'Will write into these files:'
   print, rest_outfile_aa
   print, rest_outfile_bb
   print, rest_outfile_cc
   print, rest_outfile_dd
endelse

;; two images, one for flux, one for ivar
rest_allflux_out = fltarr(nqso, nwave)
rest_allivar_out = fltarr(nqso, nwave)
print, "Memory used: ", memory(/current)/1024./1024., ' MB'

;; rest frame
for ibin=0L, nzbin-1L do begin
    zbegin = ibin*delta_z
    zend = (ibin+1)*delta_z
    print, 'redshift = ', zbegin, zend
    iz = where(qso.z gt zbegin and qso.z le zend, nz)
    if (nz gt 0) then begin
        ;; rest frame
        rest_loc = value_locate(wave, [3650./(1.+zend), 10400./(1.+zbegin)])
        rest_loglam = loglam[rest_loc[0]:rest_loc[1]+1]
        jhusdss_load_interp_spec, qso[iz], loglam=rest_loglam, $
                boss=boss, dr12=dr12, allflux=allflux, allivar=allivar
        rest_allflux_out[iz, rest_loc[0]:rest_loc[1]+1] = allflux
        rest_allivar_out[iz, rest_loc[0]:rest_loc[1]+1] = allivar
    endif
endfor

;; rest frame
;; 450-900 AA
iaa =  where(qso.z gt (3650./wave_range_aa[1]-1.-0.001) and qso.z le (10400./wave_range_aa[0]-1.+0.001))
rest_outstr_aa = {index_qso:iaa, ra:qso[iaa].ra, dec:qso[iaa].dec, z: qso[iaa].z, $
                  mjd:qso[iaa].mjd, plate:qso[iaa].plate, fiber:qso[iaa].fiber, $
                  wave:wave[rest_loc_aa[0]:rest_loc_aa[1]+1], $
                  flux:rest_allflux_out[iaa,rest_loc_aa[0]:rest_loc_aa[1]+1], $
                  ivar:rest_allivar_out[iaa,rest_loc_aa[0]:rest_loc_aa[1]+1]}
mwrfits, rest_outstr_aa, rest_outfile_aa, /create
delvar, rest_outstr_aa

;; 900-1800 AA
ibb =  where(qso.z gt (3650./wave_range_bb[1]-1.-0.001) and qso.z le (10400./wave_range_bb[0]-1.+0.001))
rest_outstr_bb = {index_qso:ibb, ra:qso[ibb].ra, dec:qso[ibb].dec, z: qso[ibb].z, $
                  mjd:qso[ibb].mjd, plate:qso[ibb].plate, fiber:qso[ibb].fiber, $
                  wave:wave[rest_loc_bb[0]:rest_loc_bb[1]+1], $
                  flux:rest_allflux_out[ibb,rest_loc_bb[0]:rest_loc_bb[1]+1], $
                  ivar:rest_allivar_out[ibb,rest_loc_bb[0]:rest_loc_bb[1]+1]}
mwrfits, rest_outstr_bb, rest_outfile_bb, /create
delvar, rest_outstr_bb

;; 1800 - 3600 AA
icc =  where(qso.z gt (3650./wave_range_cc[1]-1.-0.001) and qso.z le (10400./wave_range_cc[0]-1.+0.001))
rest_outstr_cc = {index_qso:icc, ra:qso[icc].ra, dec:qso[icc].dec, z: qso[icc].z, $
                  mjd:qso[icc].mjd, plate:qso[icc].plate, fiber:qso[icc].fiber, $
                  wave:wave[rest_loc_cc[0]:rest_loc_cc[1]+1], $
                  flux:rest_allflux_out[icc,rest_loc_cc[0]:rest_loc_cc[1]+1], $
                  ivar:rest_allivar_out[icc,rest_loc_cc[0]:rest_loc_cc[1]+1]}
mwrfits, rest_outstr_cc, rest_outfile_cc, /create
delvar, rest_outstr_cc

;; 3600 - 7200 AA
idd =  where(qso.z gt (3650./wave_range_dd[1]-1.-0.001) and qso.z le (10400./wave_range_dd[0]-1.+0.001))
rest_outstr_dd = {index_qso:idd, ra:qso[idd].ra, dec:qso[idd].dec, z: qso[idd].z, $
                  mjd:qso[idd].mjd, plate:qso[idd].plate, fiber:qso[idd].fiber, $
                  wave:wave[rest_loc_dd[0]:rest_loc_dd[1]+1], $
                  flux:rest_allflux_out[idd,rest_loc_dd[0]:rest_loc_dd[1]+1], $
                  ivar:rest_allivar_out[idd,rest_loc_dd[0]:rest_loc_dd[1]+1]}
mwrfits, rest_outstr_dd, rest_outfile_dd, /create
delvar, rest_outstr_dd

;; #######################
endif ;; rest frame
;; #######################

end
