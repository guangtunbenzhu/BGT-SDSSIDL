pro jhusdss_qsoevol_composite

;; hstfos qso
hstfos_file = '~/SDATA/Quasars/HSTFOS/hstfos_master_visual.fits'
hstfos = mrdfits(hstfos_file, 1)
nhstfos = n_elements(hstfos)

;; sdss qso
sdss_file = '~/SDATA/SDSS/QSO/HW_dr7qso_newz_absorption_info.fits'
sdss = mrdfits(sdss_file, 1)
nsdss = n_elements(sdss)

index_hstfos = where(hstfos.isgood eq 1L, nhstfos_use)
index_sdss0 = where(sdss.bal_flag eq 0L and sdss.mgii eq 0L and sdss.civ eq 0L and sdss.dla eq 0L, nsdss_use0)
index_sdss = fix(randomu(seed, nsdss_use0/10)*(nsdss_use0-1))
nsdss_use = n_elements(index_sdss)

;; wave
wave_file = '~/SDATA/SDSS/AllInOne/AIO_CommonWave.fits'
wave = (mrdfits(wave_file, 1)).wave
nwave = n_elements(wave)

;######################################
;; sdss spec
;######################################
sdss_allflux = fltarr(nsdss, nwave)
sdss_allivar = fltarr(nsdss, nwave)
print, "Memory used: ", memory(/current)/1024./1024., ' MB'

spec_aa_file = '~/SDATA/SDSS/AllInOne/AIO_QSO_SDSS_DR07_HWzzRestFrame_Wave00450_00900A.fits'
spec_aa = mrdfits(spec_aa_file, 1)
iwave_begin = where(wave eq min(spec_aa.wave), ntest)
if (ntest ne 1) then message, "something's wrong."
iwave_end = where(wave eq max(spec_aa.wave), ntest)
if (ntest ne 1) then message, "something's wrong."
sdss_allflux[spec_aa.index_qso, iwave_begin:iwave_end] = spec_aa.flux
sdss_allivar[spec_aa.index_qso, iwave_begin:iwave_end] = spec_aa.ivar
print, "Memory used: ", memory(/current)/1024./1024., ' MB'
delvar, spec_aa

spec_bb_file = '~/SDATA/SDSS/AllInOne/AIO_QSO_SDSS_DR07_HWzzRestFrame_Wave00900_01800A.fits'
spec_bb = mrdfits(spec_bb_file, 1)
iwave_begin = where(wave eq min(spec_bb.wave), ntest)
if (ntest ne 1) then message, "something's wrong."
iwave_end = where(wave eq max(spec_bb.wave), ntest)
if (ntest ne 1) then message, "something's wrong."
sdss_allflux[spec_bb.index_qso, iwave_begin:iwave_end] = spec_bb.flux
sdss_allivar[spec_bb.index_qso, iwave_begin:iwave_end] = spec_bb.ivar
print, "Memory used: ", memory(/current)/1024./1024., ' MB'
delvar, spec_bb

spec_cc_file = '~/SDATA/SDSS/AllInOne/AIO_QSO_SDSS_DR07_HWzzRestFrame_Wave01800_03600A.fits'
spec_cc = mrdfits(spec_cc_file, 1)
iwave_begin = where(wave eq min(spec_cc.wave), ntest)
if (ntest ne 1) then message, "something's wrong."
iwave_end = where(wave eq max(spec_cc.wave), ntest)
if (ntest ne 1) then message, "something's wrong."
sdss_allflux[spec_cc.index_qso, iwave_begin:iwave_end] = spec_cc.flux
sdss_allivar[spec_cc.index_qso, iwave_begin:iwave_end] = spec_cc.ivar
print, "Memory used: ", memory(/current)/1024./1024., ' MB'
delvar, spec_cc


spec_dd_file = '~/SDATA/SDSS/AllInOne/AIO_QSO_SDSS_DR07_HWzzRestFrame_Wave03600_07200A.fits'
spec_dd = mrdfits(spec_dd_file, 1)
iwave_begin = where(wave eq min(spec_dd.wave), ntest)
if (ntest ne 1) then message, "something's wrong."
iwave_end = where(wave eq max(spec_dd.wave), ntest)
if (ntest ne 1) then message, "something's wrong."
sdss_allflux[spec_dd.index_qso, iwave_begin:iwave_end] = spec_dd.flux
sdss_allivar[spec_dd.index_qso, iwave_begin:iwave_end] = spec_dd.ivar
print, "Memory used: ", memory(/current)/1024./1024., ' MB'
delvar, spec_dd

sdss_allflux_use = sdss_allflux[index_sdss,*]
sdss_allivar_use = sdss_allivar[index_sdss,*]
print, "Memory used: ", memory(/current)/1024./1024., ' MB'
delvar, sdss_allflux
delvar, sdss_allivar

;######################################
;; hstfos spec
;######################################
hstfos_allflux = fltarr(nhstfos, nwave)
hstfos_allivar = fltarr(nhstfos, nwave)
print, "Memory used: ", memory(/current)/1024./1024., ' MB'

spec_aa_file = '~/SDATA/SDSS/AllInOne/AIO_QSO_HSTFOS_NEDzRestFrame_Wave00450_00900A.fits'
spec_aa = mrdfits(spec_aa_file, 1)
iwave_begin = where(wave eq min(spec_aa.wave), ntest)
if (ntest ne 1) then message, "something's wrong."
iwave_end = where(wave eq max(spec_aa.wave), ntest)
if (ntest ne 1) then message, "something's wrong."
hstfos_allflux[spec_aa.index_qso, iwave_begin:iwave_end] = spec_aa.flux
hstfos_allivar[spec_aa.index_qso, iwave_begin:iwave_end] = spec_aa.ivar
print, "Memory used: ", memory(/current)/1024./1024., ' MB'
delvar, spec_aa

spec_bb_file = '~/SDATA/SDSS/AllInOne/AIO_QSO_HSTFOS_NEDzRestFrame_Wave00900_01800A.fits'
spec_bb = mrdfits(spec_bb_file, 1)
iwave_begin = where(wave eq min(spec_bb.wave), ntest)
if (ntest ne 1) then message, "something's wrong."
iwave_end = where(wave eq max(spec_bb.wave), ntest)
if (ntest ne 1) then message, "something's wrong."
hstfos_allflux[spec_bb.index_qso, iwave_begin:iwave_end] = spec_bb.flux
hstfos_allivar[spec_bb.index_qso, iwave_begin:iwave_end] = spec_bb.ivar
print, "Memory used: ", memory(/current)/1024./1024., ' MB'
delvar, spec_bb

spec_cc_file = '~/SDATA/SDSS/AllInOne/AIO_QSO_HSTFOS_NEDzRestFrame_Wave01800_03600A.fits'
spec_cc = mrdfits(spec_cc_file, 1)
iwave_begin = where(wave eq min(spec_cc.wave), ntest)
if (ntest ne 1) then message, "something's wrong."
iwave_end = where(wave eq max(spec_cc.wave), ntest)
if (ntest ne 1) then message, "something's wrong."
hstfos_allflux[spec_cc.index_qso, iwave_begin:iwave_end] = spec_cc.flux
hstfos_allivar[spec_cc.index_qso, iwave_begin:iwave_end] = spec_cc.ivar
print, "Memory used: ", memory(/current)/1024./1024., ' MB'
delvar, spec_cc

hstfos_allflux_use = hstfos_allflux[index_hstfos,*]
hstfos_allivar_use = hstfos_allivar[index_hstfos,*]
print, "Memory used: ", memory(/current)/1024./1024., ' MB'
delvar, hstfos_allflux
delvar, hstfos_allivar

allflux = [hstfos_allflux_use, sdss_allflux_use]
allivar = [hstfos_allivar_use, sdss_allivar_use]

delvar, hstfos_allflux_use
delvar, hstfos_allivar_use
delvar, sdss_allflux_use
delvar, sdss_allivar_use

;######################################
; Now work on what can be used
;######################################
print, "Can we make it? Let's see... "
jhusdss_composite_engine, allflux, allivar, fmean=fluxmean, fmedian=fluxmedian, fsigma=fluxsigma

djs_plot, wave, fluxmean
stop

outfile = '~/SDATA/Quasars/qso_composite.fits'
outstr = {wave:wave, fluxmean:fluxmean, fluxmedian:fluxmedian}
mwrfits, outstr, outfile, /create

end
