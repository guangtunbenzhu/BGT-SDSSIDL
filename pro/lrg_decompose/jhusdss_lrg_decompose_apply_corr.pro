;; See jhusdss_create_nmf_basis
;pro jhusdss_lrg_nmf_all, version, path=path

lrgver = 101
parentpath = '/data1/gwln2scratch/menard/gz323/SDSS/LRG/101/AllInOne/'
outfile = parentpath+'DR7_LRG_ALL_CORR.fit'
infile = jhusdss_get_path(/fitlrg)+'/'+string(lrgver,format='(i3.3)') $
          +'/'+'lrg_correction_101.fits'

;index = jhusdss_read_alllrgspec(lrgver, boss=boss, /index)
;wave = jhusdss_read_alllrgspec(lrgver, boss=boss, /wave)
index = mrdfits(parentpath+'DR7_LRG_ALL_INDEX.fit', 1)
wave = mrdfits(parentpath+'DR7_LRG_ALL_WAVE.fit', 1)
corr = mrdfits(infile, 1)

nwave = n_elements(wave.wave)
nobj = n_elements(index)
strtmp = {correction:fltarr(nwave)+1.}
out_corr = replicate(strtmp, nobj)


;   influx = allresi.residual[iuse,*]
    ;; just to use the mask
;   inivar = allspec.ivar[iuse,*]

;outwave = corr.wave
noutwave = nwave
outwave = wave.wave
out_iwave = value_locate(outwave, 5000.*(index.zlrg+1.))
out_iwave_begin = value_locate(outwave, 3800.)
out_iwave = out_iwave - out_iwave_begin

inwave = corr.wave
ninwave = n_elements(inwave)
in_iwave = value_locate(inwave, 5000.)

fluxmedian = corr.fluxmedian
ii = where(fluxmedian eq 0.)
fluxmedian[ii] = 1.

for i=0L, nobj-1L do begin
    counter, i+1, nobj
    if (index[i].zlrg gt 0.6) then continue
    wave_begin = in_iwave-out_iwave[i]
    wave_end = wave_begin+noutwave-1L-out_iwave_begin

    out_corr[i].correction[out_iwave_begin:*] = corr.fluxmedian[wave_begin:wave_end]
endfor

mwrfits, out_corr, outfile, /create
end
