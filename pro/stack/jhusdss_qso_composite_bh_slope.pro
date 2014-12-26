nmfver = 106

choice_load_data = 0
read,'load data? [1=yes, 0=no]: ',choice_load_data

if choice_load_data eq 1 then begin
   qsoshen = mrdfits('/home/menard/DATA/SDSS/QSO/dr7_bh_Jun25_2010.fits.gz', 1)
   stat = jhusdss_qsostats_readin(nmfver, boss=boss)
   allspec = jhusdss_read_allqsospec(nmfver, /flux, boss=boss)
   allcont = jhusdss_read_allqsospec(nmfver, /continuum, boss=boss)

strtmp = {logbh:-99., alpha:-99.}
bhall = replicate(strtmp, n_elements(stat))
spherematch, stat.ra, stat.dec, qsoshen.ra, qsoshen.dec, 1./3600., m1, m2
bhall[m1].logbh = qsoshen[m2].logbh_mgii_s10
bhall[m1].alpha = qsoshen[m2].alpha_mgii

endif

wavemin=1000.D0
wavemax=7000.D0
loglam = jhusdss_get_loglam(minwave=wavemin, maxwave=wavemax)
wave = 10.^loglam
noutwave = n_elements(wave)
pos_2800 = value_locate(wave, 2800.)

strtmp = {logbh:-99., alpha:-99., wave:wave, fluxmedian:fltarr(noutwave)}
composite_bh_alpha = replicate(strtmp, 25)

ninwave = n_elements(allspec.wave)

for i=0, 4 do begin
    bhmin = 8.0+0.4*i
    bhmax = bhmin+0.4
    for j=0, 4 do begin
        counter, i*5+j, 25
        alphamin = -2.5+0.5*j
        alphamax = alphamin+0.5
        ii = where(stat.spec_snr_median gt 5. and bhall.logbh gt bhmin and bhall.logbh lt bhmax $
               and bhall.alpha gt alphamin and bhall.alpha lt alphamax, nn)
        allflux = fltarr(nn, noutwave)
        allivar = fltarr(nn, noutwave)
        i2800 = value_locate(allspec.wave, 2800.*(1.+allspec.zqso[ii]))
        print, nn
        for k=0L, nn-1L do begin
            ibegin = -i2800[k]+pos_2800
            iend = ibegin+ninwave-1L
            allflux[k, ibegin:iend] = allspec.flux[ii[k],*]/allspec.norm_ratio[ii[k]]
            allivar[k, ibegin:iend] = allspec.ivar[ii[k],*]*allspec.norm_ratio[ii[k]]^2
        endfor
        jhusdss_composite_engine, allflux, allivar, fmean=fmean, fmedian=fmedian, $
           fgeomean=fgeomean, nobjuse=nobjuse, ivarweight=ivarweight, weight=weight
        composite_bh_alpha[i*5+j].fluxmedian=fmedian
        composite_bh_alpha[i*5+j].logbh=median(bhall[ii].logbh)
        composite_bh_alpha[i*5+j].alpha=median(bhall[ii].alpha)
    endfor
endfor

save, filename='composite_bh_alpha.sav', composite_bh_alpha
end
