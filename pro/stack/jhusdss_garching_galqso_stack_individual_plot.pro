
atmp = 'a'
nmfver = 106
overwrite = 1b
ivarweight = 1b
docut = 0b
qaplot = 0b
qaplot1 = 0b
loaddata = 1b
irp_min = 1
irp_max = 19
drp = 1
saveall = 1b
preprocess = 0b

sigma_smooth = 2.0
sigma_cut = 2.
snr_cut = 3.
z_cut = 0.04

;profile_slope = -1.29629
;profile_intercept = 0.636813

factor_norm1 = 0.141047
factor_norm2 = 0.0783596

read,'BOSS? [1=yes, 0=no]: ',BOSS

line_wave = [3934.79, 3969.59]

stackpath = jhusdss_get_path(/nmfqso)+'/'+string(nmfver, format='(I3.3)')+'/Stack/'
infile = stackpath + jhusdss_garching_galqso_stack_filename(nmfver, boss=boss)
spec_infile = repstr(infile, '.fits', '_spec.fits')

fit_infile = repstr(spec_infile, '.fits', '_fit.fits')
outfile = repstr(fit_infile, '.fits', '_plot.fits')

if (file_test(outfile) and ~overwrite) then begin
   splog, 'File already exists, use /overwrite to overwrite'
   return
endif else begin
   splog, 'Will write into this file: '
   print, outfile
endelse

comp0 = mrdfits(fit_infile,1)
nrp = comp0.nrp
wave = comp0.wave
nwave = n_elements(wave)
tmp = min(abs(comp0.wave-line_wave[0]), icaii)
iuse = where(wave gt xra[0] and wave lt xra[1], nuse)
iuse_nocaii = where((wave gt xra[0] and wave lt (line_wave[0]-10.) $
     or (wave gt (line_wave[0]+10.) and wave lt (line_wave[1]-10.))) $
     or (wave gt (line_wave[1]+10.) and wave lt (xra[1])), nuse_nocaii)

fix_separation = long((line_wave[1]-line_wave[0])/(line_wave[1]+line_wave[0])/alog(10.)*2E+4)

rp = fltarr(nrp)
ew_nofit_mean = fltarr(nrp)
ew_nofit_mean_2 = fltarr(nrp)
ew_nofit_median = fltarr(nrp)
ew_nofit_median_2 = fltarr(nrp)
ew_nofit_geomean = fltarr(nrp)
ew_nofit_geomean_2 = fltarr(nrp)

for i=irp_min, irp_max, drp do begin
    print, i, irp_min, irp_max
    comp = mrdfits(fit_infile, i)
  
    rp[i] = comp.rp
    ew_nofit_mean[i] = (1.-comp.fluxmean_singlet[icaii])*(wave[icaii]-wave[icaii-1])/factor_norm1
    ew_nofit_median[i] = (1.-comp.fluxmedian_singlet[icaii])*(wave[icaii]-wave[icaii-1])/factor_norm1
    ew_nofit_geomean[i] = (1.-comp.fluxgeomean_singlet[icaii])*(wave[icaii]-wave[icaii-1])/factor_norm1
    ew_nofit_mean_2[i] = (1.-comp.fluxmean_doublet[icaii])*(wave[icaii]-wave[icaii-1])/factor_norm2*6./5.
    ew_nofit_median_2[i] = (1.-comp.fluxmedian_doublet[icaii])*(wave[icaii]-wave[icaii-1])/factor_norm2*6./5.
    ew_nofit_geomean_2[i] = (1.-comp.fluxgeomean_doublet[icaii])*(wave[icaii]-wave[icaii-1])/factor_norm2*6./5.
endfor

load_dp, /b

i_indep = [1, 2, lindgen(8)*2+3]
i_indep = lindgen(nrp)

xra = [10, 1000.]
yra = [1E-5, 1E0]
djs_plot, rp[i_indep], ew_nofit_mean[i_indep], psym=4, symsize=2, /xlog, /ylog, xra=xra, yra=yra, xst=1, yst=1, thick=thick, xthick=xthick, ythick=ythick
djs_oplot, rp[i_indep], ew_nofit_median[i_indep], psym=5, symsize=2, thick=thick
djs_oplot, rp[i_indep], ew_nofit_geomean[i_indep], psym=6, symsize=2, thick=thick
djs_oplot, rp[i_indep], ew_nofit_mean_2[i_indep]*5./6., psym=4, symsize=2, thick=thick, color='red'
djs_oplot, rp[i_indep], ew_nofit_median_2[i_indep]*5./6., psym=5, symsize=2, thick=thick, color='red'
djs_oplot, rp[i_indep], ew_nofit_geomean_2[i_indep]*5./6., psym=6, symsize=2, thick=thick, color='red'

xx = findgen(1000)*0.01

profile_slope = -1.99629
profile_intercept = 1.936813
yy = profile_slope*xx+profile_intercept
djs_oplot, 10.^xx, 10.^yy

xra = [150., 1000.]
yra = [-1.E-3, 1.E-3]
djs_plot, rp, ew_nofit_mean, psym=4, symsize=2, xra=xra, yra=yra, xst=1, yst=1, thick=thick, xthick=xthick, ythick=ythick
djs_oplot, rp, ew_nofit_median, psym=5, symsize=2, thick=thick
djs_oplot, rp, ew_nofit_geomean, psym=6, symsize=2, thick=thick
djs_oplot, rp, ew_nofit_mean_2*5./6., psym=4, symsize=2, thick=thick, color='red'
djs_oplot, rp, ew_nofit_median_2*5./6., psym=5, symsize=2, thick=thick, color='red'
djs_oplot, rp, ew_nofit_geomean_2*5./6., psym=6, symsize=2, thick=thick, color='red'
djs_oplot, 10.^xx, 10.^yy
djs_oplot, !x.crange, 5.E-4*[1,1]
djs_oplot, !x.crange, -5.E-4*[1,1]

end
