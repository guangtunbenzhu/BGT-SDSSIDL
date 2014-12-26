;; This should be just a test version
;pro jhusdss_garching_galqso_match_spec_byrad_docom, lrgver, boss=boss, overwrite=overwrite, $
;       minrad=minrad, maxrad=maxrad

 lrgver = 101
 minrad = [0.003, 0.010, 0.020, 0.040, 0.080];, 0.160];, 0.320, 0.640];, 1.280, 2.560]
 maxrad = [0.010, 0.020, 0.040, 0.080, 0.160];, 0.320];, 0.640, 1.280];, 2.560, 5.120]
 zgalmin = 0.030
 zgalmax = 0.40
 rpmean = fltarr(n_elements(minrad))
 wmean = fltarr(n_elements(minrad))

 saveall = 1b
apply_cor = 0b

; maxrad = minrad*3.
print, 'radius (kpc): ', minrad*1E3, maxrad*1E3

ivarweight = 1

choice_load_data = 0
read,'load data? [1=yes, 0=no]: ',choice_load_data

;if (n_elements(lrgver) eq 0) then message, 'lrgver required'
;if (n_elements(minrad) eq 0) then minrad = 0.
;if (n_elements(maxrad) eq 0) then maxrad = 0.1
;if (minrad ge maxrad) then message, "minrad can't be larger than maxrad"

if choice_load_data eq 1 then begin

garching_path = jhusdss_get_path(/garching)
garching_file = garching_path+'/'+'gal_info_dr7_v5_2.fit.gz'
gal = mrdfits(garching_file, 1)
uniq_file = garching_path+'/'+'gal_uniq_dr7_v5_2.fits'
galuniq = mrdfits(uniq_file, 1)


;; read in SDSS qso catalog, need ra, dec
lrg = jhusdss_lrg_readin(boss=boss)
stat = jhusdss_lrgstats_readin(lrgver, boss=boss)

allspec = jhusdss_read_alllrgspec_old(lrgver, boss=boss, /flux)
allcont = jhusdss_read_alllrgspec_old(lrgver, boss=boss, /continuum)
allresi = jhusdss_read_alllrgspec_old(lrgver, boss=boss, /normresi)

endif

parentpath = '/data1/gwln2scratch/menard/gz323/SDSS/LRG/101/AllInOne/'
wave_file = 'DR7_LRG_ALL_WAVE.fit'
flux_file = 'DR7_LRG_ALL_FLUX.fit'
cont_file = 'DR7_LRG_ALL_CONTINUUM.fit'
norm_resi_file = 'DR7_LRG_ALL_NORMALIZED_RESIDUAL_101.fit'
subt_resi_file = 'DR7_LRG_ALL_SUBTRACTED_RESIDUAL_101.fit'

wave_file = parentpath+wave_file
flux_file = parentpath+flux_file
cont_file = parentpath+cont_file
norm_resi_file = parentpath+norm_resi_file
subt_resi_file = parentpath+ subt_resi_file

nwave = n_elements(allspec.wave)
nobj = n_elements(stat)

out_wave = {wave:allspec.wave}
mwrfits, out_wave, wave_file, /create

strtmp = {flux:fltarr(nwave), ivar:fltarr(nwave), norm_ratio:0.}
out_flux = replicate(strtmp, nobj)
out_flux.flux = transpose(allspec.flux)
out_flux.ivar = transpose(allspec.ivar)
out_flux.norm_ratio = allspec.norm_ratio

mwrfits, out_flux, flux_file, /create
delvar, out_flux

strtmp = {continuum:fltarr(nwave), nmf_continuum:fltarr(nwave), med_continuum:fltarr(nwave)}
out_cont = replicate(strtmp, nobj)
out_cont.nmf_continuum = transpose(allcont.nmf_continuum)
out_cont.med_continuum = transpose(allcont.med_continuum)
out_cont.continuum = out_cont.nmf_continuum*out_cont.med_continuum

mwrfits, out_cont, cont_file, /create
delvar, out_cont

strtmp = {residual:fltarr(nwave), residual_ivar:fltarr(nwave)}
out_norm_resi = replicate(strtmp, nobj)
out_norm_resi.residual = transpose(allresi.residual)
resi_ivar = allspec.ivar*allcont.nmf_continuum^2*allcont.med_continuum^2
out_norm_resi.residual_ivar = transpose(resi_ivar)

delvar, out_norm_resi




