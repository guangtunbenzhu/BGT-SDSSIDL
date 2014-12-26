;; See jhusdss_create_nmf_basis
;pro jhusdss_lrg_nmf_all, version, path=path

lrgver = 101
parentpath = '/data1/gwln2scratch/menard/gz323/SDSS/LRG/101/AllInOne/'

index = mrdfits(parentpath+'DR7_LRG_ALL_INDEX.fit', 1)
wave = mrdfits(parentpath+'DR7_LRG_ALL_WAVE.fit', 1)
flux = mrdfits(parentpath+'DR7_LRG_ALL_FLUX.fit', 1)
continuum = mrdfits(parentpath+'DR7_LRG_ALL_CONTINUUM_101.fit', 1)
norm_resi = mrdfits(parentpath+'DR7_LRG_ALL_NORMALIZED_RESIDUAL_101.fit', 1)
subt_resi = mrdfits(parentpath+'DR7_LRG_ALL_SUBTRACTED_RESIDUAL_101.fit', 1)
corr = mrdfits(parentpath+'DR7_LRG_ALL_CORRECTION_101.fit', 1)

nwave = n_elements(wave.wave)
nobj = n_elements(index)

norm_resi.residual = norm_resi.residual/corr.correction
norm_resi.residual_ivar = norm_resi.residual_ivar*corr.correction^2

mwrfits, norm_resi, parentpath+'DR7_LRG_ALL_NORMALIZED_RESIDUAL_101.fit', /create
delvar, norm_resi

continuum.nmf_continuum = continuum.nmf_continuum*corr.correction
continuum.continuum = continuum.continuum*corr.correction

subt_resi.subtracted_residual = flux.flux - continuum.continuum
mwrfits, subt_resi, parentpath+'DR7_LRG_ALL_SUBTRACTED_RESIDUAL_101.fit', /create
delvar, subt_resi

mwrfits, continuum, parentpath+'DR7_LRG_ALL_CONTINUUM_101.fit', /create
delvar, continuum

end
