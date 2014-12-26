;+
; Docementation needed!
; Need to break this into two pieces:
;    1. A generic wraper that takes in a list of objects
;    2. A piece of code that gives a lists of objects
;-
pro jhusdss_decompose_training, overwrite=overwrite, nmfver=nmfver

if (~keyword_set(nmfver)) then nmfver=102L
basis_path = jhusdss_get_path(/nmfqso)+'/'+string(nmfver, format='(I3.3)')

;; DR7
path = jhusdss_get_path(/qso)
infile =  'dr7_bh_May09_2011.fits.gz'
filename = path+'/'+infile
objs0 = mrdfits(filename, 1)
ztag = jhusdss_get_tags(/ztag)
zindex = tag_indx(objs0[0], ztag)
platetag = jhusdss_get_tags(/plate)
pindex = tag_indx(objs0[0], platetag)

;;  known MgII absorbers
;mgiifile = jhusdss_get_path(/absorber)+'/MgII/MgII_Nestor_QSOinfo.fit'
;mgii = mrdfits(mgiifile, 1)

;spherematch, objs0.ra, objs0.dec, mgii.ra, mgii.dec, 1./3600., m1, m2
;splog, 'N known absorbers', n_elements(m1)
;iabs = m1[uniq(m1, sort(m1))]
;objs = objs0[iabs]
objs = objs0

;; divide them into redshift bins
for iz=0L, 2L do begin
    basis_bins = jhusdss_decompose_bins(option=(iz+1))
    ithis = where(objs.(zindex) gt basis_bins.zmin and objs.(zindex) le basis_bins.zmax and objs.(pindex) lt 500, nthis)
    if (nthis eq 0) then continue
    thisobjs = objs[ithis]
    splog, "Perfomring NMF decomposition for ", nthis, " QSOs,"
    splog, "between redshift ", basis_bins.zmin, " and ", basis_bins.zmax

    basis_file = basis_path+'/'+basis_bins.basisfile
    if (file_test(basis_file) eq 0) then message, "Can't find the basis file."
    basis = mrdfits(basis_file, 1)
    loglam = alog10(basis.wave)
    eigen_vectors = basis.eigen_vectors

    tempwave = jhusdss_normwave_minmax(option=(iz+1))
    normminwave = tempwave[0]
    normmaxwave = tempwave[1]

    jhusdss_load_interp_spec, thisobjs, loglam=loglam, $
       ztag=ztag, fibertag=fibertag, platetag=platetag, mjdtag=mjdtag, $
       allflux=allflux, allivar=allivar, normminwave=normminwave, normmaxwave=normmaxwave, $
       normoption=1

    ;; assign 0 weight to non-positive flux.
    izero = where(allflux le 0., nzero)
    if (nzero gt 0) then allivar[izero]=0.
    ;; sanity check
    splog, "Sanity Checking ..."

    jhusdss_nmf_sanity_check, thisobjs, loglam, flux=allflux, ivar=allivar, $
       spec_index=spec_index, pix_index=pix_index
    if (n_elements(pix_index) ne (size(eigen_vectors))[2]) then eigen_vectors = eigen_vectors[*,pix_index]

    ;; no initilization
    eigen_values = randomu(seed, n_elements(thisobjs), (size(eigen_vectors))[1])+1.d-4
    jhusdss_nmf_engine_eigenvalue, allflux, eigen_vectors, weight=allivar, eigen_values=eigen_values

    nmf_continuum = eigen_vectors##eigen_values
    nmf_residual = allflux/nmf_continuum

    ;; median filtering
    mask = (allivar eq 0.)
    jhusdss_median_filter, nmf_residual, mask=mask, continuum=med_continuum, residual=residual

    ;; write out
    jhusdss_decompose_writeout, thisobjs, nmfver=nmfver, basisfile=basis_file, loglam=loglam, $
       nmf_continuum=nmf_continuum, med_continuum=med_continuum, residual=residual, $
       ivar=allivar, eigen_values=eigen_values, overwrite=overwrite

endfor

end
