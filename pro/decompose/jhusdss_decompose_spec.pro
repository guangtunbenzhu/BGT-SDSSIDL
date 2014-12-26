;+
; Decompose spetra given a list of objects
; objs includes z, plate, fiber, mjd
; Please check jhusdss_get_tags() to see what tags should be used
;
; to-do: check if file exists before running the decomposition...
;-

pro jhusdss_decompose_spec, objs, nmfver=nmfver, overwrite=overwrite, $
       qaonly=qaonly, silent=silent, boss=boss, dr12=dr12

;; Let's not deal with n(objs) > 10^4 at this moment
if (n_elements(objs) gt 10000L) then $
    message, "Please don't be too ambitious. Give me at most 10^4 objects."

if (n_elements(nmfver) eq 0) then message, "nmfver required"

basis_path = jhusdss_get_path(/nmfqso)+'/'+string(nmfver, format='(I3.3)')

;; This is by design.
;; We use different sets of basis vectors at different redshifts
;; And we have 4 sets of basis vectors

ztag = jhusdss_get_tags(/ztag)
zindex = tag_indx(objs[0], ztag)
for iz=0L, 3L do begin
    basis_bins = jhusdss_decompose_bins(option=(iz+1))
    ithis = where(objs.(zindex) gt basis_bins.zmin and $
                  objs.(zindex) le basis_bins.zmax, nthis)
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
    norm_minwave = tempwave[0]
    norm_maxwave = tempwave[1]

;; load and interpolate the spec
    jhusdss_load_interp_spec, thisobjs, loglam=loglam, $
       allflux=allflux, allivar=allivar, $
       normminwave=norm_minwave, normmaxwave=norm_maxwave, $
       normoption=1, boss=boss, dr12=dr12

    ;; assign 0 weight to non-positive flux.
    izero = where(allflux le 0., nzero)
    if (nzero gt 0) then allivar[izero]=0.
    ;; sanity check
    splog, "Sanity Checking ..."
    jhusdss_nmf_sanity_check, thisobjs, loglam, flux=allflux, ivar=allivar, $
       spec_index=spec_index, pix_index=pix_index
    if (n_elements(pix_index) ne (size(eigen_vectors))[2]) then $
       eigen_vectors = eigen_vectors[*,pix_index]

    ;; no initilization
    eigen_values = randomu(seed, n_elements(thisobjs), (size(eigen_vectors))[1])+1.d-4
    jhusdss_nmf_engine_eigenvalue, allflux, eigen_vectors, $
       weight=allivar, eigen_values=eigen_values
    nmf_continuum = eigen_vectors##eigen_values
    nmf_residual = allflux/nmf_continuum


    tmpivar = allivar*nmf_continuum^2

    ;; median filtering iteration -- Guangtun, 01/31/2012
    mask = (allivar eq 0.)
    filter_sizes=[91, 163]
    jhusdss_median_filter, nmf_residual, tmpivar, mask=mask, $
       continuum=tmp_continuum, residual=tmp_residual, filter_sizes=filter_sizes

    filter_sizes=[143, 71]
    mask = (allivar eq 0.) or (abs(tmp_continuum-nmf_residual)*sqrt(tmpivar) gt 1.5)
    jhusdss_median_filter, nmf_residual, tmpivar, mask=mask, $
       continuum=tmp_continuum, residual=tmp_residual, filter_sizes=filter_sizes

    filter_sizes=[143, 71]
    mask = (allivar eq 0.) or (abs(tmp_continuum-nmf_residual)*sqrt(tmpivar) gt 1.5)
    jhusdss_median_filter, nmf_residual, tmpivar, mask=mask, $
       continuum=med_continuum, residual=residual, filter_sizes=filter_sizes

    ;; write out
    if (~keyword_set(qaonly)) then begin
       jhusdss_decompose_writeout, thisobjs, nmfver=nmfver, $
          basisfile=basis_file, loglam=loglam, nmf_continuum=nmf_continuum, $
          med_continuum=med_continuum, residual=residual, flux=allflux, $
          ivar=allivar, eigen_values=eigen_values, overwrite=overwrite, boss=boss, dr12=dr12
    endif else begin
    ;; make the qaplots
       jhusdss_decompose_qaplot, thisobjs, loglam, nmf_continuum=nmf_continuum, $
          med_continuum=med_continuum, residual=residual, flux=allflux, ivar=allivar, $
          eigen_values=eigen_values
    endelse

endfor

end
