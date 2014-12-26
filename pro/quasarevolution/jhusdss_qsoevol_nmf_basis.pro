;+
; NAME:
;   jhusdss_create_nmf_basis
; PURPOSE:
;   Create a set of NMF basis vectors
; CALLING SEQUENCE:
;   jhusdss_create_nmf_basis, objs, loglam, nmfver, outfile=outfile, $
;      normoption=1, normminwave=normminwave, normmaxwave=normmaxwave, $
;      /reject, /qaplot
; INPUTS:
;   objs - the training set
;   loglam - wavelength grid in log scale
;   nmfver - version of NMF intended
; OPTIONAL INPUTS:
;   specpath - where the individual spectra are stored
;   nmfpath  - where the basis vectors will be stored
;   outfile  - what file the basis vectors will be stored in
;   n_dimension - how many dimension of the basis intended
;   maxiters  - maximum iterations, by default 1000.
;   maxreject - maximum iterations for rejection of outliers, by default 10.
;   normoption  - how to normalize? by default 0, no normlization. See jhusdss_load_interp_spec.pro
;   normminwave - blue end of the normalization range
;   normmaxwave - red end of the normalization range
; KEYWORDS:
;   reject    - if set, reject outliers in the training set.
;   overwrite - if set, overwrite existing file
;   qaplot    - if set, make qaplots
; OPTIONAL OUTPUTS:
;   basis - the basis that is writen in the outfile
; ROUTINES CALLED:
;   jhusdss_load_interp_spec
;   jhusdss_nmf_sanity_check
;   jhusdss_nmf_engine
;   jhusdss_nmf_basis_qaplot
; COMMENTS:
; BUGS:
; REVISION HISTORY:
;   15-Nov-2011 Guangtun Ben Zhu, JHU
;-

pro jhusdss_quaevol_nmf_basis, objs0, loglam0, nmfver, specpath=specpath, $
    nmfpath=nmfpath, outfile=outfile, ztag=ztag, fibertag=fibertag, platetag=platetag, $
    mjdtag=mjdtag, n_dimension=n_dimension, qaplot=qaplot, qapath=qapath, qafile=qafile, $
    maxiters=maxiters, reject=reject, maxreject=maxreject, nooutput=nooutput, $
    basis=basis, normminwave=normminwave, normmaxwave=normmaxwave, normoption=normoption, $
    overwrite=overwrite, init_basis=init_basis

if ((n_elements(objs0) eq 0) or (n_elements(loglam0) eq 0)) then begin
   splog, "You have to give a list of objects and a (log) wavelength vector!"
   doc_library, 'jhusdss_create_nmf_basis'
   return
endif
;; create a local copy
objs = objs0
loglam = loglam0


if (n_elements(nmfver) eq 0) then begin
   splog, "No version number of nmf given, use 000!"
   nmfver = 0L
endif

if (n_elements(nmfpath) eq 0) then nmfpath=jhusdss_get_path(/nmfqso)+'/'+string(nmfver, format='(I3.3)')
if (jhusdss_direxist(nmfpath) eq 0) then spawn, 'mkdir -p '+nmfpath
if (n_elements(outfile) eq 0) then outfile='QSO_nmf_basis.fits'
filename = nmfpath+'/'+outfile
if (file_test(filename) eq 1) then begin
   splog, "File already exists."
   if (keyword_set(overwrite)) then begin
      splog, "I am overwriting." 
   endif else begin
      splog, "I am not overwriting."
      return
   endelse
endif

if (keyword_set(qaplot)) then begin
    if (n_elements(qapath) eq 0) then qapath=nmfpath+'/QAplot'
    if (jhusdss_direxist(qapath) eq 0) then spawn, 'mkdir -p '+qapath
    if (n_elements(qafile) eq 0) then qafile=repstr(outfile,'.fits','_QAplot.ps')
endif

if (n_elements(maxreject) eq 0) then maxreject=10L

;; load spec: allflux[nspec, npix], allivar[nspec, npix]
jhusdss_load_interp_spec, objs, loglam=loglam, $
   ztag=ztag, fibertag=fibertag, platetag=platetag, mjdtag=mjdtag, $
   allflux=allflux, allivar=allivar, normminwave=normminwave, normmaxwave=normmaxwave, $
   normoption=normoption, path=specpath

;; assign 0 weight to non-positive flux.
izero = where(allflux le 0., nzero)
if (nzero gt 0) then allivar[izero]=0.
;; sanity check
splog, "Sanity Checking ..."
jhusdss_nmf_sanity_check, objs, loglam, flux=allflux, ivar=allivar

splog, "Performing NMF..."
;; ready for nmf
if (n_elements(init_basis) gt 0) then w = init_basis
jhusdss_nmf_engine, allflux, weight=allivar, n_dimension=n_dimension, $
    eigen_vectors=w, eigen_values=h, maxiters=maxiters
noriginal = (size(allflux))[1]

;; reject weird ones, sigma-clipping, 
if (keyword_set(reject)) then begin
    n_sigma_cut = 5. ;; hardcoded
    iterreject=0L

    totaloutlier = 0L
    repeat begin
       isitoutlier = bytarr((size(allflux))[1])
       ;; w[n_dimension, n_pix], h[n_spec, n_dimension], wh[n_spec, n_pix]
       wh = w##h
       htmp = h
       ;; note we did not normalize spectra when loading, so we need to normalize here temporarily
       ;; Warning, this requires normalization before hand
;      for i=0L, (size(allflux))[1]-1 do htmp[i,*] = h[i,*]/median(wh[i,*])
       for i=0L, n_dimension-1L do begin
           average = median(htmp[*,i])
           sigma = stddev(htmp[*,i])
           ii = where(abs(htmp[*,i]-average) gt sigma*n_sigma_cut, nn)
           if (nn gt 0) then isitoutlier[ii] = 1b ;;Logical OR
       endfor
       ioutliers = where(isitoutlier, comp=index_use, noutliers)
       if (noutliers gt 0) then begin
          totaloutlier = totaloutlier+noutliers
          splog, "I found "+strtrim(string(noutliers),2)+" outliers!"

          ;; chuck outliers, remember to apply it to h as well
          objs = objs[index_use]
          allflux = allflux[index_use, *]
          allivar = allivar[index_use, *]
          if (n_elements(index_use) ne (size(h))[1]) then h = h[index_use, *]

          ;; sanity check
          jhusdss_nmf_sanity_check, objs, loglam, flux=allflux, ivar=allivar, $
             spec_index=spec_index, pix_index=pix_index

          ;; reinitialization
          if (n_elements(spec_index) ne (size(h))[1]) then h = h[spec_index, *]
          if (n_elements(pix_index) ne (size(w))[2]) then w = w[*, pix_index]

          ;; ready for nmf
          jhusdss_nmf_engine, allflux, weight=allivar, n_dimension=n_dimension, $
              eigen_vectors=w, eigen_values=h, maxiters=maxiters
       endif
       iterreject++
    endrep until (noutliers eq 0 or iterreject gt maxreject)
    splog, "Total outliers found "+strtrim(string(totaloutlier),2)+" out of "+strtrim(string(noriginal),2)
    splog, "Total rejection iterations: "+strtrim(string(iterreject),2)
endif

mask = (allivar eq 0.)
nobjuse = total(~mask, 1)
basis = {eigen_vectors:w, wave:10.^loglam, nobjuse:nobjuse}

; We should only save the objid or RA/DEC or index in the original file.
zindex = tag_indx(objs[0], ztag)
pindex = tag_indx(objs[0], platetag)
findex = tag_indx(objs[0], fibertag)
mindex = tag_indx(objs[0], mjdtag)
if (n_elements(magtag) eq 0) then magtag = jhusdss_get_tags(/magtag)
magindex = tag_indx(objs[0], magtag)

trainingset = {plate:objs.(pindex), fiber:objs.(findex), mjd:objs.(mindex), z:objs.(zindex), mag:objs.(magindex), eigen_values:h} 
if (~keyword_set(nooutput)) then begin
   mwrfits, basis, filename, /create
   mwrfits, trainingset, filename
endif

if (keyword_set(qaplot)) then $
   jhusdss_nmf_basis_qaplot, basis, trainingset, qapath=qapath, qafile=qafile

end
