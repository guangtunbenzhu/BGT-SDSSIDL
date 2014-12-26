;+
; NAME:
;   jhusdss_create_composite_spectra
; PURPOSE:
;   Create composite spectra given a list of objects and weights
; CALLING SEQUENCE:
;   jhusdss_create_composite_spectra, objs, loglam, nmfver, outfile=outfile, $
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
;   < 15-Nov-2011 CREATE_NMF_EIGEN_VECTORS, Brice Menard, CITA&JHU
;   15-Nov-2011 Guangtun Ben Zhu, JHU
;-

pro jhusdss_create_composite_spectra, objs0, loglam0, compver, specpath=specpath, $
    comppath=comppath, outfile=outfile, ztag=ztag, fibertag=fibertag, platetag=platetag, $
    mjdtag=mjdtag, n_dimension=n_dimension, qaplot=qaplot, qapath=qapath, qafile=qafile, $
    maxiters=maxiters, reject=reject, maxreject=maxreject, nooutput=nooutput, $
    basis=basis, normminwave=normminwave, normmaxwave=normmaxwave, normoption=normoption, $
    overwrite=overwrite, magtag=magtag, ivarweight=ivarweight

if ((n_elements(objs0) eq 0) or (n_elements(loglam0) eq 0)) then begin
   splog, "You have to give a list of objects and a (log) wavelength vector!"
   doc_library, 'jhusdss_create_nmf_basis'
   return
endif
;; create a local copy
objs = objs0
loglam = loglam0

if (n_elements(magtag) eq 0) then magtag = jhusdss_get_tags(/magtag)
magindex = tag_indx(objs[0], magtag)

if (n_elements(compver) eq 0) then begin
   splog, "No version number of composite given, use 000!"
   compver = 0L
endif

if (n_elements(comppath) eq 0) then comppath=jhusdss_get_path(/composite)+'/'+string(compver, format='(I3.3)')
if (jhusdss_direxist(comppath) eq 0) then spawn, 'mkdir -p '+comppath
if (n_elements(outfile) eq 0) then outfile='QSO_composite.fits'
filename = comppath+'/'+outfile
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
    if (n_elements(qapath) eq 0) then qapath=comppath+'/QAplot'
    if (jhusdss_direxist(qapath) eq 0) then spawn, 'mkdir -p '+qapath
    if (n_elements(qafile) eq 0) then qafile=repstr(outfile,'.fits','_QAplot.ps')
endif

if (n_elements(maxreject) eq 0) then maxreject=10L

;; load spec: allflux[nspec, npix], allivar[nspec, npix]
jhusdss_load_interp_spec, objs, loglam=loglam, $
   ztag=ztag, fibertag=fibertag, platetag=platetag, mjdtag=mjdtag, $
   allflux=allflux, allivar=allivar, normminwave=normminwave, normmaxwave=normmaxwave, $
   normoption=normoption, path=specpath

;; sanity check
splog, "Sanity Checking ..."
jhusdss_nmf_sanity_check, objs, loglam, flux=allflux, ivar=allivar


;; weighting using lumisnoity distribution:
magnitude = objs.(magindex)
weight = jhusdss_composite_get_weight(magnitude)

jhusdss_composite_engine, allflux, allivar, weight=weight, ivarweight=ivarweight, $
   fmean=fmean, fmedian=fmedian, fgeomean=fgeomean, nobjuse=nobjuse

noriginal = n_elements(objs)
;; reject weird ones, sigma-clipping, use fmean for now
if (keyword_set(reject)) then begin
    n_sigma_cut = 50. ;; hardcoded
    iterreject=0L
    totaloutlier = 0L

    repeat begin
       isitoutlier = bytarr((size(allflux))[1])
       template = reform(fmean, 1, n_elements(fmean))
       coeffs = jhusdss_general_linfit_simple(allflux, allivar, template=template, reduced_chi2=reduced_chi2)
       ioutliers = where(reduced_chi2 gt n_sigma_cut, comp=index_use, noutliers)
       if (noutliers gt 0) then begin
          totaloutlier = totaloutlier+noutliers
          splog, "I found "+strtrim(string(noutliers),2)+" outliers!"
          objs = objs[index_use]
          allflux = allflux[index_use, *]
          allivar = allivar[index_use, *]

          ;; sanity check
          jhusdss_nmf_sanity_check, objs, loglam, flux=allflux, ivar=allivar

          ;; weighting using lumisnoity distribution:
          magnitude = objs.(magindex)
          weight = jhusdss_composite_get_weight(magnitude)
          jhusdss_composite_engine, allflux, allivar, weight=weight, ivarweight=ivarweight, $
             fmean=fmean, fmedian=fmedian, fgeomean=fgeomean, nobjuse=nobjuse
       endif
       iterreject++
    endrep until (noutliers eq 0 or iterreject gt maxreject)
    splog, "Total outliers found "+strtrim(string(totaloutlier),2)+" out of "+strtrim(string(noriginal),2)
    splog, "Total rejection iterations: "+strtrim(string(iterreject),2)
endif

stack = {nobj:n_elements(objs), wave:10.^loglam, fluxmean:fmean, fluxmedian:fmedian, $
         fluxgeomean:fgeomean, nobjuse:nobjuse}

; We should only save the objid or RA/DEC or index in the original file.
zindex = tag_indx(objs[0], ztag)
pindex = tag_indx(objs[0], platetag)
findex = tag_indx(objs[0], fibertag)
mindex = tag_indx(objs[0], mjdtag)
trainingset = {plate:objs.(pindex), fiber:objs.(findex), mjd:objs.(mindex), z:objs.(zindex), mag:objs.(magindex)}
if (~keyword_set(nooutput)) then begin
   mwrfits, stack, filename, /create
   mwrfits, trainingset, filename
endif

if (keyword_set(qaplot)) then $
   jhusdss_composite_qaplot, objs, stack, qapath=qapath, qafile=qafile

end
