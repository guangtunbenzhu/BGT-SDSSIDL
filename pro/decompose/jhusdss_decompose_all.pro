;+
; Docementation needed!
; parallelized
; use 'HW_dr7qso_newz.fits'
;-
pro jhusdss_decompose_all, nmfver, overwrite=overwrite, $
       boss=boss, dr12=dr12, iprocess=iprocess, nprocess=nprocess, noparal=noparal

;if (~keyword_set(nmfver)) then nmfver=jhusdss_get_nmf_version()
if (n_elements(nmfver) eq 0) then $
    message, "nmf version required so that you won't overwrite files unintentionally."
if (n_elements(nprocess) eq 0) then nprocess=8
if (n_elements(iprocess) eq 0) then iprocess=1
if (iprocess lt 1 or iprocess gt nprocess) then message, "iprocess should be in [1, nprocess]"

path = jhusdss_get_path(/qso)
filename =  jhusdss_dr7_qsofile()

if (keyword_set(boss)) then begin
   filename = jhusdss_boss_qsofile()
endif

if (keyword_set(dr12)) then begin
   filename = jhusdss_dr12_qsofile()
endif

infile = path+'/'+filename
splog, 'reading '+infile
objs0 = mrdfits(infile, 1)
nspec = n_elements(objs0)

;; parallelize
if (~keyword_set(noparal)) then begin
   ibegin = (iprocess-1)*nspec/nprocess
   iend = iprocess*nspec/nprocess-1
   if (iprocess eq nprocess) then iend = nspec-1
   sub_nspec = iend-ibegin+1
   print, ibegin, iend, sub_nspec
endif else begin
   ibegin = 0L
   iend = nspec-1L
endelse

;; limit the simultaneous decomposion to 10000L quasars.
dn = 10000L
nbins = sub_nspec/dn
for i=0, nbins do begin
    n_min = ibegin+i*dn
    n_max = ibegin+(((i+1)*dn-1) < (sub_nspec-1L))
    jhusdss_decompose_spec, objs0[n_min:n_max], nmfver=nmfver, overwrite=overwrite, boss=boss, dr12=dr12
endfor

end
