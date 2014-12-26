;+
; Docementation needed!
; parallelized
; use 'HW_dr7qso_newz.fits'
;-
pro jhusdss_decompose_montecarlo_all, nmfver, overwrite=overwrite, $
       boss=boss, iprocess=iprocess, nprocess=nprocess, noparal=noparal

;if (~keyword_set(nmfver)) then nmfver=jhusdss_get_nmf_version()
if (n_elements(nmfver) eq 0) then $
    message, "nmf version required so that you won't overwrite files unintentionally."
if (n_elements(nprocess) eq 0) then nprocess=8
if (n_elements(iprocess) eq 0) then iprocess=1
if (iprocess lt 1 or iprocess gt nprocess) then message, "iprocess should be in [1, nprocess]"

in_path = jhusdss_get_path(/nmfqso)+'/MC_'+string(nmfver, format='(I3.3)')
infile = in_path+'/MC_Master_absorbers.fits'
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
   sub_nspec = iend-ibegin+1
endelse

;; limit the simultaneous decomposion to 10000L quasars.
dn = 10000L
nbins = sub_nspec/dn
for i=0, nbins do begin
    n_min = ibegin+i*dn
    n_max = ibegin+(((i+1)*dn-1) < (sub_nspec-1L))
    if (n_max lt n_min) then continue
    jhusdss_decompose_montecarlo_spec, objs0[n_min:n_max], nmfver=nmfver, overwrite=overwrite, boss=boss
endfor

end
