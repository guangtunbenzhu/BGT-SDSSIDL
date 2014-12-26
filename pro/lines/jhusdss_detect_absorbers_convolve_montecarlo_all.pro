;+
; Still training
;-
pro jhusdss_detect_absorbers_convolve_montecarlo_all, nmfver, overwrite=overwrite, $
       qaonly=qaonly, boss=boss, iprocess=iprocess, nprocess=nprocess, noparal=noparal

if (n_elements(nmfver) eq 0) then message, 'nmfver required'
if (n_elements(nprocess) eq 0) then nprocess=8
if (n_elements(iprocess) eq 0) then iprocess=1
if (iprocess lt 1 or iprocess gt nprocess) then message, "iprocess should be in [1, nprocess]"

in_path = jhusdss_get_path(/nmfqso)+'/MC_'+string(nmfver, format='(I3.3)')
infile = in_path+'/MC_Master_absorbers.fits'
splog, 'reading '+infile
objs0 = mrdfits(infile, 1)
nspec = n_elements(objs0)

lines = jhusdss_abslines_all(/train)
emlines = jhusdss_emlines_all(/mask)

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

jhusdss_detect_absorbers_convolve_montecarlo_spec, objs0[ibegin:iend], nmfver, lines=lines, emlines=emlines, $
      overwrite=overwrite, boss=boss

end
