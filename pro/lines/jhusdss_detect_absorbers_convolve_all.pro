;+
; Still training
;-
pro jhusdss_detect_absorbers_convolve_all, nmfver, overwrite=overwrite, $
       qaonly=qaonly, boss=boss, dr12=dr12, iprocess=iprocess, nprocess=nprocess, noparal=noparal

if (n_elements(nmfver) eq 0) then message, 'nmfver required'
if (n_elements(nprocess) eq 0) then nprocess=8
if (n_elements(iprocess) eq 0) then iprocess=1
if (iprocess lt 1 or iprocess gt nprocess) then message, "iprocess should be in [1, nprocess]"

path = jhusdss_get_path(/qso)
if (keyword_set(dr12)) then begin
   filename = jhusdss_dr12_qsofile()
endif else begin
   if (keyword_set(boss)) then begin
      filename = jhusdss_boss_qsofile()
   endif else begin
      filename =  jhusdss_dr7_qsofile()
   endelse
endelse

infile = path+'/'+filename
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
endelse

jhusdss_detect_absorbers_convolve_spec, objs0[ibegin:iend], nmfver, lines=lines, emlines=emlines, $
      overwrite=overwrite, boss=boss, dr12=dr12

end
