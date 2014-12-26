;+
; Docementation needed!
;-
pro jhusdss_decompose_all, overwrite=overwrite, nmfver=nmfver, boss=boss

;if (~keyword_set(nmfver)) then nmfver=jhusdss_get_nmf_version()
if (n_elements(nmfver) eq 0) then $
    message, "nmf version required so that you won't overwrite files unintentionally."

path = jhusdss_get_path(/qso)
if (keyword_set(boss)) then begin
   infile = 'VAC4.fits'
endif else begin
   infile =  'dr7_bh_May09_2011.fits.gz'
endelse
filename = path+'/'+infile
splog, 'reading '+filename
objs0 = mrdfits(filename, 1)

;; could make this piece of sequential code parallel 
dn = 10000L
nbins = n_elements(objs0)/dn
for i=0, nbins do begin
    n_min = i*dn
    n_max = (((i+1)*dn-1) < (n_elements(objs0)-1L))
    jhusdss_decompose_spec, objs0[n_min:n_max], nmfver=nmfver, overwrite=overwrite, boss=boss
endfor

end
