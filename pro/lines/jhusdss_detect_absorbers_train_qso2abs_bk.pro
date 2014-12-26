;+
;-
pro jhusdss_detect_absorbers_train_qso2abs, nmfver, overwrite=overwrite 

if (n_elements(nmfver) eq 0) then nmfver = jhusdss_get_nmf_version()

;; output
path=jhusdss_get_path(/nmfqso)+'/'+string(nmfver, format='(I3.3)')+'/Absorbers'
if (jhusdss_direxist(path) eq 0) then message, "Can't find the directory."
infile = path+'/Pitts_'+jhusdss_absorbers_filename(nmfver, /mgii)
outfile = path+'/OnlyAbsorbers_Pitts_'+jhusdss_absorbers_filename(nmfver, /mgii)

if (file_test(outfile) and ~keyword_set(overwrite)) then begin
   splog, 'File already exists. Use /overwrite if you want to overwrite it.'
   return
endif else begin
   splog, 'Will write the absorber catalog into this file: '
   print, outfile
endelse

qso = mrdfits(infile, 1)
absorbers = jhusdss_absorber_cat_qso2abs(qso)
help, absorbers

mwrfits, absorbers, outfile, /create

end
