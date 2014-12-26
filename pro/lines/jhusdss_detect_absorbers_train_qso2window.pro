;+
;-
pro jhusdss_detect_absorbers_train_qso2window, nmfver, overwrite=overwrite 

if (n_elements(nmfver) eq 0) then nmfver = jhusdss_get_nmf_version()

;; input/output
path=jhusdss_get_path(/nmfqso)+'/'+string(nmfver, format='(I3.3)')+'/Absorbers'
if (jhusdss_direxist(path) eq 0) then message, "Can't find the directory."
infile = path+'/Pitts_'+jhusdss_absorbers_filename(nmfver, /mgii)

outfile = path+'/Window_Pitts_'+jhusdss_absorbers_filename(nmfver, /mgii)

if (file_test(outfile) and ~keyword_set(overwrite)) then begin
   splog, 'File already exists. Use /overwrite if you want to overwrite it.'
   return
endif else begin
   splog, 'Will write the absorber catalog into this file: '
   print, outfile
endelse

;; output from detect_absorber_all
qso = mrdfits(infile, 1)

jhusdss_detect_absorbers_qso2window_spec, qso, mgii_red=mgii_red, mgii_blue=mgii_blue, $
   feii_red=feii_red, feii_blue=feii_blue

mwrfits, mgii_red, outfile, /create
mwrfits, mgii_blue, outfile
mwrfits, feii_red, outfile
mwrfits, feii_blue, outfile

end
