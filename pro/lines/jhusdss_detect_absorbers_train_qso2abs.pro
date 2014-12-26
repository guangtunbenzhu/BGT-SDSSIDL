;+
;-
pro jhusdss_detect_absorbers_train_qso2abs, nmfver, overwrite=overwrite 

if (n_elements(nmfver) eq 0) then nmfver = jhusdss_get_nmf_version()

;; output
path=jhusdss_get_path(/nmfqso)+'/'+string(nmfver, format='(I3.3)')+'/Absorbers'
if (jhusdss_direxist(path) eq 0) then message, "Can't find the directory."
infile = path+'/Window_Pitts_'+jhusdss_absorbers_filename(nmfver, /mgii)
outfile = path+'/OnlyAbsorbers_Window_Pitts_'+jhusdss_absorbers_filename(nmfver, /mgii)

qsopath=jhusdss_get_path(/nmfqso)+'/'+string(nmfver, format='(I3.3)')+'/Decompose'
if (jhusdss_direxist(qsopath) eq 0) then message, "Can't find the directory."
statfile = qsopath+'/Pitts_QSO_decompose_NMF_'+string(nmfver, format='(I3.3)')+'stats.fits'


if (file_test(outfile) and ~keyword_set(overwrite)) then begin
   splog, 'File already exists. Use /overwrite if you want to overwrite it.'
   return
endif else begin
   splog, 'Will write the absorber catalog into this file: '
   print, outfile
endelse

stat = mrdfits(statfile, 1)

;; mgii_red
qso = mrdfits(infile, 1)
abs_mgii_red = jhusdss_absorber_cat_qso2abs(stat, qso)

;; mgii_blue
qso = mrdfits(infile, 2)
abs_mgii_blue = jhusdss_absorber_cat_qso2abs(stat, qso)

;; feii_red
qso = mrdfits(infile, 3)
abs_feii_red = jhusdss_absorber_cat_qso2abs(stat, qso)

;; feii_blue
qso = mrdfits(infile, 4)
abs_feii_blue = jhusdss_absorber_cat_qso2abs(stat, qso)

mwrfits, abs_mgii_red, outfile, /create
mwrfits, abs_mgii_blue, outfile
mwrfits, abs_feii_red, outfile
mwrfits, abs_feii_blue, outfile

end
