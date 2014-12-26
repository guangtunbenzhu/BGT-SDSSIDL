;+
; obsolete
;-
pro jhusdss_detect_absorbers_qso2window_all, nmfver, overwrite=overwrite 

if (n_elements(nmfver) eq 0) then nmfver = jhusdss_get_nmf_version()

;; output
path=jhusdss_get_path(/nmfqso)+'/'+string(nmfver, format='(I3.3)')+'/Absorbers'
if (jhusdss_direxist(path) eq 0) then message, "Can't find the directory."
infile = path+'/'+jhusdss_absorbers_filename(nmfver, /mgii)

outfile = path+'/Window_'+jhusdss_absorbers_filename(nmfver, /mgii)

;outfile1 = path+'/Blue_window_MgII_'+jhusdss_absorbers_filename(nmfver, /mgii)
;outfile2 = path+'/Red_window_MgII_'+jhusdss_absorbers_filename(nmfver, /mgii)
;outfile3 = path+'/Blue_window_FeII_'+jhusdss_absorbers_filename(nmfver, /mgii)
;outfile4 = path+'/Red_window_FeII_'+jhusdss_absorbers_filename(nmfver, /mgii)

if (file_test(outfile1) and ~keyword_set(overwrite)) then begin
   splog, 'File already exists. Use /overwrite if you want to overwrite it.'
   return
endif else begin
   splog, 'Will write the absorber catalog into this file: '
   print, outfile1
   print, outfile2
   print, outfile3 
   print, outfile4
endelse

;; output from detect_absorber_all
qso = mrdfits(infile, 1)

;; output from decompose_stats_all
;tmppath=jhusdss_get_path(/nmfqso)+'/'+string(nmfver, format='(I3.3)')+'/Decompose'
;if (jhusdss_direxist(tmppath) eq 0) then message, "Can't find the directory."
;filename = path+'/QSO_decompose_NMF_'+string(nmfver, format='(I3.3)')+'stats.fits'
;stats = mrdfits(filename, 1)

jhusdss_detect_absorbers_qso2window, qso, mgii_red=mgii_red, mgii_blue=mgii_blue, $
   feii_red=feii_red, feii_blue=feii_blue

mwrfits, mgii_red, outfile, /create
mwrfits, mgii_blue, outfile
mwrfits, feii_red, outfile
mwrfits, feii_blue, outfile

end
