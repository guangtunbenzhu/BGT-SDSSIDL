;+
;-
pro jhusdss_detect_absorbers_qso2abs, nmfver, overwrite=overwrite, boss=boss, dr12=dr12

if (n_elements(nmfver) eq 0) then message, 'nmfver required'

;; output
if (keyword_set(dr12)) then begin
   path=jhusdss_get_path(/nmfqso)+'/'+string(nmfver, format='(I3.3)')+'/Absorbers_DR12'
endif else begin
   if (keyword_set(boss)) then begin
      path=jhusdss_get_path(/nmfqso)+'/'+string(nmfver, format='(I3.3)')+'/Absorbers_BOSS'
   endif else begin
      path=jhusdss_get_path(/nmfqso)+'/'+string(nmfver, format='(I3.3)')+'/Absorbers'
   endelse
endelse

if (jhusdss_direxist(path) eq 0) then message, "Can't find the directory."
infile = path+'/Window_'+jhusdss_absorbers_filename(nmfver, /mgii, boss=boss, dr12=dr12)
outfile = path+'/OnlyAbsorbers_Window_'+jhusdss_absorbers_filename(nmfver, /mgii, boss=boss, dr12=dr12)

if (keyword_set(dr12)) then begin
   qsopath=jhusdss_get_path(/nmfqso)+'/'+string(nmfver, format='(I3.3)')+'/Decompose_DR12'
endif else begin
   if (keyword_set(boss)) then begin
      qsopath=jhusdss_get_path(/nmfqso)+'/'+string(nmfver, format='(I3.3)')+'/Decompose_BOSS'
   endif else begin
      qsopath=jhusdss_get_path(/nmfqso)+'/'+string(nmfver, format='(I3.3)')+'/Decompose'
   endelse
endelse

if (jhusdss_direxist(qsopath) eq 0) then message, "Can't find the directory."
statfile = qsopath+'/'+jhusdss_stat_filename(nmfver, boss=boss, dr12=dr12)

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
