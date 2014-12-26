;+
;-
pro jhusdss_detect_absorbers_window2pitts, nmfver, overwrite=overwrite

if (n_elements(nmfver) eq 0) then message, 'nmfver required'

;; output
path=jhusdss_get_path(/nmfqso)+'/'+string(nmfver, format='(I3.3)')+'/Absorbers'

if (jhusdss_direxist(path) eq 0) then message, "Can't find the directory."
infile = path+'/Window_'+jhusdss_absorbers_filename(nmfver, /mgii, boss=boss)
outfile = path+'/Window_Pitts_'+jhusdss_absorbers_filename(nmfver, /mgii, boss=boss)

if (file_test(outfile) and ~keyword_set(overwrite)) then begin
   splog, 'File already exists. Use /overwrite if you want to overwrite it.'
   return
endif else begin
   splog, 'Will write the matched absorber catalog into this file: '
   print, outfile
endelse

;; the Pittsburgh master catalog
qsofile = jhusdss_get_path(/absorber)+'/MgII/Master_Pitts_Catalog.fits'
qso = mrdfits(qsofile, 1)
nqso = n_elements(qso)
str_tmp = jhusdss_qso2window_blank(/absorber)

mgii_red = replicate(str_tmp, nqso)
mgii_blue = replicate(str_tmp, nqso)
feii_red = replicate(str_tmp, nqso)
feii_blue = replicate(str_tmp, nqso)

;; Master catalog
all_qsofile = jhusdss_get_path(/qso)+'/'+jhusdss_dr7_qsofile()
all_qso = mrdfits(all_qsofile, 1)

all_mgii_red = mrdfits(infile, 1)
all_mgii_blue = mrdfits(infile, 2)
all_feii_red = mrdfits(infile, 3)
all_feii_blue = mrdfits(infile, 4)

spherematch, all_qso.ra, all_qso.dec, qso.ra, qso.dec, 1./3600., m1, m2

for i=0L, n_elements(m2)-1L do begin
    counter, i+1, n_elements(m2)-1L
    for j=0L, n_tags(str_tmp)-1L do begin
	mgii_red[m2[i]].(j) = all_mgii_red[m1[i]].(j)
	mgii_blue[m2[i]].(j) = all_mgii_blue[m1[i]].(j)
	feii_red[m2[i]].(j) = all_feii_red[m1[i]].(j)
	feii_blue[m2[i]].(j) = all_feii_blue[m1[i]].(j)
    endfor
endfor

mwrfits, mgii_red, outfile, /create
mwrfits, mgii_blue, outfile
mwrfits, feii_red, outfile
mwrfits, feii_blue, outfile

end
