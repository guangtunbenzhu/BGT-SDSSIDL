;; See jhusdss_pitts_qso2abs.pro
pro jhusdss_detect_absorbers_train_match, nmfver, overwrite=overwrite
if (n_elements(nmfver) eq 0) then message, "nmfver required."

path=jhusdss_get_path(/nmfqso)+'/'+string(nmfver, format='(I3.3)')+'/Absorbers'

mfile = path+'/'+'JHU_Pitts_matched_'+string(nmfver, format='(I3.3)')+'.fits'

nomfile = path+'/'+'Pitts_nomatched_'+string(nmfver, format='(I3.3)')+'.fits'
jhu_nomfile = path+'/'+'JHU_nomatched_'+string(nmfver, format='(I3.3)')+'.fits'

if (file_test(mfile) and ~keyword_set(overwrite)) then begin
   splog, 'File already exists. Use /overwrite if you want to overwrite it.'
   return
endif else begin
   splog, 'Will write the matched absorber catalog into this file: '
   print, mfile
   print, nomfile
   print, jhu_nomfile
endelse

;; use the updated one with med_sdeviation_red
;; the Pittsburgh absorber catalog
;pittsfile = jhusdss_get_path(/absorber)+'/MgII/Master_Pitts_Catalog_Absorbers.fits' 
pittspath= jhusdss_get_path(/nmfqso)+'/'+string(nmfver, format='(I3.3)')+'/Absorbers'
pittsfile = pittspath+'/'+'Master_Pitts_Catalog_Absorbers_'+string(nmfver, format='(I3.3)')+'.fits'
pitts = mrdfits(pittsfile, 1)
npitts = n_elements(pitts)

;; the corresponding JHU catalog

jhupath= jhusdss_get_path(/nmfqso)+'/'+string(nmfver, format='(I3.3)')+'/Absorbers'
infile = jhupath+'/OnlyAbsorbers_Window_Pitts_'+jhusdss_absorbers_filename(nmfver, /mgii)
jhu = [mrdfits(infile, 1), mrdfits(infile, 3)]
;ijhu = where(jhu.criterion_mgii eq 1b, njhu)
;jhu = jhu[ijhu]
njhu = n_elements(jhu)

pitts_match = lonarr(npitts) - 1L
jhu_match = lonarr(njhu) - 1L

;; loop over all absorbers, stupid algorithm
for i=0L, npitts-1L do begin
    counter, i+1, npitts
    imatch = where(pitts[i].plate eq jhu.plate $
             and pitts[i].fiber eq jhu.fiber $
             and abs(pitts[i].zabs-jhu.zabs) le 0.003, nmatch)
    if nmatch ge 1 then begin
;      jhu_match[i] = imatch[0]
;      pitts_match[imatch[0]] = i
       pitts_match[i] = imatch[0]
       jhu_match[imatch[0]] = i
    endif
endfor

those_matched = where(pitts_match ge 0, comp=nomatch)
match_pitts = pitts[those_matched]
match_jhu = jhu[pitts_match[those_matched]]

jhu_those_matched = where(jhu_match ge 0, comp=jhu_nomatch)
nomatch_pitts = pitts[nomatch]
nomatch_jhu = jhu[jhu_nomatch]

mwrfits, match_pitts, mfile, /create
mwrfits, match_jhu, mfile
mwrfits, nomatch_pitts, nomfile, /create
mwrfits, nomatch_jhu, jhu_nomfile, /create

end
