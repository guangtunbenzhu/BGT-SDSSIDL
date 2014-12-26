;+
; NAME:
;   jhusdss_get_composite_version
; PURPOSE:
;   Return the current version of composite spectra
; CALLING SEQUENCE:
;   vers = jhusdss_get_composite_version
; COMMENTS:
;   -- Not generalized enough
;   -- No errors being output -- Do a jackknife or bootstrap to test for ensemble errors
;   -- A simplified version of gt_sdss_lrgstack.pro
; REVISION HISTORY:
;   16-Nov-2011  Guangtun Zhu, JHU
;-
function jhusdss_get_composite_version, path=path

;; directory
if (n_elements(path) eq 0) then
path = getenv("ALL_DATA")+"/SDSS/QSO/Composite"

;;Get the directories
pdir = file_search(path+'/*')
ndir = n_elements(pdir)
dirkeep = replicate(1, ndir)
pver = strarr(ndir)

;;Loop through the directories and check the versions
for ii = 0, ndir - 1 do begin
    junk = strsplit(pdir[ii], '/', /extract)
    pver[ii] = junk[n_elements(junk)-1]
    if isdigit(pver[ii]) eq 0 then dirkeep[ii] = 0 else $
    if long(pver[ii]) le 100 then dirkeep[ii] = 0
endfor

;;Take the maximum
keep = where(dirkeep eq 1)
maxval = max(long(pver[keep]), jkeep)
version = string(maxval, format='(i3.3)')

return, version

end
