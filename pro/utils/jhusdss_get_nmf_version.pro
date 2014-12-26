;+
; NAME:
;   jhusdss_get_nmf_version
; PURPOSE:
;   Return the current version of nmf eigen vectors
; CALLING SEQUENCE:
;   vers = jhusdss_get_nmf_version()
; COMMENTS:
; REVISION HISTORY:
;   27-Dec-2011  Guangtun Zhu, JHU
;-
function jhusdss_get_nmf_version, path=path

;; directory
message, "We don't use this piece of code anymore. Give the version number explicitly."

if (n_elements(path) eq 0) then $
path = jhusdss_get_path(/nmfqso)

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
