;+
; Documentation needed!
; Get a training set for line detections
;-
function jhusdss_lines_trainingset, nmfver, absinfo=absinfo, notonlyone=notonlyone

donespec = jhusdss_decompose_donespec(nmfver)
;; Pittsburgh catalog
mgiifile = jhusdss_get_path(/absorber)+'/MgII/MgII_Nestor_QSOinfo.fit'
mgii = mrdfits(mgiifile, 1)

;; no multiple systems
tmpra = mgii.ra
tmpra = mgii[UNIQ(mgii.ra, SORT(mgii.ra))].ra
;newmgii = replicate(mgii[0], n_elements(tmpra))
str = replicate({nabs:0}, n_elements(mgii))
for i=0L, n_elements(tmpra)-1 do begin
    ii = where(mgii.ra eq tmpra[i], nn)
    str[ii].nabs = nn
endfor

idone = bytarr(n_elements(mgii))
index_done = lonarr(n_elements(mgii))
for i=0L, n_elements(mgii)-1 do begin
    ii = where(donespec.plate eq mgii[i].plate and donespec.fiber eq mgii[i].fiber, nn)
    if (nn ne 1 or str[i].nabs ne 1) then continue
    idone[i] = 1b
    index_done[i] = ii
endfor

jj = where(idone, mm)
index_done = index_done[jj]

absinfo = mgii[jj]
knownabsorbers = donespec[index_done]

return, knownabsorbers

end
