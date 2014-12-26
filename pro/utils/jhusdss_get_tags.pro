;+
; Documentation needed!
;-

function jhusdss_get_tags, ztag=ztag, platetag=platetag, fibertag=fibertag, $
    mjdtag=mjdtag, magtag=magtag

;; 
if (keyword_set(ztag)) then tag='z' ;; using hewett&wild2010 catalog for DR 7 -- 03/19/2012
if (keyword_set(platetag)) then tag='plate'
if (keyword_set(fibertag)) then tag='fiber'
if (keyword_set(mjdtag)) then tag='mjd'
if (keyword_set(magtag)) then tag='mi_z2'

return, tag
end
