;+
; Documentation needed!
; Ideally, this is the only place one needs to change for path/directory
;-

function jhusdss_get_parent_path

;host = getenv('HOST')
;if (strmatch(host, '*gwln1*', /fold_case) eq 1) then begin
    parent_path = '/data1/gwln2scratch/menard/gz323/'
;endif else begin
;    parent_path = '/export/scratch1/menard/gz323/'
;endelse

return, parent_path
end
