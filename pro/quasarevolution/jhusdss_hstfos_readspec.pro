;+
; obj and channel are both scalar
;-
function jhusdss_hstfos_readspec, obj, channel, path=path, status=status

if (n_elements(path) eq 0) then path = '~/SDATA/Quasars/HSTFOS/'
subpath = repstr(repstr(strtrim(obj.name, 2), '+', 'p'), '-', 'm')
infile = path+subpath+'/'+channel+'.fits'
if (file_test(infile)) then begin
   status = 1
   return, mrdfits(infile, 1)
endif

status = -1
return, status

end
