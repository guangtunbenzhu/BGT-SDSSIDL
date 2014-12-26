;+
; NAME:
;   jhusdss_geo_mean()
; PURPOSE:
;   Calculate geometric mean
; CALLING SEQUENCE:
;   geomean = jhusdss_geo_mean(data, weight=weight, rms=rms)
; INPUTS: 
;   data - 
; OPTIONAL INPUTS:
;   weight - 
; OUTPUTS:
;   geomean - 
; OPTIONAL OUTPUTS:
;   rms - 
; COMMENTS:
;  -- Ignore all infinite/negative data points
; REVISION HISTORY:
;  16-Nov-2011, Guangtun Zhu, JHU 
;-
function jhusdss_geo_mean, data, weight=weight, rms=rms, verbose=verbose

igood = where((finite(data) ne 0) and (data gt 0.) and (weight gt 0.), ngood)
if ngood gt 1 then begin
   rms = 10^stddev(alog10(data[igood]))
   return, 10^(total(weight[igood]*alog10(data[igood])/total(weight[igood])))
endif else begin
   if (keyword_set(verbose)) then $
      print, "Warning! jhu_geo_mean finds less than 5 good values and refuses to calculate!"
   rms = 0.
   return, 0.
end


end
