;+
; Documentation needed
; min = [min_v, min_v+step, min_v+step*2, ..., max_v-delta]
; max = min+delta
; nbin = (max_v-delta-min_v)/step+1
; binsize is required
; if delta is not set, delta=binsize
; nbin is optional output
; Difference between make_bins and make_vectors:
;    1) binsize is required in make_bins and nbin is required in make_vectors
;    2) there is a delta argument in make_bins but not make_vectors
; One should add a small number to max_bin to get the desired bin number, e.g.
; bins = jhusdss_make_bins(0., 1.0+1.e-6, 0.1)
; bins = jhusdss_make_bins(0.d0, 1.0d0+1.d-14, 0.1d0)
;-
function jhusdss_average_blank

strtmp = {mean:0.0,$
          median:0.0,$
          sdev:0.d0}
return, strtmp

end

function jhusdss_compute_averages, x, y, xbin

if (n_params() ne 3) then begin
   doc_library, "jhusdss_compute_average"
   message, "You need to give the x, y and xbin."
endif

nbin = n_elements(xbin)
strtmp = jhusdss_average_blank()
outstr = replicate(strtmp, nbin)

for i=0, nbin-1L do begin
    itmp = where(x gt xbin[i].min and x le xbin[i].max, ntmp)
    if (ntmp eq 0) then continue
;   xtmp = x[itmp]
    ytmp = y[itmp]
    tmp = moment(ytmp)
    outstr[i].mean = tmp[0]
    outstr[i].sdev = sqrt(tmp[1])
    outstr[i].median = median(ytmp)
endfor

return, outstr
end
