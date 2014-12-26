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
function jhusdss_bin_blank

strtmp = {min:0.d0,$
          max:0.d0,$
          mean:0.d0,$
          mean_2d:0.d0}
;         mean_3d:0.d0}
return, strtmp

end

function jhusdss_make_bins, min_bin, max_bin, binsize, delta=delta, nbin=nbin

if (n_params() ne 3) then begin
   doc_library, "jhusdss_make_bins"
   message, "You need to specify the bin range and bin size."
endif

if (n_elements(delta) eq 0) then delta=binsize

nbin = floor((max_bin-delta-min_bin)/binsize)+1L
strtmp = jhusdss_bin_blank()
outstr = replicate(strtmp, nbin)
outstr.min = dindgen(nbin)*binsize+min_bin
outstr.max = outstr.min+delta
outstr.mean = (outstr.min+outstr.max)/2.d0
outstr.mean_2d = sqrt((outstr.min^2+outstr.max^2)/2.d0)
;outstr.mean_3d = ((outstr.min^3+outstr.max^2)/2.d0)^(1.d0/3.d0)

return, outstr
end
