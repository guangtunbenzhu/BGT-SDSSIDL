;+
; Documentation Needed!
;-
function jhusdss_composite_get_weight, mag, binsize=binsize

if (n_elements(binsize) eq 0) then begin
    splog, "You didn't specify the bin size, use 0.05 mag"
    binsize = 0.05 
endif

mag_min = min(mag)
mag_max = max(mag)

hist = histogram(mag, binsize=binsize, min=mag_min, max=mag_max)
nbin = n_elements(hist)
bin_min = findgen(nbin)*binsize+mag_min
bin_max = bin_min+binsize

bin_index = floor((mag-mag_min+1.e-6)/binsize)
weight = fltarr(n_elements(mag))
weight[*] = 1./hist[bin_index]

return, weight
end
