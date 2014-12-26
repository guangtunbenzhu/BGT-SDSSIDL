;; could crash if usemin/usemax is NaN or Infinite
function jhusdss_sample_divider, sample, nsubsample, binsize=binsize, usemin=usemin, usemax=usemax

  if (nsubsample le 1) then message, "Number of subsamples has to be larger than 1 (Use Median if it is 2!)"
  if (n_elements(usemin) eq 0) then usemin = min(sample)
  if (n_elements(usemax) eq 0) then usemax = max(sample)

  divider = fltarr(nsubsample-1)

  if (nsubsample gt 2) then begin
      hist = histogram(sample, binsize=binsize, min=usemin, max=usemax)
      cumhist = total(hist, /cum)
      cumhist = cumhist/max(cumhist)
      binvalue = findgen(n_elements(hist))*binsize+usemin+binsize/2.

       for i=0L, nsubsample-2L do begin
           tmp = min(abs(cumhist-1./float(nsubsample)*(i+1.)), itmp)
           divider[i] = binvalue[itmp]
       endfor
   endif else begin
       itmp = where(sample ge usemin and sample le usemax)
       divider = median(sample[itmp])
   endelse

   return, divider
end
