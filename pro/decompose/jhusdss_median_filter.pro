;+
; Documentation needed!
;-

;; flux[nspec, npix], weight[nspec, npix]
;; continuum[nspec, npix]
pro jhusdss_median_filter, flux, ivar, mask=mask, continuum=continuum, $
       residual=residual, filter_sizes=filter_sizes

if (n_elements(mask) eq 0) then mask = (ivar eq 0.)
;; from large to small, odd number
;; change from small to large, 01/11/2012
;if (n_elements(filter_sizes) eq 0) then filter_sizes = [189, 63, 31]
;if (n_elements(filter_sizes) eq 0) then filter_sizes = [201, 101]
;if (n_elements(filter_sizes) eq 0) then filter_sizes = [169, 81, 45]
 if (n_elements(filter_sizes) eq 0) then filter_sizes = [45, 81, 169]
n_filter_sizes = n_elements(filter_sizes)

if (size(flux, /n_dimension) ne 2) then $
   message, 'The input data matrix should be in form of [Nk(n_spec), Nn(n_pix)]'

n_spec = (size(flux))[1]
n_pix = (size(flux))[2]
continuum = fltarr(n_spec, n_pix)+1.
residual = fltarr(n_spec, n_pix)+1.


;; -- multiscale approach to filtering
;; too many loops ...
index_pix = indgen(n_pix)


for ispec=0L, n_spec-1L do begin
    last_residual_level = reform(flux[ispec,*])

    ;; median treats NaN as missing data
    imask = where(mask[ispec,*] eq 1b)

    ;; default 1
    last_continuum_level = fltarr(n_pix)+1.
    continuum_estimate = fltarr(n_pix)+1.

    for i_filter_size=0,n_filter_sizes-1 do begin
        filter_half_width = filter_sizes(i_filter_size)
        last_residual_level[imask] = !values.f_nan
        temp_for_median = [fltarr(filter_half_width/2)+!values.f_nan, last_residual_level, fltarr(filter_half_width/2)+!values.f_nan]
        ;; change from reflective to periodical
;       temp_for_median = [last_residual_level[n_pix-1-filter_half_width/2-1:n_pix-1], last_residual_level, last_residual_level[0:filter_half_width/2-1]]
        continuum_estimate = (median(temp_for_median, filter_half_width, /even))[filter_half_width/2:n_pix+filter_half_width/2-1]
        ;; division
        ;; continuum_estimate could be NaN
        last_residual_level = last_residual_level/continuum_estimate
        last_continuum_level = last_continuum_level*continuum_estimate
        inan = where(finite(last_continuum_level) eq 0, nnan)
        if (nnan gt 0) then last_continuum_level[inan] = 1.
   endfor   
   ;; return possible NaN to 1.
;  inan = where(finite(last_continuum_level) eq 0, nnan)
;  if (nnan gt 0) then last_continuum_level[inan] = 1.
   last_residual_level = reform(flux[ispec,*])/last_continuum_level
   residual[ispec,*] = last_residual_level
   continuum[ispec,*] = last_continuum_level
endfor

end
