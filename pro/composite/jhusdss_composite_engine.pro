;+
; Documentation needed
;-

;; flux[nspec, npix], ivar[nspec, npix], mask[nspec, npix], weight[nspec]
pro jhusdss_composite_engine, flux, ivar, mask=mask, weight=weight, ivarweight=ivarweight, $
         fmean=fmean, fmedian=fmedian, fgeomean=fgeomean, nobjuse=nobjuse, resistantmean=resistantmean, $
         sigma_cut=sigma_cut, fsigma=fsigma

if (size(flux, /n_dimension) ne 2) then $
   message, 'The input flux matrix should be in form of [Nk(n_spec), Nn(n_pix)]'
n_spec = (size(flux))[1]
n_pix = (size(flux))[2]
if (n_elements(sigma_cut) eq 0) then sigma_cut = 3.

if (n_elements(weight) eq 0) then weight = fltarr(n_spec)+1.
if (n_elements(mask) eq 0) then mask = (ivar le 0.)
;if (n_elements(mask) eq 0) then mask = (flux le 0. or ivar le 0.)

n_pix = (size(flux))[2]
fsigma = fltarr(n_pix)
fmean = fltarr(n_pix)
fmedian = fltarr(n_pix)
fgeomean = fltarr(n_pix)
nobjuse = fltarr(n_pix)

;; median treats NaN as missing data
;tmpflux = flux
;imask = where(mask eq 1b)
;tmpflux[imask] = !values.f_nan
;fmedian[*] = median(tmpflux, dimension=1)
;i_infinite = where(finite(fmedian) eq 0b, n_infinite)
;if (n_infinite gt 0) then fmedian[i_infinite]=0.
nobjuse[*] = total(~mask, 1)

;; mean and geometric mean
if (not keyword_set(ivarweight)) then begin
   for i=0L, n_pix-1L do begin
       iuse = where(~mask[*,i], nuse)
       if nuse gt 0 then begin
          moment_tmp = moment(flux[iuse,i], sdev=sdev_tmp)
          fsigma[i] = sdev_tmp
          resistant_mean, flux[iuse,i], sigma_cut, fmean_tmp, fmean_tmp_err, num_rejected, goodvec=igood
          fmedian[i] = median(flux[iuse[igood],i])
          fmean[i] = total(flux[iuse[igood],i]*weight[iuse[igood]])/(total(weight[iuse[igood]]) + (total(weight[iuse[igood]]) eq 0.))
          fgeomean[i] = jhusdss_geo_mean(flux[iuse[igood],i], weight=weight[iuse[igood]])
       endif
   endfor
endif else begin
   for i=0L, n_pix-1L do begin
       iuse = where(~mask[*,i], nuse)
       if nuse gt 0 then begin
          moment_tmp = moment(flux[iuse,i], sdev=sdev_tmp)
          fsigma[i] = sdev_tmp
          resistant_mean, flux[iuse,i], sigma_cut, fmean_tmp, fmean_tmp_err, num_rejected, goodvec=igood
          fmedian[i] = median(flux[iuse[igood],i])
          fmean[i] = total(flux[iuse[igood],i]*ivar[iuse[igood],i]*weight[iuse[igood]])/(total(weight[iuse[igood]]*ivar[iuse[igood],i]) + (total(weight[iuse[igood]]*ivar[iuse[igood],i]) eq 0.))
          newweight = (alog(10.)*(flux[iuse[igood],i]))^2*ivar[iuse[igood],i]*weight[iuse[igood]]
;         fgeomean[i] = jhusdss_geo_mean(flux[iuse[igood],i], weight=(ivar[iuse[igood],i]*weight[iuse[igood]]))
          fgeomean[i] = jhusdss_geo_mean(flux[iuse[igood],i], weight=newweight)
       endif
   endfor
endelse

end
