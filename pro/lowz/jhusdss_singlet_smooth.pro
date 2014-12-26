pro jhusdss_singlet_smooth, flux, ivar, sigma=sigma, $
       outflux=outflux, outivar=outivar, normalize=normalize, factor_norm=factor_norm

;; kill it -- Guangtun 03/12/2013
;common singlet_common, singlet_gaussfilter, singlet_factor_norm


if (n_elements(singlet_gaussfilter) eq 0) then begin
   ;; pixel
   if (n_elements(sigma) eq 0) then sigma = 2

   filtersize = (FIX(sigma*4)+(FIX(sigma*4+1) mod 2))*2

   singlet_gaussfilter = psf_gaussian(NPIXEL=filtersize*2, $
                  ST_DEV=sigma, NDIMEN=1, NORMALIZE=normalize)

   singlet_factor_norm = total(singlet_gaussfilter*singlet_gaussfilter)
endif

outflux = flux
smooth_err = 1./ivar
imask = where(ivar le 0., nmask)
;if (nmask gt 0) then begin
;   outflux[imask] = median(outflux)
;   smooth_err[imask] = median(smooth_err)
;endif

factor_norm = singlet_factor_norm
outflux = CONVOL(outflux, singlet_gaussfilter, $
                 /nan, /center, /edge_wrap, normalize=normalize)

;; Is this correct?
smooth_err = CONVOL(smooth_err, singlet_gaussfilter, $
                /nan, /center, /edge_wrap, normalize=normalize)

outivar = 1./smooth_err
if (nmask gt 0) then outivar[imask] = 0.

end
