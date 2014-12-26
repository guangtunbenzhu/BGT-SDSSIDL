pro jhusdss_doublet_smooth, flux, ivar, separation=separation, sigma=sigma, $
       outflux=outflux, outivar=outivar, normalize=normalize, factor_norm=factor_norm, line_ratio=line_ratio

;; common is always dangerous
;; since we don't need to call this millions of times any more
;; let's kill it -- Guangtun 03/11/2013
;; common doublet_common, doublet_gaussfilter, doublet_factor_norm

if (n_elements(doublet_gaussfilter) eq 0) then begin
   if (n_elements(separation) eq 0) then $
      splog, 'separation (in pixel) required.'
   if (n_elements(line_ratio) eq 0) then begin
      splog, 'line_ratio not given, assumed to be 2.'
      line_ratio = 2.
   endif
;; pixel
   if (n_elements(sigma) eq 0) then sigma = 2

   ;filtersize = (((FIX(sigma*16)+(FIX(sigma*16+1) mod 2))*2) < ((n_elements(flux)-separation)/2))
   filtersize = (FIX(sigma*4)+(FIX(sigma*4+1) mod 2))*2

   ;gaussfilter = fltarr(separation+filtersize)

   doublet_gaussfilter = psf_gaussian(NPIXEL=filtersize*2+separation*2, $
                  ST_DEV=sigma, NDIMEN=1, NORMALIZE=normalize)

   gaussfilter2 = psf_gaussian(NPIXEL=filtersize*2+separation*2, $
                  ST_DEV=sigma, NDIMEN=1, NORMALIZE=normalize)/line_ratio

   ;gaussfilter[0L:filtersize-1L] = gaussfilter1
   doublet_gaussfilter[separation:separation*2+filtersize*2-1L] = $
     doublet_gaussfilter[separation:separation*2+filtersize*2-1L]+gaussfilter2[0:separation+filtersize*2-1L]

   doublet_gaussfilter = doublet_gaussfilter/total(doublet_gaussfilter)

   doublet_factor_norm = total(doublet_gaussfilter*doublet_gaussfilter)
endif

outflux = flux
smooth_err = 1./ivar
imask = where(ivar le 0., nmask)
;if (nmask gt 0) then begin
;   outflux[imask] = median(outflux)
;   smooth_err[imask] = median(smooth_err)
;endif

factor_norm = doublet_factor_norm
outflux = CONVOL(outflux, doublet_gaussfilter, $
             /nan, /center, /edge_wrap, NORMALIZE=normalize)

;; Is this correct?
smooth_err = CONVOL(smooth_err, doublet_gaussfilter, $
                /nan, /center, /edge_wrap, NORMALIZE=normalize)

outivar = 1./smooth_err
if (nmask gt 0) then outivar[imask] = 0.

end
