;+
; 30-Apr-2010 Guangtun, JHU, Adopted jhusdss
; ??-???-2010 Guangtun, NYU, gt_spec_smooth
;-
pro jhusdss_lrgspec_smooth, wave, flux, ivar, $
    init_vdisp=init_vdisp, final_vdisp=final_vdisp, $
    res=res, rwave=rwave, outflux=outflux, outivar=outivar

if (not keyword_set(init_vdisp)) then init_vdisp = 150.
if (not keyword_set(final_vdisp)) then final_vdisp = 325.
if (not keyword_set(res)) then res=69.

c = 2.99792458E5 ;km/s

if (not keyword_set(rwave)) then rwave=[min(wave), max(wave)]

iwave = where(wave ge rwave[0] and wave le rwave[1], nwave)
if (nwave lt 2) then message, 'No good wavelength available'

iwave = iwave[sort(iwave)]
;; in AA
lambda = median(wave[iwave[0:nwave-1]])
dlambda = median(wave[iwave[1:nwave-1]] - wave[iwave[0:nwave-2]])

;; in km/s
sig_ori = sqrt(init_vdisp^2+res^2)
sig_final = sqrt(final_vdisp^2+res^2)

if (sig_final lt sig_ori) then begin
;  splog, 'No need to smooth'
   outflux=flux
   outivar=ivar
   return
endif

sig_cor = sqrt(sig_final^2-sig_ori^2)

;; to pixel
sigma = sig_cor/c*lambda/dlambda
filtersize=(FIX(sigma*8)+(FIX(sigma*8+1) mod 2))

gaussfilter = psf_gaussian(NPIXEL=filtersize, $
                  ST_DEV=sigma, NDIMEN=1,/NORMALIZE)
outflux = CONVOL(flux, gaussfilter, $
                 /center, /edge_truncate)

smooth_err = CONVOL(1./ivar, gaussfilter, $
                /center, /edge_truncate)
outivar = 1./smooth_err

return

end
