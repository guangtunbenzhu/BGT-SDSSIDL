;+
; NAME:
;   jhusdss_lrgfit_line
; PURPOSE:
;   fit line with a gaussian
; CALLING SEQUENCE:
;   jhusdss_lrgfit_line, lambda, flux, continuum, ivar, line [, $
;       ew=, model=, lflux= ]
; REVISION HISTORY:
;   30-Apr-2012  Guangtun, JHU, Adopted to jhusdss
;   ??-???-2010  Guangtun, NYU, sfit_line
;   31-Mar-2004  MRB, NYU, sfit_line
;-
;------------------------------------------------------------------------------
function sfit_line_func, x, p

center=p[0]
sigma=p[1]
flux=p[2]

model=fltarr(n_elements(x))
model=flux*exp(-0.5*(x-center)^2/sigma^2)/sqrt(2.*!DPI)/sigma

return, model

end
;
pro sfit_line_fit, center, sigma, lflux

common com_sfit, use_lambda, use_flux, use_err, line, maxwidth

start=fltarr(3)
start[0]=line
start[1]=2.
start[2]=0.
parinfo1={ limited:bytarr(2), limits:fltarr(2) }
parinfo=replicate(parinfo1,3)
parinfo[0].limited=1
parinfo[0].limits=line+0.5*[-maxwidth, maxwidth]
parinfo[1].limited=1
parinfo[1].limits=[maxwidth/6., maxwidth]

p=mpfitfun('sfit_line_func', use_lambda, use_flux, use_err, $
           start,ftol=1.d-10,bestnorm=chi2, parinfo=parinfo, /quiet)

lflux=p[2]
center=p[0]
sigma=p[1]

end
;
pro jhusdss_lrgfit_line, in_lambda, in_flux, in_continuum, in_ivar, in_line, $
               ew=ew, model=model, lflux=lflux, err_lflux=err_lflux, $
               err_ew=err_ew, sigma=sigma, err_sigma=err_sigma, $
               vdisp=vdisp, rwave=rwave, res=res, outclevel=outclevel

common com_sfit

if (not keyword_set(rwave)) then rwave = [min(in_lambda), max(in_lambda)]
if (not keyword_set(vdisp)) then vdisp = 150.
if (not keyword_set(res)) then res=69.

iwave = where(in_lambda ge rwave[0] and in_lambda le rwave[1], nwave)
if (nwave lt 2) then message, 'No good wavelength available'

iwave = iwave[sort(iwave)]
lambda = median(in_lambda[iwave[0:nwave-1]])
dlambda = median(in_lambda[iwave[1:nwave-1]] - in_lambda[iwave[0:nwave-2]])

sig_cor = sqrt(vdisp^2+res^2)
sigma = sig_cor/2.99792458E5*lambda/dlambda
;maxwidth=5.  ;; for PRIMUS, BC templates, weird line position ???
;maxwidth=3. ;; for dimage project, SDSS spectra
maxwidth = sigma*2.
;print, 'vdisp ', vdisp, ' maxwidth ', maxwidth
nmc=10L

line=in_line
lambda=in_lambda
continuum=in_continuum
flux=in_flux
ivar=in_ivar

;; deal with badness
ibad=where(ivar ne ivar OR flux ne flux or lambda ne lambda, nbad)
if(nbad gt 0) then begin
    ivar[ibad]=0.
    flux[ibad]=0.
    lambda[ibad]=0.
    in_continuum[ibad]=0.
endif

lflux=0.
ew=0.
err_lflux=0.
err_ew=0.
model=0.
;sigma=0.
err_sigma=0.


iuse=where(ivar gt 0. and abs(lambda-line) lt maxwidth*2.5, nuse)
if(nuse eq 0) then return

use_flux=flux[iuse]-continuum[iuse]
use_lambda=lambda[iuse]
use_err=1./sqrt(ivar[iuse])

icont=where(ivar gt 0. and abs(lambda-line) gt maxwidth*2.5 and $
            abs(lambda-line) lt maxwidth*5., ncont)
if(ncont eq 0) then return
cont_flux=flux[icont]
cont_lambda=lambda[icont]
cont_err=1./sqrt(ivar[icont])
clevel=mean(cont_flux)
outclevel = clevel


sfit_line_fit, center, sigma, lflux
ew=lflux/clevel

mc_center=fltarr(nmc)
mc_sigma=fltarr(nmc)
mc_lflux=fltarr(nmc)
mc_ew=fltarr(nmc)
use_flux_orig=use_flux
for i=0L, nmc-1L do begin
    use_flux=use_flux_orig+randomn(seed, nuse)*use_err
    sfit_line_fit, tmp_center, tmp_sigma, tmp_lflux
    
    clevel=mean(cont_flux+randomn(seed, ncont)*cont_err)

    mc_center[i]=tmp_center
    mc_sigma[i]=tmp_sigma
    mc_lflux[i]=tmp_lflux
    mc_ew[i]=tmp_lflux/clevel
endfor
err_center=sqrt(total((mc_center-center)^2)/float(nmc))
err_sigma=sqrt(total((mc_sigma-sigma)^2)/float(nmc))
err_lflux=sqrt(total((mc_lflux-lflux)^2)/float(nmc))
err_ew=sqrt(total((mc_ew-ew)^2)/float(nmc))

model=lflux*exp(-0.5*(lambda-center)^2/sigma^2)/sqrt(2.*!DPI)/sigma

end
