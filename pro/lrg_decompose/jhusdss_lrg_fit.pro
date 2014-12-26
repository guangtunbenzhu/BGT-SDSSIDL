;+
; NAME:
;   jhusdss_lrg_fit
; PURPOSE:
;   fit spectrum continuum and lines
; CALLING SEQUENCE:
;   jhusdss_lrg_fit, lambda, flux, ivar [, continuum=, lflux=, err_lflux=, $
;       ew=, err_ew= ]
; REVISION HISTORY:
;   30-Apr-2010  Guangtun, JHU, Adopted to jhusdss
;   ??-???-2010  Guangtun, NYU, mrb_sfit
;   03-Aug-2007  MRB, NYU, mrb_sfit
;-
;------------------------------------------------------------------------------
pro jhusdss_lrg_fit, in_lambda, in_flux, in_ivar, continuum=continuum, $
          lflux=lflux, err_lflux=err_lflux, ew=ew, err_ew=err_ew, $
          lines=lines, names=names, linemodel=linemodel, sigma=sigma, $
          err_sigma=err_sigma, vdisp=vdisp, clevel=clevel

lambda=in_lambda
flux=in_flux
ivar=in_ivar

continuum=jhusdss_lrgfit_continuum(lambda, flux, ivar, vdisp=vdisp)

;; see jhusdss_lrg_masklines.pro
lines= [ 4862.6778, $
          4341.6803, $
          4105.8884, $
          3727.0100, $
          2800.0000, $
          5008.1666, $
          6564.6127]
names=['Hbeta', $
       'Hgamma', $
       'Hdelta', $
       'OII', $
       'MgII', $
       'OIII5007', $
       'Halpha']

nlines=n_elements(lines)

ew=fltarr(nlines)
err_ew=fltarr(nlines)
lflux=fltarr(nlines)
clevel=fltarr(nlines)
err_lflux=fltarr(nlines)
sigma = fltarr(nlines)
err_sigma = fltarr(nlines)

linemodel=fltarr(n_elements(lambda))
for i=0L, nlines-1L do begin
    jhusdss_lrgfit_line, lambda, flux, continuum, ivar, lines[i], $
      ew=tmp_ew, lflux=tmp_lflux, err_ew=tmp_err_ew, err_lflux=tmp_err_lflux, $
      model=tmp_model, sigma=tmp_sigma, err_sigma=tmp_err_sigma, $
      vdisp=vdisp, outclevel=tmp_clevel
    if (n_elements(tmp_clevel) eq 0) then continue
    ew[i]=tmp_ew
    clevel[i]=tmp_clevel
    lflux[i]=tmp_lflux
    sigma[i]=tmp_sigma
    err_ew[i]=tmp_err_ew
    err_lflux[i]=tmp_err_lflux
    err_sigma[i]=tmp_err_sigma
    linemodel=linemodel+tmp_model
endfor

end
