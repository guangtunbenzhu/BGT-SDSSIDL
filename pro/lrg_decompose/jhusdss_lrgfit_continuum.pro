;+
; NAME:
;   jhusdss_lrgfit_continuum
; PURPOSE:
; CALLING SEQUENCE:
;   continuum = jhusdss_lrgfit_continuum(lambda, flux, ivar)
; REVISION HISTORY:
;   30-Apr-2010  Guangtun, JHU, adopted to jhusdss
;   ??-???-2010  Guangtun, NYU, sfit_continuum2
;   31-Mar-2004  MRB, NYU, sfit_continuum
;-
;------------------------------------------------------------------------------
function jhusdss_lrgfit_continuum, in_lambda, in_flux, in_ivar, coeffs=coeffs, $
    vdisp=vdisp

common sfit_common, intemplates, olambda
lambda=in_lambda
flux=in_flux
ivar=in_ivar

;; deal with badness
ibad=where(ivar ne ivar OR flux ne flux or lambda ne lambda, nbad)
if(nbad gt 0) then begin
    ivar[ibad]=0.
    flux[ibad]=0.
    lambda[ibad]=0.
endif

; lines to mask out
lines = jhusdss_lrg_masklines()

for i=0L, n_elements(lines)-1L do begin
    imask = where(abs(lambda-lines[i].line) lt 10., nmask)
    if (nmask gt 0) then ivar[imask] = 0.
endfor

;; only use <6000 AA
i_nouse = where(lambda gt 6000., n_nouse)
if (n_nouse gt 0) then ivar[i_nouse] = 0.

nt = 15
if (n_elements(olambda) eq 0) then begin
    ;; initialization
    ;; get templates
    k_reconstruct_spec, dum, loglam, /init, vname=vname, nt=knt, /nolines

    olambda=10.^loglam
    intemplates=fltarr(n_elements(olambda), nt)

;; metallicity: (3,4,5 better than 1, 2, 5????) for elliptical galaxies?
;; We just want to replace the bad pixels with continuum
;; shouldn't matter in our case since all the spectra are good
;;         0 - Z=0.0001 (m22)
;;         1 - Z=0.0004 (m32)
;;         2 - Z=0.004  (m42)
;;         3 - Z=0.008  (m52)
;;         4 - Z=0.02   (m62) (Solar, default)
;;         5 - Z=0.05   (m72)

    bc03 = k_im_read_bc03(bc03_extra=bc03_extra, metallicity=2)
    ii = where(abs(bc03.age-19.0E9)/1E9 le 0.01)
    jj = where(abs(10.^bc03_extra.logage-19.0E9)/1E9 le 0.01)
    intemplates[*,0]=primus_binspec(bc03.wave, bc03.flux[*,ii], olambda)/bc03_extra[jj].m_*1.E6

    ii = where(abs(bc03.age-17.0E9)/1E9 le 0.01)
    jj = where(abs(10.^bc03_extra.logage-17.0E9)/1E9 le 0.01)
    intemplates[*,1]=primus_binspec(bc03.wave, bc03.flux[*,ii], olambda)/bc03_extra[jj].m_*1.E6

    ii = where(abs(bc03.age-13.0E9)/1E9 le 0.01)
    jj = where(abs(10.^bc03_extra.logage-13.0E9)/1E9 le 0.01)
    intemplates[*,2]=primus_binspec(bc03.wave, bc03.flux[*,ii], olambda)/bc03_extra[jj].m_*1.E6

    bc03 = k_im_read_bc03(bc03_extra=bc03_extra, metallicity=3)
    ii = where(abs(bc03.age-18.0E9)/1E9 le 0.01)
    jj = where(abs(10.^bc03_extra.logage-18.0E9)/1E9 le 0.01)
    intemplates[*,3]=primus_binspec(bc03.wave, bc03.flux[*,ii], olambda)/bc03_extra[jj].m_*1.E6

    ii = where(abs(bc03.age-16.0E9)/1E9 le 0.01)
    jj = where(abs(10.^bc03_extra.logage-16.0E9)/1E9 le 0.01)
    intemplates[*,4]=primus_binspec(bc03.wave, bc03.flux[*,ii], olambda)/bc03_extra[jj].m_*1.E6

    ii = where(abs(bc03.age-12.0E9)/1E9 le 0.01)
    jj = where(abs(10.^bc03_extra.logage-12.0E9)/1E9 le 0.01)
    intemplates[*,5]=primus_binspec(bc03.wave, bc03.flux[*,ii], olambda)/bc03_extra[jj].m_*1.E6

    bc03 = k_im_read_bc03(bc03_extra=bc03_extra, metallicity=4)
    ii = where(abs(bc03.age-17.0E9)/1E9 le 0.01)
    jj = where(abs(10.^bc03_extra.logage-17.0E9)/1E9 le 0.01)
    intemplates[*,6]=primus_binspec(bc03.wave, bc03.flux[*,ii], olambda)/bc03_extra[jj].m_*1.E6

    ii = where(abs(bc03.age-15.0E9)/1E9 le 0.01)
    jj = where(abs(10.^bc03_extra.logage-15.0E9)/1E9 le 0.01)
    intemplates[*,7]=primus_binspec(bc03.wave, bc03.flux[*,ii], olambda)/bc03_extra[jj].m_*1.E6

    ii = where(abs(bc03.age-11.0E9)/1E9 le 0.01)
    jj = where(abs(10.^bc03_extra.logage-11.0E9)/1E9 le 0.01)
    intemplates[*,8]=primus_binspec(bc03.wave, bc03.flux[*,ii], olambda)/bc03_extra[jj].m_*1.E6

    bc03 = k_im_read_bc03(bc03_extra=bc03_extra, metallicity=5)
    ii = where(abs(bc03.age-19.0E9)/1E9 le 0.01)
    jj = where(abs(10.^bc03_extra.logage-19.0E9)/1E9 le 0.01)
    intemplates[*,9]=primus_binspec(bc03.wave, bc03.flux[*,ii], olambda)/bc03_extra[jj].m_*1.E6

    ii = where(abs(bc03.age-15.0E9)/1E9 le 0.01)
    jj = where(abs(10.^bc03_extra.logage-15.0E9)/1E9 le 0.01)
    intemplates[*,10]=primus_binspec(bc03.wave, bc03.flux[*,ii], olambda)/bc03_extra[jj].m_*1.E6

    ii = where(abs(bc03.age-10.0E9)/1E9 le 0.01)
    jj = where(abs(10.^bc03_extra.logage-10.0E9)/1E9 le 0.01)
    intemplates[*,11]=primus_binspec(bc03.wave, bc03.flux[*,ii], olambda)/bc03_extra[jj].m_*1.E6

    ii = where(abs(bc03.age-7.5E9)/1E9 le 0.01)
    jj = where(abs(10.^bc03_extra.logage-7.5E9)/1E9 le 0.01)
    intemplates[*,12]=primus_binspec(bc03.wave, bc03.flux[*,ii], olambda)/bc03_extra[jj].m_*1.E6

    ii = where(abs(bc03.age-2.5E9)/1E9 le 0.01)
    jj = where(abs(10.^bc03_extra.logage-2.5E9)/1E9 le 0.01)
    intemplates[*,13]=primus_binspec(bc03.wave, bc03.flux[*,ii], olambda)/bc03_extra[jj].m_*1.E6

    ii = where(abs(bc03.age-0.719E9)/1E9 le 0.01)
    jj = where(abs(10.^bc03_extra.logage-0.719E9)/1E9 le 0.01)
    intemplates[*,14]=primus_binspec(bc03.wave, bc03.flux[*,ii], olambda)/bc03_extra[jj].m_*1.E6

endif

nl=n_elements(olambda)
templates = intemplates
init_vdisp = 3./5000.*2.99792458E5
for i=0L, nt-1L do begin
   jhusdss_lrgspec_smooth, olambda, intemplates[*,i], fltarr(nl)+1., $
       init_vdisp=init_vdisp, final_vdisp=vdisp, outflux=outflux
   templates[*,i] = outflux
endfor

;; only keep the old and young
coeffs = fltarr(nt)
;templates = templates[*,nt:nt+1]
;; perform kcorrect fit iteratively, with successive median smooth divisions
tweak=fltarr(n_elements(flux))+1.
for i=0L, 5L do begin 
;   help,templates
    k_fit_spec, flux/tweak, ivar*tweak^2, coeffs, lambda=lambda, $
      templates=templates, olambda=olambda, /nolines, vdisp=vdisp, verbose=0
;   help,templates
;   help,templates, coeffs
    model=templates#coeffs
    continuum=primus_binspec(olambda, model, lambda)
    ratio=flux/continuum
    ibad=where(ivar le 0.,nbad)
    if(nbad gt 0) then ratio[ibad]=1.
    tweak=djs_median(ratio, width=100., boundary='reflect')
endfor

continuum=interpol(model, olambda, lambda)*tweak

return, continuum

end
