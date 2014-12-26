
;+
; load and interpolate *ALL* quasar spectra
; All in observer Frame
;-

pro jhusdss_write_alllrgspec_new, lrgver, boss=boss, overwrite=overwrite, doflux=doflux, docont=docont, donorm=donorm, dosubt=dosubt

if (n_elements(lrgver) eq 0) then message, 'lrgver required'

;; lrgpath
parent_path=jhusdss_get_parent_path()
lrgpath = parent_path+'/SDSS/LRG/'+string(lrgver, format='(I3.3)')+'/'

if (keyword_set(doflux)) then outfile = lrgpath+'/AllInOne/'+jhusdss_alllrgspec_filename(lrgver, /flux, boss=boss)
if (keyword_set(docont)) then outfile = lrgpath+'/AllInOne/'+jhusdss_alllrgspec_filename(lrgver, /continuum, boss=boss)
if (keyword_set(donorm)) then outfile = lrgpath+'/AllInOne/'+jhusdss_alllrgspec_filename(lrgver, /normresi, boss=boss)
if (keyword_set(dosubt)) then outfile = lrgpath+'/AllInOne/'+jhusdss_alllrgspec_filename(lrgver, /subtresi, boss=boss)

if (file_test(outfile) and ~keyword_set(overwrite)) then begin
   splog, 'File exists. Use Overwrite to overwrite!'
   return
endif else begin
   splog, 'Will write into these files:'
   print, outfile
endelse

lrg = jhusdss_lrg_readin(boss=boss)
nlrg = n_elements(lrg)

normwave = [5300., 5400.]

minmaxwave = jhusdss_sdsswave_minmax(boss=boss)
loglam = jhusdss_get_loglam(minwave=minmaxwave[0], maxwave=minmaxwave[1], nwave=nwave)
wave = 10.^loglam

if (keyword_set(doflux)) then $
     strtmp = {flux:fltarr(nwave), ivar:fltarr(nwave), norm_ratio:0.}

if (keyword_set(docont)) then $
     strtmp = {continuum:fltarr(nwave), nmf_continuum:fltarr(nwave), med_continuum:fltarr(nwave)}

if (keyword_set(donorm)) then $
     strtmp = {residual:fltarr(nwave), residual_ivar:fltarr(nwave)}

if (keyword_set(dosubt)) then $
     strtmp = {subtracted_residual:fltarr(nwave), subtracted_ivar:fltarr(nwave)}

outstr = replicate(strtmp, nlrg)

for ispec=0L, nlrg-1L do begin

    counter, ispec+1, nlrg

    ;; stupid -- if I saved norm_ratio in decomposition then we don't need this step
    jhusdss_load_interp_spec, lrg[ispec], loglam=loglam, $
            zuse=fltarr(1), boss=boss, allflux=allflux, allivar=allivar

    iwavemin = value_locate(wave, normwave[0]*(1.+lrg[ispec].z))
    iwavemax = value_locate(wave, normwave[1]*(1.+lrg[ispec].z))
    normflux = median(allflux[0, iwavemin:iwavemax])
    
    if (keyword_set(doflux)) then begin
        outstr[ispec].flux = allflux;/normflux
        outstr[ispec].ivar = allivar;*normflux^2/ratio[i]^2
        outstr[ispec].norm_ratio = normflux;*ratio[i]
    endif

    ;; load
    tmpspec = jhusdss_lrg_decompose_loadspec(lrg[ispec].plate, lrg[ispec].fiber, $
                          lrg[ispec].mjd, lrgver, boss=boss, error=error)
    if (error) then continue

    if (keyword_set(dosubt)) then outstr[ispec].subtracted_ivar = allivar

    tmpwave = tmpspec.wave*(1.+tmpspec.z)
    tmpflux = tmpspec.flux
    tmpivar = tmpspec.ivar
    tmpivar_residual = tmpspec.ivar*tmpspec.nmf_continuum^2*tmpspec.med_continuum^2
    tmp_nmf_continuum = tmpspec.nmf_continuum
    tmp_med_continuum = tmpspec.med_continuum
    tmp_residual= tmpspec.residual

    ;; interpolate
    tmpmask = (tmpivar le 0.)
    curr_loglam = alog10(tmpwave)

    ;; get the normalized flux
    combine1fiber, curr_loglam, tmpflux, tmpivar, newloglam=loglam, $
       newflux=finterp, newivar=iinterp, maxiter=0, $
       finalmask=tmpmask, andmask=maskinterp
    tmpratio = median(allflux[0,*]/finterp)

    if (not keyword_set(doflux)) then begin
       tmpivar_2 = fltarr(n_elements(curr_loglam))+1.
       combine1fiber, curr_loglam, tmp_nmf_continuum, tmpivar, newloglam=loglam, $
          newflux=finterp, newivar=iinterp, maxiter=0, $
          finalmask=tmpmask, andmask=maskinterp
       out_tmp_nmf_continuum = finterp*tmpratio

       combine1fiber, curr_loglam, tmp_med_continuum, tmpivar_2, newloglam=loglam, $
          newflux=finterp, newivar=iinterp, maxiter=0, $
          finalmask=tmpmask, andmask=maskinterp
       out_tmp_med_continuum = finterp

       combine1fiber, curr_loglam, tmp_residual, tmpivar_2, newloglam=loglam, $
          newflux=finterp, newivar=iinterp, maxiter=0, $
          finalmask=tmpmask, andmask=maskinterp
       out_tmp_residual = finterp
       out_tmp_residual_ivar = iinterp
    endif

    if (keyword_set(docont)) then begin
        outstr[ispec].nmf_continuum = out_tmp_nmf_continuum
        outstr[ispec].med_continuum = out_tmp_med_continuum
        outstr[ispec].continuum = out_tmp_nmf_continuum*out_tmp_med_continuum
    endif

    if (keyword_set(donorm)) then begin
        outstr[ispec].residual = out_tmp_residual
        outstr[ispec].residual_ivar = out_tmp_residual_ivar
    endif

    if (keyword_set(dosubt)) then begin
        outstr[ispec].subtracted_residual = allflux[0,*]-out_tmp_nmf_continuum*out_tmp_med_continuum
    endif
endfor

mwrfits, outstr, outfile, /create

end
