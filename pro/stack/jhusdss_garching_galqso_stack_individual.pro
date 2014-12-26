
atmp = 'a'
nmfver = 106
overwrite = 1b
ivarweight = 1b
docut = 0b
qaplot = 0b
qaplot1 = 0b
loaddata = 1b
irp_min = 1
irp_max = 21
drp = 1
saveall = 1b
preprocess = 0b

sigma_smooth = 2.0
sigma_cut = 2.
snr_cut = 3.
z_cut = 0.04

profile_slope = -1.29629
profile_intercept = 0.636813

read,'BOSS? [1=yes, 0=no]: ',BOSS

line_wave = [3934.79, 3969.59]
xra = [3800., 4300.]

stackpath = jhusdss_get_path(/nmfqso)+'/'+string(nmfver, format='(I3.3)')+'/Stack/'
infile = stackpath + jhusdss_garching_galqso_stack_filename(nmfver, boss=boss)
spec_infile = repstr(infile, '.fits', '_spec.fits')

outfile = repstr(spec_infile, '.fits', '_fit.fits')

if (file_test(outfile) and ~overwrite) then begin
   splog, 'File already exists, use /overwrite to overwrite'
   return
endif else begin
   splog, 'Will write into this file: '
   print, outfile
endelse

comp0 = mrdfits(spec_infile,1)
nrp = comp0.nrp

fix_separation = long((line_wave[1]-line_wave[0])/(line_wave[1]+line_wave[0])/alog(10.)*2E+4)

in_quadra = 0.
in_slope = 0.
in_center = line_wave[0]
in_separation = line_wave[1]-line_wave[0]
in_lflux = 0.2
in_ratio = 0.5
in_sigma = 2.0

for i=irp_min, irp_max, drp do begin
    if (keyword_set(loaddata)) then comp = mrdfits(spec_infile, i)

    wave = comp.wave
    nwave = n_elements(wave)
    tmp = min(abs(comp.wave-line_wave[0]), icaii)
    iuse = where(wave gt xra[0] and wave lt xra[1], nuse)
    iuse_nocaii = where((wave gt xra[0] and wave lt (line_wave[0]-10.) $
                     or (wave gt (line_wave[0]+10.) and wave lt (line_wave[1]-10.))) $
                     or (wave gt (line_wave[1]+10.) and wave lt (xra[1])), nuse_nocaii)
    outstr = {nrp:nrp, rp:comp.rp, rp_min:comp.rp_min, rp_max:comp.rp_max, npairs:comp.npairs, $
              wave:comp.wave, $
              singlet_residual:comp.residual, doublet_residual:comp.residual, $
              singlet_ivar:fltarr(comp.npairs, nwave), doublet_ivar:fltarr(comp.npairs, nwave), $
              fluxmean:comp.wave, fluxmedian:comp.wave, fluxgeomean:comp.wave, $
              ;; not preprocess
              fluxmean_continuum:comp.wave, fluxmedian_continuum:comp.wave, fluxgeomean_continuum:comp.wave, $
              fluxmean_singlet:comp.wave, fluxmedian_singlet:comp.wave, fluxgeomean_singlet:comp.wave, $
              fluxmean_doublet:comp.wave, fluxmedian_doublet:comp.wave, fluxgeomean_doublet:comp.wave, $
              ;; preprocess
              singlet_fluxmean:comp.wave, singlet_fluxmedian:comp.wave, singlet_fluxgeomean:comp.wave, $
              doublet_fluxmean:comp.wave, doublet_fluxmedian:comp.wave, doublet_fluxgeomean:comp.wave, $
              ew_nofit:fltarr(comp.npairs), err_ew_nofit:fltarr(comp.npairs), $
              ew_nofit_2:fltarr(comp.npairs), err_ew_nofit_2:fltarr(comp.npairs)}

    
;   tmp_kernel = [-1., 0., 1.]
    for j=0L, comp.npairs-1L do begin
        tmp = reform(comp.residual[j,*])
        tmp_ivar = reform(comp.ivar[j,*])
        if (comp.spec_snr[j] le snr_cut or comp.zgal[j] le z_cut) then begin
           comp.ivar[j,*] = 0. 
           continue
        endif
        ss = where(finite(tmp) eq 0. or finite(tmp_ivar) eq 0., tt)
        if tt gt 0 then tmp_ivar[tt] = 0.
        tmp_use = where(tmp_ivar gt 0., tmp_nuse)
        if tmp_nuse eq 0 then begin
           comp.ivar[j,*] = 0.
           continue
        endif

        if (keyword_set(docut)) then begin
        ;; first iteration
        mom_tmp = moment(tmp[tmp_use], sdev=sdev)
;       mom_tmp_ivar = moment(tmp_ivar[tmp_use], sdev=sdev_ivar)
        tmp_outlier = where(abs(tmp[tmp_use]-1.) gt 5.*sdev, tmp_noutlier)
        if (tmp_noutlier gt 0) then begin
           tmp_ivar[tmp_use[tmp_outlier]] = 0.
           tmp_ivar[tmp_use[tmp_outlier]-1] = 0.
           tmp_ivar[tmp_use[tmp_outlier]-2] = 0.
           tmp_ivar[tmp_use[tmp_outlier]+1] = 0.
           tmp_ivar[tmp_use[tmp_outlier]+2] = 0.
        endif

        ;; second iteration
        tmp_use = where(tmp_ivar gt 0., tmp_nuse)
        if tmp_nuse gt 0 then begin
           mom_tmp = moment(tmp[tmp_use], sdev=sdev)
           tmp_outlier = where(abs(tmp[tmp_use]-1.) gt 4.*sdev, tmp_noutlier)
           if (tmp_noutlier gt 0) then begin
              tmp_ivar[tmp_use[tmp_outlier]] = 0.
              tmp_ivar[tmp_use[tmp_outlier]-1] = 0.
              tmp_ivar[tmp_use[tmp_outlier]-2] = 0.
              tmp_ivar[tmp_use[tmp_outlier]+1] = 0.
              tmp_ivar[tmp_use[tmp_outlier]+2] = 0.
           endif
        endif

        ;; third iteration
        tmp_use = where(tmp_ivar gt 0., tmp_nuse)
        if tmp_nuse gt 0 then begin
           mom_tmp = moment(tmp[tmp_use], sdev=sdev)
           tmp_outlier = where(abs(tmp[tmp_use]-1.) gt 3.*sdev, tmp_noutlier)
           if (tmp_noutlier gt 0) then begin
              tmp_ivar[tmp_use[tmp_outlier]] = 0.
              tmp_ivar[tmp_use[tmp_outlier]-1] = 0.
              tmp_ivar[tmp_use[tmp_outlier]-2] = 0.
              tmp_ivar[tmp_use[tmp_outlier]+1] = 0.
              tmp_ivar[tmp_use[tmp_outlier]+2] = 0.
           endif
        endif

;       tmp_outlier = where((tmp_ivar-convol(tmp_ivar, tmp_kernel, /center, /nan)) gt 5.*sdev_ivar, tmp_noutlier) 
;       if (tmp_noutlier gt 0) then tmp_ivar[tmp_outlier] = 0.
        comp.ivar[j,*] = tmp_ivar
        endif
    endfor

    if (keyword_set(qaplot1)) then begin
    for j=0L, comp.npairs/20-1L do begin
        kmin = j*20
        kmax = kmin+20
        multiplot, [2, 20]
        for k=kmin, kmax-1 do begin
            djs_plot, comp.wave, comp.residual[k,*], xra=[3800, 4100], yra=[0., 2], xst=1, yst=1
            djs_oplot, line_wave[0]*[1,1], !y.crange, xst=1, yst=1
            djs_oplot, line_wave[1]*[1,1], !y.crange, xst=1, yst=1
            multiplot
            djs_plot, comp.wave, comp.ivar[k,*], xra=[3800, 4100], xst=1
            djs_oplot, line_wave[0]*[1,1], !y.crange, xst=1, yst=1
            djs_oplot, line_wave[1]*[1,1], !y.crange, xst=1, yst=1
            multiplot
        endfor
        read, atmp
        if (atmp eq 'q') then stop
        erase
    endfor
    endif
    !p.noerase = 0
    !p.multi = [0, 1, 1]

    if (keyword_set(preprocess)) then begin
    for j=0L, comp.npairs-1L do begin
        y = reform(comp.residual[j, iuse])
        yivar = reform(comp.ivar[j, iuse])
;       ss = where(finite(y) eq 0. or finite(yivar) eq 0., tt)
;       if (tt gt 0) then yivar[ss] = 0.

;       terror = sqrt((moment(1.-y[iuse_nocaii]))[1])

        jhusdss_singlet_smooth, 1.-y, yivar, $
           outflux=single_newy, outivar=single_newyivar, sigma=sigma_smooth, /normalize, factor_norm=factor_norm1
        jhusdss_doublet_smooth, 1.-y, yivar, separation=fix_separation, $
           outflux=double_newy, outivar=double_newyivar, sigma=sigma_smooth, /normalize, factor_norm=factor_norm2

        outstr.singlet_residual[j,iuse] = 1.-single_newy
        outstr.doublet_residual[j,iuse] = 1.-double_newy
        outstr.singlet_ivar[j,iuse] = single_newyivar
        outstr.doublet_ivar[j,iuse] = double_newyivar
        ss = where(finite(single_newyivar) eq 0 or finite(single_newy) eq 0, tt)
        if (tt gt 0) then outstr.singlet_ivar[j,iuse[ss]] = 0.
        ss = where(finite(double_newyivar) eq 0 or finite(double_newy) eq 0, tt)
        if (tt gt 0) then outstr.doublet_ivar[j,iuse[ss]] = 0.
        ;; this is wrong
        outstr.ew_nofit[j] = (1.-outstr.singlet_residual[j, icaii])*(wave[icaii]-wave[icaii-1])/factor_norm1
        outstr.err_ew_nofit[j] = 1./sqrt(outstr.singlet_ivar[j,icaii])*(wave[icaii]-wave[icaii-1])/factor_norm1
        outstr.ew_nofit_2[j] = (1.-outstr.doublet_residual[j,icaii])*(wave[icaii]-wave[icaii-1])/factor_norm2*6./5.
        outstr.err_ew_nofit_2[j] = 1./sqrt(outstr.doublet_ivar[j,icaii])*(wave[icaii]-wave[icaii-1])/factor_norm2*6./5.
    endfor
    endif

;   weight = comp.zgal^2*100.
    weight = comp.zqso*0. + 1.
;   for j=0L, comp.npairs-1L do begin
;   endfor
    jhusdss_composite_engine, comp.residual, comp.ivar, fmean=fmean, fmedian=fmedian, $
       fgeomean=fgeomean, nobjuse=nobjuse, ivarweight=ivarweight, weight=weight, sigma_cut=sigma_cut
    outstr.fluxmean = fmean
    outstr.fluxmedian = fmedian
    outstr.fluxgeomean = fgeomean

    for imean=0,2 do begin
        case imean of
           0: y = reform(fmean)
           1: y = reform(fmedian)
           2: y = reform(fgeomean)
        endcase
        med_continuum = median(y, 21, /even)
        y = (y/med_continuum)[iuse]

        in_intercept = median(1.-y[iuse_nocaii])
        terror = sqrt((moment(1.-y[iuse_nocaii]))[1])
        yivar = fltarr(nuse)+1./terror^2

        ;; get clevel
        for iter=0,1 do begin
           jhusdss_lowz_doublet_fit2, wave[iuse], 1.-y[iuse], fltarr(nuse)+1./terror^2, $
              in_quadra, in_slope, in_intercept, in_center, in_separation, in_lflux, in_ratio, in_sigma, $
              quadra=quadra, slope=slope, intercept=intercept, center=center, separation=separation, $
              lflux=lflux, ratio=ratio, sigma=sigma, $
              err_quadra=err_quadra, err_slope=err_slope, err_intercept=err_intercept, err_center=err_center, $
              err_separation=err_separation, err_lflux=err_lflux, err_ratio=err_ratio, err_sigma=err_sigma, $
              maxwidth=maxwidth

           ;; get continuum residuals
           p = [quadra, slope, intercept, center, separation, lflux, ratio, sigma]
           yfit = jhusdss_lowz_doublet_func2(wave, p)
           clevel = quadra*(wave-center)^2+slope*(wave-center)+intercept

           ;; new typical error
           terror = sqrt((moment(1.-y[iuse_nocaii]-clevel[iuse_nocaii]))[1])
        endfor

        newy = y+clevel
        jhusdss_singlet_smooth, 1.-newy, yivar, $
           outflux=single_newy, outivar=single_newyivar, sigma=sigma_smooth, /normalize, factor_norm=factor_norm1
        jhusdss_doublet_smooth, 1.-newy, yivar, separation=fix_separation, $
           outflux=double_newy, outivar=double_newyivar, sigma=sigma_smooth, /normalize, factor_norm=factor_norm2

        case imean of
           0: begin 
              outstr.fluxmean_singlet[iuse] = 1.-single_newy
              outstr.fluxmean_doublet[iuse] = 1.-double_newy
;             outstr.fluxmean_continuum[iuse] = 1.-single_newy
              end
           1: begin 
              outstr.fluxmedian_singlet[iuse] = 1.-single_newy
              outstr.fluxmedian_doublet[iuse] = 1.-double_newy
;             outstr.fluxmedian_continuum[iuse] = 1.-single_newy
              end
           2: begin 
              outstr.fluxgeomean_singlet[iuse] = 1.-single_newy
              outstr.fluxgeomean_doublet[iuse] = 1.-double_newy
;             outstr.fluxgeomean_continuum[iuse] = 1.-single_newy
              end
        endcase
    endfor

    if (keyword_set(preprocess)) then begin
    jhusdss_composite_engine, outstr.singlet_residual, outstr.singlet_ivar, fmean=fmean, fmedian=fmedian, $
       fgeomean=fgeomean, nobjuse=nobjuse, ivarweight=ivarweight, weight=weight, sigma_cut=sigma_cut
    outstr.singlet_fluxmean = fmean
    outstr.singlet_fluxmedian = fmedian
    outstr.singlet_fluxgeomean = fgeomean
    for imean=0,2 do begin
        case imean of
           0: y = reform(fmean)
           1: y = reform(fmedian)
           2: y = reform(fgeomean)
        endcase
        med_continuum = median(y, 31, /even)
        y = (y/med_continuum)[iuse]

        in_intercept = median(1.-y[iuse_nocaii])
        terror = sqrt((moment(1.-y[iuse_nocaii]))[1])
        yivar = fltarr(nuse)+1./terror^2

        ;; get clevel
        for iter=0,1 do begin
           jhusdss_lowz_doublet_fit2, wave[iuse], 1.-y[iuse], fltarr(nuse)+1./terror^2, $
              in_quadra, in_slope, in_intercept, in_center, in_separation, in_lflux, in_ratio, in_sigma, $
              quadra=quadra, slope=slope, intercept=intercept, center=center, separation=separation, $
              lflux=lflux, ratio=ratio, sigma=sigma, $
              err_quadra=err_quadra, err_slope=err_slope, err_intercept=err_intercept, err_center=err_center, $
              err_separation=err_separation, err_lflux=err_lflux, err_ratio=err_ratio, err_sigma=err_sigma, $
              maxwidth=maxwidth

           ;; get continuum residuals
           p = [quadra, slope, intercept, center, separation, lflux, ratio, sigma]
           yfit = jhusdss_lowz_doublet_func2(wave, p)
           clevel = quadra*(wave-center)^2+slope*(wave-center)+intercept

           ;; new typical error
           terror = sqrt((moment(1.-y[iuse_nocaii]-clevel[iuse_nocaii]))[1])
        endfor

        newy = y+clevel

        case imean of
           0: outstr.singlet_fluxmean[iuse] = newy
           1: outstr.singlet_fluxmedian[iuse] = newy
           2: outstr.singlet_fluxgeomean[iuse] = newy
        endcase
    endfor
    endif

    if (keyword_set(preprocess)) then begin
    jhusdss_composite_engine, outstr.doublet_residual, outstr.doublet_ivar, fmean=fmean, fmedian=fmedian, $
       fgeomean=fgeomean, nobjuse=nobjuse, ivarweight=ivarweight, weight=weight, sigma_cut=sigma_cut
    outstr.doublet_fluxmean = fmean
    outstr.doublet_fluxmedian = fmedian
    outstr.doublet_fluxgeomean = fgeomean
    for imean=0,2 do begin
        case imean of
           0: y = reform(fmean)
           1: y = reform(fmedian)
           2: y = reform(fgeomean)
        endcase
        med_continuum = median(y, 31, /even)
        y = (y/med_continuum)[iuse]

        in_intercept = median(1.-y[iuse_nocaii])
        terror = sqrt((moment(1.-y[iuse_nocaii]))[1])
        yivar = fltarr(nuse)+1./terror^2

        ;; get clevel
        for iter=0,1 do begin
           jhusdss_lowz_doublet_fit2, wave[iuse], 1.-y[iuse], fltarr(nuse)+1./terror^2, $
              in_quadra, in_slope, in_intercept, in_center, in_separation, in_lflux, in_ratio, in_sigma, $
              quadra=quadra, slope=slope, intercept=intercept, center=center, separation=separation, $
              lflux=lflux, ratio=ratio, sigma=sigma, $
              err_quadra=err_quadra, err_slope=err_slope, err_intercept=err_intercept, err_center=err_center, $
              err_separation=err_separation, err_lflux=err_lflux, err_ratio=err_ratio, err_sigma=err_sigma, $
              maxwidth=maxwidth

           ;; get continuum residuals
           p = [quadra, slope, intercept, center, separation, lflux, ratio, sigma]
           yfit = jhusdss_lowz_doublet_func2(wave, p)
           clevel = quadra*(wave-center)^2+slope*(wave-center)+intercept

           ;; new typical error
           terror = sqrt((moment(1.-y[iuse_nocaii]-clevel[iuse_nocaii]))[1])
        endfor

        newy = y+clevel

        case imean of
           0: outstr.doublet_fluxmean[iuse] = newy
           1: outstr.doublet_fluxmedian[iuse] = newy
           2: outstr.doublet_fluxgeomean[iuse] = newy
        endcase

;       print, terror
;       djs_plot, wave[iuse], y, xst=1, yst=1, xra=xra, yra=1.+[-3.*terror, 3.*terror]
;       djs_oplot, wave[iuse], newy, color='red'
;       read, atmp
    endfor
    endif

    if (saveall) then mwrfits, outstr, outfile, create=(i eq 1)
;   ;; QAplots
    if (keyword_set(qaplot)) then begin
    load_dp, /b
    print, outstr.rp, outstr.rp_min, outstr.rp_max, outstr.npairs
    mom = moment(outstr.fluxmean[iuse], sdev=sdev)
    mom_2 = moment(outstr.fluxmean_doublet[iuse], sdev=sdev_2)
;   err = sdev/sqrt(comp.npairs-1.)
;   err_2 = sdev_2/sqrt(comp.npairs-1.)
    err = sdev*3.
    err_2 = sdev_2*3
    print, mom[0], sdev, err, mom[0]/err
    print, mom_2[0], sdev_2, err_2, mom_2[0]/err_2

if (keyword_set(preprocess)) then begin
    plothist, outstr.ew_nofit, bin=sdev/20., xra=[-5.*sdev, 5.*sdev]
    read, atmp
    if (atmp eq 'q') then stop
    plothist, outstr.ew_nofit_2, bin=sdev_2/20., xra=[-5.*sdev, 5.*sdev]

    read, atmp
    if (atmp eq 'q') then stop
endif
;   djs_plot, wave, outstr.singlet_fluxmean, xst=1, yst=1, xra=[3850, 4050], yra=1.+[-3.*err, 3.*err]
    djs_plot, wave, smooth(outstr.fluxmean, 3), xst=1, yst=1, xra=[3850, 4050], yra=1.+[-2.*err, 1.*err], psym=10
    djs_oplot, line_wave[0]*[1,1], !y.crange, linestyle=1, color='green'
    djs_oplot, line_wave[1]*[1,1], !y.crange, linestyle=1, color='green'
    djs_oplot, 3990.*[1,1], !y.crange, linestyle=1
    djs_oplot, 3882.*[1,1], !y.crange, linestyle=1
    djs_oplot, !x.crange, [1,1], linestyle=1

    read, atmp
    if (atmp eq 'q') then stop
;   djs_plot, wave, outstr.singlet_fluxmean, xst=1, yst=1, xra=[3850, 4050], yra=1.+[-3.*err, 3.*err]
    djs_plot, wave, outstr.fluxgeomean_singlet, xst=1, yst=1, xra=[3850, 4050], yra=1.+[-1.*err, 1.*err], psym=10, title='geomean singlet'
    djs_oplot, wave, outstr.singlet_fluxgeomean, psym=10, color='red'
    djs_oplot, !x.crange, 1.-factor_norm1/2.*10.^(alog10(comp.rp)*profile_slope+profile_intercept)*[1,1], thick=thick
    djs_oplot, !x.crange, 1.+factor_norm1/2.*10.^(alog10(comp.rp)*profile_slope+profile_intercept)*[1,1], thick=thick
    djs_oplot, line_wave[0]*[1,1], !y.crange, linestyle=1, color='green'
    djs_oplot, line_wave[1]*[1,1], !y.crange, linestyle=1, color='green'
    djs_oplot, 3990.*[1,1], !y.crange, linestyle=1
    djs_oplot, 3882.*[1,1], !y.crange, linestyle=1
    djs_oplot, !x.crange, [1,1], linestyle=1

    read, atmp
    if (atmp eq 'q') then stop
    djs_plot, wave, outstr.fluxgeomean_doublet, xst=1, yst=1, xra=[3850, 4050], yra=1.+[-1.*err, 1.*err], psym=10, title='geomean doublet'
    djs_oplot, !x.crange, 1.-factor_norm1/2.*10.^(alog10(comp.rp)*profile_slope+profile_intercept)*[1,1], thick=thick
    djs_oplot, !x.crange, 1.+factor_norm1/2.*10.^(alog10(comp.rp)*profile_slope+profile_intercept)*[1,1], thick=thick
    djs_oplot, wave, outstr.doublet_fluxgeomean, psym=10, color='red'
    djs_oplot, line_wave[0]*[1,1], !y.crange, linestyle=1, color='green'
;   djs_oplot, line_wave[1]*[1,1], !y.crange, linestyle=1, color='green'
    djs_oplot, 3990.*[1,1], !y.crange, linestyle=1
    djs_oplot, 3882.*[1,1], !y.crange, linestyle=1
    djs_oplot, !x.crange, [1,1], linestyle=1

    read, atmp
    if (atmp eq 'q') then stop
;   djs_plot, wave, outstr.singlet_fluxmean, xst=1, yst=1, xra=[3850, 4050], yra=1.+[-3.*err, 3.*err]
    djs_plot, wave, outstr.fluxmedian_singlet, xst=1, yst=1, xra=[3850, 4050], yra=1.+[-1.*err, 1.*err], psym=10, title='median singlet'
    djs_oplot, wave, outstr.singlet_fluxmedian, psym=10, color='red'
    djs_oplot, !x.crange, 1.-factor_norm1/2.*10.^(alog10(comp.rp)*profile_slope+profile_intercept)*[1,1], thick=thick
    djs_oplot, !x.crange, 1.+factor_norm1/2.*10.^(alog10(comp.rp)*profile_slope+profile_intercept)*[1,1], thick=thick
    djs_oplot, line_wave[0]*[1,1], !y.crange, linestyle=1, color='green'
    djs_oplot, line_wave[1]*[1,1], !y.crange, linestyle=1, color='green'
    djs_oplot, 3990.*[1,1], !y.crange, linestyle=1
    djs_oplot, 3882.*[1,1], !y.crange, linestyle=1
    djs_oplot, !x.crange, [1,1], linestyle=1

    read, atmp
    if (atmp eq 'q') then stop
;   djs_plot, wave, outstr.singlet_fluxmean, xst=1, yst=1, xra=[3850, 4050], yra=1.+[-3.*err, 3.*err]
    djs_plot, wave, outstr.fluxmedian_doublet, xst=1, yst=1, xra=[3850, 4050], yra=1.+[-1.*err, 1.*err], psym=10, title='median doublet'
    djs_oplot, wave, outstr.doublet_fluxmedian, psym=10, color='red'
    djs_oplot, !x.crange, 1.-factor_norm1/2.*10.^(alog10(comp.rp)*profile_slope+profile_intercept)*[1,1], thick=thick
    djs_oplot, !x.crange, 1.+factor_norm1/2.*10.^(alog10(comp.rp)*profile_slope+profile_intercept)*[1,1], thick=thick
    djs_oplot, line_wave[0]*[1,1], !y.crange, linestyle=1, color='green'
;   djs_oplot, line_wave[1]*[1,1], !y.crange, linestyle=1, color='green'
    djs_oplot, 3990.*[1,1], !y.crange, linestyle=1
    djs_oplot, 3882.*[1,1], !y.crange, linestyle=1
    djs_oplot, !x.crange, [1,1], linestyle=1


;   read, atmp
;   djs_plot, wave, outstr.singlet_fluxmean, xst=1, yst=1, xra=[3850, 4050], yra=1.+[-3.*err, 3.*err]
;   djs_plot, wave, outstr.singlet_fluxmedian, xst=1, yst=1, xra=[3850, 4050], yra=1.+[-3.*err, 3.*err], psym=10
;   djs_oplot, line_wave[0]*[1,1], !y.crange, linestyle=1
;   djs_oplot, line_wave[1]*[1,1], !y.crange, linestyle=1

    read, atmp
    if (atmp eq 'q') then stop
;   djs_plot, wave, outstr.singlet_fluxmean, xst=1, yst=1, xra=[3850, 4050], yra=1.+[-3.*err, 3.*err]
;   djs_plot, wave, outstr.fluxmedian_doublet, xst=1, yst=1, xra=[3850, 4050], yra=1.+[-1.*err, 1.*err], psym=10
;   djs_oplot, wave, outstr.doublet_fluxmedian, psym=10, color='red'
    djs_plot, wave, outstr.fluxmean_singlet, xst=1, yst=1, xra=[3850, 4050], yra=1.+[-1.*err, 1.*err], psym=10, title='mean singlet'
    djs_oplot, wave, outstr.singlet_fluxmean, psym=10, color='red'
    djs_oplot, !x.crange, 1.-factor_norm1/2.*10.^(alog10(comp.rp)*profile_slope+profile_intercept)*[1,1], thick=thick
    djs_oplot, !x.crange, 1.+factor_norm1/2.*10.^(alog10(comp.rp)*profile_slope+profile_intercept)*[1,1], thick=thick
    djs_oplot, line_wave[0]*[1,1], !y.crange, linestyle=1, color='green'
    djs_oplot, line_wave[1]*[1,1], !y.crange, linestyle=1, color='green'
    djs_oplot, 3990.*[1,1], !y.crange, linestyle=1
    djs_oplot, 3882.*[1,1], !y.crange, linestyle=1
    djs_oplot, !x.crange, [1,1], linestyle=1

    read, atmp
    if (atmp eq 'q') then stop
;   djs_plot, wave, outstr.singlet_fluxmean, xst=1, yst=1, xra=[3850, 4050], yra=1.+[-3.*err, 3.*err]
;   djs_plot, wave, outstr.fluxmedian_doublet, xst=1, yst=1, xra=[3850, 4050], yra=1.+[-1.*err, 1.*err], psym=10
;   djs_oplot, wave, outstr.doublet_fluxmedian, psym=10, color='red'
    djs_plot, wave, outstr.fluxmean_doublet, xst=1, yst=1, xra=[3850, 4050], yra=1.+[-1.*err, 1.*err], psym=10, title='mean doublet'
    djs_oplot, wave, outstr.doublet_fluxmean, psym=10, color='red'
    djs_oplot, !x.crange, 1.-factor_norm1/2.*10.^(alog10(comp.rp)*profile_slope+profile_intercept)*[1,1], thick=thick
    djs_oplot, !x.crange, 1.+factor_norm1/2.*10.^(alog10(comp.rp)*profile_slope+profile_intercept)*[1,1], thick=thick
    djs_oplot, line_wave[0]*[1,1], !y.crange, linestyle=1, color='green'
;   djs_oplot, line_wave[1]*[1,1], !y.crange, linestyle=1, color='green'
    djs_oplot, 3990.*[1,1], !y.crange, linestyle=1
    djs_oplot, 3882.*[1,1], !y.crange, linestyle=1
    djs_oplot, !x.crange, [1,1], linestyle=1


    read, atmp
    if (atmp eq 'q') then stop
;   djs_plot, wave, outstr.doublet_fluxmean, xst=1, yst=1, xra=[3850, 4050], yra=1.+[-3.*err_2, 3.*err_2]
;   djs_plot, wave, outstr.doublet_fluxmedian, xst=1, yst=1, xra=[3850, 4050], yra=1.+[-3.*err, 3.*err], psym=10
;   djs_oplot, line_wave[0]*[1,1], !y.crange, linestyle=1
;   djs_oplot, line_wave[1]*[1,1], !y.crange, linestyle=1
;   read, atmp
    endif
endfor


end
