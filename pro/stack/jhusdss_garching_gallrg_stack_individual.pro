;pro jhusdss_lowz_caii_new

lrgver = 101
overwrite = 1b
ivarweight = 1b
read,'BOSS? [1=yes, 0=no]: ',BOSS

line_wave = [3934.79, 3969.59]
xra = [3800., 4300.]

stackpath = jhusdss_get_path(/fitlrg)+'/'+string(lrgver, format='(I3.3)')+'/Stack/'
infile = stackpath + jhusdss_garching_gallrg_stack_filename(lrgver, boss=boss)
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
sigma = 1.5
in_quadra = 0.
in_slope = 0.
in_center = line_wave[0]
in_separation = line_wave[1]-line_wave[0]
in_lflux = 0.2
in_ratio = 0.5
in_sigma = 2.0

for i=2, 2 do begin
    comp = mrdfits(spec_infile, i)

    wave = comp.wave
    nwave = n_elements(wave)
    tmp = min(abs(comp.wave-line_wave[0]), icaii)
    iuse = where(wave gt xra[0] and wave lt xra[1], nuse)
    iuse_nocaii = where((wave gt xra[0] and wave lt (line_wave[0]-10.) $
                     or (wave gt (line_wave[0]+10.) and wave lt (line_wave[1]-10.))) $
                     or (wave gt (line_wave[1]+10.) and wave lt (xra[1])), nuse_nocaii)
    outstr = {nrp:nrp, rp:comp.rp, rp_min:comp.rp_min, rp_max:comp.rp_max, npairs:comp.npairs, $
              wave:comp.wave, singlet_residual:comp.residual, doublet_residual:comp.residual, $
              singlet_ivar:fltarr(comp.npairs, nwave), doublet_ivar:fltarr(comp.npairs, nwave), $
              fluxmean:comp.wave, fluxmedian:comp.wave, fluxgeomean:comp.wave, $
              fluxmean_singlet:comp.wave, fluxmedian_singlet:comp.wave, fluxgeomean_singlet:comp.wave, $
              fluxmean_doublet:comp.wave, fluxmedian_doublet:comp.wave, fluxgeomean_doublet:comp.wave, $
              singlet_fluxmean:comp.wave, singlet_fluxmedian:comp.wave, singlet_fluxgeomean:comp.wave, $
              doublet_fluxmean:comp.wave, doublet_fluxmedian:comp.wave, doublet_fluxgeomean:comp.wave, $
              ew_nofit:fltarr(comp.npairs), err_ew_nofit:fltarr(comp.npairs), $
              ew_nofit_2:fltarr(comp.npairs), err_ew_nofit_2:fltarr(comp.npairs)}

    for j=0L, comp.npairs-1L do begin
        y = reform(comp.residual[j, iuse])
        yivar = reform(comp.ivar[j, iuse])
        ss = where(finite(y) eq 0. or finite(yivar) eq 0., tt)
        if (tt gt 0) then yivar[ss] = 0.
        if (median(yivar) le 0.) then continue
;       terror = sqrt((moment(1.-y[iuse_nocaii]))[1])

        jhusdss_singlet_smooth, 1.-y, yivar, $
           outflux=single_newy, outivar=single_newyivar, sigma=sigma, /normalize, factor_norm=factor_norm1
        jhusdss_doublet_smooth, 1.-y, yivar, separation=fix_separation, $
           outflux=double_newy, outivar=double_newyivar, sigma=sigma, /normalize, factor_norm=factor_norm2

        outstr.singlet_residual[j,iuse] = 1.-single_newy
        outstr.doublet_residual[j,iuse] = 1.-double_newy
        outstr.singlet_ivar[j,iuse] = single_newyivar
        outstr.doublet_ivar[j,iuse] = double_newyivar
        ss = where(finite(single_newyivar) eq 0 or finite(single_newy) eq 0, tt)
        if (tt gt 0) then outstr.singlet_ivar[j,iuse[ss]] = 0.
        ss = where(finite(double_newyivar) eq 0 or finite(double_newy) eq 0, tt)
        if (tt gt 0) then outstr.doublet_ivar[j,iuse[ss]] = 0.
        outstr.ew_nofit[j] = single_newy[icaii]*(wave[icaii]-wave[icaii-1])/factor_norm1
        outstr.err_ew_nofit[j] = 1./sqrt(single_newyivar[icaii])*(wave[icaii]-wave[icaii-1])/factor_norm1
        outstr.ew_nofit_2[j] = double_newy[icaii]*(wave[icaii]-wave[icaii-1])/factor_norm2*6./5.
        outstr.err_ew_nofit_2[j] = 1./sqrt(double_newyivar[icaii])*(wave[icaii]-wave[icaii-1])/factor_norm2*6./5.
    endfor

    jhusdss_composite_engine, comp.residual, comp.ivar, fmean=fmean, fmedian=fmedian, $
       fgeomean=fgeomean, nobjuse=nobjuse, ivarweight=ivarweight, weight=weight, sigma_cut=10.
    outstr.fluxmean = fmean
    outstr.fluxmedian = fmedian
    outstr.fluxgeomean = fgeomean
    for imean=0,2 do begin
        case imean of
           0: y = reform(fmean[iuse])
           1: y = reform(fmedian[iuse])
           2: y = reform(fgeomean[iuse])
        endcase
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
           outflux=single_newy, outivar=single_newyivar, sigma=sigma, /normalize, factor_norm=factor_norm1
        jhusdss_doublet_smooth, 1.-newy, yivar, separation=fix_separation, $
           outflux=double_newy, outivar=double_newyivar, sigma=sigma, /normalize, factor_norm=factor_norm2

        case imean of
           0: begin 
              outstr.fluxmean_singlet[iuse] = 1.-single_newy
              outstr.fluxmean_doublet[iuse] = 1.-double_newy
              end
           1: begin 
              outstr.fluxmedian_singlet[iuse] = 1.-single_newy
              outstr.fluxmedian_doublet[iuse] = 1.-double_newy
              end
           2: begin 
              outstr.fluxgeomean_singlet[iuse] = 1.-single_newy
              outstr.fluxgeomean_doublet[iuse] = 1.-double_newy
              end
        endcase
    endfor


    jhusdss_composite_engine, outstr.singlet_residual, outstr.singlet_ivar, fmean=fmean, fmedian=fmedian, $
       fgeomean=fgeomean, nobjuse=nobjuse, ivarweight=ivarweight, weight=weight, sigma_cut=4.
    outstr.singlet_fluxmean = fmean
    outstr.singlet_fluxmedian = fmedian
    outstr.singlet_fluxgeomean = fgeomean
    for imean=0,2 do begin
        case imean of
           0: y = reform(fmean[iuse])
           1: y = reform(fmedian[iuse])
           2: y = reform(fgeomean[iuse])
        endcase
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

    atmp = 'a'
    jhusdss_composite_engine, outstr.doublet_residual, outstr.doublet_ivar, fmean=fmean, fmedian=fmedian, $
       fgeomean=fgeomean, nobjuse=nobjuse, ivarweight=ivarweight, weight=weight, sigma_cut=4.
    outstr.doublet_fluxmean = fmean
    outstr.doublet_fluxmedian = fmedian
    outstr.doublet_fluxgeomean = fgeomean
    for imean=0,2 do begin
        case imean of
           0: y = reform(fmean[iuse])
           1: y = reform(fmedian[iuse])
           2: y = reform(fgeomean[iuse])
        endcase
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

;   ;; QAplots
    load_dp, /b
    print, outstr.rp, outstr.rp_min, outstr.rp_max, outstr.npairs
    mom = moment(outstr.ew_nofit, sdev=sdev)
    mom_2 = moment(outstr.ew_nofit_2, sdev=sdev_2)
    err = sdev/sqrt(comp.npairs-1.)
    err_2 = sdev_2/sqrt(comp.npairs-1.)
    print, mom[0], sdev, err, mom[0]/err
    print, mom_2[0], sdev_2, err_2, mom_2[0]/err_2

    atmp = 'a'
    plothist, outstr.ew_nofit, bin=sdev/20., xra=[-5.*sdev, 5.*sdev]
    read, atmp
    plothist, outstr.ew_nofit_2, bin=sdev_2/20., xra=[-5.*sdev, 5.*sdev]

    read, atmp
;   djs_plot, wave, outstr.singlet_fluxmean, xst=1, yst=1, xra=[3850, 4050], yra=1.+[-3.*err, 3.*err]
    djs_plot, wave, outstr.fluxmean, xst=1, yst=1, xra=[3850, 4050], yra=1.+[-1.*err, 1.*err], psym=10
    djs_oplot, line_wave[0]*[1,1], !y.crange, linestyle=1
    djs_oplot, line_wave[1]*[1,1], !y.crange, linestyle=1

    read, atmp
;   djs_plot, wave, outstr.singlet_fluxmean, xst=1, yst=1, xra=[3850, 4050], yra=1.+[-3.*err, 3.*err]
    djs_plot, wave, outstr.fluxmedian_singlet, xst=1, yst=1, xra=[3850, 4050], yra=1.+[-1.*err, 1.*err], psym=10
    djs_oplot, wave, outstr.singlet_fluxmedian, psym=10, color='red'
    djs_oplot, line_wave[0]*[1,1], !y.crange, linestyle=1
    djs_oplot, line_wave[1]*[1,1], !y.crange, linestyle=1

;   read, atmp
;   djs_plot, wave, outstr.singlet_fluxmean, xst=1, yst=1, xra=[3850, 4050], yra=1.+[-3.*err, 3.*err]
;   djs_plot, wave, outstr.singlet_fluxmedian, xst=1, yst=1, xra=[3850, 4050], yra=1.+[-3.*err, 3.*err], psym=10
;   djs_oplot, line_wave[0]*[1,1], !y.crange, linestyle=1
;   djs_oplot, line_wave[1]*[1,1], !y.crange, linestyle=1

    read, atmp
;   djs_plot, wave, outstr.singlet_fluxmean, xst=1, yst=1, xra=[3850, 4050], yra=1.+[-3.*err, 3.*err]
    djs_plot, wave, outstr.fluxmedian_doublet, xst=1, yst=1, xra=[3850, 4050], yra=1.+[-1.*err, 1.*err], psym=10
    djs_oplot, wave, outstr.doublet_fluxmedian, psym=10, color='red'
    djs_oplot, line_wave[0]*[1,1], !y.crange, linestyle=1
    djs_oplot, line_wave[1]*[1,1], !y.crange, linestyle=1

    read, atmp
;   djs_plot, wave, outstr.doublet_fluxmean, xst=1, yst=1, xra=[3850, 4050], yra=1.+[-3.*err_2, 3.*err_2]
;   djs_plot, wave, outstr.doublet_fluxmedian, xst=1, yst=1, xra=[3850, 4050], yra=1.+[-3.*err, 3.*err], psym=10
;   djs_oplot, line_wave[0]*[1,1], !y.crange, linestyle=1
;   djs_oplot, line_wave[1]*[1,1], !y.crange, linestyle=1
;   read, atmp
endfor

;mwrfits, outstr, outfile, /create

end
