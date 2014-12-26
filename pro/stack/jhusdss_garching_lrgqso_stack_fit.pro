;pro jhusdss_lowz_caii_new

Coarse = 0b
DoWeight = 0b
DoNaI = 0b
DoScale = 0b
Scale_factor = 0.2 ;; M* \propto M(halo)^1 \propto r(halo)^3 and r(halo) \propto M*^(1/3); Scale to 10^10.1 Msun
fid_mass = 10.3 ;; M* \propto M(halo)^1 \propto r(halo)^3 and r(halo) \propto M*^(1/3); Scale to 10^10.1 Msun
read,'Coarse? [1=yes, 0=no]: ', Coarse
read,'DoWeight? [1=yes, 0=no]: ', DoWeight
read,'DoNaI? [1=yes, 0=no]: ', DoNaI
read,'DoScale? [1=yes, 0=no]: ', DoScale

nmfver = 106
overwrite = 1b
;read,'BOSS? [1=yes, 0=no]: ',BOSS

line_wave = [2796.35, 2803.53]
xra = [2600., 3050.]

fix_separation = long((line_wave[1]-line_wave[0])/(line_wave[1]+line_wave[0])/alog(10.)*2E+4)

stackpath = jhusdss_get_path(/nmfqso)+'/'+string(nmfver, format='(I3.3)')+'/Stack/'
infile = stackpath + jhusdss_garching_lrgqso_stack_filename(nmfver, boss=boss)

if (Coarse) then begin
   infile=repstr(infile,'.fits','_coarse.fits')
endif
if (DoWeight) then begin
   infile=repstr(infile,'.fits','_weight.fits')
endif
if (DoNaI) then begin
   infile=repstr(infile,'.fits','_NaI.fits')
endif
if (DoScale) then begin
   infile=repstr(infile,'.fits','_Scale.fits')
endif

outfile = repstr(infile, '.fits', '_fit.fits')

if (file_test(outfile) and ~overwrite) then begin
   splog, 'File already exists, use /overwrite to overwrite'
   return
endif else begin
   splog, 'Will write into this file: '
   print, outfile
endelse

comp = mrdfits(infile,1)
nrp = n_elements(comp)

nmc = 200
mc_line_wave = findgen(nmc)*2.0+line_wave[0]-200.
print, mc_line_wave[0], mc_line_wave[nmc-1]

outstr = replicate({rp:0., ew:0., err_ew:0., sigma1:0., err_sigma1:0., sigma2:0., err_sigma2:0., line_ratio:0., err_line_ratio:0., $
            wave_mc:mc_line_wave, ew_mc:fltarr(nmc)-100., sigma_mc:fltarr(nmc)-100., line_ratio_mc:fltarr(nmc)-100., sdev_mc:-100., sdev_sigma_mc:-100., $
            ew_nofit:0., err_ew_nofit:0., ew_nofit_mc:fltarr(nmc)-100., sdev_nofit_mc:-100., $
            ew_nofit_2:0., err_ew_nofit_2:0., ew_nofit_2_mc:fltarr(nmc)-100., sdev_nofit_2_mc:-100., $
            ew_nofit_lr_1:0., err_ew_nofit_lr_1:0., ew_nofit_lr_1_mc:fltarr(nmc)-100., sdev_nofit_lr_1_mc:-100., $
            ew_nofit_s4:0., err_ew_nofit_s4:0., ew_nofit_s4_mc:fltarr(nmc)-100., sdev_nofit_s4_mc:-100., $
            ew_nofit_2_s4:0., err_ew_nofit_2_s4:0., ew_nofit_2_s4_mc:fltarr(nmc)-100., sdev_nofit_2_s4_mc:-100., $
            ew_nofit_lr_1_s4:0., err_ew_nofit_lr_1_s4:0., ew_nofit_lr_1_s4_mc:fltarr(nmc)-100., sdev_nofit_lr_1_s4_mc:-100., $
            ew_nofit_s5:0., err_ew_nofit_s5:0., ew_nofit_s5_mc:fltarr(nmc)-100., sdev_nofit_s5_mc:-100., $
            ew_nofit_2_s5:0., err_ew_nofit_2_s5:0., ew_nofit_2_s5_mc:fltarr(nmc)-100., sdev_nofit_2_s5_mc:-100., $
            ew_nofit_lr_1_s5:0., err_ew_nofit_lr_1_s5:0., ew_nofit_lr_1_s5_mc:fltarr(nmc)-100., sdev_nofit_lr_1_s5_mc:-100. $
            }, nrp)

;; Bordoloi et al. (2011)
;rp_bor = [40., 60., 75.]
;ew_bor = [0.20, 0.13, 0.09]
rp_bor = [30., 55., 75., 110.]
ew_bor = [0.165, 0.13, 0.18, 0.075]

for irp=1L, nrp-1L do begin
    rp = comp[irp].rp
    wave = comp[irp].wave
    tmp = min(abs(comp[irp].wave-line_wave[0]), icaii)
    tmp = min(abs(comp[irp].wave-line_wave[1]), ired)
    y = comp[irp].fluxgeomean
    med_continuum = median(y, 71, /even)
    y = (y/med_continuum)

    iuse = where(wave gt xra[0] and wave lt xra[1], nuse)
    iuse_nocaii = where((wave gt xra[0] and wave lt (line_wave[0]-10.) $
                     or (wave gt (line_wave[0]+10.) and wave lt (line_wave[1]-10.))) $
                     or (wave gt (line_wave[1]+10.) and wave lt (xra[1])), nuse_nocaii)

      ;; first fit

    in_quadra = 0.
    in_slope = 0.
    in_center = line_wave[0]
    in_separation = line_wave[1]-line_wave[0]
    in_lflux = 0.2
    in_ratio = 0.5
    in_sigma1 = 2.0
    in_sigma2 = 3.0

    in_intercept = median(1.-y[iuse_nocaii])
    terror = sqrt((moment(1.-y[iuse_nocaii]))[1])

    for iter=0,1 do begin
        jhusdss_lowz_doublet_fit3, wave[iuse], 1.-y[iuse], fltarr(nuse)+1./terror^2, $
           in_quadra, in_slope, in_intercept, in_center, in_separation, in_lflux, in_ratio, in_sigma1, in_sigma2, $
           quadra=quadra, slope=slope, intercept=intercept, center=center, separation=separation, $
           lflux=lflux, ratio=ratio, sigma1=sigma1, sigma2=sigma2, $
           err_quadra=err_quadra, err_slope=err_slope, err_intercept=err_intercept, err_center=err_center, $
           err_separation=err_separation, err_lflux=err_lflux, err_ratio=err_ratio, err_sigma1=err_sigma1, err_sigma2=err_sigma2, $
           maxwidth=maxwidth;, inratio_range=[1/3.,1/0.8]

        ;; get continuum residuals
        p = [quadra, slope, intercept, center, separation, lflux, ratio, sigma1, sigma2]
        yfit = jhusdss_lowz_doublet_func3(wave, p)
        clevel = quadra*(wave-center)^2+slope*(wave-center)+intercept

        ;; new typical error
        terror = sqrt((moment(1.-y[iuse_nocaii]-clevel[iuse_nocaii]))[1])
    endfor
    outstr[irp].rp = rp
    outstr[irp].ew = lflux
    outstr[irp].err_ew = err_lflux
    outstr[irp].sigma1 = sigma1
    outstr[irp].err_sigma1 = err_sigma1
    outstr[irp].sigma2 = sigma2
    outstr[irp].err_sigma2 = err_sigma2
    outstr[irp].line_ratio = 1./ratio
    outstr[irp].err_line_ratio = err_ratio/ratio

    newy = y+clevel

    jhusdss_singlet_smooth, yfit, fltarr(n_elements(y))+1./terror^2, $
       outflux=single_yfit, outivar=single_newyivar, sigma=2., /normalize, factor_norm=factor_norm1
    jhusdss_singlet_smooth, 1.-newy, fltarr(n_elements(y))+1./terror^2, $
       outflux=single_newy, outivar=single_newyivar, sigma=2., /normalize, factor_norm=factor_norm1
    jhusdss_doublet_smooth, 1.-newy, fltarr(n_elements(y))+1./terror^2, $
       outflux=double_newy, outivar=double_newyivar, sigma=2., /normalize, factor_norm=factor_norm2, separation=fix_separation, line_ratio=2.
    jhusdss_doublet_smooth, 1.-newy, fltarr(n_elements(y))+1./terror^2, $
       outflux=double_newy_1, outivar=double_newyivar_1, sigma=2., /normalize, factor_norm=factor_norm2_1, separation=fix_separation, line_ratio=1.

    jhusdss_singlet_smooth, 1.-newy, fltarr(n_elements(y))+1./terror^2, $
       outflux=single_newy_s4, outivar=single_newyivar_s4, sigma=3., /normalize, factor_norm=factor_norm1_s4
    jhusdss_doublet_smooth, 1.-newy, fltarr(n_elements(y))+1./terror^2, $
       outflux=double_newy_s4, outivar=double_newyivar_s4, sigma=3., /normalize, factor_norm=factor_norm2_s4, separation=fix_separation, line_ratio=2.
    jhusdss_doublet_smooth, 1.-newy, fltarr(n_elements(y))+1./terror^2, $
       outflux=double_newy_1_s4, outivar=double_newyivar_1_s4, sigma=3., /normalize, factor_norm=factor_norm2_1_s4, separation=fix_separation, line_ratio=1.

    jhusdss_singlet_smooth, 1.-newy, fltarr(n_elements(y))+1./terror^2, $
       outflux=single_newy_s5, outivar=single_newyivar_s5, sigma=4., /normalize, factor_norm=factor_norm1_s5
    jhusdss_doublet_smooth, 1.-newy, fltarr(n_elements(y))+1./terror^2, $
       outflux=double_newy_s5, outivar=double_newyivar_s5, sigma=4., /normalize, factor_norm=factor_norm2_s5, separation=fix_separation, line_ratio=2.
    jhusdss_doublet_smooth, 1.-newy, fltarr(n_elements(y))+1./terror^2, $
       outflux=double_newy_1_s5, outivar=double_newyivar_1_s5, sigma=4., /normalize, factor_norm=factor_norm2_1_s5, separation=fix_separation, line_ratio=1.

    outstr[irp].ew_nofit = single_newy[icaii]*(wave[icaii]-wave[icaii-1])/factor_norm1
    outstr[irp].err_ew_nofit = 1./sqrt(single_newyivar[icaii])*(wave[icaii]-wave[icaii-1])/factor_norm1
    outstr[irp].ew_nofit_2 = double_newy[icaii]*(wave[icaii]-wave[icaii-1])/factor_norm2*2./3.
    outstr[irp].err_ew_nofit_2 = 1./sqrt(double_newyivar[icaii])*(wave[icaii]-wave[icaii-1])/factor_norm2*2./3.
    outstr[irp].ew_nofit_lr_1 = double_newy_1[icaii]*(wave[icaii]-wave[icaii-1])/factor_norm2_1*2./4.
    outstr[irp].err_ew_nofit_lr_1 = 1./sqrt(double_newyivar_1[icaii])*(wave[icaii]-wave[icaii-1])/factor_norm2_1*2./4.

    outstr[irp].ew_nofit_s4 = single_newy_s4[icaii]*(wave[icaii]-wave[icaii-1])/factor_norm1_s4
    outstr[irp].err_ew_nofit_s4 = 1./sqrt(single_newyivar_s4[icaii])*(wave[icaii]-wave[icaii-1])/factor_norm1_s4
    outstr[irp].ew_nofit_2_s4 = double_newy_s4[icaii]*(wave[icaii]-wave[icaii-1])/factor_norm2_s4*2./3.
    outstr[irp].err_ew_nofit_2_s4 = 1./sqrt(double_newyivar_s4[icaii])*(wave[icaii]-wave[icaii-1])/factor_norm2_s4*2./3.
    outstr[irp].ew_nofit_lr_1_s4 = double_newy_1_s4[icaii]*(wave[icaii]-wave[icaii-1])/factor_norm2_1_s4*2./4.
    outstr[irp].err_ew_nofit_lr_1_s4 = 1./sqrt(double_newyivar_1_s4[icaii])*(wave[icaii]-wave[icaii-1])/factor_norm2_1_s4*2./4.

    outstr[irp].ew_nofit_s5 = single_newy_s5[icaii]*(wave[icaii]-wave[icaii-1])/factor_norm1_s5
    outstr[irp].err_ew_nofit_s5 = 1./sqrt(single_newyivar_s5[icaii])*(wave[icaii]-wave[icaii-1])/factor_norm1_s5
    outstr[irp].ew_nofit_2_s5 = double_newy_s5[icaii]*(wave[icaii]-wave[icaii-1])/factor_norm2_s5*2./3.
    outstr[irp].err_ew_nofit_2_s5 = 1./sqrt(double_newyivar_s5[icaii])*(wave[icaii]-wave[icaii-1])/factor_norm2_s5*2./3.
    outstr[irp].ew_nofit_lr_1_s5 = double_newy_1_s5[icaii]*(wave[icaii]-wave[icaii-1])/factor_norm2_1_s5*2./4.
    outstr[irp].err_ew_nofit_lr_1_s5 = 1./sqrt(double_newyivar_1_s5[icaii])*(wave[icaii]-wave[icaii-1])/factor_norm2_1_s5*2./4.

;   stop

    for imc=0L, nmc-1L do begin
        mc_center = mc_line_wave[imc]
        if (mc_center lt xra[0]+20 or mc_center gt xra[1]-20) then continue
        tmp = min(abs(comp[irp].wave-mc_line_wave[imc]), icaii)
        jhusdss_lowz_doublet_fit3, wave[iuse], 1.-newy[iuse], fltarr(nuse)+1./terror^2, $
           in_quadra, in_slope, in_intercept, mc_center, in_separation, in_lflux, in_ratio, in_sigma1, in_sigma2, $
           quadra=quadra, slope=slope, intercept=intercept, center=center, separation=separation, $
           lflux=lflux, ratio=ratio, sigma1=sigma1, sigma2=sigma2, $
           err_quadra=err_quadra, err_slope=err_slope, err_intercept=err_intercept, err_center=err_center, $
           err_separation=err_separation, err_lflux=err_lflux, err_ratio=err_ratio, err_sigma1=err_sigma1, err_sigma2=err_sigma2, $
           maxwidth=maxwidth, inratio_range=[-10, 10]

        outstr[irp].ew_mc[imc] = lflux*(1.+ratio)
        outstr[irp].line_ratio_mc[imc] = 1./ratio
        outstr[irp].sigma_mc[imc] = (sigma1+sigma2)/2.

        outstr[irp].ew_nofit_mc[imc] = single_newy[icaii]*(wave[icaii]-wave[icaii-1])/factor_norm1
        outstr[irp].ew_nofit_2_mc[imc] = double_newy[icaii]*(wave[icaii]-wave[icaii-1])/factor_norm2*2./3.
        outstr[irp].ew_nofit_lr_1_mc[imc] = double_newy_1[icaii]*(wave[icaii]-wave[icaii-1])/factor_norm2_1*2./4.

        outstr[irp].ew_nofit_s4_mc[imc] = single_newy_s4[icaii]*(wave[icaii]-wave[icaii-1])/factor_norm1_s4
        outstr[irp].ew_nofit_2_s4_mc[imc] = double_newy_s4[icaii]*(wave[icaii]-wave[icaii-1])/factor_norm2_s4*2./3.
        outstr[irp].ew_nofit_lr_1_s4_mc[imc] = double_newy_1_s4[icaii]*(wave[icaii]-wave[icaii-1])/factor_norm2_1_s4*2./4.

        outstr[irp].ew_nofit_s5_mc[imc] = single_newy_s5[icaii]*(wave[icaii]-wave[icaii-1])/factor_norm1_s5
        outstr[irp].ew_nofit_2_s5_mc[imc] = double_newy_s5[icaii]*(wave[icaii]-wave[icaii-1])/factor_norm2_s5*2./3.
        outstr[irp].ew_nofit_lr_1_s5_mc[imc] = double_newy_1_s5[icaii]*(wave[icaii]-wave[icaii-1])/factor_norm2_1_s5*2./4.
    endfor

    ;; QAplots
    load_dp, /b
    print, 'vdisp: ', outstr[irp].sigma1, outstr[irp].err_sigma1, outstr[irp].sigma2, outstr[irp].err_sigma2
    print, 'line_ratio: ', outstr[irp].line_ratio, outstr[irp].err_line_ratio
    djs_plot, wave, 1.-single_newy, psym=10, $
        xra=xra, xst=1, yra=[1.-4.*terror, 1.+4.*terror], position=[0.1, 0.5, 0.9, 0.9], xtickformat='(A1)'
    djs_oplot, wave, 1.-single_yfit, color='green'
    djs_oplot, line_wave[0]*[1,1.], !y.crange[0]+[0.5,0.9]*(!y.crange[1]-!y.crange[0]), color='red'
    djs_oplot, line_wave[1]*[1,1.], !y.crange[0]+[0.5,0.9]*(!y.crange[1]-!y.crange[0]), color='red'
    ii = where(wave gt (line_wave[0]-5.) and wave lt (line_wave[0]+5.))
    djs_oplot, wave[ii], 1.-single_newy[ii], color='red'
    ii = where(wave gt (line_wave[1]-5.) and wave lt (line_wave[1]+5.))
    djs_oplot, wave[ii], 1.-single_newy[ii], color='red'
    djs_oplot, !x.crange, [1,1], linestyle=1

    iusemc = where(outstr[irp].ew_mc gt -99., nusemc)
    mcmom = moment(outstr[irp].ew_mc[iusemc], sdev=mcsdev)
    mcmom_sigma = moment(outstr[irp].sigma_mc[iusemc], sdev=sigma_mcsdev)
    mcmom_ew_nofit = moment(outstr[irp].ew_nofit_mc[iusemc], sdev=nofit_mcsdev)
    mcmom_ew_nofit_2 = moment(outstr[irp].ew_nofit_2_mc[iusemc], sdev=nofit_2_mcsdev)
    mcmom_ew_nofit_lr_1 = moment(outstr[irp].ew_nofit_lr_1_mc[iusemc], sdev=nofit_lr_1_mcsdev)
    mcmom_ew_nofit_s4 = moment(outstr[irp].ew_nofit_s4_mc[iusemc], sdev=nofit_s4_mcsdev)
    mcmom_ew_nofit_2_s4 = moment(outstr[irp].ew_nofit_2_s4_mc[iusemc], sdev=nofit_2_s4_mcsdev)
    mcmom_ew_nofit_lr_1_s4 = moment(outstr[irp].ew_nofit_lr_1_s4_mc[iusemc], sdev=nofit_lr_1_s4_mcsdev)
    mcmom_ew_nofit_s5 = moment(outstr[irp].ew_nofit_s5_mc[iusemc], sdev=nofit_s5_mcsdev)
    mcmom_ew_nofit_2_s5 = moment(outstr[irp].ew_nofit_2_s5_mc[iusemc], sdev=nofit_2_s5_mcsdev)
    mcmom_ew_nofit_lr_1_s5 = moment(outstr[irp].ew_nofit_lr_1_s5_mc[iusemc], sdev=nofit_lr_1_s5_mcsdev)

    djs_plot, wave, 1.-double_newy, psym=10, $
        xra=xra, xst=1, yra=[1.-4.*terror, 1.+4.*terror], position=[0.1, 0.1, 0.9, 0.5], /noerase
    djs_oplot, line_wave[0]*[1,1.], !y.crange[0]+[0.5,0.9]*(!y.crange[1]-!y.crange[0]), color='green'
    ii = where(wave gt (line_wave[0]-5.) and wave lt (line_wave[0]+5.))
    djs_oplot, wave[ii], 1.-double_newy[ii], color='green'
    ii = where(wave gt (line_wave[1]-5.) and wave lt (line_wave[1]+5.))
    djs_oplot, wave[ii], 1.-double_newy[ii], color='green'
    djs_oplot, !x.crange, [1,1], linestyle=1

;   plothist, outstr[irp].ew_mc[iusemc], bin=terror, xra=[-6*mcsdev, +6*mcsdev], $
;      position=[0.3, 0.05, 0.7, 0.45], /noerase
;   itwosigma = where(outstr[irp].ew_mc gt -2.*mcsdev and outstr[irp].ew_mc lt 2.*mcsdev, ntwosigma)
;   plothist, outstr[irp].ew_mc[itwosigma], bin=terror, xra=[-6*mcsdev, +6*mcsdev], $
;      position=[0.3, 0.05, 0.7, 0.45], /fill, /noerase
;   djs_oplot, outstr[irp].ew*[1,1], !y.crange, linestyle=0, thick=15., color='red'
    print, 'rp, ew, sdev, ew_nofit, ew_nofit_2', outstr[irp].rp, outstr[irp].ew, mcsdev, outstr[irp].ew_nofit, outstr[irp].ew_nofit_2
    print, 'vdisp: ', outstr[irp].sigma1, outstr[irp].sigma2
    outstr[irp].sdev_mc = mcsdev
    outstr[irp].sdev_sigma_mc = sigma_mcsdev
    outstr[irp].sdev_nofit_mc = nofit_mcsdev
    outstr[irp].sdev_nofit_2_mc = nofit_2_mcsdev
    outstr[irp].sdev_nofit_lr_1_mc = nofit_lr_1_mcsdev

    outstr[irp].sdev_nofit_s4_mc = nofit_s4_mcsdev
    outstr[irp].sdev_nofit_2_s4_mc = nofit_2_s4_mcsdev
    outstr[irp].sdev_nofit_lr_1_s4_mc = nofit_lr_1_s4_mcsdev

    outstr[irp].sdev_nofit_s5_mc = nofit_s5_mcsdev
    outstr[irp].sdev_nofit_2_s5_mc = nofit_2_s5_mcsdev
    outstr[irp].sdev_nofit_lr_1_s5_mc = nofit_lr_1_s5_mcsdev
endfor

mwrfits, outstr, outfile, /create
xx = findgen(1000)*0.01
xx = findgen(1000)*0.01
profile_slope = -1.38
profile_intercept = alog10(0.1)
yy = profile_slope*(xx-2.)+profile_intercept

;djs_plot, outstr.rp, outstr.ew_nofit_red, psym=5, color='red', xra=[5,300], yra=[1E-3, 1E0], /xlog, /ylog, xst=1, yst=1
;djs_oplot, 10.^xx, 10.^yy/2.
;oploterror, outstr.rp, outstr.ew_nofit_red, outstr.sdev_nofit_mc, psym=5, color=djs_icolor('red'), errcolor=djs_icolor('red')
;djs_oplot, outstr.rp, outstr.sdev_nofit_mc, color='red'
;djs_oplot, outstr.rp, -outstr.sdev_nofit_mc, color='red'

;if (~Coarse) then i_indep = [0, lindgen(15)*2+1] else i_indep = [0, lindgen(5)*2+1]
if (~Coarse) then i_indep = [lindgen(16)*2+1] else i_indep = [lindgen(5)*2+1]
load_dp, /b
djs_plot, outstr.rp, outstr.ew_nofit_2, psym=4, xra=[5,10000], yra=[1E-4, 1E0], /xlog, /ylog, xst=1, yst=1
;djs_plot, outstr.rp, outstr.ew_nofit_2, psym=4, xra=[200,2000], yra=[-5E-3, 5E-3], /xlog, xst=1, yst=1
djs_oplot, 10.^xx, 10.^yy
djs_oplot, outstr.rp, outstr.sdev_nofit_2_mc, color='red'
;djs_oplot, outstr.rp, -outstr.sdev_nofit_2_mc, color='red'
djs_oplot, rp_bor, ew_bor, psym=4, color='green'
djs_plot, outstr[i_indep].rp, outstr[i_indep].ew_nofit_2, psym=4, xra=[5,50000], yra=[1E-4, 1E0], /xlog, /ylog, xst=1, yst=1
djs_oplot, rp_bor, ew_bor, psym=4, color='green'
djs_oplot, 10.^xx, 10.^yy

end
