;pro jhusdss_lowz_caii_new

;Coarse = 0b
DoWeight = 0b
;DoNaI = 0b
DoScale = 0b
Scale_factor = 0.2 ;; M* \propto M(halo)^1 \propto r(halo)^3 and r(halo) \propto M*^(1/3); Scale to 10^10.1 Msun
;read,'Coarse? [1=yes, 0=no]: ', Coarse
read,'DoWeight? [1=yes, 0=no]: ', DoWeight
;read,'DoNaI? [1=yes, 0=no]: ', DoNaI
read,'DoScale? [1=yes, 0=no]: ', DoScale

nmfver = 106
overwrite = 1b
;read,'BOSS? [1=yes, 0=no]: ',BOSS
quiescent = 0b
read,'quiescent 0/1?: ', quiescent

line_wave = [3934.79, 3969.59]
xra = [3800., 4300.]

fix_separation = long((line_wave[1]-line_wave[0])/(line_wave[1]+line_wave[0])/alog(10.)*2E+4)

stackpath = jhusdss_get_path(/nmfqso)+'/'+string(nmfver, format='(I3.3)')+'/Stack/'
infile = stackpath + jhusdss_garching_galqso_stack_filename(nmfver, boss=boss)
if (quiescent) then infile = repstr(infile, '.fits', '_quiescent.fits') else infile = repstr(infile, '.fits', '_sf.fits')

;if (Coarse) then begin
;   infile=repstr(infile,'.fits','_coarse.fits')
;endif
if (DoWeight) then begin
   infile=repstr(infile,'.fits','_weight.fits')
endif
;if (DoNaI) then begin
;   infile=repstr(infile,'.fits','_NaI.fits')
;endif
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

nmc = 250
mc_line_wave = findgen(nmc)*1.5+line_wave[0]-60.
print, mc_line_wave[0], mc_line_wave[nmc-1]
;stop

outstr = replicate({rp:0., rp_min:0., rp_max:0., ew:0., err_ew:0., sigma:0., $
            wave_mc:mc_line_wave, ew_mc:fltarr(nmc)-100., sigma_mc:fltarr(nmc)-100., sdev_mc:-100., sdev_sigma_mc:-100., $
            ew_nofit:0., err_ew_nofit:0., ew_nofit_mc:fltarr(nmc)-100., sdev_nofit_mc:-100., $
            ew_nofit_2:0., err_ew_nofit_2:0., ew_nofit_2_mc:fltarr(nmc)-100., sdev_nofit_2_mc:-100., $
            ew_nofit_red:0., err_ew_nofit_red:0. $
;           ew_nofit_red_2:0., err_ew_nofit_red_2:0., $
            }, nrp)

for irp=0L, nrp-1L do begin
    rp = comp[irp].rp
    wave = comp[irp].wave
    tmp = min(abs(comp[irp].wave-line_wave[0]), icaii)
    tmp = min(abs(comp[irp].wave-line_wave[1]), ired)
    y = comp[irp].fluxgeomean
    med_continuum = median(y, 21, /even)
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
    in_sigma = 2.0

    in_intercept = median(1.-y[iuse_nocaii])
    terror = sqrt((moment(1.-y[iuse_nocaii]))[1])

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
    outstr[irp].rp = rp
    outstr[irp].ew = lflux
    outstr[irp].err_ew = err_lflux
    outstr[irp].sigma = sigma

    newy = y+clevel

    jhusdss_singlet_smooth, 1.-newy, fltarr(n_elements(y))+1./terror^2, $
       outflux=single_newy, outivar=single_newyivar, sigma=2., /normalize, factor_norm=factor_norm1
    jhusdss_doublet_smooth, 1.-newy, fltarr(n_elements(y))+1./terror^2, $
       outflux=double_newy, outivar=double_newyivar, sigma=2., /normalize, factor_norm=factor_norm2, separation=fix_separation

    outstr[irp].ew_nofit = single_newy[icaii]*(wave[icaii]-wave[icaii-1])/factor_norm1
    outstr[irp].err_ew_nofit = 1./sqrt(single_newyivar[icaii])*(wave[icaii]-wave[icaii-1])/factor_norm1
    outstr[irp].ew_nofit_2 = double_newy[icaii]*(wave[icaii]-wave[icaii-1])/factor_norm2*2./3.
    outstr[irp].err_ew_nofit_2 = 1./sqrt(double_newyivar[icaii])*(wave[icaii]-wave[icaii-1])/factor_norm2*2./3.

    outstr[irp].ew_nofit_red = single_newy[ired]*(wave[ired]-wave[ired-1])/factor_norm1
    outstr[irp].err_ew_nofit_red = 1./sqrt(single_newyivar[ired])*(wave[ired]-wave[ired-1])/factor_norm1
;   outstr[irp].ew_nofit_red_2 = double_newy[icaii]*(wave[icaii]-wave[icaii-1])/factor_norm2*6./5.
;   outstr[irp].err_ew_nofit_red_2 = 1./sqrt(double_newyivar[icaii])*(wave[icaii]-wave[icaii-1])/factor_norm2*6./5.
;   stop

    for imc=0L, nmc-1L do begin
        mc_center = mc_line_wave[imc]
        if (mc_center lt xra[0]+20 or mc_center gt xra[1]-20) then continue
        tmp = min(abs(comp[irp].wave-mc_line_wave[imc]), icaii)
        jhusdss_lowz_doublet_fit2, wave[iuse], 1.-newy[iuse], fltarr(nuse)+1./terror^2, $
           in_quadra, in_slope, in_intercept, mc_center, in_separation, in_lflux, in_ratio, in_sigma, $
           quadra=quadra, slope=slope, intercept=intercept, center=center, separation=separation, $
           lflux=lflux, ratio=ratio, sigma=sigma, $
           err_quadra=err_quadra, err_slope=err_slope, err_intercept=err_intercept, err_center=err_center, $
           err_separation=err_separation, err_lflux=err_lflux, err_ratio=err_ratio, err_sigma=err_sigma, $
           maxwidth=maxwidth

        outstr[irp].ew_mc[imc] = lflux
        outstr[irp].sigma_mc[imc] = sigma

        outstr[irp].ew_nofit_mc[imc] = single_newy[icaii]*(wave[icaii]-wave[icaii-1])/factor_norm1
        outstr[irp].ew_nofit_2_mc[imc] = double_newy[icaii]*(wave[icaii]-wave[icaii-1])/factor_norm2*2./3.

    endfor

    ;; QAplots
    load_dp, /b
    djs_plot, wave, 1.-single_newy, psym=10, $
        xra=xra, xst=1, yra=[1.-4.*terror, 1.+4.*terror], position=[0.1, 0.5, 0.9, 0.9], xtickformat='(A1)'
    djs_oplot, line_wave[0]*[1,1.], !y.crange[0]+[0.5,0.9]*(!y.crange[1]-!y.crange[0]), color='red'
    djs_oplot, line_wave[1]*[1,1.], !y.crange[0]+[0.5,0.9]*(!y.crange[1]-!y.crange[0]), color='red'
    ii = where(wave gt (line_wave[0]-5.) and wave lt (line_wave[0]+5.))
    djs_oplot, wave[ii], 1.-single_newy[ii], color='red'
    ii = where(wave gt (line_wave[1]-5.) and wave lt (line_wave[1]+5.))
    djs_oplot, wave[ii], 1.-single_newy[ii], color='red'

    iusemc = where(outstr[irp].ew_mc gt -99., nusemc)
    mcmom = moment(outstr[irp].ew_mc[iusemc], sdev=mcsdev)
    mcmom_sigma = moment(outstr[irp].sigma_mc[iusemc], sdev=sigma_mcsdev)
    mcmom_ew_nofit = moment(outstr[irp].ew_nofit_mc[iusemc], sdev=nofit_mcsdev)
    mcmom_ew_nofit_2 = moment(outstr[irp].ew_nofit_2_mc[iusemc], sdev=nofit_2_mcsdev)

    djs_plot, wave, 1.-double_newy, psym=10, $
        xra=xra, xst=1, yra=[1.-4.*terror, 1.+4.*terror], position=[0.1, 0.1, 0.9, 0.5], /noerase
    djs_oplot, line_wave[0]*[1,1.], !y.crange[0]+[0.5,0.9]*(!y.crange[1]-!y.crange[0]), color='green'
    ii = where(wave gt (line_wave[0]-5.) and wave lt (line_wave[0]+5.))
    djs_oplot, wave[ii], 1.-double_newy[ii], color='green'
    ii = where(wave gt (line_wave[1]-5.) and wave lt (line_wave[1]+5.))
    djs_oplot, wave[ii], 1.-double_newy[ii], color='green'

;   plothist, outstr[irp].ew_mc[iusemc], bin=terror, xra=[-6*mcsdev, +6*mcsdev], $
;      position=[0.3, 0.05, 0.7, 0.45], /noerase
;   itwosigma = where(outstr[irp].ew_mc gt -2.*mcsdev and outstr[irp].ew_mc lt 2.*mcsdev, ntwosigma)
;   plothist, outstr[irp].ew_mc[itwosigma], bin=terror, xra=[-6*mcsdev, +6*mcsdev], $
;      position=[0.3, 0.05, 0.7, 0.45], /fill, /noerase
;   djs_oplot, outstr[irp].ew*[1,1], !y.crange, linestyle=0, thick=15., color='red'
    print, outstr[irp].rp, outstr[irp].ew, mcsdev, outstr[irp].sigma, outstr[irp].ew_nofit, outstr[irp].ew_nofit_2

    outstr[irp].sdev_mc = mcsdev
    outstr[irp].sdev_sigma_mc = sigma_mcsdev
    outstr[irp].sdev_nofit_mc = nofit_mcsdev
    outstr[irp].sdev_nofit_2_mc = nofit_2_mcsdev
endfor

mwrfits, outstr, outfile, /create

xx = findgen(1000)*0.01
profile_slope = -1.38
profile_intercept = alog10(0.01)
yy = profile_slope*(xx-2.)+profile_intercept
djs_plot, outstr.rp, outstr.ew_nofit_2, psym=4, xra=[5,2000], yra=[1E-4, 1E0], /xlog, /ylog, xst=1, yst=1
oploterror, outstr.rp, outstr.ew_nofit_2, outstr.sdev_nofit_2_mc, psym=5, color=djs_icolor('blue'), errcolor=djs_icolor('blue')
djs_oplot, 10.^xx, 10.^yy
djs_oplot, outstr.rp, outstr.sdev_nofit_2_mc, color='blue'
;djs_plot, outstr.rp, outstr.ew_nofit_2, psym=4, xra=[200,2000], yra=[-5E-3, 5E-3], /xlog, xst=1, yst=1
;oploterror, a_sf.rp, a.ew_nofit_2, a.sdev_nofit_2_mc, psym=5, color=djs_icolor('red'), errcolor=djs_icolor('red')
;djs_oplot, a.rp, a.sdev_nofit_2_mc, color='red'

;djs_plot, outstr.rp, outstr.ew_nofit, psym=4, xra=[5,2000], yra=[1E-4, 1E0], /xlog, /ylog, xst=1, yst=1
;;djs_plot, outstr.rp, outstr.ew_nofit, psym=4, xra=[200,2000], yra=[-5E-3, 5E-3], /xlog, xst=1, yst=1
;djs_oplot, 10.^xx, 10.^yy
;djs_oplot, outstr.rp, outstr.sdev_nofit_mc, color='red'
;djs_oplot, outstr.rp, -outstr.sdev_nofit_mc, color='red'

;djs_plot, outstr.rp, outstr.ew_nofit_red, psym=5, color='red', xra=[5,300], yra=[1E-3, 1E0], /xlog, /ylog, xst=1, yst=1
;djs_oplot, 10.^xx, 10.^yy/2.
;oploterror, outstr.rp, outstr.ew_nofit_red, outstr.sdev_nofit_mc, psym=5, color=djs_icolor('red'), errcolor=djs_icolor('red')
;djs_oplot, outstr.rp, outstr.sdev_nofit_mc, color='red'
;djs_oplot, outstr.rp, -outstr.sdev_nofit_mc, color='red'

end
