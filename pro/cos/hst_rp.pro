h100 = 0.7

cosdir = '/data1/gwln2scratch/menard/gz323/COS'
psfile = cosdir+'/QAplots/rp.ps'

choice_load_data = 0
read,'load data? [1=yes, 0=no]: ',choice_load_data

if choice_load_data eq 1 then begin
;; SDSS DR7
   gal = mrdfits('/home/gz323/DATA/BOSS/wisconsin_pca_bc03-v5_6_5.fits.gz', 1)
   gal_kcorr = mrdfits('/home/gz323/DATA/BOSS/kcorrect_granada_fsps_salp_wideform_dust-v5_6_5.fits', 1)
   gal_uniq = mrdfits('/home/gz323/DATA/BOSS/uniq_wisconsin_pca_bc03-v5_6_5.fits.gz', 1)
   qso = mrdfits('/home/menard/DATA/SDSS/QSO/dr7_bh_May09_2011_trimmed.fits', 1)

   zgal= gal[match.index_lrg].z
   zuniq = gal_uniq[match.index_lrg].choose
   umr = gal_kcorr[match.index_lrg].absmag[0]-gal_kcorr[match.index_lrg].absmag[2]
   Mr = gal_kcorr[match.index_lrg].absmag[2]+5.*alog10(h100)
   Mrivar = gal_kcorr[match.index_lrg].amivar[2]
   sn_median = gal[match.index_lrg].sn_median
   mass = gal[match.index_lrg].mstellar_median

   galex_nuv = qso[match.index_qso].galex_mag[0]
   galex_fuv = qso[match.index_qso].galex_mag[1]
   galex_sep = qso[match.index_qso].galex_offset
   zqso = qso[match.index_qso].z_hw

   hstselect = mrdfits('/home/gz323/DATA/BOSS/hstselect_final.fits',1)
   nabs = 10

   sdss = mrdfits('/home/gz323/DATA/SDSS/Garching//kcorrect.none.petro.z0.00.fits', 1)

   cosfile = '/data1/gwln2scratch/menard/gz323/COS/COS_Halo_Werk13_tbl1_simplify.txt'
   readcol, cosfile, cos_objname, cos_z, cos_rp, cos_mstellar, cos_lum, cos_flag_umr, cos_umr, cos_flag_ssfr, cos_ssfr, cos_abund, format='A, F, F, F, F, A, F, A, F, F', delimiter='|', /preserve_null

   tcosfile = '/data1/gwln2scratch/menard/gz323/COS/COS_Halo_Tumlinson13_tbl2_simplify.txt'
   readcol, tcosfile, tcos_objname, tcos_type, tcos_z, tcos_mr, tcos_mstellar, tcos_flag_ssfr, tcos_ssfr, tcos_rp, tcos_rvir, format='A, A, F, F, F, A, F, F, F', delimiter='|', /preserve_null
;  hstcos = mrdfits(cosfile, 1)
endif

xthick=8
ythick=8
cthick=4
thick=10
charthick=4
charsize=1.5
psym=6
symsize=0.2


k_print, filename=psfile, axis_char_scale=1.4, xsize=8., ysize=12.

    !p.multi = [0,1,2]
    !y.margin = 0
;   xtitle = textoidl('M_r')
;   ytitle = textoidl('M_u - M_r')
    xra = [-410, 410.]
    yra = [-410, 410.]

    djs_plot, xra, yra, /nodata, xst=4, yst=4

    ra1 = hstselect.ra_gal/360.D0*24.D0
    dec1 = hstselect.dec_gal
    ra2 = hstselect.ra_qso/360.D0*24.D0
    dec2 = hstselect.dec_qso
    posang, 1, ra1, dec1, ra2, dec2, pang
    y = hstselect.rp_mpc*1E3*cos(pang*!dpi/180.)
    x = -hstselect.rp_mpc*1E3*sin(pang*!dpi/180.)
    print, pang

    plotsym, 0, 1, /fill
    djs_oplot, x[0:nabs-1], y[0:nabs-1], psym=8, symsize=2.0, color='dark green'
    plotsym, 0, 1, /fill
    djs_oplot, x[nabs:*], y[nabs:*], psym=8, symsize=2.0, color='red'

    rad_max = [50., 100, 200, 300, 400.]
    for irad=0L, n_elements(rad_max)-1L do begin
        points = jhusdss_create_circle(0, 0, rad_max[irad])
        djs_oplot, points[0,*], points[1,*], thick=thick, linestyle=2
;       if (irad ge 2) then djs_xyouts, rad_max[irad], 0, string(rad_max[irad], format='(I3)'), charsize=charsize, charthick=charthick
        if (irad ge 1 and irad lt n_elements(rad_max)-1L) then djs_xyouts, rad_max[irad]*cos(-!dpi*0.15)*1.1, rad_max[irad]*sin(-!dpi*0.15)*1.1, string(rad_max[irad], format='(I3)'), charsize=charsize, charthick=charthick
        if (irad eq n_elements(rad_max)-1L) then djs_xyouts, rad_max[irad]*cos(-!dpi*0.15)*1.1, rad_max[irad]*sin(-!dpi*0.15)*1.1, string(rad_max[irad], format='(I3)')+' kpc', charsize=charsize, charthick=charthick
    endfor

;   xtitle = textoidl('log_{10} r_p/kpc')
    xtitle = textoidl('r_p (kpc)')
    ytitle = textoidl('\Delta N')
    xra = [1.3, 2.74]
    yra = [0.1, 5.3]
    plothist, alog10(hstselect[0:nabs-1].rp_mpc*1E3), bin=alog10(2.)/2., thick=thick, charsize=charsize, charthick=charthick, xthick=xthick, ythick=ythick, $
        xtitle=xtitle, ytitle=ytitle, fcolor=djs_icolor('dark green'), xst=1, yst=1, xra=xra, yra=yra, /fill, color=djs_icolor('dark green'), position=[0.20,0.20,0.85,0.52], xtickformat='(A1)'
    plothist, alog10(hstselect[nabs:*].rp_mpc*1E3), bin=alog10(2.)/2., thick=thick, $
        fcolor=djs_icolor('red'), /overplot, /fill, /fline, forientation=45., color=djs_icolor('red')
;   plothist, alog10(hstselect.rp_mpc*1E3), bin=0.2, thick=thick, $
;       fcolor=djs_icolor('red'), /overplot, /fill, /fline, forientation=45., color=djs_icolor('red')
;   plothist, alog10(hstselect[0:nabs-1].rp_mpc*1E3), bin=0.2, thick=thick, $
;       fcolor=djs_icolor('dark green'), /overplot, /fill, color=djs_icolor('dark green')
;   plothist, alog10(hstselect[nabs:*].rp_mpc*1E3), bin=0.2, thick=thick, charsize=charsize, charthick=charthick, $
;       xtitle=xtitle, ytitle=ytitle, fcolor=djs_icolor('red'), xst=8, yst=1, /overplot, /fill, /fline, forientation=45., color=djs_icolor('red')

    djs_xyouts, alog10(25.*0.80), -0.25, '25', charthick=charthick*1.2, charsize=charsize*1.2
    djs_xyouts, alog10(50.*0.80), -0.25, '50', charthick=charthick*1.2, charsize=charsize*1.2
    djs_xyouts, alog10(100.*0.8), -0.25, '100', charthick=charthick*1.2, charsize=charsize*1.2
    djs_xyouts, alog10(200.*0.8), -0.25, '200', charthick=charthick*1.2, charsize=charsize*1.2
    djs_xyouts, alog10(400.*0.8), -0.25, '400', charthick=charthick*1.2, charsize=charsize*1.2

    djs_xyouts, 1.35, 5.41-0.62, 'green: LRG-QSO targets with Mg II', charthick=charthick, charsize=charsize, color='dark green'
    djs_xyouts, 1.35, 5.07-0.6, 'red: LRG-QSO targets without Mg II', charthick=charthick, charsize=charsize, color='red'

k_end_print

end
