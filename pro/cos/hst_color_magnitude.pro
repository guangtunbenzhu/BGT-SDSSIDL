h100 = 0.7

cosdir = '/data1/gwln2scratch/menard/gz323/COS'
psfile = cosdir+'/QAplots/color_magnitude.ps'

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


k_print, filename=psfile, axis_char_scale=1.4, xsize=10., ysize=10.

    xtitle = textoidl('M_r')
    ytitle = textoidl('M_u - M_r')
    xra = [-16.51, -24.49]
    yra = [-0.40, 3.29]

    isdss = where(sdss.z gt 0.01 and sdss.z lt 0.10, nsdss)
    x = sdss[isdss].absmag[2]+5.*alog10(h100)
    y = sdss[isdss].absmag[0]-sdss[isdss].absmag[2]
    flevels = [0.2, 0.40, 0.78, 0.92]
    hogg_scatterplot, x, y, xthick=xthick, ythick=ythick, cthick=cthick, $
         xra=xra, yra=yra, levels=flevels, exponent=1.00, outliers=0, $
         xnpix=35L, ynpix=35L, charthick=charthick, charsize=charsize, $
         xtitle=xtitle, ytitle=ytitle, title=title, xstyle=1, ystyle=1, $
         outsymsize=0.3, position=[0.20, 0.20, 0.90, 0.90]

;   iboss = where(gal.z gt 0.4 and gal.z lt 0.75, nboss)
;   x = gal_kcorr[iboss].absmag[2]+5.*alog10(h100)
;   y = gal_kcorr[iboss].absmag[0]-gal_kcorr[iboss].absmag[2]
    iboss = where(Mr lt -21. and umr gt 1.2 and sn_median gt 3. and mass gt 7. and mass le 13., nboss)
    iran = fix(randomu(seed, 5000L)*nboss)
    x = Mr[iboss[iran]]
    y = umr[iboss[iran]]
    plotsym, 0, 1, /fill
    djs_oplot, x, y, psym=8, symsize=0.1, color='magenta'

    x = hstselect[nabs:*].absmag_gal[2]
    y = hstselect[nabs:*].absmag_gal[0] - hstselect[nabs:*].absmag_gal[2]
    plotsym, 0, 1, /fill
    djs_oplot, x, y, psym=8, symsize=1.7, color='red'


    x = hstselect[0:nabs-1].absmag_gal[2]
    y = hstselect[0:nabs-1].absmag_gal[0] - hstselect[0:nabs-1].absmag_gal[2]
;   plotsym, 0, 1, /fill
;   djs_oplot, x, y, psym=8, symsize=1.7, color='white'
;   plotsym, 0, 1, thick=5
;   djs_oplot, x, y, psym=8, symsize=1.7, color='red'
    plotsym, 0, 1, /fill
    djs_oplot, x, y, psym=8, symsize=1.7, color='dark green'

    x = tcos_Mr[0:43]
    y = cos_umr[0:43]-0.15
    plotsym, 8, 1, /fill
    djs_oplot, x, y, psym=8, symsize=1.4, color='light blue'

;   djs_xyouts, -16.8, 3.8+0.35, 'gray-scale: SDSS main sample at z\sim0.1', charsize=charsize, charthick=charthick, color='black'
;   djs_xyouts, -16.8, 3.65+0.35, 'cyan: COS-Halo galaxies at z\sim0.2', charsize=charsize, charthick=charthick, color='cyan'
;   djs_xyouts, -16.8, 3.5+0.25, 'magenta: BOSS LRG sample at z\sim0.5', charsize=charsize, charthick=charthick, color='magenta'
;   djs_xyouts, -16.8, 3.35+0.25, 'green: LRG targets with Mg II absorbers', charsize=charsize, charthick=charthick, color='dark green'
;   djs_xyouts, -16.8, 3.2+0.25, 'red: LRG targets without Mg II absorbers', charsize=charsize, charthick=charthick, color='red'

    djs_xyouts, -16.8, +0.10+0.35, 'gray-scale: SDSS main sample at z\sim0.1', charsize=charsize, charthick=charthick, color='black'
    djs_xyouts, -16.8, -0.05+0.35, 'blue: COS-Halo galaxies at z\sim0.2', charsize=charsize, charthick=charthick, color='blue'
    djs_xyouts, -16.8, -0.20+0.25, 'magenta: BOSS LRG sample at z\sim0.5', charsize=charsize, charthick=charthick, color='magenta'
    djs_xyouts, -16.8, -0.35+0.25, 'green: LRG-QSO targets with Mg II', charsize=charsize, charthick=charthick, color='dark green'
    djs_xyouts, -16.8, -0.50+0.25, 'red: LRG-QSO targets without Mg II', charsize=charsize, charthick=charthick, color='red'

k_end_print

end
