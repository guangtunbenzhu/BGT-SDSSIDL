h100 = 0.7

cosdir = '/data1/gwln2scratch/menard/gz323/COS'
psfile = cosdir+'/QAplots/zlrg.ps'

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
   nabs=10

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

    xtitle = textoidl('LRG Redshift z')
;   ytitle = textoidl('\DeltaN per \Deltaz = 0.01')
    ytitle = textoidl('\DeltaN')
    xra = [0.36, 0.66]
    yra = [0, 5]

    plothist, hstselect[0:nabs-1].z_gal, bin=0.01, xra=xra, yra=yra, $
        xtitle=xtitle, ytitle=ytitle, xthick=xthick, ythick=ythick, thick=thick, $
        charthick=charthick, charsize=charsize, color=djs_icolor('dark green'), /fill, fcolor=djs_icolor('dark green')

    plothist, hstselect[nabs:*].z_gal, bin=0.01, thick=thick, color=djs_icolor('red'), /overplot, /fill, /fline, forientation=45, fcolor=djs_icolor('red')

k_end_print

psfile = cosdir+'/QAplots/zqso.ps'
k_print, filename=psfile, axis_char_scale=1.4, xsize=10., ysize=10.

    xtitle = textoidl('Quasar Redshift z')
    ytitle = textoidl('\DeltaN')
    xra = [0.4, 1.56]
    yra = [0, 5]

    plothist, hstselect[0:nabs-1].z_qso, bin=0.05, xra=xra, yra=yra, $
        xtitle=xtitle, ytitle=ytitle, xthick=xthick, ythick=ythick, thick=thick, $
        charthick=charthick, charsize=charsize, color=djs_icolor('dark green'), /fill, fcolor=djs_icolor('dark green')

    plothist, hstselect[nabs:*].z_qso, bin=0.05, thick=thick, color=djs_icolor('red'), /overplot, /fill, /fline, forientation=45, fcolor=djs_icolor('red')

k_end_print

psfile = cosdir+'/QAplots/FUV_NUV.ps'
k_print, filename=psfile, axis_char_scale=1.4, xsize=8., ysize=12.

    !p.multi=[0,1,2]
    !y.margin=0
;   xtitle = textoidl('Galex FUV')
;   ytitle = textoidl('\DeltaN')
    xra = [16.50, 20.50]
    yra = [0, 7.49]

    plothist, hstselect[0:nabs-1].galex_mag_qso[0], bin=0.5, xra=xra, yra=yra, $
        xtitle=' ', ytitle=ytitle, xthick=xthick, ythick=ythick, thick=thick, $
        charthick=charthick, charsize=charsize, color=djs_icolor('dark green'), /fill, fcolor=djs_icolor('dark green'), xtickformat='(A1)', position=[0.2, 0.55, 0.9, 0.9]

    plothist, hstselect[nabs:*].galex_mag_qso[0], bin=0.5, thick=thick, color=djs_icolor('red'), /overplot, /fill, /fline, forientation=45, fcolor=djs_icolor('red')

    djs_xyouts, 19.5, 0.45, 'Galex NUV', charthick=charthick, charsize=charsize

    djs_xyouts, 16.6, 7.00, 'green: LRG-QSO targets with Mg II', charthick=charthick, charsize=charsize, color='dark green'
    djs_xyouts, 16.6, 6.60, 'red: LRG-QSO targets without Mg II', charthick=charthick, charsize=charsize, color='red'

    plothist, hstselect[0:nabs-1].galex_mag_qso[1], bin=0.5, xra=xra, yra=yra, $
        xtitle='Galex FUV/NUV', ytitle=ytitle, xthick=xthick, ythick=ythick, thick=thick, $
        charthick=charthick, charsize=charsize, color=djs_icolor('dark green'), /fill, fcolor=djs_icolor('dark green'), position=[0.2, 0.2, 0.9, 0.55]

    plothist, hstselect[nabs:*].galex_mag_qso[1], bin=0.5, thick=thick, color=djs_icolor('red'), /overplot, /fill, /fline, forientation=45, fcolor=djs_icolor('red')
    djs_xyouts, 19.5, 6.8, 'Galex FUV', charthick=charthick, charsize=charsize


k_end_print


psfile = cosdir+'/QAplots/FUV.ps'
k_print, filename=psfile, axis_char_scale=1.4, xsize=10., ysize=10.

    xtitle = textoidl('Galex FUV')
    ytitle = textoidl('\DeltaN')
    xra = [17.4, 20.6]
    yra = [0, 3]

    plothist, hstselect[0:nabs-1].galex_mag_qso[1], bin=0.1, xra=xra, yra=yra, $
        xtitle=xtitle, ytitle=ytitle, xthick=xthick, ythick=ythick, thick=thick, $
        charthick=charthick, charsize=charsize, color=djs_icolor('dark green'), /fill, fcolor=djs_icolor('dark green')

    plothist, hstselect[nabs:*].galex_mag_qso[1], bin=0.1, thick=thick, color=djs_icolor('red'), /overplot, /fill, /fline, forientation=45, fcolor=djs_icolor('red')

k_end_print

psfile = cosdir+'/QAplots/NUV.ps'
k_print, filename=psfile, axis_char_scale=1.4, xsize=10., ysize=10.

    xtitle = textoidl('Galex NUV')
    ytitle = textoidl('\DeltaN')
    xra = [16.8, 19.5]
    yra = [0, 4]

    plothist, hstselect[0:nabs-1].galex_mag_qso[0], bin=0.1, xra=xra, yra=yra, $
        xtitle=xtitle, ytitle=ytitle, xthick=xthick, ythick=ythick, thick=thick, $
        charthick=charthick, charsize=charsize, color=djs_icolor('dark green'), /fill, fcolor=djs_icolor('dark green')

    plothist, hstselect[nabs:*].galex_mag_qso[0], bin=0.1, thick=thick, color=djs_icolor('red'), /overplot, /fill, /fline, forientation=45, fcolor=djs_icolor('red')

k_end_print

end
