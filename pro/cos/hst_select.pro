h100 = 0.7

choice_load_data = 0
read,'load data? [1=yes, 0=no]: ',choice_load_data

if choice_load_data eq 1 then begin
   gal = mrdfits('/home/gz323/DATA/BOSS/wisconsin_pca_bc03-v5_6_5.fits.gz', 1)
   gal_kcorr = mrdfits('/home/gz323/DATA/BOSS/kcorrect_granada_fsps_salp_wideform_dust-v5_6_5.fits', 1)
   gal_uniq = mrdfits('/home/gz323/DATA/BOSS/uniq_wisconsin_pca_bc03-v5_6_5.fits.gz', 1)
   qso = mrdfits('/home/menard/DATA/SDSS/QSO/dr7_bh_May09_2011_trimmed.fits', 1)
   mgii_qso = mrdfits('/home/gz323/DATA/SDSS/QSO/NMF/107/Absorbers/ALLQSO_Trimmed_QSO_Absorbers_NMF_107_MgII_ALL.fits', 1)

   match = mrdfits('/home/gz323/DATA/BOSS/granada_fsps_salp_wideform_dust-v5_6_5_dr7_bh_match.fits', 1)

   zgal= gal[match.index_lrg].z
   zuniq = gal_uniq[match.index_lrg].choose
   umr = gal_kcorr[match.index_lrg].absmag[0]-gal_kcorr[match.index_lrg].absmag[2]
   Mr = gal_kcorr[match.index_lrg].absmag[2]+5.*alog10(h100)

   galex_nuv = qso[match.index_qso].galex_mag[0]
   galex_fuv = qso[match.index_qso].galex_mag[1]
   galex_sep = qso[match.index_qso].galex_offset
   zqso = qso[match.index_qso].z_hw
endif

ii = where(match.rp_mpc lt 0.440 $
       and umr gt 2. $
       and umr lt 2.8 $
       and Mr  lt -22.4 $
       and Mr  gt -24.2 $
       and galex_fuv lt 20.41 and galex_fuv gt 10. $
       and galex_nuv lt 19.0 and galex_nuv gt 10. $
       and galex_sep lt 5. $
       and zqso lt 2.0 $
       and zqso - zgal gt 0.1 $
       and zuniq, nn)

help, ii

;; Select Mg II absorbers
ra_select_qso = qso[match[ii].index_qso].ra
dec_select_qso = qso[match[ii].index_qso].dec
z_select_gal = gal[match[ii].index_lrg].z

;; default no Mg II absorbers
select_2 = replicate({mgii:0}, nn)

spherematch, ra_select_qso, dec_select_qso, mgii_qso.ra, mgii_qso.dec, 1./3600., m1, m2, distance12, maxmatch=0
select_2[m1].mgii = 2

jj = where (mgii_qso[m2].nabs eq 1 and abs(mgii_qso[m2].zabs[0] - z_select_gal[m1]) lt 0.004, mm)
if mm gt 0 then select_2[m1[jj]].mgii = 1

strtmp = {ra_qso:0.D, dec_qso:0.D, z_qso:0., galex_mag_qso: fltarr(2), galex_mag_err_qso:fltarr(2), $
          have_mgii:0b, z_mgii:0., rew_mgii_2796:0., rew_mgii_2796_err:0., rew_mgii_2803:0., rew_mgii_2803_err:0., $
          plate_qso:0L, fiber_qso:0L, mjd_qso:0L, plate_gal:0L, fiber_gal:0L, mjd_gal:0L, $
          ra_gal:0.D, dec_gal:0.D, z_gal:0., absmag_gal: fltarr(5), amivar_gal:fltarr(5), mstellar_gal: 0., mstellar_err_gal:0., $
          rp_deg:0., rp_mpc:0.}

ifinal_select = where(select_2.mgii ne 2, nfinal)
final_select = match[ii[ifinal_select]]
help, final_select

outstr = replicate(strtmp, nfinal)
outstr.ra_qso = qso[final_select.index_qso].ra
outstr.dec_qso = qso[final_select.index_qso].dec
outstr.z_qso = qso[final_select.index_qso].z_hw
outstr.galex_mag_qso = qso[final_select.index_qso].galex_mag
outstr.galex_mag_err_qso = qso[final_select.index_qso].galex_mag_err
outstr.plate_qso = qso[final_select.index_qso].plate
outstr.fiber_qso = qso[final_select.index_qso].fiber
outstr.mjd_qso = qso[final_select.index_qso].mjd

;;let's do spherematch again
spherematch, qso[final_select.index_qso].ra, qso[final_select.index_qso].dec, mgii_qso.ra, mgii_qso.dec, 1./3600., m1, m2, distance12, maxmatch=1

outstr[m1].have_mgii = 1b
outstr[m1].z_mgii = mgii_qso[m2].zabs[0]
outstr[m1].rew_mgii_2796 = mgii_qso[m2].rew_mgii_2796[0]
outstr[m1].rew_mgii_2796_err = mgii_qso[m2].err_rew_mgii_2796[0]
outstr[m1].rew_mgii_2803 = mgii_qso[m2].rew_mgii_2803[0]
outstr[m1].rew_mgii_2803_err = mgii_qso[m2].err_rew_mgii_2803[0]

outstr[m1].plate_qso = mgii_qso[m2].plate
outstr[m1].fiber_qso = mgii_qso[m2].fiber
outstr[m1].mjd_qso = mgii_qso[m2].mjd

outstr.plate_gal = gal[final_select.index_lrg].plate
outstr.fiber_gal = gal[final_select.index_lrg].fiberid
outstr.mjd_gal = gal[final_select.index_lrg].mjd
outstr.ra_gal = gal[final_select.index_lrg].ra
outstr.dec_gal = gal[final_select.index_lrg].dec
outstr.z_gal = gal[final_select.index_lrg].z
outstr.absmag_gal = gal_kcorr[final_select.index_lrg].absmag + 5.*alog10(h100)
outstr.amivar_gal = 1./gal_kcorr[final_select.index_lrg].amivar
outstr.mstellar_gal = gal[final_select.index_lrg].mstellar_median
outstr.mstellar_err_gal = gal[final_select.index_lrg].mstellar_err

outstr.rp_deg = final_select.rp_deg
outstr.rp_mpc = final_select.rp_mpc

mwrfits, outstr, '/home/gz323/DATA/BOSS/before_final.fits', /create

;; with Mg II absorbers
iselect = where((outstr.have_mgii and outstr.galex_mag_qso[1] lt 20.5) or (outstr.have_mgii and outstr.rp_mpc*1E3 le 50.))
nabs = n_elements(iselect)
;iselect = where((outstr.have_mgii and outstr.galex_mag_qso[0] lt 19.))
help, iselect
;rad_min = [20., 50., , 113., 160., 226.];, 270.]
;rad_max = [50., 100., 160., 226., 400.];, 405.]
rad_min = [20., 60., 100., 150., 225.]
rad_max = [60., 100., 150., 225., 450.];, 405.]
for irad = 0L, n_elements(rad_min)-1L do begin
    kk = where(~outstr.have_mgii and outstr.rp_mpc*1E3 gt rad_min[irad] and outstr.rp_mpc*1E3 le rad_max[irad], ll)
    isort = bsort(outstr[kk].galex_mag_qso[0])
    iselect = [iselect, kk[isort[0:2]]]
endfor

hstselect = outstr[iselect]
for i=0L, n_elements(hstselect)-1L do print, hstselect[i].ra_qso, hstselect[i].dec_qso, hstselect[i].have_mgii, hstselect[i].rp_mpc*1E3, hstselect[i].z_qso, hstselect[i].galex_mag_qso[1], hstselect[i].galex_mag_qso[0], hstselect[i].z_gal, (hstselect[i].z_gal+1.)*1026.

print, "Manually selected"
;imanual = [0, 1, 2, 3, 4, 5, 6, 7, 8, 10, 11, 13, 14, 16, 17, 18, 19, 20, 21, 22, 23, 24]

imanual = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 11, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22]
nabs = 8
hstselect = hstselect[imanual]

mwrfits, hstselect, '/home/gz323/DATA/BOSS/hstselect_final.fits', /create

;; print qso redshift and galex FUV mag and observed wavelength of 1030 Ang
for i=0L, n_elements(hstselect)-1L do print, hstselect[i].ra_qso, hstselect[i].dec_qso, hstselect[i].have_mgii, hstselect[i].rp_mpc*1E3, hstselect[i].z_qso, hstselect[i].galex_mag_qso[1], hstselect[i].galex_mag_qso[0], hstselect[i].z_gal, (hstselect[i].z_gal+1.)*1026.

end
