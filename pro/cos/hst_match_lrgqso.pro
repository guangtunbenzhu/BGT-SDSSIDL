
pro hst_match_lrgqso
gal = mrdfits('/home/gz323/DATA/BOSS/granada_fsps_salp_wideform_dust-v5_6_5.fits', 1)
qso = mrdfits('/home/menard/DATA/SDSS/QSO/dr7_bh_May09_2011_trimmed.fits', 1)

;; quality control later
ii = where(gal.z gt 0.4 and gal.z lt 0.75, nn)

omegaM = 0.3
omegaL = 0.7
h100 = 0.7
ch0 = 2997.92458 ;;Mpc*h^-1
gal_angdist = angdidis(gal.z, omegaM, omegaL)*ch0/h100


spherematch, gal[ii].ra, gal[ii].dec, qso.ra, qso.dec, 120./3600., m1, m2, distance12, maxmatch=0

help, ii, m1

index_gal = ii[m1]
index_qso = m2
rp_deg = distance12
rp_mpc = distance12*!dpi/180.*gal_angdist[ii[m1]]

str_tmp = {index_lrg:0L, index_qso:0L, rp_deg:0., rp_mpc:0.}
out_match = replicate(str_tmp, n_elements(index_gal))
out_match.index_lrg = index_gal
out_match.index_qso = index_qso
out_match.rp_deg = rp_deg
out_match.rp_mpc = rp_mpc


mwrfits, out_match, '/home/gz323/DATA/BOSS/granada_fsps_salp_wideform_dust-v5_6_5_dr7_bh_match.fits', /create

end
