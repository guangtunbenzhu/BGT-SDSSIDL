nmfver = 106
ivarweigth = 1b
overwrite = 1b
savespec = 0b
savemean = 1b
snrcut = 2.
sdevcut = 0.10
sigma_cut = 2.
whatthewhat = 1

   ;; foreground galaxies
   garching_path = jhusdss_get_path(/garching)
   garching_file = garching_path+'/'+'gal_info_dr7_v5_2.fit.gz'
   gal = mrdfits(garching_file, 1)
;  uniq_file = garching_path+'/'+'gal_uniq_dr7_v5_2.fits'
;  galuniq = mrdfits(uniq_file, 1)

;  sfr_file = garching_path+'/'+'gal_totsfr_dr7_v5_2.fits.gz'
;  sfr = mrdfits(sfr_file, 1)
;  mass_file = garching_path+'/'+'totlgm_dr7_v5_2.fit.gz'
;  mass = mrdfits(mass_file, 1)
;  ssfr_file = garching_path+'/'+'gal_totspecsfr_dr7_v5_2.fits.gz'
;  ssfr = mrdfits(ssfr_file, 1)
   phi_file = garching_path+'/'+'gal_phi_dr7_v5_2.fits'
   phi = mrdfits(phi_file, 1)

   ;; background quasars
;  print, 'Load SDSS DR7'
   qso = jhusdss_qso_readin(boss=boss)
;  stat = jhusdss_qsostats_readin(nmfver, boss=boss)
;  allspec = jhusdss_read_allqsospec(nmfver, /normresi, boss=boss)
;  allspec = jhusdss_read_allqsospec(nmfver, /subtresi, boss=boss)

   ;; match
   match0 = jhusdss_galqso_match_readin(boss=boss)

;; Calculate azimuthal angle
ra1 = gal[match0.index_gal].ra/360.D0*24.D0
dec1 = gal[match0.index_gal].dec
ra2 = qso[match0.index_qso].ra/360.D0*24.D0
dec2 = qso[match0.index_qso].dec

posang, 1, ra1, dec1, ra2, dec2, pang
r_phi = phi[match0.index_gal].phi_exp_deg[2]
r_ab = phi[match0.index_gal].ab_exp[2]
r_radius = phi[match0.index_gal].r_exp[2]

outstr = replicate({phi:-999., ab:-999., radius:-999., azimuthal:-999.}, n_elements(match0))

;; calculate azimuthal
;; step 1: turn position angle (of major axis) to PA (of minor axis) between -90. and 90 from North to East
new_r_phi = r_phi
ineg = where(r_phi lt 0.)
new_r_phi[ineg] = r_phi[ineg]+180.
i180 = where(r_phi gt 180.)
new_r_phi[i180] = r_phi[i180]-180.
new_r_phi = new_r_phi - 90.

;; step 2: calculate azimuthal angle (between -180. and 180.)
new_pang = pang - new_r_phi
ibig = where(new_pang gt 180.)
new_pang[ibig] = new_pang[ibig]-360.
ismall = where(new_pang lt -180.)
new_pang[ismall] = new_pang[ismall]+360.

inotgood = where(r_phi eq -1.)
new_r_phi[inotgood] = -999.
new_pang[inotgood] = -999.

outstr.phi = new_r_phi ;; original
outstr.ab = r_ab ;; original
outstr.radius = r_radius ;; original
outstr.azimuthal = new_pang

   qso_path = jhusdss_get_path(/qso)
   outfile = qso_path+'/'+jhusdss_galqso_matchfile(boss=boss)
   outfile = repstr(outfile, '.fits', '_azimuthal.fits')
 mwrfits, outstr, outfile, /create


end
