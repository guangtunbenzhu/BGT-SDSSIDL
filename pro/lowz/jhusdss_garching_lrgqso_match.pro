;+
; 200 kpc
; z=0.02, ~400''
; z=0.05, ~200''
; z=0.10, ~100''
; z=0.15, ~80''
; z=0.20, ~60''
; z=0.25, ~50''
; z=0.30, ~45''
;-
pro jhusdss_garching_lrgqso_match, sdsslrg=sdsslrg, boss=boss, overwrite=overwrite

   omegaM = 0.3
   omegaL = 0.7
   h100 = 0.7
   ch0 = 2997.92458 ;;Mpc*h^-1

;; read in SDSS qso catalog, need ra, dec
   qso_path = jhusdss_get_path(/qso)
   if (keyword_set(boss)) then begin
      qso_file = qso_path+'/'+jhusdss_boss_qsofile()
   endif else begin
      qso_file = qso_path+'/'+jhusdss_dr7_qsofile()
   endelse

;; outfile 
   outfile = qso_path+'/'+jhusdss_lrgqso_matchfile(sdsslrg=sdsslrg, boss=boss)
   if (file_test(outfile) eq 1 and ~keyword_set(overwrite)) then begin
      splog, outfile
      splog, 'File already exists, Use /overwrite to overwite'
      return
   endif else begin
      splog, 'Will write into this file:'
      splog, outfile
   endelse

   qso = mrdfits(qso_file, 1)

;; read in SDSS or BOSS lrg catalog, need ra, dec
   if (not keyword_set(sdsslrg)) then begin
      lrg_path = jhusdss_get_path(/bosslrg)
      lrg_file = lrg_path+'/'+jhusdss_boss_lrgfile()
   endif else begin
      lrg_path = jhusdss_get_path(/garching)
      lrg_file = lrg_path+'/'+jhusdss_dr7_lrgfile()
   endelse
   gal = mrdfits(lrg_file, 1)
   gal_angdist = angdidis(gal.z, omegaM, omegaL)*ch0/h100

;; redshift bins
;; at z>0.2, it's mostly LRGs
   delta_z = 0.02
   z_min = 0.02
   z_max = 0.80
   z_bin = jhusdss_make_bins(z_min, z_max, delta_z, nbin=z_nbin)

   angdist = angdidis(z_bin.min,omegaM,omegaL)*ch0/h100 ;; Mpc
   angdist[0] = angdist[1]

   rp_max = 30.000 ;; Mpc
   arc_max = rp_max/angdist/!dpi*180. ;;degree

;; match garching/SDSS_QSO catalog
   for iz = 0L, z_nbin-1L do begin
       ii = where(gal.z gt z_bin[iz].min and gal.z le z_bin[iz].max $
              and gal.ra gt 0. and gal.ra lt 360. $
              and gal.dec gt -90. and gal.dec lt 90., nn)
       if (nn gt 0) then begin
          spherematch, gal[ii].ra, gal[ii].dec, qso.ra, qso.dec, arc_max[iz], $
               m1, m2, distance12, maxmatch=0
          if (n_elements(m1) gt 1) then begin
             if (n_elements(index_gal) eq 0) then begin
                index_gal = ii[m1]
                index_qso = m2
                rp_deg = distance12
                rp_mpc = distance12*!dpi/180.*gal_angdist[ii[m1]]
             endif else begin
                index_gal = [index_gal, ii[m1]]
                index_qso = [index_qso, m2]
                rp_deg = [rp_deg, distance12]
                rp_mpc = [rp_mpc, distance12*!dpi/180.*gal_angdist[ii[m1]]]
             endelse
          endif
       endif
   endfor

   str_tmp = {index_lrg:0L, index_qso:0L, rp_deg:0., rp_mpc:0.}
   out_match = replicate(str_tmp, n_elements(index_gal))
   out_match.index_lrg = index_gal
   out_match.index_qso = index_qso
   out_match.rp_deg = rp_deg
   out_match.rp_mpc = rp_mpc

   ii = where(rp_mpc le rp_max, nn)
   mwrfits, out_match[ii], outfile, /create

end
