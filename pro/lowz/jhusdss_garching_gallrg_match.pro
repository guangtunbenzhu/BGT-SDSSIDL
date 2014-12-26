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
pro jhusdss_garching_gallrg_match, boss=boss, overwrite=overwrite

   omegaM = 0.3
   omegaL = 0.7
   h100 = 0.7
   ch0 = 2997.92458 ;;Mpc*h^-1

;; read in SDSS lrg catalog, need ra, dec
   if (keyword_set(boss)) then begin
      lrg_path = jhusdss_get_path(/bosslrg)
      lrg_file = lrg_path+'/'+jhusdss_boss_lrgfile()
   endif else begin
      lrg_path = jhusdss_get_path(/garching)
      lrg_file = lrg_path+'/'+jhusdss_dr7_lrgfile()
   endelse

;; outfile 
   outfile = lrg_path+'/'+jhusdss_gallrg_matchfile(boss=boss)
   if (file_test(outfile) eq 1 and ~keyword_set(overwrite)) then begin
      splog, 'File already exists, Use /overwrite to overwite'
      return
   endif else begin
      splog, 'Will write into this file:'
      splog, outfile
   endelse

   lrg = mrdfits(lrg_file, 1)

;; read in garching catalog, need ra, dec, z
   garching_path = jhusdss_get_path(/garching)
   garching_file = garching_path+'/'+'gal_info_dr7_v5_2.fit.gz'
   gal = mrdfits(garching_file, 1)
   gal_angdist = angdidis(gal.z, omegaM, omegaL)*ch0/h100

;; redshift bins
;; at z>0.2, it's mostly LRGs
   delta_z = 0.02
   z_min = 0.02
   z_max = 0.30
   z_bin = jhusdss_make_bins(z_min, z_max, delta_z, nbin=z_nbin)

   angdist = angdidis(z_bin.min,omegaM,omegaL)*ch0/h100 ;; Mpc
   angdist[0] = angdist[1]

   rp_max = 3.000 ;; Mpc
   arc_max = rp_max/angdist/!dpi*180. ;;degree

;; match garching/SDSS_QSO catalog
   for iz = 0L, z_nbin-1L do begin
       ii = where(gal.z gt z_bin[iz].min and gal.z le z_bin[iz].max $
              and gal.ra gt 0. and gal.ra lt 360. $
              and gal.dec gt -90. and gal.dec lt 90., nn)
       if (nn gt 0) then begin
          spherematch, gal[ii].ra, gal[ii].dec, lrg.ra, lrg.dec, arc_max[iz], $
               m1, m2, distance12, maxmatch=0
          if (n_elements(m1) gt 1) then begin
             inotself = where(distance12 gt 1./3600., nnotself)
             if (nnotself gt 1) then begin
             if (n_elements(index_gal) eq 0) then begin
                index_gal = ii[m1[inotself]]
                index_lrg = m2[inotself]
                rp_deg = distance12[inotself]
                rp_mpc = distance12[inotself]*!dpi/180.*gal_angdist[ii[m1[inotself]]]
             endif else begin
                index_gal = [index_gal, ii[m1[inotself]]]
                index_lrg = [index_lrg, m2[inotself]]
                rp_deg = [rp_deg, distance12[inotself]]
                rp_mpc = [rp_mpc, distance12[inotself]*!dpi/180.*gal_angdist[ii[m1[inotself]]]]
             endelse
             endif
          endif
       endif
   endfor

   str_tmp = {index_gal:0L, index_lrg:0L, rp_deg:0., rp_mpc:0.}
   out_match = replicate(str_tmp, n_elements(index_gal))
   out_match.index_gal = index_gal
   out_match.index_lrg = index_lrg
   out_match.rp_deg = rp_deg
   out_match.rp_mpc = rp_mpc

   ii = where(rp_mpc le rp_max, nn)
   mwrfits, out_match[ii], outfile, /create

end
