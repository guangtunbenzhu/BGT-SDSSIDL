function jhusdss_abs2qso_blank, nabsmax=nabsmax, full=full

if (n_elements(nabsmax) eq 0) then nabsmax = jhusdss_nabsmax()

if (~keyword_set(full)) then begin
   absorber = {ra:0D0, $
               dec:0D0, $
               plate:0L, $
               fiber:0L, $
               mjd:0L, $
               zqso:0., $
               err_zqso:0., $
               index_qso:0L, $

               spec_snr_median:-999., $
               med_sdeviation_red:-999., $
               med_sdeviation_blue:-999., $

               nabs:0L, $

               zabs:fltarr(nabsmax), $
               err_zabs:fltarr(nabsmax), $
               vdisp:fltarr(nabsmax), $
               err_vdisp:fltarr(nabsmax), $

               criterion_mgii:bytarr(nabsmax), $
               criterion_mgii_feii:bytarr(nabsmax), $
               criterion_feii:bytarr(nabsmax), $

               signal_mgii_2803:fltarr(nabsmax), $
               signal_mgii_2796:fltarr(nabsmax), $
               snr_mgii_2803:fltarr(nabsmax), $
               snr_mgii_2796:fltarr(nabsmax), $

               rew_mgi_2853:fltarr(nabsmax), $
               err_rew_mgi_2853:fltarr(nabsmax), $
 
               rew_mgii_2803:fltarr(nabsmax), $
               err_rew_mgii_2803:fltarr(nabsmax), $
               rew_mgii_2796:fltarr(nabsmax), $
               err_rew_mgii_2796:fltarr(nabsmax), $
 
               rew_feii_2600:fltarr(nabsmax), $
               err_rew_feii_2600:fltarr(nabsmax), $
               rew_feii_2586:fltarr(nabsmax), $
               err_rew_feii_2586:fltarr(nabsmax), $
 
               rew_feii_2383:fltarr(nabsmax), $
               err_rew_feii_2383:fltarr(nabsmax), $
               rew_feii_2374:fltarr(nabsmax), $
               err_rew_feii_2374:fltarr(nabsmax), $
 
               rew_feii_2344:fltarr(nabsmax), $
               err_rew_feii_2344:fltarr(nabsmax), $

               vdisp_mgi_2853:fltarr(nabsmax), $
               err_vdisp_mgi_2853:fltarr(nabsmax), $

               vdisp_mgii_2803:fltarr(nabsmax), $
               err_vdisp_mgii_2803:fltarr(nabsmax), $
               vdisp_mgii_2796:fltarr(nabsmax), $
               err_vdisp_mgii_2796:fltarr(nabsmax), $
 
               vdisp_feii_2600:fltarr(nabsmax), $
               err_vdisp_feii_2600:fltarr(nabsmax), $
               vdisp_feii_2586:fltarr(nabsmax), $
               err_vdisp_feii_2586:fltarr(nabsmax), $

               vdisp_feii_2383:fltarr(nabsmax), $
               err_vdisp_feii_2383:fltarr(nabsmax), $
               vdisp_feii_2374:fltarr(nabsmax), $
               err_vdisp_feii_2374:fltarr(nabsmax), $

               vdisp_feii_2344:fltarr(nabsmax), $
               err_vdisp_feii_2344:fltarr(nabsmax) $
               }
endif else begin
   absorber = {ra:0D0, $
               dec:0D0, $
               plate:0L, $
               fiber:0L, $
               mjd:0L, $
               zqso:0., $
               err_zqso:0., $
               index_qso:0L, $

               spec_snr_median:-999., $
               med_sdeviation_red:-999., $
               med_sdeviation_blue:-999., $

               nabs:0L, $

               zabs:fltarr(nabsmax), $
               err_zabs:fltarr(nabsmax), $
               vdisp:fltarr(nabsmax), $
               err_vdisp:fltarr(nabsmax), $

               criterion_mgii:bytarr(nabsmax), $
               criterion_mgii_feii:bytarr(nabsmax), $
               criterion_feii:bytarr(nabsmax), $

               rew_mgi_2853:fltarr(nabsmax), $
               err_rew_mgi_2853:fltarr(nabsmax), $

               rew_mgii_2803:fltarr(nabsmax), $
               err_rew_mgii_2803:fltarr(nabsmax), $
               rew_mgii_2796:fltarr(nabsmax), $
               err_rew_mgii_2796:fltarr(nabsmax), $

               rew_feii_2600:fltarr(nabsmax), $
               err_rew_feii_2600:fltarr(nabsmax), $
               rew_feii_2586:fltarr(nabsmax), $
               err_rew_feii_2586:fltarr(nabsmax), $

               rew_feii_2383:fltarr(nabsmax), $
               err_rew_feii_2383:fltarr(nabsmax), $
               rew_feii_2374:fltarr(nabsmax), $
               err_rew_feii_2374:fltarr(nabsmax), $

               rew_feii_2344:fltarr(nabsmax), $
               err_rew_feii_2344:fltarr(nabsmax), $

               rew_aliii_1863:fltarr(nabsmax), $
               err_rew_aliii_1863:fltarr(nabsmax), $
               rew_aliii_1855:fltarr(nabsmax), $
               err_rew_aliii_1855:fltarr(nabsmax), $

               rew_alii_1671:fltarr(nabsmax), $
               err_rew_alii_1671:fltarr(nabsmax), $

               rew_civ_1551:fltarr(nabsmax), $
               err_rew_civ_1551:fltarr(nabsmax), $
               rew_civ_1548:fltarr(nabsmax), $
               err_rew_civ_1548:fltarr(nabsmax), $

               rew_siii_1527:fltarr(nabsmax), $
               err_rew_siii_1527:fltarr(nabsmax), $

               vdisp_mgi_2853:fltarr(nabsmax), $
               err_vdisp_mgi_2853:fltarr(nabsmax), $

               vdisp_mgii_2803:fltarr(nabsmax), $
               err_vdisp_mgii_2803:fltarr(nabsmax), $
               vdisp_mgii_2796:fltarr(nabsmax), $
               err_vdisp_mgii_2796:fltarr(nabsmax), $
 
               vdisp_feii_2600:fltarr(nabsmax), $
               err_vdisp_feii_2600:fltarr(nabsmax), $
               vdisp_feii_2586:fltarr(nabsmax), $
               err_vdisp_feii_2586:fltarr(nabsmax), $

               vdisp_feii_2383:fltarr(nabsmax), $
               err_vdisp_feii_2383:fltarr(nabsmax), $
               vdisp_feii_2374:fltarr(nabsmax), $
               err_vdisp_feii_2374:fltarr(nabsmax), $

               vdisp_feii_2344:fltarr(nabsmax), $
               err_vdisp_feii_2344:fltarr(nabsmax), $

               vdisp_aliii_1863:fltarr(nabsmax), $
               err_vdisp_aliii_1863:fltarr(nabsmax), $
               vdisp_aliii_1855:fltarr(nabsmax), $
               err_vdisp_aliii_1855:fltarr(nabsmax), $
 
               vdisp_alii_1671:fltarr(nabsmax), $
               err_vdisp_alii_1671:fltarr(nabsmax), $

               vdisp_civ_1551:fltarr(nabsmax), $
               err_vdisp_civ_1551:fltarr(nabsmax), $
               vdisp_civ_1548:fltarr(nabsmax), $
               err_vdisp_civ_1548:fltarr(nabsmax), $

               vdisp_siii_1527:fltarr(nabsmax), $
               err_vdisp_siii_1527:fltarr(nabsmax) $
              }
endelse

return, absorber
end
