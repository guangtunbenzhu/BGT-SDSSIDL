function jhusdss_absorber_finalcat_blank, full=full

;if (n_elements(lines) eq 0) then lines = jhusdss_absorber_finalcat_lines()

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

               zabs:0., $
               err_zabs:0., $
               vdisp:0., $
               err_vdisp:0., $

               criterion_mgii:0b, $
               criterion_mgii_feii:0b, $
               criterion_feii:0b, $

               signal_mgii_2803:0., $
               signal_mgii_2796:0., $
               snr_mgii_2803:0., $
               snr_mgii_2796:0., $

               rew_mgi_2853:0., $
               err_rew_mgi_2853:0., $
 
               rew_mgii_2803:0., $
               err_rew_mgii_2803:0., $
               rew_mgii_2796:0., $
               err_rew_mgii_2796:0., $
 
               rew_feii_2600:0., $
               err_rew_feii_2600:0., $
               rew_feii_2586:0., $
               err_rew_feii_2586:0., $
 
               rew_feii_2383:0., $
               err_rew_feii_2383:0., $
               rew_feii_2374:0., $
               err_rew_feii_2374:0., $
 
               rew_feii_2344:0., $
               err_rew_feii_2344:0., $

               vdisp_mgi_2853:0., $
               err_vdisp_mgi_2853:0., $

               vdisp_mgii_2803:0., $
               err_vdisp_mgii_2803:0., $
               vdisp_mgii_2796:0., $
               err_vdisp_mgii_2796:0., $
 
               vdisp_feii_2600:0., $
               err_vdisp_feii_2600:0., $
               vdisp_feii_2586:0., $
               err_vdisp_feii_2586:0., $

               vdisp_feii_2383:0., $
               err_vdisp_feii_2383:0., $
               vdisp_feii_2374:0., $
               err_vdisp_feii_2374:0., $

               vdisp_feii_2344:0., $
               err_vdisp_feii_2344:0. $
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

               zabs:0., $
               err_zabs:0., $
               vdisp:0., $
               err_vdisp:0., $

               criterion_mgii:0b, $
               criterion_mgii_feii:0b, $
               criterion_feii:0b, $

               rew_mgi_2853:0., $
               err_rew_mgi_2853:0., $

               rew_mgii_2803:0., $
               err_rew_mgii_2803:0., $
               rew_mgii_2796:0., $
               err_rew_mgii_2796:0., $

               rew_feii_2600:0., $
               err_rew_feii_2600:0., $
               rew_feii_2586:0., $
               err_rew_feii_2586:0., $

               rew_feii_2383:0., $
               err_rew_feii_2383:0., $
               rew_feii_2374:0., $
               err_rew_feii_2374:0., $

               rew_feii_2344:0., $
               err_rew_feii_2344:0., $

               rew_aliii_1863:0., $
               err_rew_aliii_1863:0., $
               rew_aliii_1855:0., $
               err_rew_aliii_1855:0., $

               rew_alii_1671:0., $
               err_rew_alii_1671:0., $

               rew_civ_1551:0., $
               err_rew_civ_1551:0., $
               rew_civ_1548:0., $
               err_rew_civ_1548:0., $

               rew_siii_1527:0., $
               err_rew_siii_1527:0., $

               vdisp_mgi_2853:0., $
               err_vdisp_mgi_2853:0., $

               vdisp_mgii_2803:0., $
               err_vdisp_mgii_2803:0., $
               vdisp_mgii_2796:0., $
               err_vdisp_mgii_2796:0., $
 
               vdisp_feii_2600:0., $
               err_vdisp_feii_2600:0., $
               vdisp_feii_2586:0., $
               err_vdisp_feii_2586:0., $

               vdisp_feii_2383:0., $
               err_vdisp_feii_2383:0., $
               vdisp_feii_2374:0., $
               err_vdisp_feii_2374:0., $

               vdisp_feii_2344:0., $
               err_vdisp_feii_2344:0., $

               vdisp_aliii_1863:0., $
               err_vdisp_aliii_1863:0., $
               vdisp_aliii_1855:0., $
               err_vdisp_aliii_1855:0., $
 
               vdisp_alii_1671:0., $
               err_vdisp_alii_1671:0., $

               vdisp_civ_1551:0., $
               err_vdisp_civ_1551:0., $
               vdisp_civ_1548:0., $
               err_vdisp_civ_1548:0., $

               vdisp_siii_1527:0., $
               err_vdisp_siii_1527:0. $
              }
endelse

return, absorber
end
