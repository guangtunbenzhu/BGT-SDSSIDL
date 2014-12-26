;+
;-
function jhusdss_finalpass_blank
   return, {rew_mgi_2853:0.0,   err_rew_mgi_2853:0.0,   sigma_mgi_2853:0.0,   err_sigma_mgi_2853:0.0, $
            rew_feii_2344:0.0,  err_rew_feii_2344:0.0,  sigma_feii_2344:0.0,  err_sigma_feii_2344:0.0, $
            rew_alii_1671:0.0,  err_rew_alii_1671:0.0,  sigma_alii_1671:0.0,  err_sigma_alii_1671:0.0, $
            rew_siii_1527:0.0,  err_rew_siii_1527:0.0,  sigma_siii_1527:0.0,  err_sigma_siii_1527:0.0, $
            rew_mgii_2803:0.0,  err_rew_mgii_2803:0.0,  sigma_mgii_2803:0.0,  err_sigma_mgi_2803:0.0, $
            rew_mgii_2796:0.0,  err_rew_mgii_2796:0.0,  sigma_mgii_2796:0.0,  err_sigma_mgi_2796:0.0, $
            rew_feii_2600:0.0,  err_rew_feii_2600:0.0,  sigma_feii_2600:0.0,  err_sigma_feii_2600:0.0, $
            rew_feii_2586:0.0,  err_rew_feii_2586:0.0,  sigma_feii_2586:0.0,  err_sigma_feii_2586:0.0, $
            rew_feii_2383:0.0,  err_rew_feii_2383:0.0,  sigma_feii_2383:0.0,  err_sigma_feii_2383:0.0, $
            rew_feii_2374:0.0,  err_rew_feii_2374:0.0,  sigma_feii_2374:0.0,  err_sigma_feii_2374:0.0, $
            rew_aliii_1863:0.0, err_rew_aliii_1863:0.0, sigma_aliii_1863:0.0, err_sigma_aliii_1863:0.0, $
            rew_aliii_1855:0.0, err_rew_aliii_1855:0.0, sigma_aliii_1855:0.0, err_sigma_aliii_1855:0.0, $
            rew_civ_1551:0.0,   err_rew_civ_1551:0.0,   sigma_civ_1551:0.0,   err_sigma_civ_1551:0.0, $
            rew_civ_1548:0.0,   err_rew_civ_1548:0.0,   sigma_civ_1548:0.0,   err_sigma_civ_1548:0.0}
end

function jhusdss_finalpass_fit, spec, zabs, newzabs=newzabs

outstr = jhusdss_finalpass_blank()

lambda = spec.wave*(1.+spec.z)
flux = 1.-spec.residual
ivar = spec.ivar

;; first pass, use input
;; second pass, use the results of first pass, and fix redshift
tmpz = zabs
for ipass=0,1 do begin
    zabs_fit = -999.

    ;; center's position
    case ipass of
       0: raw_maxwidth=10.D-4
       1: raw_maxwidth=2.D-4
    endcase

;; single gaussian for MgI 2853
;; single gaussian for FeII 2344
;; single gaussian for AlII 1671
;; single gaussian for SiII 1527
   for isingle=0,3 do begin
       case isingle of
          0: begin
             thisline='MgI_2853'
             thiscenter=2852.96D0
             end
          1: begin
             thisline='FeII_2344'
             thiscenter=2344.21D0
             end
          2: begin
             thisline='AlII_1671'
             thiscenter=1670.79D0
             end
          3: begin
             thisline='SiII_1527'
             thiscenter=1526.71D0
             end
       endcase
       ilambda = where(lambda gt thiscenter/(1.+tmpz)*(1.D0-3.D-3) $
                   and lambda lt thiscenter/(1.+tmpz)*(1.D0+3.D-3) $
                   and ivar gt 0., nlambda)

       if (nlambda gt 30) then begin
           in_lambda = lambda[ilambda]
           in_flux = flux[ilambda]
           in_ivar = ivar[ilambda]
           ;; initial guess
           in_center = thisline*(1.+tmpz)
           in_line_flux = (((total(in_flux)*in_center*1.D-4*alog(10.D0)) > 0.1) < 6.0)

           maxwidth = raw_maxwidth*median(in_lambda)*alog(10D0)
           jhusdss_singlegaussian_fit, in_lambda, in_flux, in_ivar, in_center, in_line_flux, $
                center=center, err_center=err_center, lflux=lflux, err_lflux=err_lflux, $
                sigma=sigma, err_sigma=err_sigma, maxwidth=maxwidth

           if (ipass eq 1) then begin
              ;; sigma in dlambda to sigma in dpixel
              sigma = sigma/alog(10.)/in_center*1.D+4
              err_sigma = err_sigma/alog(10.)/in_center*1.D+4
              tmpexec=execute('outstr.REW_'+thisline+'=lflux/(1.+newzabs)')
              tmpexec=execute('outstr.ERR_REW_'+thisline+'=err_lflux/(1.+newzabs)')
              tmpexec=execute('outstr.SIGMA_'+thisline+'=69.03*sqrt(sigma^2-1.)')
              tmpexec=execute('outstr.ERR_SIGMA_'+thisline+'=69.03*err_sigma*sigma/sqrt(sigma^2-1.)')
           endif
           zabs_fit = [zabs_fit, (center/thisline-1.)]
       endif
   endfor

;; double gaussian for MgII  2803/2796
;; double gaussian for FeII  2600/2586
;; double gaussian for FeII  2383/2374
;; double gaussian for AlIII 1862/1854
;; double gaussian for CIV   1551/1548
   for idouble=0, 4 do begin
       case idouble of
          0: begin
             these2lines=['MgII_2803', 'MgII_2796']
             these2centers=[2803.53D0, 2796.35D0]
             end
          1: begin
             these2lines=['FeII_2600', 'FeII_2586']
             these2centers=[2600.17D0, 2586.55D0]
             end
          2: begin
             these2lines=['FeII_2383', 'FeII_2374']
             these2centers=[2382.77D0, 2374.46D0]
             end
          3: begin
             these2lines=['AlIII_1863', 'AlIII_1855']
             these2centers=[1862.79D0, 1854.72D0]
             end
          4: begin
             these2lines=['CIV_1551', 'CIV_1548']
             these2centers=[1550.78D0, 1548.20D0]
             end
       endcase

       ilambda = where(lambda gt these2centers[1]/(1.+tmpz)*(1.D0-3.D-3) $
                   and lambda lt these2centers[0]/(1.+tmpz)*(1.D0+3.D-3) $
                   and ivar gt 0., nlambda)

       if (nlambda gt 30) then begin
           in_lambda = lambda[ilambda]
           in_flux = flux[ilambda]
           in_ivar = ivar[ilambda]

           ;; initial guess
           in_center = these2lines*(1.+tmpz)
           tmpflux = (((total(in_flux)*total(in_center)/2.*1.D-4*alog(10.D0)) > 0.1) < 6.)
           ;; 2803/2796 ~ 1/1.5, 1863/1855 ~ 1/1.5, 1550/1548 ~ 1/1.5
           in_line_flux = [tmpflux/2.3*1., tmpflux/2.3*1.3]
           ;; 2600/2586 ~ 1.5, 2383/2374 ~ 1.5, 
           if (idouble eq 1 or idouble eq 2) then $
              in_line_flux = [tmpflux/2.3*1.3, tmpflux/2.3*1.]

           maxwidth = raw_maxwidth*median(in_lambda)*alog(10D0)
           jhusdss_doublegaussian_fit, in_lambda, in_flux, in_ivar, in_center, in_line_flux, $
                center=center, err_center=err_center, lflux=lflux, err_lflux=err_lflux, $
                sigma=sigma, err_sigma=err_sigma, maxwidth=maxwidth

           if (ipass eq 1) then begin
              ;; sigma in dlambda to sigma in dpixel
              sigma = sigma/alog(10.)/in_center*1.D+4
              err_sigma = err_sigma/alog(10.)/in_center*1.D+4
              for iline=0,1 do begin
                  tmpexec=execute('outstr.REW_'+these2lines[iline] $
                                + '=lflux['+string(iline,format='(i1.1)')+']/(1.+newzabs)')
                  tmpexec=execute('outstr.ERR_REW_'+these2lines[iline] $
                                + '=err_lflux['+string(iline,format='(i1.1)')+']/(1.+newzabs)')
                  tmpexec=execute('outstr.SIGMA_'+these2lines[iline] $
                                + '=69.03*sqrt(sigma['+string(iline,format='(i1.1)')+']^2-1.)')
                  tmpexec=execute('outstr.ERR_SIGMA_'+these2lines[iline] $
                                + '=69.03*err_sigma['+string(iline,format='(i1.1)') $
                                + ']*sigma['+string(iline,format='(i1.1)')+']/sqrt(sigma[' $
                                + string(iline,format='(i1.1)')+']^2-1.)')
              endfor
           endif
           zabs_fit = [zabs_fit, (center[0]/these2lines[0]-1.), (center[1]/these2lines[1]-1.)]
       endif

    endfor

;;  new redshift
    RESISTANT_Mean, zabs_fit[1:*], 3, newzabs, newzabs_sig, num_reject
    tmpz = newzabs
endfor

return, outstr
end
