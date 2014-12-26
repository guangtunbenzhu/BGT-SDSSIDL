;+
;-

pro jhusdss_finalpass_fit, spec, zabs, lines=lines, newzabs=newzabs, $
            err_newzabs=err_newzabs, allzabs=allzabs, ew=out_ew, err_ew=out_err_ew, sigma=out_sigma, $
            err_sigma=out_err_sigma, meansigma=meansigma, err_meansigma=err_meansigma

if (n_elements(lines) eq 0) then lines = jhusdss_finalpass_lines()

out_ew = fltarr(n_elements(lines))
out_err_ew = fltarr(n_elements(lines))
out_sigma = fltarr(n_elements(lines))
out_err_sigma = fltarr(n_elements(lines))
allzabs = fltarr(n_elements(lines))

lambda = spec.wave*(1.+spec.z)
flux = 1.-spec.residual
ivar = spec.ivar

;; first pass, use input
;; second pass, use the results of first pass, and fix redshift
;; It looks like the MgII redshift usually is smaller than FeII one?
;; It might be that MgII's wavelength is wrong.
;; We only use the first pass now.
tmpz = zabs
for ipass=0,0 do begin
    mgii_used = 0b
    zabs_fit = -999.
    sigma_fit = -999.

    ;; center's position
    case ipass of
       0: raw_maxwidth=5.D-4
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
             end
          1: begin
             thisline='FeII_2344'
             end
          2: begin
             thisline='AlII_1671'
             end
          3: begin
             thisline='SiII_1527'
             end
       endcase

       thisindex=where(strcmp(lines.name, thisline, /fold_case) eq 1, nindex) 
       if (nindex ne 1) then message, 'Fix this first.'
       thiscenter = lines[thisindex].wave

       ilambda = where(lambda gt thiscenter*(1.+tmpz)*(1.D0-3.D-3) $ ;; 3.D-3/alog(10.)=0.0013
                   and lambda lt thiscenter*(1.+tmpz)*(1.D0+3.D-3) $ ;; 13 pixels
                   and ivar gt 0., nlambda)

       if (nlambda ge 8) then begin
           in_lambda = lambda[ilambda]
           in_flux = flux[ilambda]
           in_ivar = ivar[ilambda]
           ;; initial guess
           in_center = thiscenter*(1.+tmpz)
           in_line_flux = (((total(in_flux)*in_center*1.D-4*alog(10.D0)) > 0.01) < 12.0)

           maxwidth = raw_maxwidth*median(in_lambda)*alog(10.D0)
           jhusdss_singlegaussian_fit, in_lambda, in_flux, in_ivar, in_center, in_line_flux, $
                center=center, err_center=err_center, lflux=lflux, err_lflux=err_lflux, $
                sigma=sigma, err_sigma=err_sigma, maxwidth=maxwidth

           if (ipass eq 0) then begin
              ;; sigma in dlambda to sigma in dpixel
              sigma = sigma/alog(10.)/in_center*1.D+4
              err_sigma = err_sigma/alog(10.)/in_center*1.D+4
              out_ew[thisindex] = lflux/(1.+tmpz)
              out_err_ew[thisindex] = err_lflux/(1.+tmpz)
              out_sigma[thisindex] = 69.03*sqrt((sigma^2-1.)>1.E-6)
              out_err_sigma[thisindex] = 69.03*err_sigma*sigma/sqrt((sigma^2-1.)>1.E-6)
              allzabs[thisindex] = center/thiscenter-1.
           endif
;          if (isingle eq 1) then begin
;             zabs_fit = [zabs_fit, (center/thiscenter-1.)]
;             sigma_fit = [sigma_fit,  69.03*sqrt((sigma^2-1.)>1.E-6)]
;          endif
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
             end
          1: begin
             these2lines=['FeII_2600', 'FeII_2586']
             end
          2: begin
             these2lines=['FeII_2383', 'FeII_2374']
             end
          3: begin
             these2lines=['AlIII_1863', 'AlIII_1855']
             end
          4: begin
             these2lines=['CIV_1551', 'CIV_1548']
             end
       endcase

       these2indices = lonarr(2)
       these2centers = dblarr(2)
       for itmp=0,1 do begin
           indextmp=where(strcmp(lines.name, these2lines[itmp], /fold_case) eq 1, nindex) 
           if (nindex ne 1) then message, 'Fix this first.'
           these2indices[itmp] = indextmp
           these2centers[itmp] = lines[indextmp].wave
       endfor

       ilambda = where(lambda gt these2centers[1]*(1.+tmpz)*(1.D0-3.D-3) $
                   and lambda lt these2centers[0]*(1.+tmpz)*(1.D0+3.D-3) $
                   and ivar gt 0., nlambda)

       if (nlambda ge 17) then begin
           in_lambda = lambda[ilambda]
           in_flux = flux[ilambda]
           in_ivar = ivar[ilambda]

           ;; initial guess
           in_center = these2centers*(1.+tmpz)
           tmpflux = (((total(in_flux)*mean(in_center)*1.D-4*alog(10.D0)) > 0.01) < 12.)
           ;; 2803/2796 ~ 1/1.5, 1863/1855 ~ 1/1.5, 1550/1548 ~ 1/1.5
           in_line_flux = [tmpflux/2.3*1., tmpflux/2.3*1.3]
           ;; 2600/2586 ~ 1.5, 2383/2374 ~ 1.5, 
           if (idouble eq 1 or idouble eq 2) then $
              in_line_flux = [tmpflux/2.3*1.3, tmpflux/2.3*1.]

           maxwidth = raw_maxwidth*median(in_lambda)*alog(10.D0)
           jhusdss_doublegaussian_fit, in_lambda, in_flux, in_ivar, in_center, in_line_flux, $
                center=center, err_center=err_center, lflux=lflux, err_lflux=err_lflux, $
                sigma=sigma, err_sigma=err_sigma, maxwidth=maxwidth

           if (ipass eq 0) then begin
              ;; sigma in dlambda to sigma in dpixel
              sigma = sigma/alog(10.)/in_center*1.D+4
              err_sigma = err_sigma/alog(10.)/in_center*1.D+4
              out_ew[these2indices] = lflux/(1.+tmpz)
              out_err_ew[these2indices] = err_lflux/(1.+tmpz)
              out_sigma[these2indices] = 69.03*sqrt((sigma^2-1.)>1.E-6)
              out_err_sigma[these2indices] = 69.03*err_sigma*sigma/sqrt((sigma^2-1.)>1.E-6)
              allzabs[these2indices] = [(center[0]/these2centers[0]-1.), (center[1]/these2centers[1]-1.)]
           endif
           if (idouble eq 0) then begin
              mgii_used = 1b
              zabs_fit = [zabs_fit, (center[0]/these2centers[0]-1.), (center[1]/these2centers[1]-1.)]
              sigma_fit = [sigma_fit,  69.03*sqrt((sigma^2-1.)>1.E-6)]
           endif
           if (~mgii_used) and (idouble eq 1 or idouble eq 2) then begin
              zabs_fit = [zabs_fit, (center[0]/these2centers[0]-1.), (center[1]/these2centers[1]-1.)]
              sigma_fit = [sigma_fit,  69.03*sqrt((sigma^2-1.)>1.E-6)]
           endif
       endif

    endfor

;;  new redshift
;   print, zabs_fit[1:*]
    RESISTANT_Mean, zabs_fit[1:*], 5, newzabs, err_newzabs, num_reject
    tmpz = newzabs
endfor
RESISTANT_Mean, sigma_fit[1:*], 5, meansigma, err_meansigma, num_reject
;print, newzabs, err_newzabs
;print, meansigma, err_meansigma

end
