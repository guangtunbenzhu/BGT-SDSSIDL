;; temporary, should use gaussian fit of ew 
function jhusdss_montecarlo_absorbers_pool, nmfver

absorber0 = jhusdss_absorber_readin(nmfver)

imgii = where(absorber0.snr_mgii_2796 gt 6.  $
          and absorber0.snr_mgii_2803 gt 2.  $
          and absorber0.signal_mgii_2803 lt absorber0.rew_mgii_2803*1.6+0.1 $
          and absorber0.signal_mgii_2803 gt -2.0/36.*(absorber0.rew_mgii_2803-6)^2+2.0 $
          and absorber0.signal_mgii_2796 lt absorber0.rew_mgii_2796*1.6+0.1 $
          and absorber0.signal_mgii_2796 gt -2.0/36.*(absorber0.rew_mgii_2796-6)^2+2.0, nmgii)

absorber = absorber0[imgii]

;; upweight strong systems
ii = where(absorber.rew_mgii_2796 gt 2.5 and absorber.rew_mgii_2796 le 3., nn)
absorber = [absorber, absorber[ii]]
ii = where(absorber.rew_mgii_2796 gt 3. and absorber.rew_mgii_2796 le 3.5, nn)
absorber = [absorber, absorber[ii], absorber[ii]] 
ii = where(absorber.rew_mgii_2796 gt 3.5 and absorber.rew_mgii_2796 le 4., nn)
absorber = [absorber, absorber[ii], absorber[ii], absorber[ii], absorber[ii]]
ii = where(absorber.rew_mgii_2796 gt 4. and absorber.rew_mgii_2796 le 5., nn)
absorber = [absorber, absorber[ii], absorber[ii], absorber[ii], absorber[ii], absorber[ii], absorber[ii], $
                      absorber[ii], absorber[ii], absorber[ii]] 
ii = where(absorber.rew_mgii_2796 gt 5., nn)
absorber = [absorber, absorber[ii], absorber[ii], absorber[ii], absorber[ii], absorber[ii], absorber[ii], $
                      absorber[ii], absorber[ii], absorber[ii], absorber[ii], absorber[ii], absorber[ii], $
                      absorber[ii], absorber[ii], absorber[ii], absorber[ii], absorber[ii], absorber[ii], $
                      absorber[ii], absorber[ii], absorber[ii], absorber[ii], absorber[ii], absorber[ii]]
ii = where(absorber.rew_mgii_2796 lt 0.35, nn)
absorber = [absorber, absorber[ii], absorber[ii], absorber[ii], absorber[ii], absorber[ii], absorber[ii], $
                      absorber[ii], absorber[ii], absorber[ii]]

nmgii = n_elements(absorber)
;; I could have just used struct_selecttags, which would be slower
strtmp = {rew_mgii_2803:0., rew_mgii_2796:0., signal_mgii_2803:0., signal_mgii_2796:0., vdisp:0.}
outstr = replicate(strtmp, nmgii)

for i=0L, nmgii-1L do begin
    outstr[i].rew_mgii_2803 = absorber[i].rew_mgii_2803
    outstr[i].rew_mgii_2796 = absorber[i].rew_mgii_2796
    outstr[i].signal_mgii_2803 = absorber[i].signal_mgii_2803
    outstr[i].signal_mgii_2796 = absorber[i].signal_mgii_2796
    outstr[i].vdisp = absorber[i].vdisp
endfor

return, outstr
end
