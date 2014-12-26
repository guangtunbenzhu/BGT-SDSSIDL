pro jhusdss_qaplot_completeness_makeimage, nmfver
    mc = jhusdss_montecarlo_readin(nmfver)
    fall = float(mc.ndetected)/(float(mc.ncovered+(mc.ncovered eq 0)))
    fexpand = congrid(fall, 38*20, 776, /interp) 
    ;[0.2->5.96, 4000->9200]
    outfile = jhusdss_get_path(/nmfqso)+'/'+string(nmfver, format='(i3.3)')+'/MonteCarlo/Completeness_Image_'$
              + string(nmfver, format='(i3.3)')+'.fits'
;   mwrfits, fexpand[0:330*2, 0:371*2], outfile, /create
    ;; w[33] = 5.97
    ;; z[760] = 2.26
    mwrfits, fexpand[0:33*20, 0:760], outfile, /create
end
