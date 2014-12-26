pro jhusdss_absorbers_composite_allabs, nmfver, overwrite=overwrite

    compath = jhusdss_get_path(/nmfqso)+'/'+string(nmfver, format='(I3.3)')+'/Composite'
    filename = 'Absorbers_Composite_Allabs_0.8AA.fits'
    outfile = compath + '/' + filename
    if (file_test(outfile) and ~keyword_set(overwrite)) then begin 
       splog, 'File already exists. Use /overwrite if you want to overwrite it.'
       return
    endif

    absorber = jhusdss_absorber_readin(nmfver)
    tmp = replicate({boss:0}, n_elements(absorber))
    absorber = struct_addtags(absorber, tmp)

    absorber_boss = jhusdss_absorber_readin(nmfver, /boss)
    tmp = replicate({boss:1}, n_elements(absorber_boss))
    absorber_boss = struct_addtags(absorber_boss, tmp)

    allabsorber = [absorber, absorber_boss]

    ;; zbins
    zmin = 0.30
    zmax = 2.79
    iabs = where(allabsorber.zabs gt zmin and allabsorber.zabs lt zmax and allabsorber.rew_mgii_2796 gt 0.8, nabs)
    help, allabsorber, nabs

    composite = jhusdss_absorbers_composite_engine(allabsorber[iabs], nmfver=nmfver)

    mwrfits, composite, outfile, /create

end
