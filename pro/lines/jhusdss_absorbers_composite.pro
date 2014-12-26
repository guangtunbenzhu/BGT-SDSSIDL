pro jhusdss_absorbers_composite, nmfver, boss=boss, overwrite=overwrite, deltaz=dz

    absorber = jhusdss_absorber_readin(nmfver, /byabs, /trim, boss=boss)

    ;; zbins
    if (n_elements(dz) eq 0) then dz = 0.02
    zmin = 0.43
    zmax = 2.29
    nzbin = floor((zmax-zmin)/dz)
    zbin_min = findgen(nzbin)*dz+zmin
    zbin_max = findgen(nzbin)*dz+zmin+dz

    compath = jhusdss_get_path(/nmfqso)+'/'+string(nmfver, format='(I3.3)')+'/Composite'
    for iz=0L, nzbin-1L do begin
        filename = 'Absorbers_Composite_z'+string(zbin_min[iz], format='(f4.2)')+'_'$
                 + string(zbin_max[iz], format='(f4.2)')+'.fits'
        outfile = compath + '/' + filename
        if (file_test(outfile) and ~keyword_set(overwrite)) then begin 
           splog, 'File already exists. Use /overwrite if you want to overwrite it.'
           return
        endif
        iabs = where(absorber.zabs gt zbin_min[iz] and absorber.zabs le zbin_max[iz], nabs)
        if (nabs eq 0) then continue

        composite = jhusdss_absorbers_composite_engine(absorber[iabs], nmfver=nmfver)

        mwrfits, composite, outfile, /create
    endfor

end
