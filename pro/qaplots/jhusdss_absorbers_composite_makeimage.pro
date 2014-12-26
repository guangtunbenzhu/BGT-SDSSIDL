;;  See jhusdss_absorbers_composite, nmfver, boss=boss, overwrite=overwrite
pro jhusdss_absorbers_composite_makeimage, nmfver, boss=boss, overwrite=overwrite, deltaz=dz

    ;; zbins
    if (n_elements(dz) eq 0) then dz = 0.02
    compath = jhusdss_get_path(/nmfqso)+'/'+string(nmfver, format='(I3.3)')+'/Composite'
    outfile=compath+'/'+'Composite_Image_dz'+string(dz, format='(f4.2)')+'.fits'
    if (file_test(outfile) and ~keyword_set(overwrite)) then begin
       splog, 'File exists. Not overwriting.'
       return
    endif

    zmin = 0.43
    zmax = 2.29
    nzbin = floor((zmax-zmin)/dz)
    zbin_min = findgen(nzbin)*dz+zmin
    zbin_max = findgen(nzbin)*dz+zmin+dz

;   nz_all = 4000L ; == 2.0/0.4/0.00125
    dz_all = 0.0005
    nexpand = fix(dz/dz_all)
    nz_all = floor((zmax-zmin)/dz_all)
    zbin_all = findgen(nz_all)*dz_all+zmin+0.5*dz_all

    ;; observer-frame wavelength
    minwave = 3750D0
    maxwave = 9850D0
    loglam = jhusdss_get_loglam(minwave=minwave, maxwave=maxwave)
    iwave = where(loglam gt alog10(3800.) and loglam lt alog10(9200), nwave)
    outimage = fltarr(nwave, nz_all)

    ;; read in all then interpolate
    for iz=0L, nzbin-1L do begin
        filename = 'Absorbers_Composite_z'+string(zbin_min[iz], format='(f4.2)')+'_'$
                 + string(zbin_max[iz], format='(f4.2)')+'.fits'
        infile = compath + '/' + filename
        if (file_test(infile) eq 0) then begin
           splog, 'File doesnot exist. Please Check!'
           return
        endif
        comp = mrdfits(infile, 1)
        if (n_elements(inimage) eq 0) then begin
           inwave = comp.wave
           inimage = fltarr(n_elements(comp.wave), nzbin)+1.
        endif
        inimage[*,iz] = smooth(comp.fluxmedian, 5)
    endfor

    inimage_expand = congrid(inimage, n_elements(inwave), nz_all, /interp)
    ;; take care of the edge using nearest neighbor
    for iz=0L, nzbin-2L do begin
        istart = (where(inwave gt 9200D0/(1.+zbin_min[iz+1])))[0]
        izstart = ((iz*nexpand)>0)
        izend = (((iz+1)*nexpand-1L)<(nz_all-1))
        for j=izstart, izend do inimage_expand[istart:n_elements(inwave)-1, j] = inimage[istart:n_elements(inwave)-1, iz]
    endfor

;   save, filename='temp.sav', inimage_expand, inimage_expand0, inimage
;   stop

    for jz=0L, nz_all-1L do begin
        counter, jz+1, nz_all
        zz = zbin_all[jz]
        wave = inwave*(1.+zz)
        flux = inimage_expand[*,jz]
        ivar = fltarr(n_elements(wave))+1.
        curr_loglam = alog10(wave)
        combine1fiber, curr_loglam, flux, ivar, newloglam=loglam, newflux=newflux
;       outimage[*,jz] = (smooth(newflux, 5))[iwave]
        outimage[*,jz] = newflux[iwave]
    endfor

    mwrfits, outimage, outfile, /create
    mwrfits, inimage_expand, outfile
    mwrfits, inimage, outfile
    
end
