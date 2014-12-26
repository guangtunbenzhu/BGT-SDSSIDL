pro jhusdss_detect_absorbers_convolve_spec, spec0, nmfver, lines=lines, emlines=emlines, $
       overwrite=overwrite, boss=boss, path=path

nspec = n_elements(spec0)
if (n_elements(lines) eq 0) then lines = jhusdss_abslines_all(/train)
if (n_elements(emlines) eq 0) then emlines = jhusdss_emlines_all(/mask)
if (n_elements(waveminmax) eq 0) then waveminmax = jhusdss_sdsswave_minmax()

;; path
if (~keyword_set(path)) then begin
   if (keyword_set(boss)) then begin
       path=jhusdss_get_path(/nmfqso)+'/'+$
          string(nmfver, format='(I3.3)')+'/Convolve_BOSS'
   endif else begin
       path=jhusdss_get_path(/nmfqso)+'/'+$
          string(nmfver, format='(I3.3)')+'/Convolve'
   endelse
endif

if (jhusdss_direxist(path) eq 0) then spawn, 'mkdir -p '+path

;; Make a  common redshift grid
zmin = waveminmax[0]/2796D0-1.
zmax = waveminmax[1]/2803D0-1.
;dz = 0.001 ;; This is just magic, isn't it?
dz = jhusdss_concolve_dz() ;; change to 0.0005, 02/01/2012
nz = floor((zmax-zmin)/dz)

zgrid = findgen(nz)*dz+zmin

for ispec=0L, nspec-1L do begin
    spec = spec0[ispec]

    outfile = jhusdss_detect_absorbers_convolve_filename(spec.plate, spec.fiber)
    subpath = path+'/'+string(spec.plate, format='(i4.4)')
    if (jhusdss_direxist(subpath) eq 0) then spawn, 'mkdir -p '+subpath
    filename = subpath+'/'+outfile

    if (file_test(filename) and (not keyword_set(overwrite))) then begin
       splog, filename+' file exists, not overwriting ...'
       continue
    endif

    ;; initialization
    outstr = {ra:spec.ra, dec:spec.dec, plate:spec.plate, fiber:spec.fiber, mjd:spec.mjd, $
              z:spec.z, lines:lines.name, zgrid:zgrid, $
              signal:fltarr(n_elements(lines), nz), ivar:fltarr(n_elements(lines), nz), $
              ewall:fltarr(n_elements(lines), nz),  err_ewall:fltarr(n_elements(lines), nz)}

    ;; izmin and izmax for this spec

    ;; start 3000 km/s away from spec.z
    izmax = value_locate(zgrid, spec.z-0.01) < (n_elements(zgrid)-1)

    ;; start from 1220 AA. -- 01/25/2012, Guangtun
    tmp_wavemin = 1220.*(1.+spec.z)
    tmp_zmin = tmp_wavemin/2796.D0-1.
    izmin = value_locate(zgrid, tmp_zmin) > 0
    if ((izmax-izmin) lt 5) then continue

    wave = spec.wave*(1.+spec.z)

    ;; kill linemask -- 01/25/2012, Guangtun
    ;; mask out emissing lines (only CIV, and MgII as above)
;   linemask = jhusdss_linemask(wave, emlines, spec.z)
;   ivar = spec.ivar*(~linemask)

    ivar = spec.ivar
    iwave = where(wave gt (1.+zgrid[izmin])*2796.D0  and wave le (1.+zgrid[izmax])*2803.D0, nwave)
    wave = wave[iwave]
    residual = spec.residual[iwave]
    ivar = ivar[iwave]*(spec.nmf_continuum[iwave]*spec.med_continuum[iwave])^2
    
    ;; wavelength grid size
    dwave = jhusdss_dwave(wave)

    for iz=izmin, izmax do begin
        outstr.signal[*,iz] = jhusdss_detect_absorbers_convolve(wave, dwave, residual, ivar, zgrid[iz], lines, ew=ew, errew=errew, newivar=tmp_newivar)
        outstr.ivar[*, iz] = tmp_newivar
        outstr.ewall[*, iz] = ew
        outstr.err_ewall[*, iz] = errew
    endfor

    mwrfits, outstr, filename, /create
endfor

end
