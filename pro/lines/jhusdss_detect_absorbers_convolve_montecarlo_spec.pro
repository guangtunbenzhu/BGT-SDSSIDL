pro jhusdss_detect_absorbers_convolve_montecarlo_spec, objs0, nmfver, $
       lines=lines, emlines=emlines, $
       overwrite=overwrite, boss=boss

nspec = n_elements(objs0)
if (n_elements(nmfver) eq 0) then $
    message, "nmf version required so that you won't overwrite files unintentionally."

if (n_elements(lines) eq 0) then lines = jhusdss_abslines_all(/train)
if (n_elements(emlines) eq 0) then emlines = jhusdss_emlines_all(/mask)
if (n_elements(waveminmax) eq 0) then waveminmax = jhusdss_sdsswave_minmax(boss=boss)

;; path
if (~keyword_set(path)) then begin
   if (keyword_set(boss)) then begin
       path=jhusdss_get_path(/nmfqso)+'/MC_'+$
          string(nmfver, format='(I3.3)')+'/Convolve_BOSS'
   endif else begin
       path=jhusdss_get_path(/nmfqso)+'/MC_'+$
          string(nmfver, format='(I3.3)')+'/Convolve'
   endelse
endif

if (jhusdss_direxist(path) eq 0) then spawn, 'mkdir -p '+path

;; Make a  common redshift grid
zmin = waveminmax[0]/2796D0-1.
zmax = waveminmax[1]/2803D0-1.
;dz = 0.001 ;; This is just magic, isn't it?
dz = jhusdss_convolve_dz() ;; change to 0.0005, 02/01/2012
nz = floor((zmax-zmin)/dz)

zgrid = findgen(nz)*dz+zmin

for ispec=0L, nspec-1L do begin
    counter, ispec+1, nspec
    obj = objs0[ispec]

    outfile = jhusdss_detect_absorbers_convolve_filename(obj.plate, obj.fiber)
    subpath = path+'/'+string(obj.plate, format='(i4.4)')
    if (jhusdss_direxist(subpath) eq 0) then spawn, 'mkdir -p '+subpath
    filename = subpath+'/'+outfile

    if (file_test(filename) and (not keyword_set(overwrite))) then begin
       splog, filename+' file exists, not overwriting ...'
       continue
    endif

    ;; read in the spectrum:
    spec = jhusdss_decompose_montecarlo_loadspec(obj.plate, obj.fiber, nmfver, boss=boss, error=error)
    if (error eq 1b) then begin
       splog, "Decomposed spectrum not found"
       continue
    endif
    ;; initialization
    outstr = {ra:spec.ra, dec:spec.dec, plate:spec.plate, fiber:spec.fiber, mjd:spec.mjd, $
              z:spec.z, lines:lines.name, zgrid:zgrid, $
              signal:fltarr(n_elements(lines), nz), ivar:fltarr(n_elements(lines), nz), $
              ewall:fltarr(n_elements(lines), nz),  err_ewall:fltarr(n_elements(lines), nz)}

    ;; izmin and izmax for this spec

    ;; start 3000 km/s away from spec.z
    ;; start 12000 km/s REDSHIFTED from spec.z to include quasar-associated absorbers -- 03/20/2012
    izmax = value_locate(zgrid, spec.z+0.04) < (n_elements(zgrid)-1)

    ;; start from 1200 AA. -- 03/20/2012, Guangtun
    tmp_wavemin = 1200.*(1.+spec.z)
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
