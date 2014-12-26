;; convolve_all should be of similar structure!
pro jhusdss_montecarlo_spec, qsos, absorbers, nmfver, $
       overwrite=overwrite, boss=boss, path=path, noresidual=noresidual

if (n_elements(nmfver) eq 0) then $
    message, "nmf version required so that you won't overwrite files unintentionally."

;; path
if (~keyword_set(path)) then begin
   if (keyword_set(boss)) then begin
       path=jhusdss_get_path(/nmfqso)+'/'+$
          string(nmfver, format='(I3.3)')+'/MonteCarlo_BOSS'
   endif else begin
       path=jhusdss_get_path(/nmfqso)+'/'+$
          string(nmfver, format='(I3.3)')+'/MonteCarlo'
   endelse
endif

nqso = n_elements(qsos)
nparent = n_elements(absorbers)

;; residual
average_residual_file = jhusdss_get_path(/nmfqso)+'/'+ string(nmfver, format='(I3.3)')+$
                        '/Composite/Residual_Composite_Convolved_'+ string(nmfver, format='(I3.3)')+'.fits'
average_residual =  mrdfits(average_residual_file, 1)
;;

for iqso=0L, nqso-1L do begin

    counter, iqso+1L, nqso
    ;; prepare for write out
    outfile = jhusdss_montecarlo_filename(qsos[iqso].plate, qsos[iqso].fiber)
    subpath = path+'/'+string(qsos[iqso].plate, format='(i4.4)')
    if (jhusdss_direxist(subpath) eq 0) then spawn, 'mkdir -p '+subpath
    filename = subpath+'/'+outfile

    if (file_test(filename) and (not keyword_set(overwrite))) then begin
       splog, filename+' file exists, not overwriting ...'
       continue
    endif
 
    ;; read in convolved spectra
    spec = jhusdss_convolve_loadspec(qsos[iqso].plate, qsos[iqso].fiber, nmfver, error=error)
    if error then begin
       splog, "Can't find the convolved spectrum."
       continue
    endif

;   orispec = jhusdss_decompose_loadspec(qsos[iqso].plate, qsos[iqso].fiber, nmfver)

    ;; mc structure
    if (total(spec.zgrid - average_residual.zgrid) ne 0.) then message, 'Fix this first'
    nz = n_elements(spec.zgrid)
    mc = jhusdss_montecarlo_blank(nz)
    mc.plate = spec.plate
    mc.fiber = spec.fiber
    mc.zgrid = spec.zgrid

    ;; covered redshift
    icovered = where(spec.ivar[1,*] gt 0. and spec.ivar[0,*] gt 0., ncovered)
    if ncovered eq 0b then begin
       splog, 'No useful pixels.'
       continue 
    endif

    mc.isitcovered[icovered] = 1b

    ;; assign ew and signal
    iran = floor(randomu(seed, ncovered)*nparent)
    mc.rew_mgii_2796[icovered] = absorbers[iran].rew_mgii_2796
    mc.rew_mgii_2803[icovered] = absorbers[iran].rew_mgii_2803
;   mc.signal_mgii_2796[icovered] = absorbers[iran].signal_mgii_2796 + randomn(seed, ncovered)/sqrt(spec.ivar[1,icovered])
;   mc.signal_mgii_2803[icovered] = absorbers[iran].signal_mgii_2803 + randomn(seed, ncovered)/sqrt(spec.ivar[0,icovered])
    mc.signal_mgii_2796[icovered] = absorbers[iran].signal_mgii_2796 + randomn(seed, ncovered)/sqrt(spec.ivar[1,icovered]) $
                                  + average_residual.signal[1,icovered]
    mc.signal_mgii_2803[icovered] = absorbers[iran].signal_mgii_2803 + randomn(seed, ncovered)/sqrt(spec.ivar[0,icovered]) $
                                  + average_residual.signal[0,icovered]

    ;; detected?
    idetected = where(mc.signal_mgii_2796[icovered]*sqrt(spec.ivar[1,icovered]) gt 4. and mc.signal_mgii_2803[icovered]*sqrt(spec.ivar[0,icovered]) gt 2., ndetected)
    if ndetected eq 0b then begin
       splog, 'No detections.'
       continue 
    endif
    mc.isitdetected[icovered[idetected]] = 1b

    ;; write out mc sample
    mwrfits, mc, filename, /create
;   print, filename
endfor

end
