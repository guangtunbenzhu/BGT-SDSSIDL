;; Signal - EW: EW = 6.-sqrt((4.5-(absorber0.signal_mgii_2796<4.5))*36./4.5)

pro jhusdss_montecarlo_wmin_spec, qsos, nmfver, $
       overwrite=overwrite, boss=boss, path=path, noresidual=noresidual

if (n_elements(nmfver) eq 0) then $
    message, "nmf version required so that you won't overwrite files unintentionally."

;; path
if (~keyword_set(path)) then begin
   if (keyword_set(boss)) then begin
       path=jhusdss_get_path(/nmfqso)+'/'+$
          string(nmfver, format='(I3.3)')+'/Wmin_BOSS'
   endif else begin
       path=jhusdss_get_path(/nmfqso)+'/'+$
          string(nmfver, format='(I3.3)')+'/Wmin'
   endelse
endif

nqso = n_elements(qsos)

;; residual
average_residual_file = jhusdss_get_path(/nmfqso)+'/'+ string(nmfver, format='(I3.3)')+$
                        '/Composite/Residual_Composite_Convolved_'+ string(nmfver, format='(I3.3)')+'.fits'
average_residual =  mrdfits(average_residual_file, 1)
;;

for iqso=0L, nqso-1L do begin

    counter, iqso+1L, nqso
    ;; prepare for write out
    outfile = jhusdss_montecarlo_wmin_filename(qsos[iqso].plate, qsos[iqso].fiber)
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

    ;; mc structure
    if (total(spec.zgrid - average_residual.zgrid) ne 0.) then message, 'Fix this first'
    nz = n_elements(spec.zgrid)
    mc = jhusdss_montecarlo_wmin_blank(nz)
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

    ;; get signal_min
    mc.signalmin_mgii_2796[icovered] = randomn(seed, ncovered)/sqrt(spec.ivar[1,icovered]) + 4.*1./sqrt(spec.ivar[1,icovered]) $
                                        - average_residual.signal[1,icovered]
    mc.signalmin_mgii_2803[icovered] = randomn(seed, ncovered)/sqrt(spec.ivar[0,icovered]) + 2.*1./sqrt(spec.ivar[0,icovered]) $
                                        - average_residual.signal[0,icovered]
    mc.rewmin_mgii_2796[icovered] = jhusdss_montecarlo_signal_wmin(mc.signalmin_mgii_2796[icovered])
    mc.rewmin_mgii_2803[icovered] = jhusdss_montecarlo_signal_wmin(mc.signalmin_mgii_2803[icovered])
    
    ;; write out mc sample
    mwrfits, mc, filename, /create
;   print, filename
endfor

end
