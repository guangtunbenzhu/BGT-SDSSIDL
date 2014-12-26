pro jhusdss_montecarlo_wmin_gigantic, nmfver, $
       overwrite=overwrite, boss=boss, path=path

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

;; qsos
qsos = jhusdss_absorber_readin(nmfver, /byqso)
nqso = n_elements(qsos)

outfile = path+'/'+jhusdss_montecarlo_wmin_gigantic_filename(nmfver)

for i=0L, nqso-1L do begin
    counter, i+1L, nqso
    ;; prepare for write out
    infile = jhusdss_montecarlo_wmin_filename(qsos[i].plate, qsos[i].fiber)
    subpath = path+'/'+string(qsos[i].plate, format='(i4.4)')
    filename = subpath+'/'+infile
    wmin = mrdfits(filename, 1, status=error, /silent)
    if (i eq 0L) then begin
       zgrid = {zgrid:wmin.zgrid}
       strtmp = {plate:0L, fiber:0L, isitcovered:bytarr(n_elements(wmin.zgrid)), rewmin_mgii_2796:fltarr(n_elements(wmin.zgrid)), rewmin_mgii_2803:fltarr(n_elements(wmin.zgrid))}
       outstr = replicate(strtmp, nqso)
    endif
    if (error ne 0) then continue

    outstr[i].plate = wmin.plate
    outstr[i].fiber = wmin.fiber
    outstr[i].isitcovered = wmin.isitcovered
    outstr[i].rewmin_mgii_2796 = wmin.rewmin_mgii_2796
    outstr[i].rewmin_mgii_2803 = wmin.rewmin_mgii_2803
endfor

mwrfits, zgrid, outfile, /create
mwrfits, outstr, outfile
stop

end
