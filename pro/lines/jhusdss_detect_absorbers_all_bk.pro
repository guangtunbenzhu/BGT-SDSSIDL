;+
; Still training
;-
pro jhusdss_detect_absorbers_all, nmfver=nmfver, overwrite=overwrite, qaonly=qaonly, boss=boss

if (n_elements(nmfver) eq 0) then nmfver = jhusdss_get_nmf_version()

;; output
path=jhusdss_get_path(/nmfqso)+'/'+string(nmfver, format='(I3.3)')+'/Absorbers'
if (jhusdss_direxist(path) eq 0) then message, "Can't find the directory."
filename = path+'/'+jhusdss_absorbers_filename(nmfver, /mgii, boss=boss)
if (~keyword_set(qaonly)) then begin
   if (file_test(filename) and ~keyword_set(overwrite)) then begin
      splog, 'File already exists. Use /overwrite if you want to overwrite it.'
      return
   endif else begin
      splog, 'Will write the absorber catalog into this file: '
      print, filename
   endelse
endif

dsswave = jhusdss_sdsswave_minmax()
lines = jhusdss_abslines_all(/train)
emlines = jhusdss_emlines_all(/mask)

donespec = jhusdss_decompose_donespec(nmfver, boss=boss)
nspec = n_elements(donespec)
nabsmax = 4
absorbers = replicate({ra:0D0, dec:0D0, plate:0L, fiber:0L, mjd:0L, zqso:0., $
                       zabs:fltarr(nabsmax), snr:fltarr(n_elements(lines), nabsmax), $
                       ew:fltarr(n_elements(lines), nabsmax), err_ew: fltarr(n_elements(lines), nabsmax), $
                       lines:lines.name}, nspec)
nabs = 0L
nabsall = 0L
for itest=20000L, nspec-1L do begin
   counter, itest+1, nspec
   spec = jhusdss_decompose_loadspec(donespec[itest].plate, donespec[itest].fiber, nmfver, boss=boss)
;  if (spec.z lt 2.) then continue

   absorbers0 = jhusdss_detect_absorbers_engine(spec, lines=lines, emlines=emlines, nabsmax=nabsmax)
   if (absorbers0.zabs[0] gt 0.) then begin
      nabs++
      nabsall++
      absorbers[itest].ra = spec.ra
      absorbers[itest].dec = spec.dec
      absorbers[itest].plate = spec.plate
      absorbers[itest].fiber= spec.fiber
      absorbers[itest].mjd = spec.mjd
      absorbers[itest].zqso  = spec.z
      absorbers[itest].zabs = absorbers0.zabs
      absorbers[itest].snr = absorbers0.snr
      absorbers[itest].ew = absorbers0.ew
      absorbers[itest].err_ew = absorbers0.err_ew
      if (keyword_set(qaonly)) then begin
         for k=1,3 do if (absorbers0.zabs[k] gt 0.) then nabsall++
         jhusdss_detect_absorbers_qaplot, spec, absorbers0, lines
            print, ' Found '+strtrim(string(nabs),2)+' ('+strtrim(string(nabsall),2)+') absorbers in '+$
                   strtrim(string(itest+1),2)+' QSOs.'
            a = 'a'
            print, ' Write out the qaplot (otherwise continue)? Y/y?'
            read, a
            if (a eq 'Y' or a eq 'y') then jhusdss_detect_absorbers_qaplot, spec, absorbers0, lines, $
                path=path, /writeout, boss=boss
            a = 'a'
            print, ' Stop (otherwise continue)? Y/y?'
            read, a
            if (a eq 'Y' or a eq 'y') then stop
;        endif
      endif
   endif

endfor

ii = where(absorbers.zabs[0] gt 0., nn)
if (nn gt 0 and ~keyword_set(qaonly)) then begin
   mwrfits, absorbers[ii], filename, /create
endif else begin
   splog, "I didn't find any absorber"
endelse

end
