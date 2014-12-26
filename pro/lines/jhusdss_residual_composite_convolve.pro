pro  jhusdss_residual_composite_convolve, nmfver, overwrite=overwrite, boss=boss

if (n_elements(nmfver) eq 0) then message, 'nmfver required'

;; output 
if (~keyword_set(boss)) then begin
   path=jhusdss_get_path(/nmfqso)+'/'+string(nmfver, format='(I3.3)')+'/Composite'
endif else begin
   path=jhusdss_get_path(/nmfqso)+'/'+string(nmfver, format='(I3.3)')+'/Composite_BOSS'
endelse
if (jhusdss_direxist(path) eq 0) then message, "Can't find the directory."

if (~keyword_set(boss)) then begin
   outfile = path+'/Residual_Composite_Convolved_'+string(nmfver, format='(i3.3)') + '.fits'
endif else begin
   outfile = path+'/Residual_Composite_Convolved_'+string(nmfver, format='(i3.3)') + '_BOSS.fits'
endelse 

if (~keyword_set(boss)) then begin
   infile = path+'/Residual_Composite_'+string(nmfver, format='(i3.3)') + '.fits'
endif else begin
   infile = path+'/Residual_Composite_'+string(nmfver, format='(i3.3)') + '_BOSS.fits'
endelse 

;; see jhusdss_detect_absorbers_convolve_spec
if (n_elements(lines) eq 0) then lines = jhusdss_abslines_all(/train)
if (n_elements(waveminmax) eq 0) then waveminmax = jhusdss_sdsswave_minmax()

zmin = waveminmax[0]/2796D0-1.
zmax = waveminmax[1]/2803D0-1.
;dz = 0.001 ;; This is just magic, isn't it?
dz = jhusdss_convolve_dz() ;; change to 0.0005, 02/01/2012
nz = floor((zmax-zmin)/dz)

zgrid = findgen(nz)*dz+zmin

;; initialization
outstr = {ra:0D0, dec:0D0, plate:0L, fiber:0L, mjd:0L, $
          z:0., lines:lines.name, zgrid:zgrid, $
          signal:fltarr(n_elements(lines), nz), ivar:fltarr(n_elements(lines), nz), $
          ewall:fltarr(n_elements(lines), nz),  err_ewall:fltarr(n_elements(lines), nz)}

spec = mrdfits(infile, 1)
wave = spec.wave
residual = spec.fluxmedian
ivar = fltarr(n_elements(residual))+1.
dwave = jhusdss_dwave(wave)

for iz=0L, nz-1L do begin
        outstr.signal[*,iz] = jhusdss_detect_absorbers_convolve(wave, dwave, residual, ivar, zgrid[iz], lines, ew=ew, errew=errew, newivar=tmp_newivar)
        outstr.ivar[*, iz] = tmp_newivar
        outstr.ewall[*, iz] = ew
        outstr.err_ewall[*, iz] = errew
endfor

mwrfits, outstr, outfile, /create

end
