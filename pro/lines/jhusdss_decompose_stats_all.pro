;+
;-
;; To create a qso parameter file. We will trim those with 'bad' statistics: 
;; weird eigenvalues, large skewness, large width, low s/n 
pro jhusdss_decompose_stats_all, nmfver, boss=boss, dr12=dr12

if (n_elements(nmfver) eq 0) then message, 'nmfver required'

;; output
if (keyword_set(dr12)) then begin
   path=jhusdss_get_path(/nmfqso)+'/'+string(nmfver, format='(I3.3)')+'/Decompose_DR12'
endif else begin
   if (keyword_set(boss)) then begin
       path=jhusdss_get_path(/nmfqso)+'/'+string(nmfver, format='(I3.3)')+'/Decompose_BOSS'
   endif else begin
       path=jhusdss_get_path(/nmfqso)+'/'+string(nmfver, format='(I3.3)')+'/Decompose'
   endelse
endelse

if (jhusdss_direxist(path) eq 0) then message, "Can't find the directory."
outfile = path+'/'+jhusdss_stat_filename(nmfver, boss=boss, dr12=dr12)

qsopath = jhusdss_get_path(/qso)

if (keyword_set(dr12)) then begin
   infile = jhusdss_dr12_qsofile()
endif else begin
   if (keyword_set(boss)) then begin
      infile = jhusdss_boss_qsofile()
   endif else begin
      infile =  jhusdss_dr7_qsofile()
   endelse
endelse

qsofile = qsopath+'/'+infile
splog, 'reading '+qsofile
qso = mrdfits(qsofile, 1)
nqso = n_elements(qso)

; 11-29-2014, add eigen values and basisfile
str_tmp = {ra:0d0, dec:0d0, plate:0L, fiber:0L, mjd:0L, zqso:0., err_zqso:0., $
           spec_snr_median:0., isitdecomposed:1b, isitstated_red:0b, isitstated_blue:0b, isitconvolved:1b, $
           med_mean_red:-999., med_sdeviation_red:-999., med_skewness_red:-999., $
           med_mean_blue:-999., med_sdeviation_blue:-999., med_skewness_blue:-999., $
           eigen_values: fltarr(12), basisfile:''}

outstats = replicate(str_tmp, nqso)

for i=0L, nqso-1L do begin

   counter, i+1, nqso

   ;; load convolved spectra
   spec = jhusdss_decompose_loadspec(qso[i].plate, qso[i].fiber, nmfver, boss=boss, dr12=dr12, error=error)
   if error then begin
      splog, "Can't find the decomposed spectrum."
      outstats[i].isitdecomposed = 0b
      continue
   endif

   ;; from QSO
   outstats[i].ra = spec.ra
   outstats[i].dec = spec.dec
   outstats[i].plate = spec.plate
   outstats[i].fiber= spec.fiber
   outstats[i].mjd = spec.mjd
   outstats[i].zqso  = spec.z
   if ((~keyword_set(boss)) and (~keyword_set(dr12))) then outstats[i].err_zqso  = qso[i].zerr

   ii = where(spec.ivar gt 0., nn)
   if (nn gt 0.) then outstats[i].spec_snr_median = median(spec.flux[ii]*sqrt(spec.ivar[ii]))

   ;; red window
   ii = where(spec.ivar gt 0. and spec.wave*(1.+spec.z) gt 1550.*(1.+spec.z+0.02), nn)
   if (nn gt 10) then begin
      outstats[i].isitstated_red = 1b

      ;; medium moments
      residual = spec.med_continuum[ii]-1.
      moms = moment(residual, /nan)
      outstats[i].med_mean_red = moms[0]
      outstats[i].med_sdeviation_red = sqrt(moms[1])
      outstats[i].med_skewness_red = moms[2]
   endif

   ;; blue window
   ii = where(spec.ivar gt 0. and spec.wave*(1.+spec.z) gt 1250.*(1.+spec.z+0.02) and spec.wave*(1.+spec.z) lt 1550.*(1.+spec.z-0.06), nn)
   if (nn gt 10) then begin
      outstats[i].isitstated_blue = 1b

      ;; medium moments
      residual = spec.med_continuum[ii]-1.
      moms = moment(residual, /nan)
      outstats[i].med_mean_blue= moms[0]
      outstats[i].med_sdeviation_blue = sqrt(moms[1])
      outstats[i].med_skewness_blue = moms[2]
   endif

   convolved_spec = jhusdss_convolve_loadspec(qso[i].plate, qso[i].fiber, nmfver, boss=boss, dr12=dr12, error=error)
   if error then begin
      splog, "Can't find the convolved spectrum."
      outstats[i].isitconvolved= 0b
      continue
   endif

   outstats[i].basisfile = spec.basisfile
   outstats[i].eigen_values = spec.eigen_values

endfor

mwrfits, outstats, outfile, /create

end
