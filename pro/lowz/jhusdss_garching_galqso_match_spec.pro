;; obsolete
;; This should be just a test version
pro jhusdss_garching_galqso_match_spec, nmfver, boss=boss, overwrite=overwrite, $
       iprocess=iprocess, nprocess=nprocess, noparal=noparal

if (n_elements(nmfver) eq 0) then message, 'nmfver required'
if (n_elements(nprocess) eq 0) then nprocess = 8
if (n_elements(iprocess) eq 0) then iprocess = 1
if (iprocess gt nprocess) then message, "iprocess can't be larger than nprocess"

if (keyword_set(boss)) then begin
   lowzpath=jhusdss_get_path(/nmfqso)+'/'+string(nmfver, format='(I3.3)')+'/Lowz_BOSS'
endif else begin
   lowzpath=jhusdss_get_path(/nmfqso)+'/'+string(nmfver, format='(I3.3)')+'/Lowz'
endelse

if (~keyword_set(noparal)) then $
   outfile = lowzpath+'/Pro_'+string(iprocess, format='(i1.1)')+'_'+jhusdss_lowz_spec_filename(nmfver) $
else $
   outfile = lowzpath+'/'+jhusdss_lowz_spec_filename(nmfver)

if (file_test(outfile) and ~keyword_set(overwrite)) then begin
   splog, 'File already exists. Use /overwrite if you want to overwrite it.'
   return
endif else begin
   splog, 'Will write the composite into this file: '
   print, outfile
endelse

garching_path = jhusdss_get_path(/garching)
garching_file = garching_path+'/'+'gal_info_dr7_v5_2.fit.gz'
gal = mrdfits(garching_file, 1)
sfr_file = garching_path+'/'+'gal_totsfr_dr7_v5_2.fits.gz'
sfr = mrdfits(sfr_file, 1)

;; read in SDSS qso catalog, need ra, dec
qso_path = jhusdss_get_path(/qso)
if (keyword_set(boss)) then begin
   qso_file = qso_path+'/'+jhusdss_boss_qsofile()
endif else begin
   qso_file = qso_path+'/'+jhusdss_dr7_qsofile()
endelse
qso = mrdfits(qso_file, 1)

match = jhusdss_galqso_match_readin(boss=boss)
nmatch = n_elements(match)

;; parallelize
if (~keyword_set(noparal)) then begin
   ibegin = (iprocess-1)*nmatch/nprocess
   iend = iprocess*nmatch/nprocess-1
   if (iprocess eq nprocess) then iend = nmatch-1
   sub_nmatch = iend-ibegin+1
   print, ibegin, iend, sub_nmatch
endif else begin
   ibegin = 0L
   iend = nmatch-1L
endelse

minwave = 3850D0
maxwave = 4050D0
loglam = jhusdss_get_loglam(minwave=minwave, maxwave=maxwave)
nwave = n_elements(loglam)
allflux = fltarr(sub_nmatch, nwave)
allivar = fltarr(sub_nmatch, nwave)

for i=ibegin, iend do begin
    j = i-ibegin
    counter, i+1, nmatch
    plate = qso[match[i].index_qso].plate
    fiber = qso[match[i].index_qso].fiber
    mjd = qso[match[i].index_qso].mjd
    zqso = qso[match[i].index_qso].z
    zgal = gal[match[i].index_gal].z

    spec = jhusdss_decompose_loadspec(plate, fiber, nmfver, boss=boss, error=error)
    if (error) then continue
    residual = jhusdss_redecompose_spec(spec, zgal, new_med_continuum=new_med_continuum)
    wave = spec.wave*(1.+spec.z)/(1.+zgal)
    ivar = spec.ivar*spec.nmf_continuum^2*new_med_continuum^2
    mask = (ivar le 0.)
    curr_loglam = alog10(wave)

    iuse = where(wave gt (minwave-50D0) and wave lt (maxwave+50D0), nuse)
    if (nuse le 5) then continue
    combine1fiber, curr_loglam[iuse], residual[iuse], ivar[iuse], newloglam=loglam, $
           newflux=finterp, newivar=iinterp, maxiter=0, $
           finalmask=mask[iuse], andmask=maskinterp

    allflux[j, *] = finterp
    allivar[j, *] = iinterp
endfor

mwrfits, allflux, outfile, /create
mwrfits, allivar, outfile
end
