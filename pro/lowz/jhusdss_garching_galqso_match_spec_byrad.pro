;; Obsolete
;; This should be just a test version
pro jhusdss_garching_galqso_match_spec_byrad, nmfver, boss=boss, overwrite=overwrite, $
       minrad=minrad, maxrad=maxrad

if (n_elements(nmfver) eq 0) then message, 'nmfver required'
if (n_elements(minrad) eq 0) then minrad = 0.
if (n_elements(maxrad) eq 0) then maxrad = 0.1
if (minrad ge maxrad) then message, "minrad can't be larger than maxrad"

if (keyword_set(boss)) then begin
   lowzpath=jhusdss_get_path(/nmfqso)+'/'+string(nmfver, format='(I3.3)')+'/Lowz_BOSS'
endif else begin
   lowzpath=jhusdss_get_path(/nmfqso)+'/'+string(nmfver, format='(I3.3)')+'/Lowz'
endelse

outfile = lowzpath+'/Rad'+string(minrad,format='(f4.2)')+'_'+string(maxrad,format='(f4.2)')+'_' $
        + jhusdss_lowz_spec_filename(nmfver)

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
qso = jhusdss_qso_readin(boss=boss)
allspec = jhusdss_read_allqsospec(nmfver, boss=boss)

match = jhusdss_galqso_match_readin(boss=boss)

isub = where(match.rp_mpc gt minrad $
         and match.rp_mpc le maxrad, nmatch)
if nmatch eq 0 then message, "Can't find any pair within the annulus"
match = match[isub]
nmatch = n_elements(match)

minwave = 3800D0
maxwave = 4200D0
loglam = jhusdss_get_loglam(minwave=minwave, maxwave=maxwave)
nwave = n_elements(loglam)

outstr = {wave:10.^loglam, allflux:fltarr(nmatch, nwave), allivar:fltarr(nmatch, nwave)}

nallwave = n_elements(allspec.wave)

for i=0L, nmatch-1L do begin
    counter, i+1, nmatch
    plate = qso[match[i].index_qso].plate
    fiber = qso[match[i].index_qso].fiber
    mjd = qso[match[i].index_qso].mjd
    zqso = qso[match[i].index_qso].z
    zgal = gal[match[i].index_gal].z

    tmpqsoindex = match[i].index_qso
    tmpflux = reform(allspec.flux[tmpqsoindex,*])
    tmpnmfcontinuum = reform(allspec.nmf_continuum[tmpqsoindex,*])
    tmpivar = reform(allspec.ivar[tmpqsoindex,*])

    thisspec = {z:0., wave:allspec.wave, ivar:tmpivar, flux:tmpflux, nmf_continuum:tmpnmfcontinuum}

    residual = jhusdss_redecompose_spec(thisspec, zgal, new_med_continuum=new_med_continuum)

;   spec = jhusdss_decompose_loadspec(plate, fiber, nmfver, boss=boss, error=error)
;   if (error) then continue
;   residual = jhusdss_redecompose_spec(spec, zgal, new_med_continuum=new_med_continuum)
;   wave = spec.wave*(1.+spec.z)/(1.+zgal)
;   ivar = spec.ivar*spec.nmf_continuum^2*new_med_continuum^2
;   residual = allspec.flux[match[i].index_qso,*]

    ivar = allspec.ivar[match[i].index_qso,*]*allspec.nmf_continuum[match[i].index_qso,*]^2*allspec.med_continuum[match[i].index_qso,*]^2
    wave = allspec.wave/(1.+zgal)

    mask = (ivar le 0.)
    curr_loglam = alog10(wave)

    iuse = where(wave gt (minwave-50D0) and wave lt (maxwave+50D0), nuse)
    if (nuse le 5) then continue
    combine1fiber, curr_loglam[iuse], residual[iuse], ivar[iuse], newloglam=loglam, $
           newflux=finterp, newivar=iinterp, maxiter=0, $
           finalmask=mask[iuse], andmask=maskinterp

    outstr.allflux[i, *] = finterp
    outstr.allivar[i, *] = iinterp
endfor

mwrfits, match, outfile, /create
mwrfits, outstr, outfile
end
