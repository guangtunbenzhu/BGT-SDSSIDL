;+
;; Need to rm -rf MC_nmfver every time we run
;; rm -rf /home/gz323/DATA/SDSS/QSO/NMF//MC_107/*.fits 
;-
pro jhusdss_real_montecarlo_all, nmfver, overwrite=overwrite

if (n_elements(nmfver) eq 0) then $
   message, 'nmfver required'

;; Let's start with 1000
nabs = 10000L
strtmp = {rew_mgii_2796:0., rew_mgii_2803:0., vdisp:0., zabs:0.}
mc_absorbers = replicate(strtmp, nabs)

;; qso pool
qsos = jhusdss_absorber_readin(nmfver, /byqso)
imc_qso = floor(randomu(seed, nabs)*n_elements(qsos))
mc_qsos = qsos[imc_qso]

;; absorber pool
absorbers = jhusdss_montecarlo_absorbers_pool(nmfver)

;; log10w
log10w_min = alog10(0.2)
log10w_max = alog10(6.0)
log10w = randomu(seed, nabs)*(log10w_max-log10w_min)+log10w_min
mc_absorbers.rew_mgii_2796 = 10.^log10w
for i=0L, nabs-1L do begin
    tmp = min(abs(mc_absorbers[i].rew_mgii_2796 - absorbers.rew_mgii_2796), imin)
    mc_absorbers[i].rew_mgii_2803 = absorbers[imin].rew_mgii_2803*mc_absorbers[i].rew_mgii_2796/absorbers[imin].rew_mgii_2796
    mc_absorbers[i].vdisp = absorbers[imin].vdisp
endfor

;; z
z_min = 0.36
z_max = 2.3
z = randomu(seed, nabs)*(z_max-z_min)+z_min
mc_absorbers.zabs = z

out_path = jhusdss_get_path(/nmfqso)+'/MC_'+string(nmfver, format='(I3.3)')
outfile = out_path+'/MC_Master_absorbers.fits'
if (file_test(outfile) and ~keyword_set(overwrite)) then begin
    message, 'File already exists. Use /overwrite if you want to overwrite it.'
endif else begin
    mwrfits, mc_qsos, outfile, /create
    mwrfits, mc_absorbers, outfile
endelse

jhusdss_real_montecarlo_spec, nmfver, mc_qsos, mc_absorbers, overwrite=overwrite

end
