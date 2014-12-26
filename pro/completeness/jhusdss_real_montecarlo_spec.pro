;+
;; Given a random QSO and a random absorber, insert the absorber into the QSO
;; Not for BOSS
;; nqsos = nabsorbers
;-
pro jhusdss_real_montecarlo_spec, nmfver, qsos, absorbers, overwrite=overwrite

in_path = jhusdss_get_path(/spectra)
out_path = jhusdss_get_path(/nmfqso)+'/MC_'+string(nmfver, format='(I3.3)')+'/Flux'

for i=0L, n_elements(qsos)-1L do begin
    counter, i+1, n_elements(qsos)

    outfile = out_path+'/MC_flux_'+string(qsos[i].plate, format='(i4.4)')+'_'+string(qsos[i].fiber, format='(i3.3)')+'_'+string(qsos[i].mjd, format='(i5.5)')+'.fits'
    if (file_test(outfile) and ~keyword_set(overwrite)) then begin
       splog, 'File already exists. Use /overwrite if you want to overwrite it.'
       continue
    endif

    if (file_test(outfile)) then begin
       splog, 'Overwriting existing file!!!'
       tmp_input = mrdfits(outfile, 1, /silent)
       wave = tmp_input.wave
       flux = tmp_input.flux
    endif else begin
    ;;  get wavelength and flux
       readspec, qsos[i].plate, qsos[i].fiber, mjd=qsos[i].mjd, path=in_path, $
           wave=wave, flux=flux, invvar=ivar, andmask=andmask, ormask=ormask
    endelse

    ;; generate the absorber
    ;; line profile fit
    ;; 2803.53, 2796.35, 2600.17, 2586.65, 2382.77, 2344.21,

    ew2796 = absorbers[i].rew_mgii_2796
    ew2803 = absorbers[i].rew_mgii_2803
    ew2600 = absorbers[i].rew_mgii_2796*0.25
    ew2586 = absorbers[i].rew_mgii_2796*0.125
    ew2383 = absorbers[i].rew_mgii_2796*0.25
    ew2344 = absorbers[i].rew_mgii_2796*0.20
    sigma = sqrt(((absorbers[i].vdisp/69.03)^2+1.))*1.D-4*2800.*alog(10.)

    par = [2796.35, 2803.53, ew2796, ew2803, sigma]*(1.+absorbers[i].zabs)
    e_tao_lambda= doublegaussian_line_func(wave, par)
;   e_tao_lambda= ((1.-doublegaussian_line_func(wave, par))>0.)

    par = [2600.17, 2586.65, ew2600, ew2586, sigma]*(1.+absorbers[i].zabs)
    e_tao_lambda= e_tao_lambda+doublegaussian_line_func(wave, par)
;   e_tao_lambda= e_tao_lambda*((1.-doublegaussian_line_func(wave, par))>0.)

    par = [2382.77, 2344.21, ew2383, ew2344, sigma]*(1.+absorbers[i].zabs)
    e_tao_lambda= e_tao_lambda+doublegaussian_line_func(wave, par)
;   e_tao_lambda= e_tao_lambda*((1.-doublegaussian_line_func(wave, par))>0.)
    e_tao_lambda = 1.-e_tao_lambda

    ;; insert absorber
    output = {wave:wave, flux:flux*e_tao_lambda}

    ;; write out
    mwrfits, output, outfile, /create
endfor

end
