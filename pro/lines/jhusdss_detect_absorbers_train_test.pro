;+
;-
pro jhusdss_detect_absorbers_train_test, nmfver, overwrite=overwrite 

;; output
path=jhusdss_get_path(/nmfqso)+'/'+string(nmfver, format='(I3.3)')+'/Absorbers'
if (jhusdss_direxist(path) eq 0) then message, "Can't find the directory."
infile = path+'/OnlyAbsorbers_Window_Pitts_'+jhusdss_absorbers_filename(nmfver, /mgii)

;; mgii_red
mgii_red = mrdfits(infile, 1)
mgii_blue = mrdfits(infile, 2)
feii_red = mrdfits(infile, 3)
feii_blue = mrdfits(infile, 4)

ii = where(mgii_blue.med_sdeviation_blue lt 0.05 and mgii_blue.rew_mgii_2803 gt 3., nn)

print, nn
for i=0L, nn-1L do begin
    spec = jhusdss_decompose_loadspec(mgii_blue[ii[i]].plate, mgii_blue[ii[i]].fiber, nmfver)
    djs_plot, spec.wave*(1.+spec.z), smooth(spec.flux, 5), xra=[2000., 3100.]*(1.+mgii_blue[ii[i]].zabs)
    djs_oplot, [2803., 2803.]*(1.+mgii_blue[ii[i]].zabs), !y.crange, color='green'
    djs_oplot, [2796., 2796.]*(1.+mgii_blue[ii[i]].zabs), !y.crange, color='green'
    djs_oplot, [2600., 2600.]*(1.+mgii_blue[ii[i]].zabs), !y.crange, color='green'
    djs_oplot, [2586., 2586.]*(1.+mgii_blue[ii[i]].zabs), !y.crange, color='green'
    djs_oplot, [2383., 2383.]*(1.+mgii_blue[ii[i]].zabs), !y.crange, color='green'
    djs_oplot, [2344., 2344.]*(1.+mgii_blue[ii[i]].zabs), !y.crange, color='green'
    print, mgii_blue[ii[i]].plate, mgii_blue[ii[i]].fiber, mgii_blue[ii[i]].zabs
    print, mgii_blue[ii[i]].med_sdeviation_blue, mgii_blue[ii[i]].med_sdeviation_red
    a = 'a'
    read, a
    if a eq 'q' then stop
endfor
;; mgii_blue
;; qso = mrdfits(infile, 2)
;; abs_mgii_blue = jhusdss_absorber_cat_qso2abs(stat, qso)

;; feii_red
;; qso = mrdfits(infile, 3)
;; abs_feii_red = jhusdss_absorber_cat_qso2abs(stat, qso)

;; feii_blue
;; qso = mrdfits(infile, 4)
;; abs_feii_blue = jhusdss_absorber_cat_qso2abs(stat, qso)

;; help, absorbers

;; mwrfits, abs_mgii_red, outfile, /create
;; mwrfits, abs_mgii_blue, outfile
;; mwrfits, abs_feii_red, outfile
;; mwrfits, abs_feii_blue, outfile

end
