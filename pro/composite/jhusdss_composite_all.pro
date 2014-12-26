;+
; Documentation Needed
;-
pro jhusdss_composite_all, version=version, path=path, infile=infile, zallmin=zallmin, $
       zallmax=zallmax, zstep=zstep, delta_z=delta_z, normwaveoption=normwaveoption, $
       delta_mag=delta_mag

;; the minimum size of the sample to create the composite spectra from
n_training_min = 100L

;; version number
if (n_elements(version) eq 0) then version=001L
splog, 'version: ', version

;; sdss wavelength coverage in the observer frame
sdsswave = jhusdss_sdsswave_minmax()

;; the default quasar catalog
if (n_elements(path) eq 0) then path = getenv('RAW_DATA')+'/SDSS/QSO'
if (n_elements(infile) eq 0) then infile =  'dr7_bh_May09_2011.fits.gz'
filename = path+'/'+infile
objs0 = mrdfits(filename, 1)

;; redshift tag
ztag = jhusdss_get_tags(/ztag)
zindex = tag_indx(objs0[0], ztag)

;; magnitude tag
magtag = jhusdss_get_tags(/magtag)
magindex = tag_indx(objs0[0], magtag)

;; redshift bins
if(n_elements(zallmin) eq 0) then zallmin = 0.3
if(n_elements(zallmax) eq 0) then zallmax = 1.0
if(n_elements(zstep) eq 0) then zstep = 0.1
if(n_elements(delta_z) eq 0) then delta_z = 0.3
zbins = jhusdss_make_bins(zallmin, zallmax+1.e-6, zstep, delta=delta_z, nbin=nz)
splog, 'z_min: ', zbins.min
splog, 'z_max: ', zbins.max

;; luminosity (magnitude) bins
if(n_elements(magallmin) eq 0) then magallmin = -28.5
if(n_elements(magallmax) eq 0) then magallmax = -22.5
if(n_elements(magstep) eq 0) then magstep = 0.1
if(n_elements(delta_mag) eq 0) then delta_mag = 0.3
magbins = jhusdss_make_bins(magallmin, magallmax+1.e-6, magstep, delta=delta_mag, nbin=nmag)
splog, 'mag_min: ', magbins.min
splog, 'mag_max: ', magbins.max

;; Make a common wave vector
loglam = jhusdss_get_loglam()
wave = 10.^loglam
nwave = n_elements(wave)
    
;; normalization wavelength range
tempwave = jhusdss_normwave_minmax(option=normwaveoption)
normminwave = tempwave[0]
normmaxwave = tempwave[1]

;; Loop over redshift/luminosity bins
for i=0L, nz-1L do begin
    tmp_zbin = zbins[i]

    for j=0L, nmag-1L do begin
        tmp_magbin = magbins[j]

        ;; create automatic file names
        outfile = jhusdss_composite_basis_name(tmp_zbin, tmp_magbin, normminwave=normminwave)

        ;; find the list of objects in each bin
        iuse = where(objs0.(zindex) gt tmp_zbin.min and objs0.(zindex) le tmp_zbin.max $
                 and objs0.(magindex) gt tmp_magbin.min and objs0.(magindex) le tmp_magbin.max, nuse)
        if (nuse lt n_training_min) then begin
           splog, "The number of quasars is less than ", n_training_min
           splog, "I am skipping this redshift/luminosity bin."
           continue
        endif
        newobjs = objs0[iuse]

        ;; wavelength grid
        iwave = where(loglam gt alog10(sdsswave[0]/(1.+tmp_zbin.max)) $
                  and loglam lt alog10(sdsswave[1]/(1.+tmp_zbin.min)), nwave)
        if (nwave eq 0) then begin
           splog, "No useful wavelength coverage in this redshift/luminosity bin."
           splog, "Check your wavelength grid."
           Continue
        endif
        newloglam = loglam[iwave]
        
        ;; create composite_spectra
        jhusdss_create_composite_spectra, newobjs, newloglam, version, outfile=outfile, $
           normoption=1, normminwave=normminwave, normmaxwave=normmaxwave, $
           /ivarweight, /reject, /qaplot
    endfor
endfor

end
