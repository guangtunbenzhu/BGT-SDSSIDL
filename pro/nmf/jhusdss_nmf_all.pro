;+
; Documentation needed!
; The wrapper
; No need to use so many bins eventually
; we can normalize the spectra using mean flux within a certain range but then re-normzalize them
; to another certain range using mean conversion factor obtained with quasars with common coverage.
;-
pro jhusdss_nmf_all, version, path=path, infile=infile, zallmin=zallmin, $
       zallmax=zallmax, zstep=zstep, delta_z=delta_z, normwaveoption=normwaveoption, $
       n_training=n_training

;; the minimum size of the sample to create the basis vectors from
n_training_min = 100L

if (n_elements(n_training) eq 0) then n_training = 6000L

;; version number
if (n_elements(version) eq 0) then message, 'nmf version is required'
splog, 'version: ', version

;; sdss wavelength coverage in the observer frame
sdsswave = jhusdss_sdsswave_minmax()

;; the default quasar catalog
if (n_elements(path) eq 0) then path = jhusdss_get_path(/qso)
if (n_elements(infile) eq 0) then infile =  'dr7_bh_May09_2011.fits.gz'
;if (n_elements(parfile) eq 0) then parfile =  'default.par'
filename = path+'/'+infile
objs0 = mrdfits(filename, 1)

;; redshift tag
ztag = jhusdss_get_tags(/ztag)
zindex = tag_indx(objs0[0], ztag)

;; redshift bins
if(n_elements(zallmin) eq 0) then zallmin = 0.3
if(n_elements(zallmax) eq 0) then zallmax = 1.0
if(n_elements(zstep) eq 0) then zstep = 0.1
if(n_elements(delta_z) eq 0) then delta_z = 0.5

zbins = jhusdss_make_bins(zallmin, zallmax+1.e-6, zstep, delta=delta_z, nbin=nz)
splog, 'z_min: ', zbins.min
splog, 'z_max: ', zbins.max

;; Make a common wave vector
loglam = jhusdss_get_loglam()
wave = 10.^loglam
nwave = n_elements(wave)

;; normalization wavelength range
tempwave = jhusdss_normwave_minmax(option=normwaveoption)
normminwave = tempwave[0]
normmaxwave = tempwave[1]

;; Loop over redshift bins
for i=0L, nz-1 do begin
    tmp_zbin = zbins[i]

    ;; file names
    outfile = jhusdss_nmf_basis_name(tmp_zbin.min, tmp_zbin.max, normminwave=normminwave)

    ;; Training Set
    iuse = where(objs0.(zindex) gt tmp_zbin.min and objs0.(zindex) le tmp_zbin.max, nuse)
    ;; the random function will give ~10% repeated index ...
    ;; to-do: create the list of n_training objects iteratively
    if nuse gt n_training then begin
       iran = floor(randomu(seed, n_training)*float(nuse))
       iran = iran[uniq(iran, sort(iran))]
       iuse = iuse[iran]
    endif

    nuse = n_elements(iuse)
    if nuse lt n_training_min then begin
       splog, "The number of usable spectra is smaller than "+strtrim(string(n_training_min),2)+". Not performing nmf."
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

    jhusdss_create_nmf_basis, newobjs, newloglam, version, outfile=outfile, $
       normoption=1, normminwave=normminwave, normmaxwave=normmaxwave, $
       init_basis=init_basis, /reject, /qaplot

endfor

end
