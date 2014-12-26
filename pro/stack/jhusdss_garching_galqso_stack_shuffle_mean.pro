;pro jhusdss_lowz_caii_new

Coarse = 0b
DoWeight = 0b
DoNaI = 0b
DoScale = 0b
Scale_factor = 0.2 ;; M* \propto M(halo)^1 \propto r(halo)^3 and r(halo) \propto M*^(1/3); Scale to 10^10.1 Msun
fid_mass = 10.3 ;; M* \propto M(halo)^1 \propto r(halo)^3 and r(halo) \propto M*^(1/3); Scale to 10^10.1 Msun
read,'Coarse? [1=yes, 0=no]: ', Coarse
read,'DoWeight? [1=yes, 0=no]: ', DoWeight
read,'DoNaI? [1=yes, 0=no]: ', DoNaI
read,'DoScale? [1=yes, 0=no]: ', DoScale

nmfver = 106
overwrite = 1b
;read,'BOSS? [1=yes, 0=no]: ',BOSS

line_wave = [3934.79, 3969.59]
xra = [3800., 4300.]

fix_separation = long((line_wave[1]-line_wave[0])/(line_wave[1]+line_wave[0])/alog(10.)*2E+4)

stackpath = jhusdss_get_path(/nmfqso)+'/'+string(nmfver, format='(I3.3)')+'/Stack/'
infile0 = stackpath + jhusdss_garching_galqso_stack_filename(nmfver, boss=boss)

nshuffle = 12
for ishuffle=1L, nshuffle do begin
    infile = repstr(infile0, '.fits', '_shuffle_'+string(ishuffle, format='(i2.2)')+'.fits')
    comp = mrdfits(infile, 1)
    nrp = n_elements(comp)
    if ishuffle eq 1 then begin
       wave = comp[0].wave
       nwave = n_elements(wave)
       fluxmedian = fltarr(nwave, nrp, nshuffle)
       fluxmean = fltarr(nwave, nrp, nshuffle)
       fluxgeomean = fltarr(nwave, nrp, nshuffle)
    endif
    fluxmedian[*,*,ishuffle-1] = comp.fluxmedian
    fluxmean[*,*,ishuffle-1] = comp.fluxmean
    fluxgeomean[*,*,ishuffle-1] = comp.fluxgeomean
endfor

outfile = repstr(infile, '.fits', '_shuffle_mean.fits')

if (file_test(outfile) and ~overwrite) then begin
   splog, 'File already exists, use /overwrite to overwrite'
   return
endif else begin
   splog, 'Will write into this file: '
   print, outfile
endelse

outstr = replicate({wave:wave, fluxmedian:fltarr(nwave), fluxmean:fltarr(nwave), fluxgeomean:fltarr(nwave)}, nrp)

outstr.fluxmedian = median(fluxmedian, dimension=3)
outstr.fluxmean = median(fluxmean, dimension=3)
outstr.fluxgeomean = median(fluxgeomean, dimension=3)

;for irp = 0L, nrp-1L do begin
;endfor
mwrfits, outstr, outfile, /create

end
