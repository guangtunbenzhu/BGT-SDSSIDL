;+
; Documentation needed
;- 

function jhusdss_detect_absorbers_engine, spec0, lines=lines, emlines=emlines, nabsmax=nabsmax

nspec = n_elements(spec0)

sdsswave = jhusdss_sdsswave_minmax()
wavemax = sdsswave[1]

;; This set of absorption lines should not be changed:
;; name = ['MgII_2803', 'MgII_2796', 'FeII_2600',  'FeII_2587', 'FeII_2383', 'FeII_2344', 'CIV_1551']
;; threshold = [2.0, 2.0,  1.0,  1.0,  1.0,  1.0, 2.0 (doesn't matter)]
if (n_elements(lines) eq 0) then lines = jhusdss_abslines_all(/train)
if (n_elements(emlines) eq 0) then emlines = jhusdss_emlines_all(/mask)
if (n_elements(nabsmax) eq 0) then nabsmax = 4
zmaxtmp = sdsswave[1]/2803D0-1.-0.01

;; for each spectra, at most we found 4 absorbers
absorbers = replicate({zabs:fltarr(nabsmax), snr:fltarr(n_elements(lines), nabsmax), $
                       ew:fltarr(n_elements(lines), nabsmax), $
                       err_ew:fltarr(n_elements(lines), nabsmax), lines:lines.name}, nspec)

dz = 0.001 ;; the magic number
for ispec=0L, nspec-1L do begin
   spec = spec0[ispec]

   zmax = (spec.z-0.01) < zmaxtmp
   wavemin = ((1250.*(1.+spec.z)) > sdsswave[0])
   zmin = wavemin/2600D0-1+0.01
   if (zmax lt zmin+0.01) then continue

   nz = floor((zmax-zmin)/dz)
   zz = findgen(nz)*dz+zmin
   ff = fltarr(nz)
   snr = fltarr(n_elements(lines), nz)
   ewall = fltarr(n_elements(lines), nz)
   err_ewall = fltarr(n_elements(lines), nz)
   magic = fltarr(nz)
   magic2 = fltarr(nz)

   wave = spec.wave*(1.+spec.z)

   ;; mask out emission lines
   linemask = jhusdss_linemask(wave, emlines, spec.z)
   ivar = spec.ivar*(~linemask)

   iwave = where(wave gt wavemin and wave le wavemax, nwave)
   wave = wave[iwave]
   residual = spec.residual[iwave]
   ivar = ivar[iwave]*(spec.nmf_continuum[iwave]*spec.med_continuum[iwave])^2
   iuse = where(ivar gt 0., nvar)
   if (nvar le nwave/2*1) then continue

   dwave = jhusdss_dwave(wave)

   for i=0L, nz-1L do begin
       snr[*,i] = jhusdss_detect_absorbers_convolve(wave, dwave, residual, ivar, zz[i], lines, ew=ew, errew=errew)
       ewall[*,i] = ew
       err_ewall[*,i] = errew
       magic[i] = total(snr[0:5,i] ge lines[0:5].threshold)
       magic2[i] = total(snr[0:5,i] ge lines[0:5].threshold*2)
   endfor

   ii = where((magic ge 3 and snr[0,*] gt lines[0].threshold $
          and snr[1,*] gt lines[1].threshold) $
          or  (magic ge 5) or (magic2 ge 4), nn)
   if (nn eq 0) then continue

   ;; get rid of continuous redshifts (Ideally should use snr ...)
   if (nn gt 1) then begin
      dii = ii[1:nn-1]-ii[0:nn-2]
      jj = where(dii eq 1, comp=kk, mm)
      if (mm gt 0) then ii[jj+1] = -1
      jj = where(ii ne -1)
      ii = ii[jj]
   endif

   nn = n_elements(ii)
   zuse = zz[ii]
   snralluse = snr[*,ii]
   ewalluse = ewall[*,ii]
   err_ewalluse = err_ewall[*,ii]
   snruse = sqrt(total(snr[0:3,ii]^2, 1))
   snruse2 = sqrt(total(snr[2:5,ii]^2, 1))
;  snruse2 = sqrt(snr[1,ii]^2+snr[3,ii]^2+snr[5,ii]^2)
   isort = reverse(bsort(snruse))

   ;; resolve line confusion
   ;; Hate this chunk of code!
   istart = 1
   while (n_elements(isort) gt istart) do begin
      iuse = bytarr(n_elements(isort)) + 1b
      ;; don't want continuous redshift
      ;; 3000 km/s Dz=0.01
      if ((abs(zuse[isort[istart]]-lines[1].wave*(1.+zuse[isort[istart-1]])/2803.+1.)) lt 0.01) or $
         ((abs(zuse[isort[istart-1]]-lines[1].wave*(1.+zuse[isort[istart]])/2803.+1.)) lt 0.01) then begin
            iuse[istart]=0b
      endif else begin
         ;; if the second is at lower redshift than the first, then compare snruse1, likely the second is wrong
         if (min(abs(zuse[isort[istart]]-lines[2:6].wave*(1.+zuse[isort[istart-1]])/2803.+1.)) lt 0.01) then begin
            iuse[istart]=0b
         ;; if the first is at lower redshift than the second, then compare snruse2, likely the first is wrong
         endif else begin
            if min(abs(zuse[isort[istart-1]]-lines[2:6].wave*(1.+zuse[isort[istart]])/2803.+1.)) lt 0.01 then begin
               if snruse2[isort[istart-1]] le snruse2[isort[istart]] $
                  then iuse[istart-1]=0b else iuse[istart]=0b
            endif
         endelse
      endelse
      kk = where(iuse)
      isort = isort[kk]
      istart++
   endwhile

   nmax = n_elements(isort) < nabsmax

   absorbers[ispec].zabs[0:nmax-1] = zuse[isort[0:nmax-1]]
   absorbers[ispec].snr[*, 0:nmax-1] = snralluse[*, isort[0:nmax-1]]
   absorbers[ispec].ew[*, 0:nmax-1] = ewalluse[*, isort[0:nmax-1]]
   absorbers[ispec].err_ew[*, 0:nmax-1] = err_ewalluse[*, isort[0:nmax-1]]
endfor

return, absorbers

end
