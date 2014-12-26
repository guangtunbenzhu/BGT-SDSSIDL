;+
; Documentation needed
; on convolved spectra
;- 

function jhusdss_detect_absorbers_newengine, spec, orispec, lines=lines, nabsmax=nabsmax, nodofinal=nodofinal;, skymask=skymask

;; use first six lines: 2803, 2796, 2600, 2587, 2374, 2344 
if (n_elements(spec) ne 1) then message, "I can only deal with one spectrum."
if (n_elements(lines) eq 0) then lines = jhusdss_train_lines()
if (n_elements(nabsmax) eq 0) then nabsmax = jhusdss_nabsmax()

absorbers_tmp = jhusdss_absorber_blank(nabsmax=nabsmax)
final_tmp = jhusdss_finalpass_blank(nabsmax=nabsmax)
absorbers = struct_addtags(absorbers_tmp, final_tmp)

nz = n_elements(spec.zgrid)
snr = fltarr(7, nz)
for i=0L, 6L do snr[i,*] = spec.signal[i,*]*sqrt(spec.ivar[i,*])


;; mask out CIV and CIII
lines_tmp = lines.wave
;[2803.53,     2796.35,     2600.17,      2586.65,     2382.77,     $
;             2344.21,     1550.78]

for iline=0L, 6L do begin
    jj = where(lines_tmp[iline]*(1.+spec.zgrid) lt 1550.*(1.+spec.z+0.02) $
           and lines_tmp[iline]*(1.+spec.zgrid) gt 1550.*(1.+spec.z-0.06), mm)
    if (mm gt 0) then snr[iline, jj]=0.
;; Kill this mask, deal with it in the end
;   jj = where(lines_tmp[iline]*(1.+spec.zgrid) lt 1909.*(1.+spec.z+0.01) $ 
;          and lines_tmp[iline]*(1.+spec.zgrid) gt 1909.*(1.+spec.z-0.03), mm)
;   if (mm gt 0) then snr[iline, jj]=0.
endfor

;; magic numbers
magic = bytarr(4,nz)
;; Criterion 1: MgII 2803&2796 > 3 and Fe II 2600&2587 Not available
magic[0,*] = (snr[0,*] gt 2. $
         and snr[1,*] gt 4.)

;; Criterion 2: MgII 2803&2796 > 2 and Fe II 2600|2587 > 2
magic[1,*] = (snr[0,*] gt 2. $
         and snr[1,*] gt 4. $
         and (fix(snr[2,*] gt 2.0)+fix(snr[3,*] gt 2.0)+fix(snr[4,*] gt 2.0)+fix(snr[5,*] gt 2.0)) ge 1)

;; Criterion 3: MgII 2803|2796 > 2 and All FeII > 1.
magic[2,*] = ((fix(snr[2,*] gt 2.0)+fix(snr[3,*] gt 2.0)+fix(snr[4,*] gt 2.0)+fix(snr[5,*] gt 2.0)) ge 4)

;; Criterion 4: All FeII > 2.
;magic[3,*] = ((fix(snr[2,*] gt 1.5)+fix(snr[3,*] gt 1.5)+fix(snr[4,*] gt 1.5)+fix(snr[5,*] gt 1.5)) ge 3)

;; All candidate candidates
i_candidate = where(magic[0,*] or magic[1,*] or magic[2,*] or magic[3,*], n_candidate)

if n_candidate gt 0 then begin

;; group continuous redshifts
if n_candidate gt 1 then begin
   z_candidate= spec.zgrid[i_candidate]
   dz_candidate = z_candidate[1:n_candidate-1] - z_candidate[0:n_candidate-2]
   ;; Assuming it's very rare to have two absorbers within 900 km/s
   i_partition = where(dz_candidate gt 0.003, n_partition)
   if n_partition eq 0 then begin
      snrusetmp_candidate = (snr[0:3,i_candidate] < 20.)
      snrusetmp_candidate = (snrusetmp_candidate > 0.)
      snruse_candidate = sqrt(total(snrusetmp_candidate^2, 1))
      snrusetmp_max = max(snruse_candidate, imax)
      new_i_candidate = i_candidate[imax]
;     new_i_candidate = median(i_candidate)
   endif else begin
      new_i_candidate = lindgen(n_partition+1)
      ;; the zeroth
      snrusetmp_candidate = (snr[0:3,i_candidate[0:i_partition[0]]] < 20.)
      snrusetmp_candidate = (snrusetmp_candidate > 0.)
      snruse_candidate = sqrt(total(snrusetmp_candidate^2, 1))
      snrusetmp_max = max(snruse_candidate, imax)
      new_i_candidate[0] = i_candidate[0+imax]
;     new_i_candidate[0] = median(i_candidate[0:i_partition[0]])

      if n_partition gt 1 then begin
         for j=1L, n_partition-1L do begin
             snrusetmp_candidate = (snr[0:3,i_candidate[i_partition[j-1]+1:i_partition[j]]] < 20.)
             snrusetmp_candidate = (snrusetmp_candidate > 0.)
             snruse_candidate = sqrt(total(snrusetmp_candidate^2, 1))
             snrusetmp_max = max(snruse_candidate, imax)
             new_i_candidate[j] = i_candidate[i_partition[j-1]+1+imax]
;            new_i_candidate[j] = median(i_candidate[i_partition[j-1]+1:i_partition[j]])
         endfor
      endif

      snrusetmp_candidate = (snr[0:3,i_candidate[i_partition[n_partition-1]+1:n_candidate-1]] < 20.)
      snrusetmp_candidate = (snrusetmp_candidate > 0.)
      snruse_candidate = sqrt(total(snrusetmp_candidate^2, 1))
      snrusetmp_max = max(snruse_candidate, imax)
      new_i_candidate[n_partition] = i_candidate[i_partition[n_partition-1]+1+imax]
;     new_i_candidate[n_partition] = median(i_candidate[i_partition[n_partition-1]+1:n_candidate-1])
   endelse
endif else begin
   new_i_candidate = i_candidate
endelse
new_n_candidate = n_elements(new_i_candidate)

;; sort with SNR of MgII
snrusetmp = (snr[0:3,new_i_candidate] < 20.)
snrusetmp = (snrusetmp > 0.)
snruse = sqrt(total(snrusetmp^2, 1))
isort = reverse(bsort(snruse))

new_i_candidate = new_i_candidate[isort]
snruse2tmp = (snr[2:5,new_i_candidate] < 20.)
snruse2tmp = (snruse2tmp > 0.)
snruse2 = sqrt(total(snruse2tmp^2, 1))

iuse = bytarr(new_n_candidate) + 1b
new_z_candidate = spec.zgrid[new_i_candidate]

;print, new_z_candidate

;; mask out 1550. CIV

;; first pass
for i=0L, new_n_candidate-1L do begin
    tmpz = new_z_candidate[i]
;   print, tmpz

    ;; 0.02 beyond 1250. is useless
    if (1220.*(1.+orispec.z+0.02) gt 2803.53*(tmpz+1.)) then iuse[i] = 0b
    ;; 0.02 near qso is likely associated with qsos, useful for other studies, but not ours
    if (orispec.z+0.04 lt tmpz) then iuse[i] = 0b

;   print, new_z_candidate[i], magic[*,new_i_candidate[i]], iuse[i] 

    ;; double Gaussian fit to eliminate bad ones
    ;; mostly affect regsion at 1910\pm0.01 and <1550+0.02
    if ((iuse[i]) and (magic[0,new_i_candidate[i]] eq 1b or magic[1, new_i_candidate[i]] eq 1b) $
      and (magic[2, new_i_candidate[i]] ne 1b)) then begin
       lambda = orispec.wave*(1.+orispec.z)
       ilambda = where((lambda gt 2796.35D0*(1.+tmpz)*(1.D0-3.D-3)) and $
                       (lambda lt 2803.53D0*(1.+tmpz)*(1.D0+3.D-3)) and $
                       (orispec.ivar gt 0.), nlambda)
       if nlambda lt 17 then begin
          iuse[i]=0b
       endif else begin
          in_lambda = lambda[ilambda]
          in_flux = 1.-orispec.residual[ilambda]
          in_ivar = orispec.ivar[ilambda]
          in_center = [2796.35, 2803.53]*(1.+tmpz)

          tmpflux = (((total(in_flux)*mean(in_center)*1.D-4*alog(10.D0)) > 0.01) < 12.)
          in_line_flux = [tmpflux/2.3*1.3, tmpflux/2.3*1.]

;         in_line_flux = [((spec.ewall[1, new_i_candidate[i]]>spec.ewall[0, new_i_candidate[i]])>0.01), $
;                     (((spec.ewall[1, new_i_candidate[i]]>spec.ewall[0, new_i_candidate[i]])/1.3)>0.01)]*(1.+tmpz)

          jhusdss_doublegaussian_fit, in_lambda, in_flux, in_ivar, in_center, in_line_flux, $
             center=center, err_center=err_center, lflux=lflux, err_lflux=err_lflux, $
             sigma=sigma, err_sigma=err_sigma

;      if ((1550.*(1.+orispec.z)/2803.+0.02) gt (tmpz+1.)) $
;           or (abs(1910.*(1.+orispec.z)/2803.-(1.+tmpz)) le 0.02) then begin
;         dcenter_max = 2.0D-4
;         sigma_max = 5.0D-4
;      endif else begin
;         dcenter_max = 5.0D-4
;         sigma_max = 5.0D-4
;       endelse

          dcenter_max = 3.0D-4
          real_sigma = sqrt(((sigma/alog(10.)/2800./(1.+tmpz)*1.D+4)^2-1.)>1.E-6)
          real_ew =lflux[0]/(1.+tmpz)
          real_dcenter = abs((center[1]-center[0])/2800./(1.+tmpz) - 2.56426D-3)
          real_ratio = lflux[0]/lflux[1]

          if ((real_ratio gt 4. or real_ratio lt 0.1) $
           or (real_dcenter gt dcenter_max)) then begin
;          or (real_sigma gt real_ew*40./69.+1.5)) then begin
;         if (abs((center[1]-center[0])/2800./(1.+tmpz) - 2.56426D-3) gt dcenter_max $
;            or sigma/alog(10.)/2800./(1.+tmpz) gt sigma_max) then begin
             iuse[i] = 0b
          endif
       endelse
    endif
endfor

;; second pass
;; eliminate line confusions.
;; stupid algorithm, but for new_n_candidate < ~4, it should be Okay
;print, iuse
if (new_n_candidate gt 1) then begin
for i=0L, new_n_candidate-2L do begin
    if iuse[i] then begin
       for j=i+1L, new_n_candidate-1L do begin
           if iuse[j] then begin
              delta_z = new_z_candidate[i] - new_z_candidate[j]
;             print, new_z_candidate[i], new_z_candidate[j], delta_z
              ;; 2803 & 2796 Mg II confusion
              if ((abs(delta_z - (new_z_candidate[i]+1.-2803.53*(1.+new_z_candidate[i])/2796.35)) le 0.004) $
               or (abs(delta_z - (-new_z_candidate[j]-1.+2803.53*(1.+new_z_candidate[j])/2796.35)) le 0.004)) $
                  then iuse[j]=0b
              ;; 2853 Mg I confusion
              if ((abs(delta_z - (new_z_candidate[i]+1.-2852.96*(1.+new_z_candidate[i])/2796.35)) le 0.004) $
                 and (magic[2, new_i_candidate[j]] ne 1b))then iuse[j]=0b
              ;; Fe II & C IV confusion
              possible_dz = [new_z_candidate[i]+1.-lines[2:6].wave*(1.+new_z_candidate[i])/2803.35, $ 
                            -new_z_candidate[j]-1.+lines[2:6].wave*(1.+new_z_candidate[j])/2803.35, $
                             new_z_candidate[i]+1.-lines[2:6].wave*(1.+new_z_candidate[i])/2796.53, $ 
                            -new_z_candidate[j]-1.+lines[2:6].wave*(1.+new_z_candidate[j])/2796.53]
              if (min(abs(delta_z - possible_dz)) le 0.003) then begin
                 if (new_z_candidate[i] gt new_z_candidate[j]) then begin
                    if (magic[2, new_i_candidate[j]] ne 1b) then iuse[j] = 0b
;                   iuse[j] = 0b
                 endif else begin
                    if (magic[1, new_i_candidate[j]] ne 1b and $
                        magic[0, new_i_candidate[j]] ne 1b) then iuse[j]=0b
                 endelse
              endif
           endif
       endfor
    endif
endfor
endif

istillalive = where(iuse, nstillalive)
if nstillalive gt 0 then begin
   ifinal = new_i_candidate[istillalive]
   zfinal = new_z_candidate[istillalive]

   nmax = n_elements(ifinal) < nabsmax

   absorbers.nabs = nmax
   absorbers.snr[*,0:nmax-1] = snr[*,ifinal[0:nmax-1]]
   absorbers.signal[*,0:nmax-1] = spec.signal[*,ifinal[0:nmax-1]]
   absorbers.ivar[*,0:nmax-1] = spec.ivar[*,ifinal[0:nmax-1]]
   absorbers.magic[*,0:nmax-1] = magic[*,ifinal[0:nmax-1]]
   absorbers.zabs_firstpass[0:nmax-1] = zfinal[0:nmax-1]
   absorbers.ew_firstpass[*,0:nmax-1] = spec.ewall[*,ifinal[0:nmax-1]]
   absorbers.err_ew_firstpass[*,0:nmax-1] = spec.err_ewall[*,ifinal[0:nmax-1]]

   if (~keyword_set(nodofinal)) then begin
      final0 = jhusdss_detect_absorbers_finalpass(orispec, zfinal[0:nmax-1], nabsmax=nabsmax)
      absorbers.zabs = final0.zabs
      absorbers.err_zabs = final0.err_zabs
      absorbers.zabs_all = final0.zabs_all
      absorbers.ew = final0.ew
      absorbers.err_ew = final0.err_ew
      absorbers.vdisp = final0.vdisp
      absorbers.err_vdisp = final0.err_vdisp
      absorbers.vdisp_all = final0.vdisp_all
      absorbers.err_vdisp_all = final0.err_vdisp_all
   endif
endif

endif

return, absorbers

end
