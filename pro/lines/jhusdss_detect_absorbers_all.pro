;+
;-
pro jhusdss_detect_absorbers_all, nmfver, overwrite=overwrite, iprocess=iprocess, $
    nprocess=nprocess, noparal=noparal, nodofinal=nodofinal, boss=boss, dr12=dr12

if (n_elements(nmfver) eq 0) then nmfver = jhusdss_get_nmf_version()
if (n_elements(nprocess) eq 0) then nprocess = 8
if (n_elements(iprocess) eq 0) then iprocess = 1
if (iprocess gt nprocess) then message, "iprocess can't be larger than nprocess"

;; output
if (keyword_set(dr12)) then begin
   path=jhusdss_get_path(/nmfqso)+'/'+string(nmfver, format='(I3.3)')+'/Absorbers_DR12'
endif else begin
   if (keyword_set(boss)) then begin
      path=jhusdss_get_path(/nmfqso)+'/'+string(nmfver, format='(I3.3)')+'/Absorbers_BOSS'
   endif else begin
      path=jhusdss_get_path(/nmfqso)+'/'+string(nmfver, format='(I3.3)')+'/Absorbers'
   endelse
endelse

if (jhusdss_direxist(path) eq 0) then message, "Can't find the directory."

if (~keyword_set(noparal)) then $
   outfile = path+'/Pro_'+string(iprocess, format='(i1.1)')+'_'+jhusdss_absorbers_filename(nmfver, /mgii, boss=boss, dr12=dr12) $
else $
   outfile = path+'/'+jhusdss_absorbers_filename(nmfver, /mgii, boss=boss, dr12=dr12)

if (file_test(outfile) and ~keyword_set(overwrite)) then begin
   splog, 'File already exists. Use /overwrite if you want to overwrite it.'
   return
endif else begin
   splog, 'Will write the absorber catalog into this file: '
   print, outfile
endelse

lines = jhusdss_train_lines()

;; Yue Shen's catalog
qsopath = jhusdss_get_path(/qso)
if (keyword_set(dr12)) then begin
   infile = jhusdss_dr12_qsofile()
endif else begin
   if (keyword_set(boss)) then begin
      infile = jhusdss_boss_qsofile()
   endif else begin
      infile =  jhusdss_dr7_qsofile()
   endelse
endelse

mgiifile = qsopath+'/'+infile
splog, 'reading '+mgiifile
mgii = mrdfits(mgiifile, 1)
nspec = n_elements(mgii)

nabsmax = jhusdss_nabsmax()

absorbers_tmp = jhusdss_absorber_blank(nabsmax=nabsmax)
final_tmp = jhusdss_finalpass_blank(nabsmax=nabsmax)
absorbers_tmp = struct_addtags(absorbers_tmp, final_tmp)

str_tmp = {spec_snr_median:0., isitconvolved:1b}
absorbers_tmp = struct_addtags(absorbers_tmp, str_tmp)

;; parallelize
if (~keyword_set(noparal)) then begin
   ibegin = (iprocess-1)*nspec/nprocess
   iend = iprocess*nspec/nprocess-1
   if (iprocess eq nprocess) then iend = nspec-1
   sub_nspec = iend-ibegin+1
   print, ibegin, iend, sub_nspec
   absorbers = replicate(absorbers_tmp, sub_nspec)
endif else begin
   ibegin = 0L
   iend = nspec-1L
   absorbers = replicate(absorbers_tmp, nspec)
endelse


for i=ibegin, iend do begin
   j=i-ibegin
   counter, i+1, n_elements(mgii)

   ;; load convolved spectra
   spec = jhusdss_convolve_loadspec(mgii[i].plate, mgii[i].fiber, nmfver, boss=boss, dr12=dr12, error=error)
   if error then begin
      splog, "Can't find the convolved spectrum."
      absorbers[j].isitconvolved = 0b
      continue
   endif
   orispec = jhusdss_decompose_loadspec(mgii[i].plate, mgii[i].fiber, nmfver, boss=boss, dr12=dr12, error=error)

   absorbers0 = jhusdss_detect_absorbers_newengine(spec, orispec, lines=lines, nabsmax=nabsmax, nodofinal=nodofinal)

   ;; from QSO
   absorbers[j].ra = spec.ra
   absorbers[j].dec = spec.dec
   absorbers[j].plate = spec.plate
   absorbers[j].fiber= spec.fiber
   absorbers[j].mjd = spec.mjd
   ;; we use z_hw (spec.z == mgii[iz].z_hw, see jhusdss_decompose_writeout.pro)
   ;; now we the redshift 'z' in the HWDR7 catalog
   absorbers[j].zqso  = mgii[i].z
   if ((~keyword_set(boss)) and (~keyword_set(dr12))) then absorbers[j].err_zqso  = mgii[i].zerr
;  absorbers[j].spec_snr_median = mgii[i].spec_snr_median
   jj = where(orispec.ivar gt 0., nn)
   if (nn gt 0) then absorbers[j].spec_snr_median = median(orispec.flux[jj]*sqrt(orispec.ivar[jj]))

   ;; from first pass
   absorbers[j].nabs = absorbers0.nabs
   absorbers[j].snr = absorbers0.snr
   absorbers[j].signal = absorbers0.signal
   absorbers[j].ivar = absorbers0.ivar
   absorbers[j].magic = absorbers0.magic

   ;; from first pass
   absorbers[j].zabs_firstpass = absorbers0.zabs_firstpass
   absorbers[j].ew_firstpass = absorbers0.ew_firstpass
   absorbers[j].err_ew_firstpass = absorbers0.err_ew_firstpass

   ;; from final pass
   absorbers[j].zabs = absorbers0.zabs         ;; mean redshift
   absorbers[j].err_zabs = absorbers0.err_zabs ;; ensemple error of mean redshift
   absorbers[j].zabs_all = absorbers0.zabs_all ;; redshift for each line
   absorbers[j].ew = absorbers0.ew             ;; ew from Gaussian fit
   absorbers[j].err_ew = absorbers0.err_ew     ;; error of ew from Gaussian fit
   absorbers[j].vdisp = absorbers0.vdisp       ;; mean vdisp, likely have not physical meaning
   absorbers[j].err_vdisp = absorbers0.err_vdisp ;; ensemple error of mean vdisp
   absorbers[j].vdisp_all = absorbers0.vdisp_all ;; vdisp of each singlet/doublet
   absorbers[j].err_vdisp_all = absorbers0.err_vdisp_all ;; error of vdisp from Gaussian fit 
endfor

mwrfits, absorbers, outfile, /create

end
