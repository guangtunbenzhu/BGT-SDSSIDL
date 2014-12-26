;+
; Guangtun Ben Zhu, JHU
;-
function jhusdss_absorber_readin, nmfver, byqso=byqso, $
   train=train, boss=boss, dr12=dr12, notrim=notrim, feii=feii, keepasso=keepasso, $
   keepciii=keepciii, keeplowzlowsnr=keeplowzlowsnr, addall=addall

if (n_elements(nmfver) eq 0) then message, "nmfver required."

;; output
if (keyword_set(dr12)) then begin
   path=jhusdss_get_path(/nmfqso)+'/'+string(nmfver, format='(I3.3)')+'/Absorbers_DR12'
   statpath=jhusdss_get_path(/nmfqso)+'/'+string(nmfver, format='(I3.3)')+'/Decompose_DR12'
endif else begin
   if (keyword_set(boss)) then begin
      path=jhusdss_get_path(/nmfqso)+'/'+string(nmfver, format='(I3.3)')+'/Absorbers_BOSS'
      statpath=jhusdss_get_path(/nmfqso)+'/'+string(nmfver, format='(I3.3)')+'/Decompose_BOSS'
   endif else begin
      path=jhusdss_get_path(/nmfqso)+'/'+string(nmfver, format='(I3.3)')+'/Absorbers'
      statpath=jhusdss_get_path(/nmfqso)+'/'+string(nmfver, format='(I3.3)')+'/Decompose'
   endelse
endelse

if (jhusdss_direxist(path) eq 0) then message, "Can't find the directory."
if (jhusdss_direxist(statpath) eq 0) then message, "Can't find the directory."

if (~keyword_set(train)) then begin
   filename = path+'/Window_'+jhusdss_absorbers_filename(nmfver, /mgii, boss=boss, dr12=dr12)
   abs_filename = path+'/OnlyAbsorbers_Window_'+jhusdss_absorbers_filename(nmfver, /mgii, boss=boss, dr12=dr12)
   stat_filename = statpath+'/'+jhusdss_stat_filename(nmfver, boss=boss, dr12=dr12)
endif else begin
   filename = path+'/Window_Pitts_'+jhusdss_absorbers_filename(nmfver, /mgii, boss=boss, dr12=dr12)
   abs_filename = path+'/OnlyAbsorbers_Window_Pitts_'+jhusdss_absorbers_filename(nmfver, /mgii, boss=boss, dr12=dr12)
   stat_filename = statpath+'/Pitts_'+jhusdss_stat_filename(nmfver, boss=boss, dr12=dr12)
endelse

if (keyword_set(byqso)) then begin

   if (~keyword_set(feii)) then begin
      qso0 = mrdfits(filename, 1)
   endif else begin
      qso0 = mrdfits(filename, 3)
   endelse

   stat0 = mrdfits(stat_filename, 1)
   if (~keyword_set(notrim)) then begin
       index = jhusdss_absorber_trim(stat0, /qso)
       qso0 = qso0[index]
       stat0 = stat0[index]
       help, index
   endif
   outqso0 = struct_addtags(qso0, stat0)
   return, outqso0
endif 

if (~keyword_set(feii)) then begin
   absorbers0 = mrdfits(abs_filename, 1) 
endif else begin
   absorbers0 = mrdfits(abs_filename, 3) 
endelse

if (keyword_set(addall)) then begin
   absorbers0 = [mrdfits(abs_filename,1), mrdfits(abs_filename,2), mrdfits(abs_filename,3), mrdfits(abs_filename,4)]
endif

if (~keyword_set(notrim)) then begin
   index = jhusdss_absorber_trim(absorbers0, keepciii=keepciii, keepasso=keepasso, keeplowzlowsnr=keeplowzlowsnr, addall=addall)
   absorbers0 = absorbers0[index]
endif

;; Note N_abs is before trimming
return, absorbers0

end
