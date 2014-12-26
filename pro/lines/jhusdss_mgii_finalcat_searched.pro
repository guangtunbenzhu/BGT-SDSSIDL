pro jhusdss_mgii_finalcat_searched, nmfver, boss=boss, dr12=dr12, overwrite=overwrite, addall=addall

if (n_elements(nmfver) eq 0) then message, 'nmfver required'

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

    qsofile = path+'/ALLQSO_Trimmed_QSO_'+jhusdss_absorbers_filename(nmfver, boss=boss, dr12=dr12, /mgii)
    if (keyword_set(addall)) then  qsofile = repstr(qsofile, '.fits', '_ALL.fits')
    if (~file_test(qsofile)) then message, 'Need '+qsofile+' to fix N_abs'

    outfile = path+'/Trimmed_OnlyAbsorbers_Window_'+jhusdss_absorbers_filename(nmfver, /mgii, boss=boss, dr12=dr12)
    outfile = repstr(outfile, '.fits', '_Searched_Quasars.fits')
    if (keyword_set(addall)) then outfile = repstr(outfile, '.fits', '_ALL.fits')

    if (file_test(outfile) and ~keyword_set(overwrite)) then begin
       splog, "File already exists, use /overwrite to overwrite."
       return
    endif else begin
       splog, "Will write into this file: "
       splog, outfile
    endelse

;;  by absorber
    a = jhusdss_absorber_readin(nmfver, boss=boss, dr12=dr12, /byqso)

;;  Fix N_abs
    print, "use "+qsofile+" to fix N_abs"
    qso = mrdfits(qsofile, 1)
    spherematch, a.ra, a.dec, qso.ra, qso.dec, 1./3600., m1, m2, maxmatch=0
    a[m1].nabs = qso[m2].nabs

    select_tags = ['RA', 'DEC', 'PLATE', 'MJD', 'FIBER', 'ZQSO', 'NABS']
    outstr = struct_selecttags(a, select_tags=select_tags)
    mwrfits, outstr, outfile, /create
    
;;  by qso
;;  This needs to be re-treated because the quasar catalog includes 'bad' absorbers

end
