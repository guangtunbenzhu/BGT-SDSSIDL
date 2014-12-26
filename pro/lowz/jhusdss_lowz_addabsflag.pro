pro jhusdss_lowz_addabsflag, nmfver, minrad=minrad, maxrad=maxrad, boss=boss, overwrite=overwrite

if (n_elements(nmfver) eq 0) then message, 'nmfver required'
if (n_elements(minrad) eq 0) then minrad=0.0
if (n_elements(maxrad) eq 0) then maxrad=0.8

if (keyword_set(boss)) then begin
   lowzpath=jhusdss_get_path(/nmfqso)+'/'+string(nmfver, format='(I3.3)')+'/Lowz_BOSS'
endif else begin
   lowzpath=jhusdss_get_path(/nmfqso)+'/'+string(nmfver, format='(I3.3)')+'/Lowz'
endelse

outfile = lowzpath+'/Flag_Rad'+string(minrad,format='(f4.2)')+'_'+string(maxrad,format='(f4.2)')+'_' $
       + jhusdss_lowz_spec_filename(nmfver)

;; match file
spec = jhusdss_galqso_match_spec_readin(nmfver, minrad, maxrad, boss=boss, match=match)
;; mgii catalog
qso = jhusdss_absorber_readin(nmfver, /byqso, /notrim)
qso_feii = jhusdss_absorber_readin(nmfver, /byqso, /feii, /notrim)
;; garching
garching_path = jhusdss_get_path(/garching)
garching_file = garching_path+'/'+'gal_info_dr7_v5_2.fit.gz'
gal = mrdfits(garching_file, 1)

flag = replicate({flag:0b}, n_elements(match))
lines = jhusdss_finalpass_lines()

for ife=0L, 1L do begin
    case ife of
         0: zabs = qso.zabs
         1: zabs = qso_feii.zabs
    endcase
    for iline=0L, n_elements(lines)-1L do begin
        for jabs=0L, 9L do begin
            lineswave_observer = lines[iline].wave*(1.+zabs[jabs, match.index_qso])
            caii_blue = 3900.*(1.+gal[match.index_gal].z)
            caii_red = 4000.*(1.+gal[match.index_gal].z)
            iflag = where(lineswave_observer gt caii_blue and lineswave_observer lt caii_red, nflag)
            if nflag gt 0 then flag[iflag].flag = 1b
        endfor
    endfor
endfor

mwrfits, flag, outfile, /create

end
