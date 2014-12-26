listfile = getenv('ALL_DATA')+'/Masha/Weirong/spectralist.dat'
readcol, listfile, atemp, format='A'

nmfver = 106L
qso = jhusdss_qso_readin(boss=boss)
stats = jhusdss_qsostats_readin(nmfver, boss=boss)
list = replicate(qso[0], n_elements(atemp))
stats_list = replicate(stats[0], n_elements(atemp))

for i=0L, n_elements(atemp)-1L do begin
    list[i].MJD =  long(strmid(atemp[i], 7, 5))
    list[i].PLATE =  long(strmid(atemp[i], 13, 4))
    list[i].FIBER =  long(strmid(atemp[i], 18, 3))
    list[i].name = 'N/A'
    list[i].ra = -999.D0
    list[i].dec = -999.D0
    list[i].z = -999.
    list[i].zerr = -999.
    list[i].flag_first = -1
    list[i].zother = -999.
    list[i].zcode = -1
    list[i].flag_spec = -1
    list[i].flag_cat = -1

    stats_list[i].MJD =  long(strmid(atemp[i], 7, 5))
    stats_list[i].PLATE =  long(strmid(atemp[i], 13, 4))
    stats_list[i].FIBER =  long(strmid(atemp[i], 18, 3))
    stats_list[i].ra = -999.D0
    stats_list[i].dec = -999.D0
    stats_list[i].zqso = -999.
    stats_list[i].err_zqso = -999.
    stats_list[i].SPEC_SNR_MEDIAN = -999.
    stats_list[i].ISITDECOMPOSED = -1.
    stats_list[i].ISITSTATED_RED = -1
    stats_list[i].ISITSTATED_BLUE = -1
    stats_list[i].ISITCONVOLVED = -1
    stats_list[i].MED_MEAN_RED = -999.
    stats_list[i].MED_SDEVIATION_RED = -999.
    stats_list[i].MED_SKEWNESS_RED = -999.
    stats_list[i].MED_MEAN_BLUE = -999.
    stats_list[i].MED_SDEVIATION_BLUE = -999.
    stats_list[i].MED_SKEWNESS_BLUE = -999.
endfor

n_nomatch = 0L

for i=0L, n_elements(atemp)-1L do begin
    iobj = where(qso.mjd eq list[i].MJD and qso.plate eq list[i].PLATE and qso.fiber eq list[i].FIBER, nobj)
    if nobj eq 0 then begin
       print, "No Match"
       print, atemp[i]
       n_nomatch = n_nomatch+1L
    endif else begin
       if nobj gt 1 then stop, "Something's wrong"
       list[i].name = qso[iobj].name
       list[i].ra = qso[iobj].ra
       list[i].dec = qso[iobj].dec
       list[i].z = qso[iobj].z
       list[i].zerr = qso[iobj].zerr
       list[i].flag_first = qso[iobj].flag_first
       list[i].zother = qso[iobj].zother
       list[i].zcode = qso[iobj].zcode
       list[i].flag_spec = qso[iobj].flag_spec
       list[i].flag_cat = qso[iobj].flag_cat

       stats_list[i].ra = stats[iobj].ra
       stats_list[i].dec = stats[iobj].dec
       stats_list[i].zqso = stats[iobj].zqso
       stats_list[i].err_zqso = stats[iobj].err_zqso
       stats_list[i].SPEC_SNR_MEDIAN = stats[iobj].SPEC_SNR_MEDIAN
       stats_list[i].ISITDECOMPOSED = stats[iobj].ISITDECOMPOSED
       stats_list[i].ISITSTATED_RED = stats[iobj].ISITSTATED_RED
       stats_list[i].ISITSTATED_BLUE = stats[iobj].ISITSTATED_BLUE
       stats_list[i].ISITCONVOLVED = stats[iobj].ISITCONVOLVED
       stats_list[i].MED_MEAN_RED = stats[iobj].MED_MEAN_RED
       stats_list[i].MED_SDEVIATION_RED = stats[iobj].MED_SDEVIATION_RED
       stats_list[i].MED_SKEWNESS_RED = stats[iobj].MED_SKEWNESS_RED
       stats_list[i].MED_MEAN_BLUE = stats[iobj].MED_MEAN_BLUE
       stats_list[i].MED_SDEVIATION_BLUE = stats[iobj].MED_SDEVIATION_BLUE
       stats_list[i].MED_SKEWNESS_BLUE = stats[iobj].MED_SKEWNESS_BLUE
    endelse
endfor
print, n_nomatch

outfile = getenv('ALL_DATA')+'/Masha/Weirong/spectralist.fits'
mwrfits, list, outfile, /create

stats_outfile = getenv('ALL_DATA')+'/Masha/Weirong/spectralist_stats.fits'
mwrfits, stats_list, stats_outfile, /create
end
