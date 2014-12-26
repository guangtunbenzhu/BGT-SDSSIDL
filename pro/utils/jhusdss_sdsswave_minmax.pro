function jhusdss_sdsswave_minmax, boss=boss, dr12=dr12

if (keyword_set(boss) or keyword_set(dr12)) then begin
    return, [3600., 10400.]
endif else begin
    return, [3700., 9200.]
endelse

end
