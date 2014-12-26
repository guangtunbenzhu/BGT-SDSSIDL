function jhusdss_lrg_readin, boss=boss

if (keyword_set(boss)) then begin
   lrgpath = jhusdss_get_path(/bosslrg)
   filename = jhusdss_boss_lrgfile()
endif else begin
   lrgpath = jhusdss_get_path(/garching)
   filename =  jhusdss_dr7_lrgfile()
endelse

infile = lrgpath+'/'+filename

return, mrdfits(infile, 1)

end
