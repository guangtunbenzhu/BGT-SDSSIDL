function jhusdss_gallrg_match_readin, boss=boss

   if (keyword_set(boss)) then begin
      lrg_path = jhusdss_get_path(/bosslrg)
      infile = lrg_path+'/'+jhusdss_gallrg_matchfile(boss=boss)
   endif else begin
      lrg_path = jhusdss_get_path(/garching)
      infile = lrg_path+'/'+jhusdss_gallrg_matchfile(boss=boss)
   endelse

   return, mrdfits(infile, 1)

end
