function jhusdss_galqso_match_readin, boss=boss

   qso_path = jhusdss_get_path(/qso)
   infile = qso_path+'/'+jhusdss_galqso_matchfile(boss=boss)
   return, mrdfits(infile, 1)

end
