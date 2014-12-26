function jhusdss_lrgqso_match_readin, sdsslrg=sdsslrg, boss=boss

   qso_path = jhusdss_get_path(/qso)
   infile = qso_path+'/'+jhusdss_lrgqso_matchfile(sdsslrg=sdsslrg, boss=boss)
   return, mrdfits(infile, 1)

end
