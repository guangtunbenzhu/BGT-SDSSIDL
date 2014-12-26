function jhusdss_galqso_match_readin, boss=boss, azimuthal=azimuthal, nbckde=nbckde

   qso_path = jhusdss_get_path(/qso)
   infile = qso_path+'/'+jhusdss_galqso_matchfile(boss=boss, nbckde=nbckde)
   if (keyword_set(azimuthal)) then infile = repstr(infile, '.fits', '_azimuthal.fits')
   return, mrdfits(infile, 1)

end
