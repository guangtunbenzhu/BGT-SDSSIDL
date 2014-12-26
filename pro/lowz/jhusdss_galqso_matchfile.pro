function jhusdss_galqso_matchfile, boss=boss, nbckde=nbckde
    if (keyword_set(boss)) then return, 'GAL_QSO_Match_BOSS.fits'
    if (keyword_set(nbckde)) then return, 'GAL_QSO_Match_NBCKDE.fits'
    return, 'GAL_QSO_Match.fits'

end
