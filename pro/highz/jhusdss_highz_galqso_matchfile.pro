function jhusdss_highz_galqso_matchfile, boss=boss
    if (keyword_set(boss)) then return, 'LRG_QSO_Match_BOSS.fits'
    return, 'LRG_QSO_Match.fits'

end
