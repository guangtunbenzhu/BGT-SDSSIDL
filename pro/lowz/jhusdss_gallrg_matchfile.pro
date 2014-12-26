function jhusdss_gallrg_matchfile, boss=boss
    if (keyword_set(boss)) then return, 'GAL_LRG_Match_BOSS.fits'
    return, 'GAL_LRG_Match.fits'

end
