function jhusdss_lrgqso_matchfile, sdsslrg=sdsslrg, boss=boss

    if (keyword_set(sdsslrg)) then prefix='DR7_' else prefix='BOSS_'
    if (keyword_set(boss)) then filename=prefix+'LRG_QSO_Match_BOSS.fits' else filename=prefix+'LRG_QSO_Match.fits'
    return, filename

end
