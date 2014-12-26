jhusdss_garching_galqso_match
    Galaxies (DR7, gal_info_dr7_v5_2.fit.gz) - Quasars (DR7, HW_dr7qso_newz.fits)
    Galaxies (DR7, gal_info_dr7_v5_2.fit.gz) - Quasars (DR9, VAC5.fits)
    Galaxies (DR7, gal_info_dr7_v5_2.fit.gz) - Quasars (DR6, NBCKDE_z0.5_good.fits)

    jhusdss_galqso_match_readin: jhusdss_get_path(/qso)+jhusdss_galqso_matchfile
    /home/menard/DATA/SDSS/QSO/GAL_QSO_Match.fits (DR7-DR7)
    /home/menard/DATA/SDSS/QSO/GAL_QSO_Match_BOSS.fits (DR7-DR9)
    /home/menard/DATA/SDSS/QSO/GAL_QSO_Match_NBCKDE.fits (DR7-DR6 photo)

jhusdss_garching_gallrg_match
    Galaxies (DR7, gal_info_dr7_v5_2.fit.gz) - LRGs (DR7, gal_lrg_index_dr7_v5_2.fits)
    Galaxies (DR7, gal_info_dr7_v5_2.fit.gz) - LRGs (DR9, PCA_all_dr9.fits)

    jhusdss_gallrg_match_readin: jhusdss_get_path(/garching)+jhusdss_gallrg_matchfile
    /export/scratch1/menard/gz323/SDSS/Garching/GAL_LRG_Match.fits (DR7-DR7)
    /export/scratch1/menard/gz323/SDSS/BOSS/GAL_LRG_Match_BOSS.fits (DR7-DR9)

    Need to apply an empirical correction:

jhusdss_garching_lrgqso_match
    LRGs (DR7, gal_lrg_index_dr7_v5_2.fits) - Quasars (DR7, HW_dr7qso_newz.fits)
    LRGs (DR7, gal_lrg_index_dr7_v5_2.fits) - Quasars (DR9, VAC5.fits)
#   LRGs (DR9, PCA_all_dr9.fits) - Quasars (DR7, HW_dr7qso_newz.fits)
#   02/18/2013
    LRGs (DR9, wisconsin_pca_bc03-v5_6_0_z0.4.fits.gz) - Quasars (DR7, HW_dr7qso_newz.fits)
    LRGs (DR9, PCA_all_dr9.fits) - Quasars (DR9, VAC5.fits)

    jhusdss_lrgqso_match_readin: jhusdss_get_path(/garching)+jhusdss_gallrg_matchfile
    /export/scratch1/menard/gz323/SDSS/QSO/DR7_LRG_QSO_Match.fits (DR7-DR7)
    /export/scratch1/menard/gz323/SDSS/QSO/DR7_LRG_QSO_Match_BOSS.fits (DR7-DR9)
    /export/scratch1/menard/gz323/SDSS/QSO/BOSS_LRG_QSO_Match.fits (DR9-DR7)
    /export/scratch1/menard/gz323/SDSS/QSO/BOSS_LRG_QSO_Match_BOSS.fits (DR9-DR9)

###
Wait
jhusdss_garching_lrglrg_match
    LRGs (DR9, PCA_all_dr9.fits) - LRGs (DR9, PCA_all_dr9.fits)

    jhusdss_lrglrg_match_readin: jhusdss_get_path(/garching)+jhusdss_gallrg_matchfile
    /export/scratch1/menard/gz323/SDSS/Garching/GAL_LRG_Match.fits (DR7-DR7)
    /export/scratch1/menard/gz323/SDSS/Garching/GAL_LRG_Match_BOSS.fits (DR7-DR9)


