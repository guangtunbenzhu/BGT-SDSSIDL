pro jhusdss_add_hw_absorptioninfo

qso = jhusdss_qso_readin()
n_qso = n_elements(qso)

balfile = '~/SDATA/SDSS/QSO/HW_dr7qso_newz_photo_more.fits'
bal = mrdfits(balfile, 1)

mgiifile = '~/DATA/SDSS/QSO/NMF/107/Absorbers/ALLQSO_Trimmed_QSO_Absorbers_NMF_107_MgII_ALL.fits'
mgii = mrdfits(mgiifile, 1)

civfile = '~/SDATA/SDSS/CIV/Master_Cooksey_CIV_Catalog.fits'
civ = mrdfits(civfile, 1)

dlafile = '~/SDATA/SDSS/DLA/Noterdaeme_DLA.fits'
dla = mrdfits(dlafile, 1)

llsfile1 = '~/SDATA/SDSS/DLA/Prochaska_LLS_Intervening.fits'
lls1 = mrdfits(llsfile1, 1)
llsfile2 = '~/SDATA/SDSS/DLA/Prochaska_LLS_Proximate.fits'
lls2 = mrdfits(llsfile2, 1)

outstr = replicate({ra:0D0, dec:0D0, $
          BAL_flag:0L, $
          MgII:0B, N_MgII:0L, REW_MgII_Weak:0., $
          CIV:0B, N_CIV:0L, REW_CIV_Weak:0., $
          DLA:0B, LLS_Intervening:0B, LLS_Proximate:0B}, n_qso)

;; BAL
outstr.ra = qso.ra
outstr.dec = qso.dec
outstr.bal_flag = bal.bal_flag

;; C IV
spherematch, qso.ra, qso.dec, civ.ra, civ.dec, 1./3600., m1, m2
outstr[m1].civ = 1B
outstr[m1].n_civ = civ[m2].nabs
for i=0L, n_elements(m1)-1L do outstr[m1[i]].rew_civ_weak = min(civ[m2[i]].ew[0,0:civ[m2[i]].nabs-1])

;; Mg II
spherematch, qso.ra, qso.dec, mgii.ra, mgii.dec, 1./3600., m1, m2
outstr[m1].mgii = 1B
outstr[m1].n_mgii = mgii[m2].nabs
for i=0L, n_elements(m1)-1L do outstr[m1[i]].rew_mgii_weak = min(mgii[m2[i]].rew_mgii_2796[0:mgii[m2[i]].nabs-1])

;; DLA
for i=0L, n_elements(dla)-1L do begin
    imatch= where(qso.plate eq dla[i].plate and qso.fiber eq dla[i].fiber, nmatch)
    if nmatch eq 0 or nmatch gt 1 then begin
       message, "something's wrong, nmatch="+string(nmatch)
    endif
    outstr[imatch].dla = 1B
endfor

;; LLS
for i=0L, n_elements(lls1)-1L do begin
    imatch= where(qso.plate eq lls1[i].plate and qso.fiber eq lls1[i].fiber, nmatch)
    if nmatch eq 0 or nmatch gt 1 then begin
       message, "something's wrong, nmatch="+string(nmatch)
    endif
    outstr[imatch].lls_intervening = 1B
endfor

;; LLS
for i=0L, n_elements(lls2)-1L do begin
    imatch= where(qso.plate eq lls2[i].plate and qso.fiber eq lls2[i].fiber, nmatch)
    if nmatch eq 0 or nmatch gt 1 then begin
       message, "something's wrong, nmatch="+string(nmatch)
    endif
    outstr[imatch].lls_proximate = 1B
endfor
    
outfile = '~/SDATA/SDSS/QSO/HW_dr7qso_newz_absorption_info.fits'
mwrfits, outstr, outfile, /create

end
