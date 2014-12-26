;+
; Documentation needed!
; Ideally, this is the only place one needs to change for path/directory
;-

function jhusdss_get_path, spectra=spectra, qso=qso, composite=composite, $
   nmfqso=nmfqso, absorber=absorber, boss_spectra=boss_spectra, garching=garching, $
   bosslrg=bosslrg, fitlrg=fitlrg, nmflrg=nmflrg, sky=sky, dr12_spectra=dr12_spectra

parent_path = jhusdss_get_parent_path()

;; where the spectra are
if (keyword_set(spectra)) then path=getenv('RAW_DATA')+'/SDSS/Spectra/'
if (keyword_set(boss_spectra)) then path=getenv('RAW_DATA')+'/SDSS3/BOSS/groups/boss/spectro/redux/v5_4_45/'
if (keyword_set(dr12_spectra)) then path=getenv('ALL_DATA')+'/SDSS3/v5_7_0/'

;; where the QSO catalogs are
if (keyword_set(qso)) then path=getenv('RAW_DATA')+'/SDSS/QSO/'

;; where the known absorber catalogs are
if (keyword_set(absorber)) then path=getenv('RAW_DATA')+'/SDSS/Absorber/'

;; where to store the composite spectra
if (keyword_set(composite)) then path=getenv('ALL_DATA')+'/SDSS/QSO/Composite/'

;; where to store NMF basis vectors and decomposed spectra
if (keyword_set(nmfqso)) then path=getenv('ALL_DATA')+'/SDSS/QSO/NMF/'
if (keyword_set(nmflrg)) then path=getenv('ALL_DATA')+'/SDSS/LRG/NMF/'

;; garching
if (keyword_set(garching)) then path = parent_path+'/SDSS/Garching/'

;; BOSS lrg
if (keyword_set(bosslrg)) then path=parent_path+'/SDSS/BOSS/'

;; where to store lrg basis vectors and decomposed spectra
if (keyword_set(fitlrg)) then path=getenv('ALL_DATA')+'/SDSS/LRG/'

;; sky
if (keyword_set(sky)) then path=parent_path+'/SDSS/SKY/'

;; test
if (jhusdss_direxist(path) eq 0) then message, "Well, I can't find your data!"

return, path
end
