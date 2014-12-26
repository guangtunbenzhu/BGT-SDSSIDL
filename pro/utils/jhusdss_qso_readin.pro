function jhusdss_qso_readin, boss=boss, dr12=dr12, nbckde=nbckde

qsopath = jhusdss_get_path(/qso)
filename =  jhusdss_dr7_qsofile()
if (keyword_set(nbckde)) then begin
   filename = jhusdss_nbckde_qsofile()
endif 
if (keyword_set(boss)) then begin
   filename = jhusdss_boss_qsofile()
endif 
if (keyword_set(dr12)) then begin
   filename = jhusdss_dr12_qsofile()
endif 

infile = qsopath+'/'+filename

return, mrdfits(infile, 1)

end
