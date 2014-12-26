pro jhusdss_grb_mc, overwrite=overwrite

nmfver = 107
path = jhusdss_get_path(/nmfqso)+'/'+string(nmfver, format='(I3.3)')+'/Absorbers'
outfile = path+'/ALLQSO_GRB_MC.fits'
if (file_test(outfile) and ~keyword_set(overwrite)) then begin
   splog, 'File already exists. Use /overwrite if you want to overwrite it.'
   return
endif else begin
   splog, 'Will write the absorber catalog into this file: '
   print, outfile
endelse

qsofile = path+'/ALLQSO_Trimmed_QSO_FEII_'+jhusdss_absorbers_filename(nmfver, boss=boss, /mgii)
qso = mrdfits(qsofile, 1)

grbfile = getenv('ALL_DATA')+'/GRB/'+'gzfunc.dat'
readcol, grbfile, z, gz_f, gz_i, format='(f,f,f)'
nz = n_elements(z)
deltaz = median(z[1:nz-1] - z[0:nz-2])
ngrb_f = round(gz_f/deltaz)
ngrb_i = round(gz_i/deltaz)

zmgii_min = (qso.zqso*1550./2800.+0.02) > (3800./2800.-1.)
zmgii_max = (qso.zqso-0.04) < (9200./2800.-1.)

nabsmax = 10
nmc = 10000L
nabs_mc_f = lonarr(nz, nmc)
nabs_mc_i = lonarr(nz, nmc)

zrange = where(z gt 0.4 and z lt 2.0, nzrange)
min_zrange = min(zrange)
max_zrange = max(zrange)

for iz=min_zrange, max_zrange do begin
    counter, iz, nz
    this_z = z[iz]
    this_ngrb_f = ngrb_f[iz]
    this_ngrb_i = ngrb_i[iz]
    if (this_ngrb_f eq 0) then continue

    ;; select all quasars with Zmin<Z<Zmax and Wmin(Z) > 1
    iqso = where(this_z gt zmgii_min and this_z lt zmgii_max, nqso)

    ;; randomly select ngrb
    iran_f = floor(randomu(seed, this_ngrb_f, nmc)*nqso)
    zabs_f = fltarr(nabsmax, this_ngrb_f, nmc)
    for igrb_f=0L, this_ngrb_f-1L do zabs_f[*, igrb_f, *] = qso[iqso[iran_f[igrb_f, *]]].zabs
    for imc=0L, nmc-1L do begin
        iabs_f = where(zabs_f[*,*,imc] gt this_z and zabs_f[*,*,imc] le this_z+deltaz, nabs_f)
        nabs_mc_f[iz, imc] = nabs_f
    endfor

    if (this_ngrb_i gt 0) then begin
       iran_i = floor(randomu(seed, this_ngrb_i, nmc)*nqso)
       zabs_i = fltarr(nabsmax, this_ngrb_i, nmc)
       for igrb_i=0L, this_ngrb_i-1L do zabs_i[*, igrb_i, *] = qso[iqso[iran_i[igrb_i, *]]].zabs
       for imc=0L, nmc-1L do begin
           iabs_i = where(zabs_i[*,*,imc] gt this_z and zabs_i[*,*,imc] le this_z+deltaz, nabs_i)
           nabs_mc_i[iz, imc] = nabs_i
       endfor
    endif

endfor

strtmp = {z:z, ngrb_f:ngrb_f, ngrb_i:ngrb_i, nabs_f:nabs_mc_f, nabs_i:nabs_mc_i}
mwrfits, strtmp, outfile, /create

end
