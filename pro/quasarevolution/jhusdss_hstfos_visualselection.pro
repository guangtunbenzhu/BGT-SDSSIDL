pro jhusdss_hstfos_visualselection

qsofile = '~/SDATA/Quasars/HSTFOS/table1_master.fits'
qso = mrdfits(qsofile, 1)
nqso = n_elements(qso)

observer_file = '~/SDATA/SDSS/AllInOne/AIO_QSO_HSTFOS_ObserverFrame_Wave01100_03350A.fits'
spec = mrdfits(observer_file, 1)

outfile = '~/SDATA/Quasars/HSTFOS/hstfos_master_visual.fits'
outstr = replicate({name:'Name', ra:0.D0, dec:0.D0, z:0., h130:0L, h190:0L, h270:0L, isgood:0L}, nqso)
outstr.name = qso.name
outstr.ra = qso.ra2000
outstr.dec = qso.dec2000
outstr.z = qso.z
outstr.h130 = qso.h130
outstr.h190 = qso.h190
outstr.h270 = qso.h270

isgood = 'n'
for i=0L, nqso-1L do begin
    djs_plot, spec.wave, smooth(spec.flux[i,*], 5), xra=[1150, 3300], xst=1, yst=1

    djs_oplot, [2796., 2796.], !y.crange, color='green', thick=2, linestyle=2
    djs_oplot, [2803., 2803.], !y.crange, color='green', thick=2, linestyle=2
    djs_oplot, [1260., 1260.], !y.crange, color='green', thick=2, linestyle=2
    djs_oplot, [1393., 1393.], !y.crange, color='green', thick=2, linestyle=2

    djs_oplot, [2796., 2796.]*(1.+qso[i].z), !y.crange, color='red'
    djs_oplot, [2803., 2803.]*(1.+qso[i].z), !y.crange, color='magenta'
    djs_oplot, [1550., 1550.]*(1.+qso[i].z), !y.crange, color='blue'
    djs_oplot, [1216., 1216.]*(1.+qso[i].z), !y.crange, color='green'
    djs_oplot, [1032., 1032.]*(1.+qso[i].z), !y.crange, color='yellow'

    ; print, 'absorption, (y/n)'
    ; read, absorption
    ; if (absorption eq 'y') then outstr[i].absorption = 1L
    print, 'isgood, (y/n):'
    read, isgood 
    if (isgood eq 'y') then outstr[i].isgood= 1L
endfor

mwrfits, outstr, outfile, /create

end
