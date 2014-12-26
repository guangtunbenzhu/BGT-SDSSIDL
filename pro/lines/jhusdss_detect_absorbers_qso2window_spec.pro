;+
;-
pro jhusdss_detect_absorbers_qso2window_spec, qso, mgii_red=mgii_red, mgii_blue=mgii_blue, $
            feii_red=feii_red, feii_blue=feii_blue

    nqso = n_elements(qso)

    str_tmp = jhusdss_qso2window_blank(/absorber)
    mgii_red = replicate(str_tmp, nqso)
    mgii_blue = replicate(str_tmp, nqso)
    feii_red = replicate(str_tmp, nqso)
    feii_blue = replicate(str_tmp, nqso)
    for i=0L, nqso-1L do begin
        counter, i+1, nqso
        nabs = qso[i].nabs
        nmgii_red = -1L
        nmgii_blue = -1L
        nfeii_red = -1L
        nfeii_blue = -1L
        for j=0L, nabs-1L do begin
            ;; criterion_mgii and in red window
            string_tmp = ''
            ;; hack 09-19-September 2014
            ;; keep mgii_red and feii_red, but change the other two to the shorter wavelength
            if (qso[i].magic[0,j] eq 1b) and (qso[i].zabs[j] gt 1550.*(1.+qso[i].zqso)/2796.35-1.) then string_tmp = 'mgii_red'
            if (qso[i].magic[0,j] eq 1b) and (qso[i].zabs[j] le 1550.*(1.+qso[i].zqso)/2796.35-1.) then string_tmp = 'mgii_blue'
;           if (qso[i].magic[0,j] eq 1b) and (qso[i].zabs[j] lt 1550.*(1.+qso[i].zqso)/2803.53-1.) then string_tmp = 'mgii_blue'
            if (qso[i].magic[0,j] eq 0b) and (qso[i].zabs[j] gt 1550.*(1.+qso[i].zqso)/2586.65-1.) then string_tmp = 'feii_red'
            if (qso[i].magic[0,j] eq 0b) and (qso[i].zabs[j] le 1550.*(1.+qso[i].zqso)/2586.65-1.) then string_tmp = 'feii_blue'
;           if (qso[i].magic[0,j] eq 0b) and (qso[i].zabs[j] lt 1550.*(1.+qso[i].zqso)/2600.17-1.) then string_tmp = 'feii_blue'

            if (string_tmp eq '') then message, 'Fix this first!'
;           stop
            temp = execute('n'+string_tmp+'++')
            temp = execute(string_tmp+'[i].nabs=n'+string_tmp+'+1')
            temp = execute(string_tmp+'[i].snr[*, n'+string_tmp+']=qso[i].snr[*,j]')
            temp = execute(string_tmp+'[i].signal[*, n'+string_tmp+']=qso[i].signal[*,j]')
            temp = execute(string_tmp+'[i].criterion_mgii[n'+string_tmp+']=qso[i].magic[0,j]')
            temp = execute(string_tmp+'[i].criterion_mgii_feii[n'+string_tmp+']=qso[i].magic[1,j]')
            temp = execute(string_tmp+'[i].criterion_feii[n'+string_tmp+']=qso[i].magic[2,j]')
            temp = execute(string_tmp+'[i].zabs[n'+string_tmp+']=qso[i].zabs[j]')
            temp = execute(string_tmp+'[i].err_zabs[n'+string_tmp+']=qso[i].err_zabs[j]')
            temp = execute(string_tmp+'[i].zabs_all[*, n'+string_tmp+']=qso[i].zabs_all[*, j]')
            temp = execute(string_tmp+'[i].ew[*, n'+string_tmp+']=qso[i].ew[*, j]')
            temp = execute(string_tmp+'[i].err_ew[*, n'+string_tmp+']=qso[i].err_ew[*, j]')
            temp = execute(string_tmp+'[i].vdisp[n'+string_tmp+']=qso[i].vdisp[j]')
            temp = execute(string_tmp+'[i].err_vdisp[n'+string_tmp+']=qso[i].err_vdisp[j]')
            temp = execute(string_tmp+'[i].vdisp_all[*, n'+string_tmp+']=qso[i].vdisp_all[*, j]')
            temp = execute(string_tmp+'[i].err_vdisp_all[*, n'+string_tmp+']=qso[i].err_vdisp_all[*, j]')
       endfor
    endfor

end
