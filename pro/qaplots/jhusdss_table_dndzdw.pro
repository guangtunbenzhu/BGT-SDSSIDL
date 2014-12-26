pro jhusdss_table_dndzdw, nmfver

if (n_elements(nmfver) eq 0) then begin
   splog, 'nmfver required.'
   return
endif

dndz = jhusdss_montecarlo_dndzdw_readin(nmfver)

path=jhusdss_get_path(/nmfqso)+'/'+string(nmfver, format='(I3.3)')
qapath = path+'/QAplots'

;; make the table
texfile = qapath+'/'+'dndzdw_z_'+string(nmfver, format='(I3.3)')+'.tex'

openw, lun, texfile, /get_lun
printf, lun, '  & & & & & & & & & & & & \\'
text=' & '
for j=0L, n_elements(dndz.z)-2L do $
    text=text+' $'+string(dndz.z_min[j], format='(f4.2)')+'-'+string(dndz.z_max[j], format='(f4.2)')+'$ & '
printf, lun, text+' $'+string(dndz.z_min[j], format='(f4.2)')+'-'+string(dndz.z_max[j], format='(f4.2)')+'$ \\'

text=' & '
for j=0L, n_elements(dndz.z)-2L do $
    text=text+' $'+string(dndz.zbin_median[j], format='(f4.2)')+'$ & '
printf, lun, text+' $'+string(dndz.zbin_median[j], format='(f4.2)')+'$ \\'

for i=0L, n_elements(dndz.w)-1L do begin 
    text = ' $'+string(dndz.w_min[i], format='(f4.2)')+'-'+string(dndz.w_max[i],format='(f4.2)')+'$ & '
    for j=0L, n_elements(dndz.z)-2L do $
        text = text+' $'+string(alog10(dndz.phi[i,j]), format='(f5.2)')+'^{+'+string(alog10(dndz.phi_upper_limit[i,j])-alog10(dndz.phi[i,j]), format='(f5.2)') $
            + '}_{-'+string(alog10(dndz.phi[i,j])-alog10(dndz.phi_lower_limit[i,j]), format='(f5.2)')+'}'+'$ & '
    printf, lun, text+' $'+string(alog10(dndz.phi[i,j]), format='(f5.2)')+'^{+'+string(alog10(dndz.phi_upper_limit[i,j])-alog10(dndz.phi[i,j]), format='(f5.2)') $
        + '}_{-'+string(alog10(dndz.phi[i,j])-alog10(dndz.phi_lower_limit[i,j]), format='(f5.2)')+'}'+'$ \\'
endfor

free_lun, lun

;; make the table
texfile = qapath+'/'+'Cum_dndzdw_z_'+string(nmfver, format='(I3.3)')+'.tex'

openw, lun, texfile, /get_lun

for i=0L, n_elements(dndz.w)-1L do begin 
    text = ' $'+string(dndz.w_min[i], format='(f4.2)')+'$ & '
    iw = value_locate(dndz.w_cum_min, dndz.w_min[i])
    for j=0L, n_elements(dndz.z)-2L do $
        text = text+' $'+string(alog10(dndz.cum_phi[iw,j]), format='(f5.2)')+'^{+'+string(alog10(dndz.cum_phi_upper_limit[iw,j])-alog10(dndz.cum_phi[iw,j]), format='(f5.2)') $
            + '}_{-'+string(alog10(dndz.cum_phi[iw,j])-alog10(dndz.cum_phi_lower_limit[iw,j]), format='(f5.2)')+'}'+'$ & '
    printf, lun, text+' $'+string(alog10(dndz.cum_phi[iw,j]), format='(f5.2)')+'^{+'+string(alog10(dndz.cum_phi_upper_limit[iw,j])-alog10(dndz.cum_phi[iw,j]), format='(f5.2)') $
        + '}_{-'+string(alog10(dndz.cum_phi[iw,j])-alog10(dndz.cum_phi_lower_limit[iw,j]), format='(f5.2)')+'}'+'$ \\'
endfor

free_lun, lun

;; make the table
texfile = qapath+'/'+'Wstar_z_'+string(nmfver, format='(I3.3)')+'.tex'
openw, lun, texfile, /get_lun
text = '$N^*$ & '
for i=0L, n_elements(dndz.n_star)-2L do $
    text =  text+' $ & $'+string(dndz.n_star[i], format='(f4.2)')+'\pm'+string(dndz.n_star_err[i], format='(f4.2)')
    printf, lun, text+' $ & $'+string(dndz.n_star[i], format='(f4.2)')+'\pm+'+string(dndz.n_star_err[i], format='(f4.2)')+'$ \\' 
text = '$W^*$ & '
for i=0L, n_elements(dndz.w_star)-2L do $
    text =  text+' $ & $'+string(dndz.w_star[i], format='(f4.2)')+'\pm'+string(dndz.w_star_err[i], format='(f4.2)')
    printf, lun, text+' $ & $'+string(dndz.w_star[i], format='(f4.2)')+'\pm+'+string(dndz.w_star_err[i], format='(f4.2)')+'$ \\' 

free_lun, lun

;; make the table
texfile = qapath+'/'+'Wstar_z_column_'+string(nmfver, format='(I3.3)')+'.tex'
openw, lun, texfile, /get_lun
for i=0L, n_elements(dndz.n_star)-1L do begin
    text = '$'+ string(dndz.z_min[i], format='(f4.2)')+'-'+string(dndz.z_max[i], format='(f4.2)') $
         + '$ & $' + string(dndz.z[i], format='(f4.2)') $
         + '$ & $' + string(dndz.n_star[i], format='(f4.2)')+'\pm'+string(dndz.n_star_err[i], format='(f4.2)') $
         + '$ & $' + string(dndz.w_star[i], format='(f4.2)')+'\pm'+string(dndz.w_star_err[i], format='(f4.2)') $
         + '$ \\'
    printf, lun, text
endfor
free_lun, lun

;; make the table
texfile = qapath+'/'+'Evolution_parameter_'+string(nmfver, format='(I3.3)')+'.tex'
openw, lun, texfile, /get_lun
;printf, lun, '    & $g_0/W_0$ & $\alpha_g/\alpha_W$ & $\beta_g/\beta_W$ & $\gamma_g/\gamma_W$ \\'
printf, lun, '$g$ & $'+string(dndz.f0, format='(f4.2)')+'\pm'+string(dndz.f0_err, format='(f4.2)') $
             +'$ & $'+string(dndz.f0_alpha, format='(f4.2)')+'\pm'+string(dndz.f0_alpha_err, format='(f4.2)') $
             +'$ & $'+string(dndz.f0_bbeta, format='(f4.2)')+'\pm'+string(dndz.f0_bbeta_err, format='(f4.2)') $
             +'$ & $'+string(dndz.f0_ggamma, format='(f4.2)')+'\pm'+string(dndz.f0_ggamma_err, format='(f4.2)') +'$ \\$'
printf, lun, '$W$ & $'+string(dndz.w0, format='(f4.2)')+'\pm'+string(dndz.w0_err, format='(f4.2)') $
             +'$ & $'+string(dndz.w0_alpha, format='(f4.2)')+'\pm'+string(dndz.w0_alpha_err, format='(f4.2)') $
             +'$ & $'+string(dndz.w0_bbeta, format='(f4.2)')+'\pm'+string(dndz.w0_bbeta_err, format='(f4.2)') $
             +'$ & $'+string(dndz.w0_ggamma, format='(f4.2)')+'\pm'+string(dndz.w0_ggamma_err, format='(f4.2)') +'$ \\$'

free_lun, lun


end
