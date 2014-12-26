;pro brice_composite, nmfver, boss=boss, overwrite=overwrite, deltaz=dz

@jhusdss_composite_engine

nmfver = 005
;absorber_all = jhusdss_absorber_readin(nmfver, /byabs, /trim, boss=boss)
absorber_all = jhusdss_absorber_readin(nmfver);, /trim, boss=boss)
absorber = absorber_all
;absorber = absorber_all(where(absorber_all.zabs lt 0.5))
absorber = absorber(SORT(absorber.zabs))
n_abs = n_elements(absorber)


;; number of absorbers to group
n_in_group = 50
n_z_bin = n_abs / n_in_group



;; observer-frame wavelength
minwave = 3800.
maxwave = 9200.
loglam = jhusdss_get_loglam(minwave=minwave, maxwave=maxwave)
n_lambda = n_elements(loglam)

inimage = fltarr(n_lambda, n_z_bin)+1.


print,n_lambda,n_z_bin

print,'n_abs=',n_abs
print,'n_in_group=',n_in_group
print,'n_z_bin=',n_z_bin

abs_counter = 0L
i_group = 0L
WHILE i_group lt n_z_bin-1 do begin
        
    ind_min = abs_counter
    ind_max = abs_counter + n_in_group - 1
    
    print,ind_min,ind_max,abs_counter,absorber[ind_max].zabs
    z_median = MEDIAN(absorber[ind_min:ind_max].zabs)
    comp = jhusdss_absorbers_composite_engine(absorber[ind_min:ind_max], nmfver=nmfver)
    

    abs_counter = abs_counter + n_in_group * 1.
    i_group++

    y = INTERPOL(comp.fluxmedian,comp.wave*(1.+z_median),10^loglam)
    inimage[*,i_group] = smooth(y,5)

;    contour,(inimage>0.8)<1.2,/fill,nlevels=10
;    plot,smooth(comp.fluxmedian, 5)

endwhile

mwrfits,inimage,'my_image.fits',/create

stop

;contour,(inimage>0.8)<1.2,/fill,nlevels=50,$
;  xr=[minwave,maxwave],xstyle=1

end

