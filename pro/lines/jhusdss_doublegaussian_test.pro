function doublegaussian_line_func, x, p

center1=p[0]
center2=p[1]
flux1=p[2]
flux2=p[3]
sigma=p[4]

model = fltarr(n_elements(x))
model = flux1*exp(-0.5*(x-center1)^2/sigma^2)/sqrt(2.*!DPI)/sigma $
      + flux2*exp(-0.5*(x-center2)^2/sigma^2)/sqrt(2.*!DPI)/sigma

return, model
end

pro jhusdss_doublegaussian_test, plate, fiber, zabs, nmfver=nmfver

    if (n_elements(plate) eq 0) then  plate = 2035
    if (n_elements(fiber) eq 0) then  fiber = 353
    if (n_elements(zabs) eq 0) then  zabs = 0.859225
    if (n_elements(nmfver) eq 0) then  nmfver = 106
    a = jhusdss_decompose_loadspec(plate, fiber, nmfver)
    b = jhusdss_convolve_loadspec(plate, fiber, nmfver)
    newabs = jhusdss_detect_absorbers_newengine(b, a)
    print, 'zqso= ', a.z
    if (newabs.nabs gt 0) then begin
       print, 'zabs= ', newabs.zabs[0:newabs.nabs-1]
       print, 'zabs_firstpass= ', newabs.zabs_firstpass[0:newabs.nabs-1]
;   print, reform(newabs.ew_fit[4,0:newabs.nabs-1])
;   print, reform(newabs.ew[0,0:newabs.nabs-1])
       print, 'w(2796)_fit= ', reform(newabs.ew[5,0:newabs.nabs-1])
       print, 'w(2796)_sum= ', reform(newabs.ew_firstpass[1,0:newabs.nabs-1])
    endif
;   print, reform(newabs.magic[0,*])
;   print, reform(newabs.magic[1,*])
;   print, reform(newabs.magic[2,*])
;   print, newabs.sigma_fit_mean
;   print, reform(newabs.ew_fit[4, *])
;   print, reform(newabs.ew_fit[5, *])
;   print, reform(newabs.ew_fit[5, *]/(newabs.ew_fit[4,*]+(newabs.ew_fit[4,*] eq 0.)))
;   print, newabs.zabs
;   print, a.z, 1216.*(1.+a.z), 1550.*(1.+a.z)

    lambda = a.wave*(1.+a.z)
    tmpivar = a.ivar*a.nmf_continuum^2*a.med_continuum^2
;   ilambda = where(lambda gt 2796.*(1.+zabs)-20 and lambda lt 2803.*(1.+zabs)+20.)
    ilambda = where(lambda gt (2796.35-8.)*(1.+zabs) and lambda lt (2803.53+8.)*(1.+zabs))
    
    djs_plot, lambda, smooth(a.flux, 5), xra=[2000., 3500]*(1.+zabs), xstyle=1, yra=[0.0, 1.5], ystyle=1
;   djs_plot, lambda, a.flux, xra=[3700., 9200], xstyle=1, yra=[0.0, 1.5], ystyle=1
    djs_oplot, lambda, smooth(a.residual,5), color='red'
    djs_oplot, lambda, 1.-1./sqrt(tmpivar), color='gray'
    djs_oplot, lambda, 1.+1./sqrt(tmpivar), color='gray'

    djs_oplot, [2803,2803]*(1.+zabs), !y.crange, color='green'
    djs_oplot, [2796,2796]*(1.+zabs), !y.crange, color='green'
    djs_oplot, [2600,2600]*(1.+zabs), !y.crange, color='green'
    djs_oplot, [2586,2586]*(1.+zabs), !y.crange, color='green'
    djs_oplot, [2374,2374]*(1.+zabs), !y.crange, color='green'
    djs_oplot, [2344,2344]*(1.+zabs), !y.crange, color='green'
    stop

    zz = min(abs(b.zgrid-zabs), imin)
    print, b.ewall[*,imin]
    print, b.signal[*,imin]*sqrt(b.ivar[*,imin])
    in_lambda = lambda[ilambda]
    in_flux = 1.-a.residual[ilambda]
    in_ivar = a.ivar[ilambda]
    in_line = [2796.35, 2803.53]*(1.+zabs)
;   in_lflux = [1.45, 1.01]
    in_lflux = [((b.ewall[1, imin]>b.ewall[0, imin])>0.), $
                (((b.ewall[1, imin]>b.ewall[0, imin])/1.3)>0.)]

    jhusdss_doublegaussian_fit, in_lambda, in_flux, in_ivar, in_line, in_lflux, $
               center=center, err_center=err_center, $
               lflux=lflux, err_lflux=err_lflux, sigma=sigma, err_sigma=err_sigma

    print, in_line
    print, center
    print, (center[1]-center[0])/2800./(1.+zabs)
    print, center[0]/2796.53-1., center[1]/2803.35-1.
;   print, in_lflux
    print, lflux/(1.+zabs)
    print, sigma/alog(10.)/2800./(1.+zabs)

;   print, ew
;   print, sigma*alog(10.)*3.0D+5

    x = in_lambda
    p = [center, lflux, sigma]
    y = doublegaussian_line_func(x,p)
    djs_plot, in_lambda, in_flux, psym=4
    djs_oplot, in_lambda, y, color='red'
    djs_oplot, in_lambda, in_ivar
;   print, newabs.magic[*, 0:newabs.nabs-1]
    print, median(a.flux*sqrt(a.ivar))
    stop


end
