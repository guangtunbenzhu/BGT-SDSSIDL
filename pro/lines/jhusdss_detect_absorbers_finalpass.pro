;+
;-
function jhusdss_detect_absorbers_finalpass, spec, zabs, nabsmax=nabsmax

if (n_elements(nabsmax) eq 0) then nabsmax = jhusdss_nabsmax()
outstr = jhusdss_finalpass_blank(nabsmax=nabsmax)
nzabs = (n_elements(zabs) < nabsmax)

for i=0, nzabs-1 do begin
    jhusdss_finalpass_fit, spec, zabs[i], newzabs=newzabs, $
       err_newzabs=err_newzabs, allzabs=allzabs, $
       ew=ew, err_ew=err_ew, sigma=sigma, err_sigma=err_sigma, $
       meansigma=meansigma, err_meansigma=err_meansigma
;   print, zabs[i], newzabs
    outstr.zabs[i] = newzabs
    outstr.err_zabs[i] = err_newzabs
    outstr.zabs_all[*, i] = allzabs
    outstr.ew[*, i] = ew
    outstr.err_ew[*, i] = err_ew
    outstr.vdisp[i] = meansigma
    outstr.err_vdisp[i] = err_meansigma
    outstr.vdisp_all[*, i] = sigma
    outstr.err_vdisp_all[*, i] = err_sigma
endfor

return, outstr
end
