function jhusdss_redecompose_spec, spec, zabs, abslines=abslines, new_med_continuum=new_med_continuum

   if (n_elements(spec) gt 1) then begin 
      splog, 'I can only deal with one spectrum at a time for now.'
      return, -1
   endif
   if (n_elements(abslines) eq 0) then abslines=jhusdss_redecompose_lines()

   tmpmask = jhusdss_redecompose_mask(spec.wave*(1.+spec.z), abslines, zabs)
   mask  = ((tmpmask) or (spec.ivar le 0.))

   ;; redo median filtering
   nmf_residual = spec.flux/spec.nmf_continuum
   tmpivar = spec.ivar*spec.nmf_continuum^2
   nmf_residual = reform(nmf_residual, 1, n_elements(nmf_residual))
   tmpivar = reform(tmpivar, 1, n_elements(nmf_residual))
   mask = reform(mask, 1, n_elements(nmf_residual))

   newmask = mask
   filter_sizes=[91, 163]
   jhusdss_median_filter, nmf_residual, tmpivar, mask=newmask, $
      continuum=tmp_med_continuum, residual=tmp_residual, filter_sizes=filter_sizes

   newmask = ((mask) or (abs(tmp_med_continuum-nmf_residual)*sqrt(tmpivar) gt 1.5))
   filter_sizes=[143, 71]
   jhusdss_median_filter, nmf_residual, tmpivar, mask=newmask, $
      continuum=tmp_med_continuum, residual=tmp_residual, filter_sizes=filter_sizes

   newmask = ((mask) or (abs(tmp_med_continuum-nmf_residual)*sqrt(tmpivar) gt 1.5))
   filter_sizes=[143, 71]
   jhusdss_median_filter, nmf_residual, tmpivar, mask=newmask, $
      continuum=tmp_med_continuum, residual=tmp_residual, filter_sizes=filter_sizes


   new_med_continuum = reform(tmp_med_continuum)
   residual = reform(tmp_residual)
   return, residual
end
