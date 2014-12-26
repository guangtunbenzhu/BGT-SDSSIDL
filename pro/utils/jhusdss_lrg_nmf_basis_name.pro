function jhusdss_lrg_nmf_basis_name, nmfver

   if (n_elements(nmfver) eq 0) then message, 'nmfver required.'
   return, 'NMF_LRG_BASISSET_'+string(nmfver, format='(i3.3)')+'.fits'
end
