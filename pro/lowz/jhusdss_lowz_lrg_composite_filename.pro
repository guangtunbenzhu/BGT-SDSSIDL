function jhusdss_lowz_lrg_composite_filename, nmfver

filename = 'Lowz_LRG_Composite_'+string(nmfver, format='(i3.3)')+'.fits'
return, filename

end
