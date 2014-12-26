function jhusdss_lowz_lrgqso_composite_filename, nmfver

filename = 'Lowz_LRGQSO_Composite_'+string(nmfver, format='(i3.3)')+'.fits'
return, filename

end
