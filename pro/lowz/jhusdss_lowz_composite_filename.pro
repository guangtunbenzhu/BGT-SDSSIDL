function jhusdss_lowz_composite_filename, nmfver

filename = 'Lowz_Composite_'+string(nmfver, format='(i3.3)')+'.fits'
return, filename

end
