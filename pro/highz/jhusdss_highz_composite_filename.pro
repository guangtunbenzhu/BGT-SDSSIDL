function jhusdss_highz_composite_filename, nmfver

filename = 'Highz_Composite_'+string(nmfver, format='(i3.3)')+'.fits'
return, filename

end
