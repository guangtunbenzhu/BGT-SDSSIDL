;+
; Documentation Needed
;-
function jhusdss_decompose_bins, option=option

case option of
   1: strtmp = {zmin:0.1, zmax:0.6, basisfile:'QSO_NMF_basis_z000_100_norm4150.fits'}
   2: strtmp = {zmin:0.6, zmax:1.0, basisfile:'QSO_NMF_basis_z040_179_norm3020.fits'}
   3: strtmp = {zmin:1.0, zmax:2.5, basisfile:'QSO_NMF_basis_z080_280_norm2150.fits'}
   4: strtmp = {zmin:2.5, zmax:4.7, basisfile:'QSO_NMF_basis_z200_479_norm1420.fits'}
   else: strtmp = {zmin:0.1, zmax:0.6, basisfile:'QSO_NMF_basis_z000_100_norm4150.fits'}
endcase

return, strtmp

end
