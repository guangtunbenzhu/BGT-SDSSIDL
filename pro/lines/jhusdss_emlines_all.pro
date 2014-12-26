function  jhusdss_emlines_all, mask=mask

if (~keyword_set(mask)) then begin
   name = ['FeII_2964', 'MgII_2800', 'FeII_2325',  'FeIII_2077',  'CIII_1906', $
           'FeII_1664', 'CIV_1546',  'SiIV_1398',  'OI_1305',     'NV_1240',   $
           'Lya_1216']
   wave = [2964.28,     2800.26,     2324.58,     2076.62,      1905.97, $
           1664.74,     1546.15,     1398.33,     1305.42,      1239.85, $
           1216.25]
   ;; Tao (optical depth), scaled
   flux = [2.02,        14.72,        2.01,       1.58,         15.94, $
           0.48,        25.29,        8.92,       1.99,         2.46,  $
           100.]
   sigma = [22.92,      34.95,        22.23,      16.99,        23.58, $
            5.50,       14.33,        12.50,      5.42,         2.71,  $
            19.46]
endif else begin
   name = ['CIV_1546']
   wave = [1546.15]
   flux = [100.]
   sigma = [19.46]
endelse

lines = replicate({name:'Lines', wave:0.D, sigma:3.D, flux:0.D}, n_elements(wave))
lines.name = name
lines.sigma = sigma
lines.wave = wave
lines.flux = flux

return, lines
end
