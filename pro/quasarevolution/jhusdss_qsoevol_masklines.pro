;+
; vacuum wavelength
;-

function jhusdss_qsoevol_masklines

   name = ['NaI_5898',  'NaI_5892',  'CaII_3970',  'CaII_3935', $ ;4
           'MgI_2853',  'MgII_2803', 'MgII_2796',  'FeII_2600',  'FeII_2586',  $ ;5
           'FeII_2383', 'FeII_2374', 'FeII_2344',  'AlIII_1863', 'AlIII_1855', $ ;5
           'AlII_1671', 'FeII_1608', 'CIV_1551',   'CIV_1548',   'SiII_1527',  'SiIV_1403',  $ ;6
           'SiIV_1394', 'NiII_1370', 'CII_1336',   'SiII_1304', $ ;4
           'SiII_1260', 'Lyal_1216', 'SiIII_1206', 'SiII_1193',  'SiIII_1190', 'PII_1153', $ ;6
           'FeII_1145', 'OVI_1038',  'OVI_1032',   'Lybe_1026', $ ; 4
           'Lyga_972'] ; 1
   wave = [5897.56,     5891.58,     3969.59,     3934.77, $ ;4
           2852.96,     2803.53,     2796.35,     2600.17,      2586.55, $ ;5
           2382.77,     2374.46,     2344.21,     1862.79,      1854.72, $ ;5
           1670.79,     1608.45,     1550.78,     1548.20,      1526.71,      1402.77, $ ;6
           1393.76,     1370.13,     1335.66,     1304.37, $ ;4
           1260.42,     1215.67,     1206.50,     1193.29,      1190.42,      1152.82, $ ;6
           1144.94,     1037.62,     1031.93,     1025.72, $ ;4
           972.54] ;1
   ;; Tao (optical depth), scaled, need to be revisited
   flux = [0.1,         0.1,         0.1,         0.1,  $ ; 4
           0.97,        6.99,        10.00,       5.23,         3.01, $ ;5
           3.98,        1.55,        3.98,        0.97,         1.55, $ ;5
           3.98,        0.1,         3.98,        5.23,         3.98,         3.01, $ ;6
           3.98,        3.98,        3.98,        3.98, $ ; 4
           3.98,        3.98,        3.98,        3.98,         3.98,         3.98, $ ;6
           3.98,        3.98,        3.98,        3.98, $ ; 4
           3.98] ; 1

lines = replicate({name:'Lines', wave:0.D, sigma:3.D, flux:0.D, threshold:1.D}, n_elements(wave))
lines.name = name
lines.wave = wave
lines.flux = flux
if (n_elements(threshold) gt 0) then lines.threshold = threshold

return, lines
end
