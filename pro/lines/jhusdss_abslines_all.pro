function  jhusdss_abslines_all, mgii=mgii, feii=feii, train=train, forplot=forplot, ew=ew

if (~keyword_set(mgii) and ~keyword_set(feii) and ~keyword_set(train) $
 and ~keyword_set(forplot) and ~keyword_set(ew)) then begin
   name = ['CaII_3970', 'CaII_3935', $
           'MgI_2853',  'MgII_2803', 'MgII_2796', 'FeII_2600',  'FeII_2586',  $
           'FeII_2383', 'FeII_2374', 'FeII_2344', 'AlIII_1863', 'AlIII_1855', $
           'AlII_1671', 'FeII_1608', 'CIV_1551',  'CIV_1548',  'SiII_1527',  'SiIV_1403',  $
           'SiIV_1394', 'NiII_1370']
   wave = [3969.59,     3934.77, $
           2852.96,     2803.53,     2796.35,     2600.17,      2586.55, $
           2382.77,     2374.46,     2344.21,     1862.79,      1854.72, $
           1670.79,     1608.45,     1550.78,     1548.20,     1526.71,      1402.77, $
           1393.76,     1370.13]
   ;; Tao (optical depth), scaled
   flux = [0.1,         0.1, $
           0.97,        6.99,        10.00,       5.23,         3.01, $
           3.98,        1.55,        3.98,        0.97,         1.55, $
           3.98,        0.1,         3.98,        5.23,        3.98,         3.01, $
           3.98,        3.98]
endif

if (keyword_set(forplot)) then begin
   name = ['CaII_3970', 'CaII_3935', $
           'MgI_2853',  'MgII_2803', 'MgII_2796', 'FeII_2600',  'FeII_2586',  $
           'FeII_2383', 'FeII_2374', 'FeII_2344', 'AlIII_1863', 'AlIII_1855', $
           'AlII_1671', 'FeII_1608', 'CIV_1551',  'CIV_1548',  'SiII_1527',  'SiIV_1403',  $
           'SiIV_1394', 'CII_1335', 'SiII_1304', 'SiII_1260']
   wave = [3969.59,     3934.77, $
           2852.96,     2803.53,     2796.35,     2600.17,      2586.55, $
           2382.77,     2374.46,     2344.21,     1862.79,      1854.72, $
           1670.79,     1608.45,     1555.78,     1544.20,     1526.71,      1404.77, $
           1391.76,     1335.71,     1304.37,     1260.42]
   ;; Tao (optical depth), scaled
   flux = [0.1,         0.1, $
           0.97,        6.99,        10.00,       5.23,         3.01, $
           3.98,        1.55,        3.98,        0.97,         1.55, $
           3.98,        0.1,         3.98,        5.23,        3.98,         3.01, $
           3.98,        5.23,        5.23,        5.23]
endif


if (keyword_set(mgii)) then begin
   name = ['MgII_2803', 'MgII_2796']
   wave = [2803.53,     2796.35]
   flux = [6.99,        10.00]
endif

if (keyword_set(feii)) then begin
   name = ['FeII_2600',  'FeII_2586',  'FeII_2383', 'FeII_2374', 'FeII_2344']
   wave = [2600.17,      2586.55,      2382.77,     2374.46,     2344.21]
   flux = [5.23,         3.01,         3.98,        1.55,        3.98]
endif

;; the first six are used to find absorbers,
;; the ones (only CIV for now) after are used to resolve line confusion
if (keyword_set(train)) then begin
   name = ['MgII_2803', 'MgII_2796', 'FeII_2600',  'FeII_2586', 'FeII_2383', $
           'FeII_2344', 'CIV_1551']
   wave = [2803.53,     2796.35,     2600.17,      2586.55,     2382.77,     $
           2344.21,     1550.78]
   flux = [6.99,        10.00,       5.23,         3.01,        3.98,        $
           3.98,        10.00]
   threshold = [2.0,     2.0,         1.0,         1.0,         1.0, $
                1.0,     2.0]
;  name = ['MgII_2803', 'MgII_2796', 'FeII_2600',  'FeII_2383', 'FeII_2344', 'FeII_2586', 'FeII_2374']
;  wave = [2803.53,     2796.35,     2600.17,      2382.77,     2344.21,      2586.55,      2374.46]
;  flux = [6.99,        10.00,       5.23,         33.98,       3.98,         3.01,         1.55]
;  threshold = [2.0,     2.0,         1.0,         1.0,         1.0,          1.2,          0.8]
endif

if (keyword_set(ew)) then begin
   name = ['CaII_3970', 'CaII_3935', $
           'MgI_2853',  'MgII_2803', 'MgII_2796', 'FeII_2600',  'FeII_2586',  $
           'FeII_2383', 'FeII_2374', 'FeII_2344', 'AlIII_1863', 'AlIII_1855', $
           'AlII_1671', 'FeII_1608', 'CIV_1551',  'CIV_1548',  'SiII_1527',  'SiIV_1403',  $
           'SiIV_1394', 'NiII_1370']
   wave = [3969.59,     3934.77, $
           2852.96,     2803.53,     2796.35,     2600.17,      2586.55, $
           2382.77,     2374.46,     2344.21,     1862.79,      1854.72, $
           1670.79,     1608.45,     1550.78,     1548.20,     1526.71,      1402.77, $
           1393.76,     1370.13]
   ;; Tao (optical depth), scaled
   flux = [0.1,         0.1, $
           0.97,        6.99,        10.00,       5.23,         3.01, $
           3.98,        1.55,        3.98,        0.97,         1.55, $
           3.98,        0.1,         3.98,        5.23,        3.98,         3.01, $
           3.98,        3.98]

   name = ['MgII_2803', 'MgII_2796', 'FeII_2600',  'FeII_2586', 'FeII_2383', $
           'FeII_2344', 'CIV_1551']
   wave = [2803.53,     2796.35,     2600.17,      2586.55,     2382.77,     $
           2344.21,     1550.78]
   flux = [6.99,        10.00,       5.23,         3.01,        3.98,        $
           3.98,        10.00]
   threshold = [2.0,     2.0,         1.0,         1.0,         1.0, $
                1.0,     2.0]
endif

lines = replicate({name:'Lines', wave:0.D, sigma:3.D, flux:0.D, threshold:1.D}, n_elements(wave))
lines.name = name
lines.wave = wave
lines.flux = flux
if (n_elements(threshold) gt 0) then lines.threshold = threshold

return, lines
end
