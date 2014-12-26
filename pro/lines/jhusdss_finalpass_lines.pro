function  jhusdss_finalpass_lines

   ;; singlets first, then doublets
   name = ['MgI_2853', 'FeII_2344', 'AlII_1671', 'SiII_1527', $
            'MgII_2803', 'MgII_2796', 'FeII_2600', 'FeII_2586', $
            'FeII_2383', 'FeII_2374', 'AlIII_1863', 'AlIII_1855', $
            'CIV_1551', 'CIV_1548']
   wave = [2852.96,     2344.21,     1670.79,     1526.71,     $
           2803.53,     2796.35,     2600.17,     2586.55,     $
           2382.77,     2374.46,     1862.79,     1854.72,     $
           1550.78,     1548.20]

lines = replicate({name:'Lines', wave:0.D}, n_elements(wave))
lines.name = name
lines.wave = wave

return, lines
end
