function  jhusdss_train_lines

   name = ['MgII_2803', 'MgII_2796', 'FeII_2600',  'FeII_2586', 'FeII_2383', $
           'FeII_2344', 'CIV_1551']
   wave = [2803.53,     2796.35,     2600.17,      2586.65,     2382.77,     $
           2344.21,     1550.78]

lines = replicate({name:'Lines', wave:0.D}, n_elements(wave))
lines.name = name
lines.wave = wave

return, lines
end
