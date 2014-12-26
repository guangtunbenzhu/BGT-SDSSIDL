;+
; Mask out all the emission lines
; 30-Apr-2010 Guangtun, JHU, Adopted to jhusdss
; ??-???-2010 Guangtun, NYU, es0_lines
;-
function jhusdss_lrg_masklines

   strtmp = {line:0., name:' '}

lines= [6564.6127, $
        5008.1666, $
        3727.0898, $
        4862.6778, $
        4960.2140, $
        4341.6803, $
        6585.1583, $
        3869.7867, $
        4105.8884, $
        6718.1642, $
        7137.6370, $
        3971.1933, $
        2800.0, $
        6732.5382, $
        6549.7689]
names=['Halpha', $
       'OIII5007', $
       'OII3727', $
       'Hbeta', $
       'OIII4959', $
       'Hgamma', $
       'NII6584', $
       'NeIII3869', $
       'Hdelta', $
       'SII6717', $
       'ArIII7137', $
       'Hepsilon', $
       'MgII', $
       'SII6732', $
       'NII6549']
  outstr = replicate(strtmp, n_elements(lines))
  outstr.line = lines
  outstr.name = names

  return, outstr
end
