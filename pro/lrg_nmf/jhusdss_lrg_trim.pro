function jhusdss_lrg_trim, lrg, sncut=sncut

if (n_elements(sncut) eq 0) then sncut = 3.0

    itrim = where((lrg.sn_median gt sncut $
              and  lrg.mass gt 10.7 and lrg.ssfr lt -0.40*(lrg.mass-11.)-11.4 $
              and  lrg.mass lt 12.1 and lrg.ssfr gt -13.) $
               or (lrg.sn_median gt sncut and lrg.z gt 0.32999 $
              and  lrg.z lt 0.6), ntrim)

return, itrim

end
