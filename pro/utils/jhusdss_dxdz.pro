;; 
;; in units of c/H0
;;
function jhusdss_dxdz, z, omegam=omegam, omegal=omegal
   if (n_elements(omegam) eq 0) then omegam=0.3
   if (n_elements(omegal) eq 0) then omegal=0.7
   return, (1+z)^2/sqrt(omegal+omegam*(1+z)^3)
end
