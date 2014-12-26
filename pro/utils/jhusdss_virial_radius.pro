function jhusdss_virial_radius, mass, redshift=redshift, Omega_m=Omega_m, Omega_L=Omega_L

  if (n_elements(redshift) eq 0) then redshift=fltarr(n_elements(mass))+0.1
  if (n_elements(Omega_m) eq 0) then Omega_m=0.3
  if (n_elements(Omega_L) eq 0) then Omega_L=0.7

  H0 = 70. ;;km s-1 Mpc-1
  G0 = 4.3D-9 ;; Mpc km2 s-2 MSun-1
  rho_crit = 3.*H0^2/8./!dpi/G0
  ;; This is not correct!
  xx = Omega_m*(1.+redshift)^3/(Omega_m*(1.+redshift)^3+Omega_L)
  factor = 3./4./!dpi/rho_crit/(18.*!dpi*!dpi+82.*xx-39.*xx^2)*(1.+xx)

  return, (factor*mass)^(1./3.)

end
