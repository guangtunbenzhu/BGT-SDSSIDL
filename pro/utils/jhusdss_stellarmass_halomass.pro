;; Behroozi 2010, 0<z<1
function jhusdss_stellarmass_halomass, mass, redshift=redshift, Omega_m=Omega_m, Omega_L=Omega_L

  if (n_elements(redshift) eq 0) then redshift=fltarr(n_elements(mass))+0.1
  if (n_elements(Omega_m) eq 0) then Omega_m=0.3
  if (n_elements(Omega_L) eq 0) then Omega_L=0.7

  a = 1./(1.+redshift)
  M00 = 10.72
  M0a = 0.55
  M10 = 12.35
  M1a = 0.28
  beta_0 = 0.44
  beta_a = 0.18
  delta_0 = 0.57
  delta_a = 0.17
  gamma_0 = 1.56
  gamma_a = 2.51

  log10_M1 = M10 + M1a*(a-1.)
  log10_M0 = M00 + M0a*(a-1.)
  betaA = beta_0 + beta_a*(a-1.)
  deltaA = delta_0 + delta_a*(a-1.)
  gammaA = gamma_0 + gamma_a*(a-1.)

  M_to_M0 = mass/10.^log10_M0
  halomass = log10_M1+betaA*(alog10(mass)-log10_M0)+(M_to_M0^deltaA)/(1.+M_to_M0^(-gammaA))-0.5

  return, halomass
end
