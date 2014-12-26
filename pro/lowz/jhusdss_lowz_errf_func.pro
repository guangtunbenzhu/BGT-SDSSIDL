function jhusdss_lowz_errf_func, x, p
   return, -0.5*erf((x-p[0])/p[1])+1.5
end
