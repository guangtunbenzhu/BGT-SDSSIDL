function jhusdss_dndzdw_func, xw, xz, p
n_star = p[0] 
w_star_strong = p[1] 
alpha_strong = p[2] 
bbeta_strong = p[3] 
n_star_weak = p[4] 
w_star_weak = p[5] 
alpha_weak = p[6] 
bbeta_weak = p[7] 
 
model = n_star_strong*(1.+xz)^alpha_strong*exp(-xw/w_star_strong/(1.+xz)^bbeta_strong) $ 
      + n_star_weak*(1.+xz)^alpha_weak*exp(-xw/w_star_weak/(1.+xz)^bbeta_weak) 
return, model 
end 
