function jhusdss_dndzdw_noz_func2, xw, p

n_star_strong = p[0] 
w_star_strong = p[1] 
n_star_weak = p[2] 
w_star_weak = p[3] 
 
model = n_star_strong*exp(-xw/w_star_strong) $ 
      + n_star_weak*exp(-xw/w_star_weak) 
return, model 
end 
