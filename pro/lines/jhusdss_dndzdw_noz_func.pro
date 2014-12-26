function jhusdss_dndzdw_noz_func, xw, p

n_star_strong = p[0] 
w_star_strong = p[1] 
 
model = n_star_strong*exp(-xw/w_star_strong)
return, model 
end 
