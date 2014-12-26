function jhusdss_dndzdw_func3, xw, xz, p
f0 = p[0]
f0_alpha = p[1]
f0_bbeta = p[2]
f0_ggamma = p[3]
w0 = p[4]
w0_alpha = p[5]
w0_bbeta = p[6]
w0_ggamma = p[7]

;model = f0*(f0_alpha+xz)/(1.+(xz/f0_bbeta)^f0_ggamma)*exp(-xw/w0/(w0_alpha+xz)*(1.+(xz/w0_bbeta)^w0_ggamma))
;model = f0*(1.+xz)/(1.+(xz/f0_bbeta)^f0_ggamma)*exp(-xw/w0/(1.+xz)^w0_alpha*(1.+(xz/w0_bbeta)^w0_ggamma))
;model = f0*(1.+xz)^f0_alpha/(1.+(xz/f0_bbeta)^f0_ggamma)*exp(-xw/w0/(1.+xz)^w0_alpha*(1.+(xz/w0_bbeta)^w0_ggamma))
model = f0*(1.+xz)^f0_alpha/(1.+(xz/f0_bbeta)^f0_ggamma)*exp(-(xw)/w0/(1.+xz)^w0_alpha*(1.+(xz/w0_bbeta)^w0_ggamma))
;model = f0*(1.+xz)^f0_alpha/(1.+(xz/f0_bbeta))^f0_ggamma*exp(-xw/w0/(1.+xz)^w0_alpha*(1.+(xz/w0_bbeta))^w0_ggamma)
return, model
end
