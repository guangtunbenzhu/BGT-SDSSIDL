;+
;w = jhusdss_rew(nl, flu, lambda, bb=bb)
; 
;-
function jhusdss_rew, nl, flu, lambda
return, nl*flu*lambda^2/1.13D20
end
