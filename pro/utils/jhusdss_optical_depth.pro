;+
;tao = jhusdss_optical_depth(nl, flu, lambda, bb=bb)
;-
function jhusdss_optical_depth, nl, flu, lambda, bb=bb
if (n_elements(bb) eq 0) then bb=10.
return, 0.7580*(nl/1D13)*(flu/0.4164)*(lambda/1215.7)*(10./bb)
end
