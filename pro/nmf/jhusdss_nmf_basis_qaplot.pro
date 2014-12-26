;+
; Documentation Needed!
;-
pro jhusdss_nmf_basis_qaplot, basis, trainingset, qapath=qapath, qafile=qafile

if (n_elements(qapath) eq 0) then qapath='.'
if (n_elements(qafile) eq 0) then qafile='QSO_NMF_basis_qaplot.ps'
n_dimension = (size(basis.eigen_vectors))[1]
thick=3
charthick=3
charsize=1.5
k_print, filename=qapath+'/'+qafile, axis_char_scale=1.3, xsize=8, ysize=12
   !p.multi = [0, 1, n_dimension]
   !y.margin = 0
   xtitle = textoidl('Wavelength (\AA)')
   for i=0L, n_dimension-1L do $
       djs_plot, basis.wave, basis.eigen_vectors[i, *], xtickformat='(A1)', thick=thick, $
           xthick=thick, ythick=thick, charthick=charthick, charsize=charsize, xra=minmax(basis.wave)
   djs_axis, xaxis=0, xtitle=xtitle, charthick=charthick, charsize=charsize

   xtitle = textoidl('a (coefficient)')
   for i=0L, n_dimension-1L do $
       plothist, trainingset.eigen_values[*,i], bin=0.02, peak=1, xtickformat='(A1)', thick=thick, $
          xthick=xthick, ythick=ythick, charthick=charthick, charsize=charsize, xra=[0,1]
   djs_axis, xaxis=0, xtitle=xtitle, charthick=charthick, charsize=charsize
k_end_print

end
