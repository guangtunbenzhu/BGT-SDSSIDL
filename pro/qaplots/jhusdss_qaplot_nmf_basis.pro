;+
; Documentation Needed!
;-
pro jhusdss_qaplot_nmf_basis, nmfver

path=jhusdss_get_path(/nmfqso)+'/'+string(nmfver, format='(I3.3)')
if (n_elements(qafile) eq 0) then $
    qafile='QSO_NMF_basis_qaplot_'+string(nmfver, format='(I3.3)')+'.ps'

psfile = path+'/QAplots/'+qafile
; 11-30-2014
psfile = repstr(psfile, '.ps', '_new.ps')

infile = path+'/'+['QSO_NMF_basis_z000_100_norm4150.fits', 'QSO_NMF_basis_z040_179_norm3020.fits', 'QSO_NMF_basis_z080_280_norm2150.fits', 'QSO_NMF_basis_z200_479_norm1420.fits']

basis = mrdfits(infile[1], 1)
n_dimension = (size(basis.eigen_vectors))[1]
thick=6
charthick=3.0
charsize=2.5

;ivector = [8, 4, 2, 1, 3, 0, 5, 6, 7, 9, 10, 11]
ivector = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11] ;8, 4, 2, 1, 3, 0, 5, 6, 7, 9, 10, 11]
nvector = n_elements(ivector)

k_print, filename=psfile, axis_char_scale=1.3, xsize=10, ysize=10
   !p.multi = [0, 1, nvector]
   !y.margin = 0
   !x.margin = 0
   xra = [1400., 5500.]
   yra = [0., 7.]
   xtitle = textoidl('Rest-frame \lambda (\AA)')
   title = textoidl('NMF basis set at 0.4<z<1.8')
   nwave = n_elements(basis.wave)
   iwave = floor(findgen(nwave/5)*5)
   for i=0L, nvector-1L do begin
       djs_plot, basis.wave[iwave], basis.eigen_vectors[ivector[i], iwave], xtickformat='(A1)', $
           ytickformat='(A1)', thick=thick, xthick=thick, ythick=thick, $
           charthick=charthick, charsize=charsize, xra=xra, yra=yra, ytickinterval=10.
       if (keyword_set(labels)) then begin
       if i eq 0L then begin
          djs_xyouts, !x.crange[0]+0.3*(!x.crange[1]-!x.crange[0]),$
                      !y.crange[0]+1.08*(!y.crange[1]-!y.crange[0]),$
                      title, charthick=3.0, charsize=1.5
          djs_xyouts, 1550.*(1.+0.04), !y.crange[0]+0.73*(!y.crange[1]-!y.crange[0]), $
              'C_{ }IV', charsize=1.2, charthick=3.0, color='dark green'
          djs_xyouts, 1909.*(1.-0.04), !y.crange[0]+0.65*(!y.crange[1]-!y.crange[0]), $
              'C_{ }III]', charsize=1.2, charthick=3.0, color='dark green'
          djs_xyouts, 2800.*(1.-0.04), !y.crange[0]+0.70*(!y.crange[1]-!y.crange[0]), $
              'Mg_{ }II', charsize=1.2, charthick=3.0, color='dark green'
       endif
       if i eq 1L then begin
          djs_xyouts, 1550.*(1.+0.02), !y.crange[0]+0.73*(!y.crange[1]-!y.crange[0]), $
              'C_{ }IV', charsize=1.2, charthick=3.0, color='dark green'
          djs_xyouts, 1909.*(1.+0.02), !y.crange[0]+0.50*(!y.crange[1]-!y.crange[0]), $
              'C_{ }III]', charsize=1.2, charthick=3.0, color='dark green'
          djs_xyouts, 2800.*(1.-0.04), !y.crange[0]+0.25*(!y.crange[1]-!y.crange[0]), $
              'Mg_{ }II', charsize=1.2, charthick=3.0, color='dark green'
       endif
       if i eq 2L then begin
          djs_xyouts, 1550.*(1.+0.04), !y.crange[0]+0.75*(!y.crange[1]-!y.crange[0]), $
              'C_{ }IV', charsize=1.2, charthick=3.0, color='dark green'
          djs_xyouts, 1909.*(1.-0.04), !y.crange[0]+0.20*(!y.crange[1]-!y.crange[0]), $
              'C_{ }III]', charsize=1.2, charthick=3.0, color='dark green'
          djs_xyouts, 2800.*(1.-0.04), !y.crange[0]+0.20*(!y.crange[1]-!y.crange[0]), $
              'Mg_{ }II', charsize=1.2, charthick=3.0, color='dark green'
       endif


       if i eq 3L then begin
          djs_xyouts, 5007.*(1.+0.01), !y.crange[0]+0.73*(!y.crange[1]-!y.crange[0]), $
              '[O_{ }III]', charsize=1.2, charthick=3.0, color='dark green'
          djs_xyouts, 4959.*(1.-0.06), !y.crange[0]+0.60*(!y.crange[1]-!y.crange[0]), $
              '[O_{ }III]', charsize=1.3, charthick=3.0, color='dark green'
          djs_xyouts, 3870.*(1.-0.02), !y.crange[0]+0.37*(!y.crange[1]-!y.crange[0]), $
              '[Ne_{ }III]', charsize=1.1, charthick=3.0, color='dark green'
          djs_xyouts, 3727.*(1.-0.05), !y.crange[0]+0.42*(!y.crange[1]-!y.crange[0]), $
              '[O_{ }II]', charsize=1.1, charthick=3.0, color='dark green'
          djs_xyouts, 3426.*(1.-0.06), !y.crange[0]+0.32*(!y.crange[1]-!y.crange[0]), $
              '[Ne_{ }V]', charsize=1.1, charthick=3.0, color='dark green'
       endif
       if i eq 4L then begin
          djs_xyouts, 4863.*(1.-0.02), !y.crange[0]+0.52*(!y.crange[1]-!y.crange[0]), $
              'H\beta', charsize=1.2, charthick=3.0, color='dark green'
          djs_xyouts, 4341.*(1.-0.02), !y.crange[0]+0.32*(!y.crange[1]-!y.crange[0]), $
              'H\gamma', charsize=1.2, charthick=3.0, color='dark green'
          djs_xyouts, 4102.*(1.-0.02), !y.crange[0]+0.32*(!y.crange[1]-!y.crange[0]), $
              'H\delta', charsize=1.2, charthick=3.0, color='dark green'
       endif

       endif

       if i eq 8 then $
          djs_xyouts, !x.crange[0]-0.04*(!x.crange[1]-!x.crange[0]), $
                      !y.crange[0]+0.6*(!y.crange[1]-!y.crange[0]), $
                      'Flux (Arbitray Unit)', orientation=90, $
                      charsize=2, charthick=charthick
   endfor
   djs_axis, xaxis=0, xtitle=xtitle, charthick=charthick, charsize=charsize
k_end_print

end
