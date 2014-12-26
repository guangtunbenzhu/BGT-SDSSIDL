;pro jhusdss_grb_plot

nmfver = 107
path = jhusdss_get_path(/nmfqso)+'/'+string(nmfver, format='(I3.3)')+'/Absorbers'
infile = path+'/ALLQSO_GRB_MC.fits'

choice_load_data = 0
read,'load data? [1=yes]: ',choice_load_data
if choice_load_data eq 1 then begin

  a = mrdfits(infile, 1)
  nabs_i = total(a.nabs_i, 1)
  nabs_f = total(a.nabs_f, 1)

  grbfile = getenv('ALL_DATA')+'/GRB/'+'gzfunc.dat'
  readcol, grbfile, z, gz_f, gz_i, format='(f,f,f)'
  xmodel = z
  nz = n_elements(z)
  deltaz = median(z[1:nz-1] - z[0:nz-2])

  dndz = jhusdss_montecarlo_dndzdw_readin(106)
  w_min_tmp = 1.0
  ymodel = dndz.f0*(1.+xmodel)^(dndz.f0_alpha+dndz.w0_alpha)/(1.+(xmodel/dndz.f0_bbeta)^dndz.f0_ggamma)*dndz.w0/(1.+(xmodel/dndz.w0_bbeta)^dndz.w0_ggamma)*exp(-w_min_tmp/dndz.w0*(1.+(xmodel/dndz.w0_bbeta)^dndz.w0_ggamma)/(1.+xmodel)^dndz.w0_alpha)
  iz = where(xmodel gt 0.4 and xmodel lt 2.0, nz)
  nf_model = total((ymodel*gz_f)[iz])
  ni_model = total((ymodel*gz_i)[iz])
endif

k_print, filename='grb_mc.ps', axis_char_scale=1.4, xsize=14, ysize=8
  !p.multi=[0,2,1]
  !x.margin=0
  plothist, nabs_f, xra=[0.1,30], xst=1, yra=[0,1200], $
     thick=8, xthick=6, ythick=6, xtitle='N(absorbers)', $
     ytitle='Number of MC realization', charsize=1.3, charthick=3.0
  djs_oplot, [1,1]*nf_model, !y.crange, thick=10, color='red'
  djs_xyouts, !x.crange[0]+0.6*(!x.crange[1]-!x.crange[0]), $
      !y.crange[0]+0.9*(!y.crange[1]-!y.crange[0]), $
      'Final Sample', charsize=1.5, charthick=3

  plothist, nabs_i, xra=[0.1,30], xst=1, yra=[0,1200], $
     thick=8, xthick=6, ythick=6, xtitle='N(absorbers)', $
     ytickformat='(A1)', charsize=1.3, charthick=3.0
  djs_oplot, [1,1]*ni_model, !y.crange, thick=10, color='red'
  djs_xyouts, !x.crange[0]+0.55*(!x.crange[1]-!x.crange[0]), $
      !y.crange[0]+0.9*(!y.crange[1]-!y.crange[0]), $
      'Independent Sample', charsize=1.5, charthick=3
k_end_print

end
