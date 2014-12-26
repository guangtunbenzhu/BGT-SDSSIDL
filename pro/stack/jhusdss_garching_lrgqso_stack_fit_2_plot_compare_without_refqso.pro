
;pro jhusdss_lowz_caii_new

Coarse = 0b
DoWeight = 0b
DoNaI = 0b
DoScale = 0b
Scale_factor = 0.2 ;; M* \propto M(halo)^1 \propto r(halo)^3 and r(halo) \propto M*^(1/3); Scale to 10^10.1 Msun
fid_mass = 10.3 ;; M* \propto M(halo)^1 \propto r(halo)^3 and r(halo) \propto M*^(1/3); Scale to 10^10.1 Msun
read,'Coarse? [1=yes, 0=no]: ', Coarse
read,'DoWeight? [1=yes, 0=no]: ', DoWeight
read,'DoNaI? [1=yes, 0=no]: ', DoNaI
read,'DoScale? [1=yes, 0=no]: ', DoScale

nmfver = 106
overwrite = 1b
;read,'BOSS? [1=yes, 0=no]: ',BOSS

line_wave = [2796.35, 2803.53]
xra = [2600., 3050.]
;xra = [2650., 2950.]

fix_separation = long((line_wave[1]-line_wave[0])/(line_wave[1]+line_wave[0])/alog(10.)*2E+4)

stackpath = jhusdss_get_path(/nmfqso)+'/'+string(nmfver, format='(I3.3)')+'/Stack/'
infile = stackpath + jhusdss_garching_lrgqso_stack_filename(nmfver, boss=boss)
shuffle_infile = repstr(infile, '.fits', '_shuffle_mean.fits')

if (Coarse) then begin
   infile=repstr(infile,'.fits','_coarse.fits')
endif
if (DoWeight) then begin
   infile=repstr(infile,'.fits','_weight.fits')
endif
if (DoNaI) then begin
   infile=repstr(infile,'.fits','_NaI.fits')
endif
if (DoScale) then begin
   infile=repstr(infile,'.fits','_Scale.fits')
endif

psfile = repstr(infile, '.fits', '_fit_jhusdss_garching_lrgqso_stack_fit_2.ps')

comp = mrdfits(infile,1)
comp_shuffle = mrdfits(shuffle_infile,1)
nrp = n_elements(comp)

i_indep = [1, 7, 13, 17, 25, 29]

thick=8
xthick=8
ythick=8
charsize=1.4
charthick=3

xra = [2650., 2950.]
xtitle='\lambda (\AA)' 
;xtitle1='(1-R_{\lambda,d})/\sigma(R_{\lambda,d})' 
xtitle1='(1-<R_c>)/\sigma(<R_c>)' 
;xtitle1='(1-!3'+string("303B)+'!XF_\lambda)/\sigma_F' 
;ytitle='Normalized Flux F_\lambda'
ytitle='Normalized Flux <R>';_\lambda'
;title='Single Gaussian Line Profile Measurement'
xra1=[-4.9,6.1]
yra1=[0,1.1]

xx = 10.^(findgen(1000)*0.0024+1.)

plotsym, 0, /fill

k_print, filename=psfile, axis_char_scale=1.3, xsize=6, ysize=3

nrp=6

;pos = [0.15, 0.9-0.8/nrp*(0+1), 0.95, 0.9-0.8/nrp*0]
;pos1 = [0.70, 0.9-0.8/nrp*(0+1), 0.95, 0.9-0.8/nrp*0]

;djs_plot, xra, [1., 1.], psym=10, xra=xra, yra=[1.-0.001, 1.+0.001], $
;    position=pos, /nodata, xtickformat='(A1)', ytickformat='(A1)', $ 
;    thick=2, xst=5, yst=5, xthick=xthick, ythick=ythick
pos = [0.175, 0.20, 0.95, 0.9]

for i=3L, nrp-1 do begin
    rp = comp[i_indep[i]].rp
    wave = comp[i_indep[i]].wave
    tmp = min(abs(comp[i_indep[i]].wave-line_wave[0]), icaii)
    y_ori = comp[i_indep[i]].fluxgeomean
    y_ref = comp_shuffle[i_indep[i]].fluxgeomean
    y = comp[i_indep[i]].fluxgeomean/comp_shuffle[i_indep[i]].fluxgeomean

    ; page 1: Original stack
    print, rp
    djs_plot, wave, smooth(y_ori, 5), xra=xra, xst=1, yra=[1.-0.015, 1.+0.005], $
        thick=5, xthick=xthick, ythick=ythick, color='blue', yminor=2, ytickinterval=interval, $
        xtitle=xtitle, ytitle=ytitle, pos=pos
    ; page 2: Original stack
    djs_plot, wave, smooth(y_ori, 5), xra=xra, xst=1, yra=[1.-0.015, 1.+0.005], $
        thick=5, xthick=xthick, ythick=ythick, color='blue', yminor=2, ytickinterval=interval, $
        xtitle=xtitle, ytitle=ytitle, pos=pos
    djs_oplot, wave, smooth(y_ref, 5), color='red', thick=5
    ; page 3: Original stack
    djs_plot, wave, smooth(y, 5), xra=xra, xst=1, yra=[1.-0.005, 1.+0.005], $
        thick=5, xthick=xthick, ythick=ythick, color='black', yminor=2, ytickinterval=interval, $
        xtitle=xtitle, ytitle=ytitle, pos=pos
    djs_oplot, xra, [1.,1], color='gray', linestyle=2, thick=5
endfor

k_end_print

end
