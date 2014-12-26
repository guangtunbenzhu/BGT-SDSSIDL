;; This should in Python

pro jhusdss_montecarlo_completeness_plot, nmfver, boss=boss
 
;; path
if (~keyword_set(path)) then begin
   if (keyword_set(boss)) then begin
       path=jhusdss_get_path(/nmfqso)+'/'+$
          string(nmfver, format='(I3.3)')+'/MonteCarlo_BOSS'
   endif else begin
       path=jhusdss_get_path(/nmfqso)+'/'+$
          string(nmfver, format='(I3.3)')+'/MonteCarlo'
   endelse
endif

filename = jhusdss_montecarlo_completeness_filename(nmfver)
infile = path+'/'+filename

comp = mrdfits(infile, 1)

charsize = 1.5
charthick = 2
thick=7

qapath = path+'/QAplot'
psfile = qapth+'/'+repstr(filename, '.fits', '.ps')

k_print, filename=psfile, axis_char_scale=1.3, xsize=8, ysize=8

rimage = comp.ndetected/comp.ncovered
nw_rgb_make, rimage, rimage,rimage, colors=rcolors, name='tmp.jpg'
k_end_print

end
