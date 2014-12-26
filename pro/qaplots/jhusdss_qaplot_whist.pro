FUNCTION tick_exponent, axis, index, number

     ; A special case.
     IF number EQ 0 THEN RETURN, '0'

     ; Assuming multiples of 10 with format.
     ex = String(number, Format='(e8.0)')
     pt = StrPos(ex, '.')

     first = StrMid(ex, 0, pt)
     sign = StrMid(ex, pt+2, 1)
     thisExponent = StrMid(ex, pt+3)

     ; Shave off leading zero in exponent
     WHILE StrMid(thisExponent, 0, 1) EQ '0' DO thisExponent = StrMid(thisExponent, 1)

     ; Fix for sign and missing zero problem.
     IF (Long(thisExponent) EQ 0) THEN BEGIN
        sign = ''
        thisExponent = '0'
     ENDIF

     ; Make the exponent a superscript.
     IF sign EQ '-' THEN BEGIN
;       RETURN, first + 'x10!U' + sign + thisExponent + '!N'
        RETURN, '10!U' + sign + thisExponent + '!N'
     ENDIF ELSE BEGIN
;       RETURN, first + 'x10!U' + thisExponent + '!N'
        RETURN, '10!U' + thisExponent + '!N'
     ENDELSE

END

pro jhusdss_qaplot_whist, nmfver, train=train

absorbers = jhusdss_absorber_readin(nmfver, train=train, /byabs, /trim)
;; read in monte carlo completeness
mc = jhusdss_montecarlo_readin(nmfver)
;; post-processing mc
f_log10w_z = float(mc.ndetected)/float((mc.ncovered+(mc.ncovered eq 0)))
nlog10w = n_elements(mc.log10w)

;; for strong absorbers
istrong = where(mc.log10w_min gt alog10(4.), nstrong)
ftmp = total(mc.ndetected[istrong,*], 1)/total(mc.ncovered[istrong,*], 1)
for i=0L, nstrong-1L do f_log10w_z[istrong[i], *] = ftmp

;; for weak absorbers we extrapolate
iweak = where(mc.log10w_max le alog10(0.30), nweak)
ftmp = total(mc.ndetected[iweak,*], 1)/total(mc.ncovered[iweak,*], 1)
for i=0L, nweak-1L do $
    for j=0L, n_elements(f_log10w_z[iweak[i],*])-1L do $
        f_log10w_z[iweak[i], j] = $
        interpol(f_log10w_z[iweak[nweak-1]+1:*, j], mc.log10w[iweak[nweak-1]+1:*], $
                 mc.log10w[iweak[i]])

;; get f for each absorber
log10w_abs = alog10(absorbers.rew_mgii_2796)
w_abs = absorbers.rew_mgii_2796
z_abs = absorbers.zabs
;; bilinear interpolation
;; virtual subscript
nabs = n_elements(absorbers)
v_log10w_abs = (log10w_abs-min(mc.log10w_min))/mc.dlog10w
v_z_abs = (z_abs-min(mc.z_min))/mc.dz
f_abs = fltarr(nabs)
for i=0L, nabs-1L do f_abs[i] = bilinear(f_log10w_z, v_log10w_abs[i], v_z_abs[i])

;; init
thick=8
charsize=1.4
charthick=3
;; path
path=jhusdss_get_path(/nmfqso)+'/'+string(nmfver, format='(I3.3)')
qapath = path+'/QAplots'
if (jhusdss_direxist(qapath) eq 0) then spawn, 'mkdir -p '+qapath

;; w(2796) hist
psfile = qapath+'/'+'w2796_hist_corr_'+string(nmfver, format='(I3.3)')+'.ps'
k_print, filename=psfile, axis_char_scale=1.5, xsize=8, ysize=8

  xra = [-0.2, 7]
  yra = [0.8, 60000]
; xtitle = textoidl('log_{10}(W_0^{2796}) (\AA)')
  xtitle = textoidl('W^{\lambda2796}_0 (\AA)')
  ytitle = textoidl('\DeltaN (\DeltaW^{\lambda2796}_0=0.1 \AA)')
  x = absorbers.rew_mgii_2796
  ii = where(x gt 0.1)

  plothist, x[ii], xra=xra, yra=yra, bin=0.10, $
      thick=thick, xthick=thick, ythick=thick, $
      xtitle=xtitle, ytitle=ytitle, $
      charsize=charsize, charthick=charthick, $
      xstyle=1, ystyle=1, /ylog, ytickformat='tick_exponent'

delta_w = 0.10
w_min = 0.10
w_max = 7.0
w_bin = jhusdss_make_bins(w_min, w_max, delta_w, nbin=w_nbin)
w_median = fltarr(w_nbin)
n_all = fltarr(w_nbin)
for iw=0L, w_nbin-1L do begin
    this_w_min = w_bin[iw].min
    this_w_max = w_bin[iw].max
    tt = where(w_abs gt this_w_min and w_abs lt this_w_max, ss)
    w_median[iw] = w_bin[iw].mean
    if (ss gt 0) then begin
       n_all[iw] = total(1./(f_abs[tt]+(f_abs[tt] le 0.)))
    endif
endfor

  djs_oplot, [0., w_median], [0., n_all], $
      psym=10, color='red', thick=thick+2, linestyle=0

  items = ['Observed Distribution', 'Completeness-corrected']
  colors = [djs_icolor('black'), djs_icolor('red')]
  legend, items, textcolor=colors, color=colors, linestyle=[0,0], $
    charsize=1.5, charthick=charthick, thick=thick, $
    box=0, /right, /top, pspacing=1.6

k_end_print


stop
end
