;pro jhusdss_mgii_lineratio

;read in Churchill+1999
readcol, '/home/menard/DATA/SDSS/Absorber/MgII/Churchill_weak.dat', cc_zabs, cc_qso, cc_wv, cc_err_wv, cc_W2796, cc_err_W2796, cc_dr, cc_err_dr, cc_z_w_dr, format='F,A,F,F,F,F,F,F,F'
cc_wtot = cc_w2796*(1.+1./cc_dr)
jj = where(cc_dr gt 0.1 and cc_dr lt 5)
cc_wtot = cc_wtot[jj]
cc_dr = cc_dr[jj]

;read in Narayanan+2007
readcol, '/home/menard/DATA/SDSS/Absorber/MgII/Narayanan_weak.dat', an_qso, an_zabs, an_W2796, an_err_W2796, an_w2803, an_err_w2803, an_dr, an_err_dr, an_z_w_dr, format='A,F,F,F,F,F,F,F,F'
an_wtot = an_w2796*(1.+1./an_dr)
jj = where(an_dr gt 0.1 and an_dr lt 5)
an_wtot = an_wtot[jj]
an_dr = an_dr[jj]

;read in JHUSDSS
absorbers = jhusdss_absorber_readin(107)
zm_wtot = absorbers.rew_mgii_2796+absorbers.rew_mgii_2803
zm_dr = absorbers.rew_mgii_2796/absorbers.rew_mgii_2803
jj = where(zm_dr gt 0.1 and zm_dr lt 5)
zm_wtot = zm_wtot[jj]
zm_dr = zm_dr[jj]

all_wtot = [cc_wtot, an_wtot, zm_wtot]
all_dr = [cc_dr, an_dr, zm_dr]

w_min = 0.02
w_max = 6
dlogw = 0.15

w_bins = jhusdss_make_bins(alog10(w_min), alog10(w_max), dlogw)
n_bins = n_elements(w_bins)
w_bins[0].min = 0.
w_mean = fltarr(n_bins)
dr_mean = fltarr(n_bins)
dr_sdev = fltarr(n_bins)

for i=0L, n_elements(w_bins)-1L do begin
    w_left = 10.^w_bins[i].min
    w_right = 10.^w_bins[i].max
    ii = where(all_wtot gt w_left and all_wtot le w_right, nn)
;   izm = where(zm_wtot gt w_left and zm_wtot le w_right, nzm)
;   if ncc gt 0 and nzm gt 0 then begin
;       all_wtot = [cc_wtot[icc], zm_wtot[izm]]
;       all_dr = [cc_dr[icc], zm_dr[izm]]
;   endif else begin
;      if ncc gt 0 then
;         all_wtot = cc_wtot[icc]
;         all_dr = cc_dr[icc]
;      endif
;      if nzm gt 0 then
;         all_wtot = zm_wtot[izm]
;         all_dr = zm_dr[izm]
;      endif
;   endelse
    if (nn gt 0) then begin
       w_mean[i] = median(all_wtot[ii])
       dr_mean[i] = median(all_dr[ii])
       tmp = moment(all_dr[ii], sdev=sdev)
       dr_sdev[i] = sdev
    endif
endfor

djs_plot, w_mean, dr_mean, psym=4, /xlog, xra=[0.01, 4], xst=1, yra=[0,3], yst=1
oploterror, w_mean, dr_mean, dr_sdev, psym=4
outstr = {w_mean:w_mean, dr_mean:dr_mean, dr_sdev:dr_sdev}
mwrfits, outstr, '/home/menard/DATA/SDSS/Absorber/MgII/mgii_lineratio.fits', /create

end
