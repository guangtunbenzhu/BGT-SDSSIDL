function jhusdss_montecarlo_wmin_blank, nz

if (~keyword_set(nz)) then begin
    waveminmax = jhusdss_sdsswave_minmax()
    ;; Make a  common redshift grid
    zmin = waveminmax[0]/2796D0-1.
    zmax = waveminmax[1]/2803D0-1.
    dz = 0.0005 ;; change to 0.0005, 02/01/2012
    nz = floor((zmax-zmin)/dz)
endif

strtmp = {plate:0L, fiber:0L, zgrid:fltarr(nz), $
          isitcovered:bytarr(nz), $
          rewmin_mgii_2796:fltarr(nz), rewmin_mgii_2803:fltarr(nz), $
          signalmin_mgii_2796:fltarr(nz), signalmin_mgii_2803:fltarr(nz)}

return, strtmp
end
