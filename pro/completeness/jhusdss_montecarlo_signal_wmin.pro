function jhusdss_montecarlo_signal_wmin, signal
   return, 6.-sqrt((5.0-(signal<5.0))*36./5.0)
end
