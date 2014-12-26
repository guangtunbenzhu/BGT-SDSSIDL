pro jhusdss_lowz_deadalive, mass, ssfr, index_dead=index_dead, index_alive=index_alive, index_use=index_use

    a = -0.3
    b = -10.9
    index_dead = where(mass gt 7. and mass lt 13. and ssfr gt -13. and ssfr lt -7 $
              and ssfr lt a*(mass-10.5)+b, ndead)
    index_alive = where(mass gt 7. and mass lt 13. and ssfr gt -13. and ssfr lt -7 $
              and ssfr ge a*(mass-10.5)+b, nalive)
    index_use = where(mass gt 7. and mass lt 13. and ssfr gt -13. and ssfr lt -7, nall)
end
