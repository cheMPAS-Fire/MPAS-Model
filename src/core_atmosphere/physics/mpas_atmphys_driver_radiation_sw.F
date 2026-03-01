! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!=================================================================================================================
 module mpas_atmphys_driver_radiation_sw
 use mpas_kind_types
 use mpas_pool_routines
 use mpas_timer,only: mpas_timer_start,mpas_timer_stop

 use mpas_atmphys_constants
 use mpas_atmphys_manager, only: gmt,curr_julday,julday,year
 use mpas_atmphys_camrad_init
 use mpas_atmphys_rrtmg_swinit
 use mpas_atmphys_vars
 
!wrf physics:
 use module_mp_thompson_aerosols
 use module_ra_rrtmg_sw_aerosols
 use module_ra_cam
 use module_ra_rrtmg_sw
 use mpas_smoke_wrapper 


 implicit none
 private
 public:: allocate_radiation_sw,   &
          deallocate_radiation_sw, &
          driver_radiation_sw,     &
          init_radiation_sw,       &
          radconst

!MPAS driver for parameterization of shortwave radiation codes.
!Laura D. Fowler (send comments to laura@ucar.edu).
!2013-05-01.
!
! subroutines in mpas_atmphys_driver_radiation_sw:
! ------------------------------------------------
! allocate_radiation_sw  : allocate local arrays for parameterization of sw radiation codes.
! deallocate_radiation_sw: deallocate local arrays for parameterization of sw radiation codes.
! init_radiation_sw      : initialization of individual sw radiation codes.
! driver_radiation_sw    : main driver (called from subroutine physics_driver).
! radiation_sw_from_MPAS : initialize local arrays.
! radiation_sw_to_MPAS   : copy local arrays to MPAS arrays.
! radconst               : calculate solar declination,...
!
! WRF physics called from driver_radiation_sw:
! --------------------------------------------
! * module_ra_cam        : CAM short wave radiation code.
! * module_ra_rrtmg_sw   : RRTMG short wave radiation code.
!
! add-ons and modifications to sourcecode:
! ----------------------------------------
! * removed the pre-processor option "do_hydrostatic_pressure" before call to subroutines
!   rrtmg_sw and camrad.
!   Laura D. Fowler (laura@ucar.edu) / 2013-05-29.
! * added structure diag in the call to subroutine init_radiation_sw and call to subroutine
!   camradinit for initialization of variable mxaerl.
!   Laura D. Fowler (laura@ucar.edu) / 2013-07-01.
! * modified the call to subroutine rrtmg_swrad to include the option of using the same ozone
!   climatology as the one used in the CAM radiation codes.
!   Laura D. Fowler (laura@ucar.edu) / 2013-07-17.
! * in call to subroutine rrtmg_swrad, replaced the variable g (that originally pointed to
!   gravity) with gravity, for simplicity.
!   Laura D. Fowler (laura@ucar.edu) / 2014-03-21.
! * throughout the sourcecode, replaced all "var_struct" defined arrays by local pointers.
!   Laura D. Fowler (laura@ucar.edu) / 2014-04-22.
! * modified sourcecode to use pools.
!   Laura D. Fowler (laura@ucar.edu) / 2014-05-15.
! * cleaned-up the call to rrtmg_swrad after cleaning up subroutine rrtmg_swrad in module_ra_rrtmg_sw.F.
!   Laura D. Fowler (laura@ucar.edu) / 2016-07-05.
! * added the cloud radii for cloud water, cloud ice, and snow calculated in the Thompson cloud microphysics
!   scheme in the call to subroutine rrtmg_swrad.
!   Laura D. Fowler (laura@ucar.edu) / 2016-07-07.
! * removed qr_p, and qg_p in the call to rrtmg_swrad since not used in the calculation of the cloud optical
!   properties.
!   Laura D. Fowler (laura@ucar.edu) / 2016-07-08.
! * in the call to rrtmg_swrad, substituted the variables qv_p, qc_p, qi_p, and qs_p with qvrad_p, qcrad_p,
!   qirad_p, and qsrad_p initialized in subroutine cloudiness_from_MPAS.
!   Laura D. Fowler (laura@ucar.edu) / 2016-07-09.
! * in subroutines radiation_sw_from_MPAS, revised the initialization of re_cloud, re_ice, re_snow, to
!   handle the case when the cloud microphysics parameterization is turned off, i.e. config_microp_scheme='off'.
!   Laura D. Fowler (laura@ucar.edu) / 2017-02-10.
! * since we removed the local variable radt_sw_scheme from mpas_atmphys_vars.F, now defines radt_sw_scheme
!   as a pointer to config_radt_sw_scheme.
!   Laura D. Fowler (laura@ucar.edu) / 2017-02-16.
! * added the variables swddir,swddni,swddif for use in the updated version of the Noah LSM.
!   Laura D. Fowler (laura@ucar.edu) / 2023-04-21.
! * removed the variables f_qv and f_qg in the call to subroutine camrad.
!   Laura D. Fowler (laura@ucar.edu) / 2024-02-13.
! * in subroutine radiation_sw_from_MPAS, added the calculation of the optical properties of "water-friendly" and
!   "ice-friendly" aerosols from the Thompson cloud microphysics scheme for use in the RRTMG short-wave radiation
!   code.
!   Laura D. Fowler (laura@ucar.edu) / 2024-05-16.
! * in subroutine driver_radiation_sw, modified the argument list in the call to subroutine rrtmg_sw to include
!   the optical properties of "water-friendly" and "ice-friendly" aerosols from the Thompson cloud microphysics
!   scheme.
!   Laura D. Fowler (laura@ucar.edu) / 2024-05-16.


 contains


!=================================================================================================================
 subroutine allocate_radiation_sw(configs,xtime_s)
!=================================================================================================================

!input arguments:
 type(mpas_pool_type),intent(in):: configs
 real(kind=RKIND),intent(in):: xtime_s

!local pointers:
 character(len=StrKIND),pointer:: mp_scheme,     &
                                  radt_sw_scheme
 logical,pointer:: config_tempo_aerosolaware
 integer,pointer:: aero_dir_rad_fdb
!-----------------------------------------------------------------------------------------------------------------

 call mpas_pool_get_config(configs,'config_microp_scheme' ,mp_scheme     )
 call mpas_pool_get_config(configs,'config_radt_sw_scheme',radt_sw_scheme)
 call mpas_pool_get_config(configs,'config_tempo_aerosolaware',config_tempo_aerosolaware)
 call mpas_pool_get_config(configs,'aero_dir_rad_fdb',aero_dir_rad_fdb)
      
 if(.not.allocated(f_ice)        ) allocate(f_ice(ims:ime,kms:kme,jms:jme)        )
 if(.not.allocated(f_rain)       ) allocate(f_rain(ims:ime,kms:kme,jms:jme)       )

 if(.not.allocated(xlat_p)       ) allocate(xlat_p(ims:ime,jms:jme)               )
 if(.not.allocated(xlon_p)       ) allocate(xlon_p(ims:ime,jms:jme)               )

 if(.not.allocated(sfc_albedo_p) ) allocate(sfc_albedo_p(ims:ime,jms:jme)         )
 if(.not.allocated(snow_p)       ) allocate(snow_p(ims:ime,jms:jme)               )
 if(.not.allocated(tsk_p)        ) allocate(tsk_p(ims:ime,jms:jme)                )
 if(.not.allocated(xice_p)       ) allocate(xice_p(ims:ime,jms:jme)               )
 if(.not.allocated(xland_p)      ) allocate(xland_p(ims:ime,jms:jme)              )

 if(.not.allocated(coszr_p)      ) allocate(coszr_p(ims:ime,jms:jme)              )
 if(.not.allocated(gsw_p)        ) allocate(gsw_p(ims:ime,jms:jme)                )
 if(.not.allocated(swcf_p)       ) allocate(swcf_p(ims:ime,jms:jme)               )
 if(.not.allocated(swdnb_p)      ) allocate(swdnb_p(ims:ime,jms:jme)              )
 if(.not.allocated(swdnbc_p)     ) allocate(swdnbc_p(ims:ime,jms:jme)             )
 if(.not.allocated(swdnbcln_p)   ) allocate(swdnbcln_p(ims:ime,jms:jme)           ) !Clean sky testing
 if(.not.allocated(swdnt_p)      ) allocate(swdnt_p(ims:ime,jms:jme)              )
 if(.not.allocated(swdntc_p)     ) allocate(swdntc_p(ims:ime,jms:jme)             )
 if(.not.allocated(swdntcln_p)   ) allocate(swdntcln_p(ims:ime,jms:jme)           ) !Clean sky testing
 if(.not.allocated(swupb_p)      ) allocate(swupb_p(ims:ime,jms:jme)              )
 if(.not.allocated(swupbc_p)     ) allocate(swupbc_p(ims:ime,jms:jme)             )
 if(.not.allocated(swupbcln_p)   ) allocate(swupbcln_p(ims:ime,jms:jme)           ) !Clean sky testing
 if(.not.allocated(swupt_p)      ) allocate(swupt_p(ims:ime,jms:jme)              )
 if(.not.allocated(swuptc_p)     ) allocate(swuptc_p(ims:ime,jms:jme)             )
 if(.not.allocated(swuptcln_p)   ) allocate(swuptcln_p(ims:ime,jms:jme)           ) !Clean sky testing
 
 if(.not.allocated(rthratensw_p) ) allocate(rthratensw_p(ims:ime,kms:kme,jms:jme) )
! if(.not.allocated(rthratenswc_p) ) allocate(rthratenswc_p(ims:ime,kms:kme,jms:jme) )

 radiation_sw_select: select case (trim(radt_sw_scheme))
    case("rrtmg_sw")
       if(.not.allocated(recloud_p)    ) allocate(recloud_p(ims:ime,kms:kme,jms:jme)       )
       if(.not.allocated(reice_p)      ) allocate(reice_p(ims:ime,kms:kme,jms:jme)         )
       if(.not.allocated(resnow_p)     ) allocate(resnow_p(ims:ime,kms:kme,jms:jme)        )

       if(.not.allocated(alswvisdir_p) ) allocate(alswvisdir_p(ims:ime,jms:jme)            )
       if(.not.allocated(alswvisdif_p) ) allocate(alswvisdif_p(ims:ime,jms:jme)            )
       if(.not.allocated(alswnirdir_p) ) allocate(alswnirdir_p(ims:ime,jms:jme)            )
       if(.not.allocated(alswnirdif_p) ) allocate(alswnirdif_p(ims:ime,jms:jme)            )
       if(.not.allocated(swvisdir_p)   ) allocate(swvisdir_p(ims:ime,jms:jme)              )
       if(.not.allocated(swvisdif_p)   ) allocate(swvisdif_p(ims:ime,jms:jme)              )
       if(.not.allocated(swnirdir_p)   ) allocate(swnirdir_p(ims:ime,jms:jme)              )
       if(.not.allocated(swnirdif_p)   ) allocate(swnirdif_p(ims:ime,jms:jme)              )

       if(.not.allocated(swddir_p)     ) allocate(swddir_p(ims:ime,jms:jme)                )
       if(.not.allocated(swddni_p)     ) allocate(swddni_p(ims:ime,jms:jme)                )
       if(.not.allocated(swddif_p)     ) allocate(swddif_p(ims:ime,jms:jme)                )

       if(.not.allocated(swdnflx_p)    ) allocate(swdnflx_p(ims:ime,kms:kme+1,jms:jme)     )
       if(.not.allocated(swdnflxc_p)   ) allocate(swdnflxc_p(ims:ime,kms:kme+1,jms:jme)    )
       if(.not.allocated(swdnflxcln_p) ) allocate(swdnflxcln_p(ims:ime,kms:kme+1,jms:jme)  ) !Clean sky testing
       if(.not.allocated(swupflx_p)    ) allocate(swupflx_p(ims:ime,kms:kme+1,jms:jme)     )
       if(.not.allocated(swupflxc_p)   ) allocate(swupflxc_p(ims:ime,kms:kme+1,jms:jme)    )
       if(.not.allocated(swupflxcln_p) ) allocate(swupflxcln_p(ims:ime,kms:kme+1,jms:jme)  ) !Clean sky testing

       if(.not.allocated(pin_p)        ) allocate(pin_p(num_oznlevels)                     )
       if(.not.allocated(o3clim_p)     ) allocate(o3clim_p(ims:ime,1:num_oznlevels,jms:jme))

       if(.not.allocated(tauaer_p)     ) allocate(tauaer_p(ims:ime,kms:kme,jms:jme,nbndsw) )
       if(.not.allocated(ssaaer_p)     ) allocate(ssaaer_p(ims:ime,kms:kme,jms:jme,nbndsw) )
       if(.not.allocated(asyaer_p)     ) allocate(asyaer_p(ims:ime,kms:kme,jms:jme,nbndsw) )
       if(.not.allocated(taod5503d_sim_p)) allocate(taod5503d_sim_p(ims:ime,kms:kme,jms:jme)   )
       if(.not.allocated(taod5502d_sim_p)) allocate(taod5502d_sim_p(ims:ime,kms:kme)   )
       if (aero_dir_rad_fdb > 0 .and. aero_dir_rad_fdb /= 2) then
         if(.not.allocated(ht_p)       ) allocate(ht_p(ims:ime,jms:jme)       )
         if(.not.allocated(taer_type_p)) allocate(taer_type_p(ims:ime,jms:jme))
       endif

       aerosol_select: select case(mp_scheme)
          case("mp_thompson_aerosols")
!             if(.not.allocated(ht_p)       ) allocate(ht_p(ims:ime,jms:jme)       )
!             if(.not.allocated(taer_type_p)) allocate(taer_type_p(ims:ime,jms:jme))
             if(.not.allocated(taod5502d_p)) allocate(taod5502d_p(ims:ime,jms:jme))
             if(.not.allocated(taod5503d_p)) allocate(taod5503d_p(ims:ime,kms:kme,jms:jme))
          case("mp_tempo")
             if (config_tempo_aerosolaware) then
!                if(.not.allocated(ht_p)       ) allocate(ht_p(ims:ime,jms:jme)       )
!                if(.not.allocated(taer_type_p)) allocate(taer_type_p(ims:ime,jms:jme))
                if(.not.allocated(taod5502d_p)) allocate(taod5502d_p(ims:ime,jms:jme))
                if(.not.allocated(taod5503d_p)) allocate(taod5503d_p(ims:ime,kms:kme,jms:jme))
             endif
          case default
       end select aerosol_select

    case("cam_sw")
       if(.not.allocated(glw_p)        ) allocate(glw_p(ims:ime,jms:jme)                )
       if(.not.allocated(lwcf_p)       ) allocate(lwcf_p(ims:ime,jms:jme)               )
       if(.not.allocated(lwdnb_p)      ) allocate(lwdnb_p(ims:ime,jms:jme)              )
       if(.not.allocated(lwdnbc_p)     ) allocate(lwdnbc_p(ims:ime,jms:jme)             )
       if(.not.allocated(lwdnt_p)      ) allocate(lwdnt_p(ims:ime,jms:jme)              )
       if(.not.allocated(lwdntc_p)     ) allocate(lwdntc_p(ims:ime,jms:jme)             )
       if(.not.allocated(lwupb_p)      ) allocate(lwupb_p(ims:ime,jms:jme)              )
       if(.not.allocated(lwupbc_p)     ) allocate(lwupbc_p(ims:ime,jms:jme)             )
       if(.not.allocated(lwupt_p)      ) allocate(lwupt_p(ims:ime,jms:jme)              )
       if(.not.allocated(lwuptc_p)     ) allocate(lwuptc_p(ims:ime,jms:jme)             )
       if(.not.allocated(olrtoa_p)     ) allocate(olrtoa_p(ims:ime,jms:jme)             )
       if(.not.allocated(sfc_emiss_p)  ) allocate(sfc_emiss_p(ims:ime,jms:jme)          )
       if(.not.allocated(rthratenlw_p) ) allocate(rthratenlw_p(ims:ime,kms:kme,jms:jme) )

       if(.not.allocated(cemiss_p)     ) allocate(cemiss_p(ims:ime,kms:kme,jms:jme)     )
       if(.not.allocated(taucldc_p)    ) allocate(taucldc_p(ims:ime,kms:kme,jms:jme)    )
       if(.not.allocated(taucldi_p)    ) allocate(taucldi_p(ims:ime,kms:kme,jms:jme)    )

       if(.not.allocated(pin_p)        ) allocate(pin_p(num_oznlevels)                  )
       if(.not.allocated(ozmixm_p) ) &
          allocate(ozmixm_p(ims:ime,1:num_oznlevels,jms:jme,num_months) )
       
       if(.not.allocated(m_hybi_p)     ) allocate(m_hybi_p(num_aerlevels)               )
       if(.not.allocated(m_psn_p)      ) allocate(m_psn_p(ims:ime,jms:jme)              )
       if(.not.allocated(m_psp_p)      ) allocate(m_psp_p(ims:ime,jms:jme)              )
       if(.not.allocated(aerosolcn_p)  ) &
          allocate(aerosolcn_p(ims:ime,1:num_aerlevels,jms:jme,num_aerosols) )
       if(.not.allocated(aerosolcp_p)  ) &
          allocate(aerosolcp_p(ims:ime,1:num_aerlevels,jms:jme,num_aerosols) )

       !allocate these arrays on the first time step, only:
       if(xtime_s .lt. 1.e-12) then

          if(.not.allocated(emstot_p) ) allocate(emstot_p(ims:ime,kms:kme,jms:jme) )
          if(.not.allocated(abstot_p) ) &
             allocate(abstot_p(ims:ime,kms:kme,cam_abs_dim2,jms:jme) )
          if(.not.allocated(absnxt_p) ) &
             allocate(absnxt_p(ims:ime,kms:kme,cam_abs_dim1,jms:jme) )

       endif

    case default

 end select radiation_sw_select

 end subroutine allocate_radiation_sw

!=================================================================================================================
 subroutine deallocate_radiation_sw(configs)
!=================================================================================================================

!input arguments:
 type(mpas_pool_type),intent(in):: configs

!local pointers:
 character(len=StrKIND),pointer:: mp_scheme,     &
                                  radt_sw_scheme
 logical,pointer:: config_tempo_aerosolaware
 integer,pointer:: aero_dir_rad_fdb

!-----------------------------------------------------------------------------------------------------------------

 call mpas_pool_get_config(configs,'config_microp_scheme' ,mp_scheme     )
 call mpas_pool_get_config(configs,'config_radt_sw_scheme',radt_sw_scheme)
 call mpas_pool_get_config(configs,'config_tempo_aerosolaware',config_tempo_aerosolaware)     
 call mpas_pool_get_config(configs,'aero_dir_rad_fdb',aero_dir_rad_fdb)

 if(allocated(f_ice)        ) deallocate(f_ice        )
 if(allocated(f_rain)       ) deallocate(f_rain       )
 if(allocated(xlat_p)       ) deallocate(xlat_p       )
 if(allocated(xlon_p)       ) deallocate(xlon_p       )
 if(allocated(sfc_albedo_p) ) deallocate(sfc_albedo_p )
 if(allocated(snow_p)       ) deallocate(snow_p       )
 if(allocated(tsk_p)        ) deallocate(tsk_p        )
 if(allocated(xice_p)       ) deallocate(xice_p       )
 if(allocated(xland_p)      ) deallocate(xland_p      )
 if(allocated(coszr_p)      ) deallocate(coszr_p      )
 if(allocated(gsw_p)        ) deallocate(gsw_p        )
 if(allocated(swcf_p)       ) deallocate(swcf_p       )
 if(allocated(swdnb_p)      ) deallocate(swdnb_p      )
 if(allocated(swdnbc_p)     ) deallocate(swdnbc_p     )
 if(allocated(swdnbcln_p)   ) deallocate(swdnbcln_p   )!Clean sky
 if(allocated(swdnt_p)      ) deallocate(swdnt_p      )
 if(allocated(swdntc_p)     ) deallocate(swdntc_p     )
 if(allocated(swdntcln_p)   ) deallocate(swdntcln_p   )!Clean sky
 if(allocated(swupb_p)      ) deallocate(swupb_p      )
 if(allocated(swupbc_p)     ) deallocate(swupbc_p     )
 if(allocated(swupbcln_p)   ) deallocate(swupbcln_p   )!Clean sky
 if(allocated(swupt_p)      ) deallocate(swupt_p      )
 if(allocated(swuptc_p)     ) deallocate(swuptc_p     )
 if(allocated(swuptcln_p)   ) deallocate(swuptcln_p   )!Clean sky
 
 if(allocated(rthratensw_p) ) deallocate(rthratensw_p )
! if(allocated(rthratenswc_p) ) deallocate(rthratenswc_p )

 radiation_sw_select: select case (trim(radt_sw_scheme))
    case("rrtmg_sw")
       if(allocated(recloud_p)    ) deallocate(recloud_p    )
       if(allocated(reice_p)      ) deallocate(reice_p      )
       if(allocated(resnow_p)     ) deallocate(resnow_p     )

       if(allocated(alswvisdir_p) ) deallocate(alswvisdir_p )
       if(allocated(alswvisdif_p) ) deallocate(alswvisdif_p )
       if(allocated(alswnirdir_p) ) deallocate(alswnirdir_p )
       if(allocated(alswnirdif_p) ) deallocate(alswnirdif_p )
       if(allocated(swvisdir_p)   ) deallocate(swvisdir_p   )
       if(allocated(swvisdif_p)   ) deallocate(swvisdif_p   )
       if(allocated(swnirdir_p)   ) deallocate(swnirdir_p   )
       if(allocated(swnirdif_p)   ) deallocate(swnirdif_p   )

       if(allocated(swddir_p)     ) deallocate(swddir_p     )
       if(allocated(swddni_p)     ) deallocate(swddni_p     )
       if(allocated(swddif_p)     ) deallocate(swddif_p     )

       if(allocated(swdnflx_p)    ) deallocate(swdnflx_p    )
       if(allocated(swdnflxc_p)   ) deallocate(swdnflxc_p   )
       if(allocated(swdnflxcln_p) ) deallocate(swdnflxcln_p )
       if(allocated(swupflx_p)    ) deallocate(swupflx_p    )
       if(allocated(swupflxc_p)   ) deallocate(swupflxc_p   )
       if(allocated(swupflxcln_p) ) deallocate(swupflxcln_p )

       if(allocated(pin_p)        ) deallocate(pin_p        )
       if(allocated(o3clim_p)     ) deallocate(o3clim_p     )

       if(allocated(tauaer_p)     ) deallocate(tauaer_p     )
       if(allocated(ssaaer_p)     ) deallocate(ssaaer_p     )
       if(allocated(asyaer_p)     ) deallocate(asyaer_p     )
       if(allocated(taod5503d_sim_p)) deallocate(taod5503d_sim_p)
       if(allocated(taod5502d_sim_p)) deallocate(taod5502d_sim_p)
       if(aero_dir_rad_fdb > 0 .and. aero_dir_rad_fdb /= 2) then
         if(allocated(ht_p)       ) deallocate(ht_p       )
         if(allocated(taer_type_p)) deallocate(taer_type_p) 
       endif


       aerosol_select: select case(mp_scheme)
          case("mp_thompson","mp_thompson_aerosols")
!             if(allocated(ht_p)       ) deallocate(ht_p       )
!             if(allocated(taer_type_p)) deallocate(taer_type_p)
             if(allocated(taod5502d_p)) deallocate(taod5502d_p)
             if(allocated(taod5503d_p)) deallocate(taod5503d_p)
          case("mp_tempo")
             if	(config_tempo_aerosolaware) then
 !               if(allocated(ht_p)       ) deallocate(ht_p       )
 !               if(allocated(taer_type_p)) deallocate(taer_type_p)
                if(allocated(taod5502d_p)) deallocate(taod5502d_p)
                if(allocated(taod5503d_p)) deallocate(taod5503d_p)
             endif
          case default
       end select aerosol_select

    case("cam_sw")
       if(allocated(pin_p)        ) deallocate(pin_p        )
       if(allocated(m_hybi_p)     ) deallocate(m_hybi_p     )

       if(allocated(xlat_p)       ) deallocate(xlat_p       )
       if(allocated(xlon_p)       ) deallocate(xlon_p       )

       if(allocated(glw_p)        ) deallocate(glw_p        )
       if(allocated(lwcf_p)       ) deallocate(lwcf_p       )
       if(allocated(lwdnb_p)      ) deallocate(lwdnb_p      )
       if(allocated(lwdnbc_p)     ) deallocate(lwdnbc_p     )
       if(allocated(lwdnt_p)      ) deallocate(lwdnt_p      )
       if(allocated(lwdntc_p)     ) deallocate(lwdntc_p     )
       if(allocated(lwupb_p)      ) deallocate(lwupb_p      )
       if(allocated(lwupbc_p)     ) deallocate(lwupbc_p     )
       if(allocated(lwupt_p)      ) deallocate(lwupt_p      )
       if(allocated(lwuptc_p)     ) deallocate(lwuptc_p     )
       if(allocated(olrtoa_p)     ) deallocate(olrtoa_p     )
       if(allocated(sfc_emiss_p)  ) deallocate(sfc_emiss_p  )
       if(allocated(rthratenlw_p) ) deallocate(rthratenlw_p )

       if(allocated(cemiss_p)     ) deallocate(cemiss_p     )
       if(allocated(ozmixm_p)     ) deallocate(ozmixm_p     )
       if(allocated(taucldc_p)    ) deallocate(taucldc_p    )
       if(allocated(taucldi_p)    ) deallocate(taucldi_p    )

       if(allocated(m_psn_p)      ) deallocate(m_psn_p      )
       if(allocated(m_psp_p)      ) deallocate(m_psp_p      )
       if(allocated(aerosolcn_p)  ) deallocate(aerosolcn_p  )
       if(allocated(aerosolcp_p)  ) deallocate(aerosolcp_p  )

    case default

 end select radiation_sw_select

 end subroutine deallocate_radiation_sw

!=================================================================================================================
 subroutine radiation_sw_from_MPAS(configs,mesh,state,time_lev,diag_physics,atm_input, &
                                   sfc_input,xtime_s,its,ite)
!=================================================================================================================
! Initializes local arrays from MPAS pools for use in shortwave radiation schemes.
! Handles:
!   (1) Cloud effective radii for RRTMG (indirect effect via microphysics)
!   (2) Aerosol optical properties for RRTMG direct radiative feedback
!       controlled by two independent namelist options:
!         aer_opt          : RRTMG-level aerosol pathway (0=none, 1=ECMWF, 2=direct AOD, 3=Thompson WFA/IFA)
!         aero_dir_rad_fdb : MPAS driver-level direct feedback source
!                            (0=none, 1=simple AOD, 2=Mie/SMOKE, 3=Thompson internal, 4=Thompson+Mie combined)
!
! Code has been written by Minsu Choi CIRES/NOAA GSL(Minsu.Choi@noaa.gov)/ 2026-02-16
! Atmospheric Composition Branch at NOAA GSL
!=================================================================================================================

!input arguments:
 type(mpas_pool_type),intent(in):: mesh
 type(mpas_pool_type),intent(in):: configs
 type(mpas_pool_type),intent(in):: state
 type(mpas_pool_type),intent(in):: diag_physics
 type(mpas_pool_type),intent(in):: atm_input
 type(mpas_pool_type),intent(in):: sfc_input

 integer,intent(in):: its,ite
 integer,intent(in):: time_lev

 real(kind=RKIND),intent(in):: xtime_s

!local variables:
 integer:: i,j,k,n
 real(kind=RKIND):: tau_thompson, tau_mie, ssa_thompson, ssa_mie, asy_thompson, asy_mie
 real(kind=RKIND):: tau_total, ssa_num, asy_num

!local pointers:
 logical,pointer:: config_o3climatology
 logical,pointer:: config_microp_re
 logical,pointer:: config_tempo_aerosolaware
 integer,pointer:: aero_dir_rad_fdb
 integer,pointer:: aer_opt
 character(len=StrKIND),pointer:: radt_sw_scheme
 character(len=StrKIND),pointer:: microp_scheme

 real(kind=RKIND),dimension(:),pointer    :: latCell,lonCell
 real(kind=RKIND),dimension(:),pointer    :: skintemp,snow,xice,xland
 real(kind=RKIND),dimension(:),pointer    :: m_ps,pin
 real(kind=RKIND),dimension(:),pointer    :: sfc_albedo,sfc_emiss
 real(kind=RKIND),dimension(:),pointer    :: taod5502d
 real(kind=RKIND),dimension(:,:),pointer  :: zgrid
 real(kind=RKIND),dimension(:,:),pointer  :: cldfrac,m_hybi,o3clim
 real(kind=RKIND),dimension(:,:),pointer  :: re_cloud,re_ice,re_snow
 real(kind=RKIND),dimension(:,:),pointer  :: taod5503d
 real(kind=RKIND),dimension(:,:,:),pointer:: aerosols,ozmixm
 real(kind=RKIND),dimension(:,:,:),pointer:: tauaer_sw, ssaaer_sw, asyaer_sw
 real(kind=RKIND),dimension(:,:),pointer  :: aod3d_simple_pool_ptr

 ! Temporary arrays for aero_dir_rad_fdb=4 (Thompson + Mie combination)
 ! This arrays only used in this driver. 
 real(kind=RKIND), allocatable:: tauaer_thompson_tmp(:,:,:,:)
 real(kind=RKIND), allocatable:: ssaaer_thompson_tmp(:,:,:,:)
 real(kind=RKIND), allocatable:: asyaer_thompson_tmp(:,:,:,:)

!-----------------------------------------------------------------------------------------------------------------

 call mpas_pool_get_config(configs,'config_microp_re'        ,config_microp_re        )
 call mpas_pool_get_config(configs,'config_o3climatology'    ,config_o3climatology    )
 call mpas_pool_get_config(configs,'config_radt_sw_scheme'   ,radt_sw_scheme          )
 call mpas_pool_get_config(configs,'config_microp_scheme'    ,microp_scheme           )
 call mpas_pool_get_config(configs,'config_tempo_aerosolaware',config_tempo_aerosolaware)
 call mpas_pool_get_config(configs,'aero_dir_rad_fdb'        ,aero_dir_rad_fdb        )
 call mpas_pool_get_config(configs,'aer_opt'                 ,aer_opt                 )

 call mpas_pool_get_array(mesh,'latCell',latCell)
 call mpas_pool_get_array(mesh,'lonCell',lonCell)

 call mpas_pool_get_array(state,'aerosols',aerosols,time_lev)

 call mpas_pool_get_array(sfc_input,'skintemp',skintemp)
 call mpas_pool_get_array(sfc_input,'snow'    ,snow    )
 call mpas_pool_get_array(sfc_input,'xice'    ,xice    )
 call mpas_pool_get_array(sfc_input,'xland'   ,xland   )

 call mpas_pool_get_array(atm_input,'pin'  ,pin  )
 call mpas_pool_get_array(atm_input,'ozmixm',ozmixm)

 call mpas_pool_get_array(diag_physics,'sfc_albedo',sfc_albedo)
 call mpas_pool_get_array(diag_physics,'sfc_emiss' ,sfc_emiss )
 call mpas_pool_get_array(diag_physics,'cldfrac'   ,cldfrac   )
 call mpas_pool_get_array(diag_physics,'o3clim'    ,o3clim    )
 call mpas_pool_get_array(diag_physics,'m_hybi'    ,m_hybi    )
 call mpas_pool_get_array(diag_physics,'m_ps'      ,m_ps      )

 call mpas_pool_get_array(diag_physics,'tauaer_sw',tauaer_sw,time_lev)
 call mpas_pool_get_array(diag_physics,'ssaaer_sw',ssaaer_sw,time_lev)
 call mpas_pool_get_array(diag_physics,'asyaer_sw',asyaer_sw,time_lev)
 call mpas_pool_get_array(diag_physics,'aod3d_simple',aod3d_simple_pool_ptr)

!-----------------------------------------------------------------------------------------------------------------
! Section 1: Copy basic surface and atmospheric fields
!-----------------------------------------------------------------------------------------------------------------
 do j = jts,jte
 do i = its,ite
    xlat_p(i,j)       = latCell(i) / degrad
    xlon_p(i,j)       = lonCell(i) / degrad
    sfc_albedo_p(i,j) = sfc_albedo(i)
    snow_p(i,j)       = snow(i)
    tsk_p(i,j)        = skintemp(i)
    xice_p(i,j)       = xice(i)
    xland_p(i,j)      = xland(i)
 enddo
 enddo

 do j = jts,jte
 do k = kts,kte
 do i = its,ite
    cldfrac_p(i,k,j) = cldfrac(k,i)
    f_ice(i,k,j)     = 0.0_RKIND
    f_rain(i,k,j)    = 0.0_RKIND
 enddo
 enddo
 enddo

!-----------------------------------------------------------------------------------------------------------------
! Section 2: Initialize output SW radiation arrays to zero
!-----------------------------------------------------------------------------------------------------------------
 do j = jts,jte
 do i = its,ite
    coszr_p(i,j)    = 0.0_RKIND
    gsw_p(i,j)      = 0.0_RKIND
    swcf_p(i,j)     = 0.0_RKIND
    swdnb_p(i,j)    = 0.0_RKIND  ;  swdnbc_p(i,j)   = 0.0_RKIND  ;  swdnbcln_p(i,j) = 0.0_RKIND
    swdnt_p(i,j)    = 0.0_RKIND  ;  swdntc_p(i,j)   = 0.0_RKIND  ;  swdntcln_p(i,j) = 0.0_RKIND
    swupb_p(i,j)    = 0.0_RKIND  ;  swupbc_p(i,j)   = 0.0_RKIND  ;  swupbcln_p(i,j) = 0.0_RKIND
    swupt_p(i,j)    = 0.0_RKIND  ;  swuptc_p(i,j)   = 0.0_RKIND  ;  swuptcln_p(i,j) = 0.0_RKIND
    swddir_p(i,j)   = 0.0_RKIND  ;  swddni_p(i,j)   = 0.0_RKIND  ;  swddif_p(i,j)   = 0.0_RKIND
 enddo
 do k = kts,kte
 do i = its,ite
    rthratensw_p(i,k,j) = 0.0_RKIND
 enddo
 enddo
 enddo

!-----------------------------------------------------------------------------------------------------------------
! Section 3: RRTMG-specific initialization
!-----------------------------------------------------------------------------------------------------------------
 radiation_sw_select: select case (trim(radt_sw_scheme))

    case("rrtmg_sw")

       !--- Initialize aerosol optical property arrays to no-aerosol defaults
       !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j,k,n)
       do n = 1,nbndsw
       do j = jts,jte
       do k = kts,kte
       do i = its,ite
          tauaer_p(i,k,j,n) = 0.0_RKIND
          ssaaer_p(i,k,j,n) = 1.0_RKIND
          asyaer_p(i,k,j,n) = 0.0_RKIND
       enddo
       enddo
       enddo
       enddo
       !$OMP END PARALLEL DO

       !--- Initialize flux arrays
       do j = jts,jte
       do k = kts,kte+2
       do i = its,ite
          swdnflx_p(i,k,j)    = 0.0_RKIND  ;  swdnflxc_p(i,k,j)   = 0.0_RKIND  ;  swdnflxcln_p(i,k,j) = 0.0_RKIND
          swupflx_p(i,k,j)    = 0.0_RKIND  ;  swupflxc_p(i,k,j)   = 0.0_RKIND  ;  swupflxcln_p(i,k,j) = 0.0_RKIND
       enddo
       enddo
       enddo

       !----------------------------------------------------------------------
       ! Section 3a: Cloud effective radii (indirect effect via microphysics)
       ! Independent of aero_dir_rad_fdb
       !----------------------------------------------------------------------
       microp_re_select: select case(microp_scheme)
          case("mp_thompson","mp_thompson_aerosols","mp_wsm6","mp_nssl2m","mp_tempo")
             if(config_microp_re) then
                call mpas_pool_get_array(diag_physics,'re_cloud',re_cloud)
                call mpas_pool_get_array(diag_physics,'re_ice'  ,re_ice  )
                call mpas_pool_get_array(diag_physics,'re_snow' ,re_snow )
                do j = jts,jte
                do k = kts,kte
                do i = its,ite
                   recloud_p(i,k,j) = re_cloud(k,i)
                   reice_p(i,k,j)   = re_ice(k,i)
                   resnow_p(i,k,j)  = re_snow(k,i)
                enddo
                enddo
                enddo
             else
                do j = jts,jte
                do k = kts,kte
                do i = its,ite
                   recloud_p(i,k,j) = 0._RKIND
                   reice_p(i,k,j)   = 0._RKIND
                   resnow_p(i,k,j)  = 0._RKIND
                enddo
                enddo
                enddo
             endif
          case default
       end select microp_re_select

       !----------------------------------------------------------------------
       ! Section 3b: Aerosol direct radiative feedback
       ! Controlled by aero_dir_rad_fdb (MPAS driver-level)
       ! Options:
       !   0 = no direct feedback (arrays remain at no-aerosol defaults)
       !   1 = simple AOD parameterization (diagnostic RH-based)
       !   2 = Mie calculation from SMOKE driver
       !   3 = Thompson internal pathway (WFA/IFA climatology)
       !   4 = Thompson WFA/IFA + Mie combined (optical weighting)
       !----------------------------------------------------------------------
       direct_feedback_select: select case(aero_dir_rad_fdb)

          !--- Case 0: No direct feedback
          !    Arrays already initialized to no-aerosol defaults above
          case(0)
             ! nothing to do

          !--- Case 1: Simple AOD parameterization
          !--- For code developer, to understand the simplest approach
          !--- We recommend go to /physics_noaa/SMOKE/mpas_smoke_diagnostics.F90
          case(1)
             call mpas_pool_get_array(mesh,'zgrid',zgrid)

             !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j,k)
             do j = jts,jte
             do i = its,ite
                ht_p(i,j) = zgrid(1,i)
                if(xland_p(i,j) == 1._RKIND) then
                   taer_type_p(i,j) = 1
                else
                   taer_type_p(i,j) = 3
                endif
                taod5502d_sim_p(i,j) = 0._RKIND
                do k = kts,kte
                   taod5503d_sim_p(i,k,j) = aod3d_simple_pool_ptr(k,i)
                   taod5502d_sim_p(i,j)   = taod5502d_sim_p(i,j) + taod5503d_sim_p(i,k,j)
                enddo
             enddo
             enddo
             !$OMP END PARALLEL DO

             call calc_aerosol_rrtmg_sw( &
                      ht       = ht_p,           dz8w     = dz_p,           &
                      p        = pres_hyd_p,      t3d      = t_p,            &
                      qv3d     = qv_p,                                       &
                      tauaer   = tauaer_p,        ssaaer   = ssaaer_p,       &
                      asyaer   = asyaer_p,                                   &
                      aod5503d = taod5503d_sim_p, aod5502d = taod5502d_sim_p,&
                      aer_type = taer_type_p,                                &
                      aer_aod550_opt = 1, aer_angexp_opt = 3,                &
                      aer_ssa_opt    = 3, aer_asy_opt    = 3,                &
                      aer_aod550_val = 0., aer_angexp_val = 0.,              &
                      aer_ssa_val    = 0., aer_asy_val    = 0.,              &
                      ims=ims, ime=ime, jms=jms, jme=jme, kms=kms, kme=kme, &
                      its=its, ite=ite, jts=jts, jte=jte, kts=kts, kte=kte  )

          !--- Case 2: Mie calculation from SMOKE driver
          case(2)
             !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j,k,n)
             do n = 1,nbndsw
             do j = jts,jte
             do k = kts,kte
             do i = its,ite
                tauaer_p(i,k,j,n) = tauaer_sw(n,k,i)
                ssaaer_p(i,k,j,n) = ssaaer_sw(n,k,i)
                asyaer_p(i,k,j,n) = asyaer_sw(n,k,i)
             enddo
             enddo
             enddo
             enddo
             !$OMP END PARALLEL DO

          !--- Case 3: Thompson internal pathway (WFA/IFA from climatology only)
          !    Valid for mp_thompson_aerosols and mp_tempo with aerosolaware=.true.
          case(3)
             call mpas_pool_get_array(mesh,'zgrid',zgrid)
             call mpas_pool_get_array(diag_physics,'taod5502d',taod5502d)
             call mpas_pool_get_array(diag_physics,'taod5503d',taod5503d)

             if(microp_scheme == 'mp_thompson_aerosols' .or. &
                (microp_scheme == 'mp_tempo' .and. config_tempo_aerosolaware)) then

                do j = jts,jte
                do i = its,ite
                   ht_p(i,j) = zgrid(1,i)
                   if(xland_p(i,j)==1._RKIND) then
                      taer_type_p(i,j) = 1
                   elseif(xland_p(i,j)==2._RKIND) then
                      taer_type_p(i,j) = 3
                   endif
                enddo
                enddo

                call gt_aod( &
                        p_phy = pres_hyd_p , dz8w = dz_p  , t_phy     = t_p        , qvapor = qv_p , &
                        nwfa  = nwfa_p     , nifa  = nifa_p, taod5503d = taod5503d_p,                 &
                        ims=ims, ime=ime, jms=jms, jme=jme, kms=kms, kme=kme,                        &
                        its=its, ite=ite, jts=jts, jte=jte, kts=kts, kte=kte                         )

                do j = jts,jte
                do i = its,ite
                   taod5502d_p(i,j) = 0._RKIND
                   do k = kts,kte
                      taod5502d_p(i,j) = taod5502d_p(i,j) + taod5503d_p(i,k,j)
                      taod5503d(k,i)   = taod5503d_p(i,k,j)
                   enddo
                   taod5502d(i) = taod5502d_p(i,j)
                enddo
                enddo

                call calc_aerosol_rrtmg_sw( &
                              ht       = ht_p        , dz8w     = dz_p        , &
                              p        = pres_hyd_p  , t3d      = t_p         , &
                              qv3d     = qv_p        , tauaer   = tauaer_p    , &
                              ssaaer   = ssaaer_p    , asyaer   = asyaer_p    , &
                              aod5502d = taod5502d_p , aod5503d = taod5503d_p , &
                              aer_type = taer_type_p ,                                               &
                              aer_aod550_opt = taer_aod550_opt, aer_angexp_opt = taer_angexp_opt,    &
                              aer_ssa_opt    = taer_ssa_opt   , aer_asy_opt    = taer_asy_opt   ,    &
                              aer_aod550_val = aer_aod550_val , aer_angexp_val = aer_angexp_val ,    &
                              aer_ssa_val    = aer_ssa_val    , aer_asy_val    = aer_asy_val    ,    &
                              ims=ims, ime=ime, jms=jms, jme=jme, kms=kms, kme=kme,                 &
                              its=its, ite=ite, jts=jts, jte=jte, kts=kts, kte=kte                  )
             else
                call mpas_log_write('WARNING: aero_dir_rad_fdb=3 requires mp_thompson_aerosols or ' // &
                                    'mp_tempo with config_tempo_aerosolaware=.true. ' //               &
                                    'Falling back to no direct feedback (aero_dir_rad_fdb=0).',        &
                                    messageType=MPAS_LOG_WARN)
             endif

          !--- Case 4: Thompson WFA/IFA + Mie combined (optical weighting)
          !    tauaer = tauaer_thompson + tauaer_mie
          !    ssaaer = (ssa_t*tau_t + ssa_m*tau_m) / tau_total
          !    asyaer = (asy_t*ssa_t*tau_t + asy_m*ssa_m*tau_m) / (ssa_t*tau_t + ssa_m*tau_m)
          case(4)
             if(microp_scheme == 'mp_thompson_aerosols' .or. &
                (microp_scheme == 'mp_tempo' .and. config_tempo_aerosolaware)) then

                call mpas_pool_get_array(mesh,'zgrid',zgrid)
                call mpas_pool_get_array(diag_physics,'taod5502d',taod5502d)
                call mpas_pool_get_array(diag_physics,'taod5503d',taod5503d)

                !--- Step 1: Compute Thompson optical properties into temporary arrays
                ! This array is local array
                allocate(tauaer_thompson_tmp(ims:ime,kms:kme,jms:jme,nbndsw))
                allocate(ssaaer_thompson_tmp(ims:ime,kms:kme,jms:jme,nbndsw))
                allocate(asyaer_thompson_tmp(ims:ime,kms:kme,jms:jme,nbndsw))

                do j = jts,jte
                do i = its,ite
                   ht_p(i,j) = zgrid(1,i)
                   if(xland_p(i,j)==1._RKIND) then
                      taer_type_p(i,j) = 1
                   elseif(xland_p(i,j)==2._RKIND) then
                      taer_type_p(i,j) = 3
                   endif
                enddo
                enddo

                call gt_aod( &
                        p_phy = pres_hyd_p , dz8w = dz_p  , t_phy     = t_p        , qvapor = qv_p , &
                        nwfa  = nwfa_p     , nifa  = nifa_p, taod5503d = taod5503d_p,                 &
                        ims=ims, ime=ime, jms=jms, jme=jme, kms=kms, kme=kme,                        &
                        its=its, ite=ite, jts=jts, jte=jte, kts=kts, kte=kte                         )

                do j = jts,jte
                do i = its,ite
                   taod5502d_p(i,j) = 0._RKIND
                   do k = kts,kte
                      taod5502d_p(i,j) = taod5502d_p(i,j) + taod5503d_p(i,k,j)
                      taod5503d(k,i)   = taod5503d_p(i,k,j)
                   enddo
                   taod5502d(i) = taod5502d_p(i,j)
                enddo
                enddo

                call calc_aerosol_rrtmg_sw( &
                              ht       = ht_p        , dz8w     = dz_p        , &
                              p        = pres_hyd_p  , t3d      = t_p         , &
                              qv3d     = qv_p        , tauaer   = tauaer_thompson_tmp , &
                              ssaaer   = ssaaer_thompson_tmp    , asyaer   = asyaer_thompson_tmp    , &
                              aod5502d = taod5502d_p , aod5503d = taod5503d_p , &
                              aer_type = taer_type_p ,                                               &
                              aer_aod550_opt = taer_aod550_opt, aer_angexp_opt = taer_angexp_opt,    &
                              aer_ssa_opt    = taer_ssa_opt   , aer_asy_opt    = taer_asy_opt   ,    &
                              aer_aod550_val = aer_aod550_val , aer_angexp_val = aer_angexp_val ,    &
                              aer_ssa_val    = aer_ssa_val    , aer_asy_val    = aer_asy_val    ,    &
                              ims=ims, ime=ime, jms=jms, jme=jme, kms=kms, kme=kme,                 &
                              its=its, ite=ite, jts=jts, jte=jte, kts=kts, kte=kte                  )

                !--- Step 2: Combine Thompson + Mie 
                !--- Weight approach
                !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j,k,n, &
                !$OMP   tau_thompson,tau_mie,ssa_thompson,ssa_mie,asy_thompson,asy_mie, &
                !$OMP   tau_total,ssa_num,asy_num)
                do n = 1,nbndsw
                do j = jts,jte
                do k = kts,kte
                do i = its,ite
                   tau_thompson = tauaer_thompson_tmp(i,k,j,n)
                   tau_mie      = tauaer_sw(n,k,i)
                   ssa_thompson = ssaaer_thompson_tmp(i,k,j,n)
                   ssa_mie      = ssaaer_sw(n,k,i)
                   asy_thompson = asyaer_thompson_tmp(i,k,j,n)
                   asy_mie      = asyaer_sw(n,k,i)

                   tau_total = tau_thompson + tau_mie

                   if(tau_total > 0.0_RKIND) then
                      ssa_num = ssa_thompson * tau_thompson + ssa_mie * tau_mie
                      tauaer_p(i,k,j,n) = tau_total
                      ssaaer_p(i,k,j,n) = ssa_num / tau_total

                      asy_num = asy_thompson * ssa_thompson * tau_thompson + &
                                asy_mie      * ssa_mie      * tau_mie
                      if(ssa_num > 0.0_RKIND) then
                         asyaer_p(i,k,j,n) = asy_num / ssa_num
                      else
                         asyaer_p(i,k,j,n) = 0.0_RKIND
                      endif
                   else
                      ! No aerosol in this cell: no-aerosol (defaults)
                      tauaer_p(i,k,j,n) = 0.0_RKIND
                      ssaaer_p(i,k,j,n) = 1.0_RKIND
                      asyaer_p(i,k,j,n) = 0.0_RKIND
                   endif
                enddo
                enddo
                enddo
                enddo
                !$OMP END PARALLEL DO

                deallocate(tauaer_thompson_tmp)
                deallocate(ssaaer_thompson_tmp)
                deallocate(asyaer_thompson_tmp)

             else
                call mpas_log_write('WARNING: aero_dir_rad_fdb=4 requires mp_thompson_aerosols or ' // &
                                    'mp_tempo with config_tempo_aerosolaware=.true. ' //               &
                                    'Falling back to Mie-only direct feedback (aero_dir_rad_fdb=2).',  &
                                    messageType=MPAS_LOG_WARN)
                !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j,k,n)
                do n = 1,nbndsw
                do j = jts,jte
                do k = kts,kte
                do i = its,ite
                   tauaer_p(i,k,j,n) = tauaer_sw(n,k,i)
                   ssaaer_p(i,k,j,n) = ssaaer_sw(n,k,i)
                   asyaer_p(i,k,j,n) = asyaer_sw(n,k,i)
                enddo
                enddo
                enddo
                enddo
                !$OMP END PARALLEL DO
             endif

          case default
             call mpas_log_write('WARNING: Unknown aero_dir_rad_fdb value. ' // &
                                 'Defaulting to no direct feedback.',            &
                                 messageType=MPAS_LOG_WARN)

       end select direct_feedback_select

       !----------------------------------------------------------------------
       ! Section 3c: Ozone volume mixing ratio (No changes from NOAA GSL)
       !----------------------------------------------------------------------
       if(config_o3climatology) then
          do k = 1, num_oznLevels
             pin_p(k) = pin(k)
          enddo
          do j = jts,jte
          do k = 1, num_oznLevels
          do i = its,ite
             o3clim_p(i,k,j) = o3clim(k,i)
          enddo
          enddo
          enddo
       else
          do k = 1, num_oznLevels
             pin_p(k) = 0.0_RKIND
          enddo
          do j = jts,jte
          do k = 1, num_oznLevels
          do i = its,ite
             o3clim_p(i,k,j) = 0.0_RKIND
          enddo
          enddo
          enddo
       endif

!-----------------------------------------------------------------------------------------------------------------
! Section 4: CAM SW initialization (No changes from NOAA GSL)
!-----------------------------------------------------------------------------------------------------------------
    case("cam_sw")
       do j = jts,jte
       do i = its,ite
          sfc_emiss_p(i,j) = sfc_emiss(i)
          olrtoa_p(i,j)    = 0.0_RKIND  ;  glw_p(i,j)    = 0.0_RKIND
          lwcf_p(i,j)      = 0.0_RKIND  ;  lwdnb_p(i,j)  = 0.0_RKIND
          lwdnbc_p(i,j)    = 0.0_RKIND  ;  lwdnt_p(i,j)  = 0.0_RKIND
          lwdntc_p(i,j)    = 0.0_RKIND  ;  lwupb_p(i,j)  = 0.0_RKIND
          lwupbc_p(i,j)    = 0.0_RKIND  ;  lwupt_p(i,j)  = 0.0_RKIND
          lwuptc_p(i,j)    = 0.0_RKIND
       enddo
       do k = kts,kte
       do i = its,ite
          rthratenlw_p(i,k,j) = 0.0_RKIND
          cemiss_p(i,k,j)     = 0.0_RKIND
          taucldc_p(i,k,j)    = 0.0_RKIND
          taucldi_p(i,k,j)    = 0.0_RKIND
       enddo
       enddo
       enddo

       if(xtime_s .lt. 1.e-12) then
          do j = jts,jte
          do n = 1,cam_abs_dim1
          do k = kts,kte
          do i = its,ite
             absnxt_p(i,k,n,j) = 0.0_RKIND
          enddo
          enddo
          enddo
          do n = 1,cam_abs_dim2
          do k = kts,kte+1
          do i = its,ite
             abstot_p(i,k,n,j) = 0.0_RKIND
          enddo
          enddo
          enddo
          do k = kts,kte+1
          do i = its,ite
             emstot_p(i,k,j) = 0.0_RKIND
          enddo
          enddo
          enddo
       endif

       do k = 1, num_oznlevels
          pin_p(k) = pin(k)
       enddo
       do n = 1, num_months
       do j = jts,jte
       do k = 1, num_oznlevels
       do i = its,ite
          ozmixm_p(i,k,j,n) = ozmixm(n,k,i)
       enddo
       enddo
       enddo
       enddo

       do k = 1, num_aerlevels
          m_hybi_p(k) = m_hybi(k,1)
       enddo
       do j = jts,jte
       do i = its,ite
          m_psp_p(i,j) = m_ps(i)
          m_psn_p(i,j) = m_ps(i)
       enddo
       enddo
       do n = 1,num_aerosols
       do j = jts,jte
       do k = 1, num_aerlevels
       do i = its,ite
          aerosolcp_p(i,k,j,n) = aerosols(n,k,i)
          aerosolcn_p(i,k,j,n) = aerosols(n,k,i)
       enddo
       enddo
       enddo
       enddo

    case default

 end select radiation_sw_select

 end subroutine radiation_sw_from_MPAS
!=================================================================================================================
 subroutine radiation_sw_to_MPAS(diag_physics,tend_physics,its,ite)
!=================================================================================================================

!input arguments:
 type(mpas_pool_type),intent(inout):: diag_physics
 type(mpas_pool_type),intent(inout):: tend_physics

 integer,intent(in):: its,ite

!local variables:
 integer:: i,j,k,n

!local pointers:
 real(kind=RKIND),dimension(:),pointer  :: coszr,gsw,swcf,swdnb,swdnbc,swdnt,swdntc, &
                                           swupb,swupbc,swupt,swuptc,swddir,swddni,  &
                                           swddif
 real(kind=RKIND),dimension(:),pointer  :: swdnbcln, swdntcln, swuptcln, swupbcln                                   
 real(kind=RKIND),dimension(:,:),pointer:: rthratensw !,rthratenswc

!-----------------------------------------------------------------------------------------------------------------
!call mpas_log_write(' ')
!call mpas_log_write('--- enter subroutine radiation_sw_to_MPAS:')

 call mpas_pool_get_array(diag_physics,'coszr'     ,coszr     )
 call mpas_pool_get_array(diag_physics,'gsw'       ,gsw       )
 call mpas_pool_get_array(diag_physics,'swcf'      ,swcf      )
 call mpas_pool_get_array(diag_physics,'swdnb'     ,swdnb     )
 call mpas_pool_get_array(diag_physics,'swdnbc'    ,swdnbc    )
 call mpas_pool_get_array(diag_physics,'swdnbcln'  ,swdnbcln  ) ! Clean Diagnostic
 call mpas_pool_get_array(diag_physics,'swdnt'     ,swdnt     )
 call mpas_pool_get_array(diag_physics,'swdntc'    ,swdntc    )
 call mpas_pool_get_array(diag_physics,'swdntcln'  ,swdntcln  ) ! Clean Diagnostic
 call mpas_pool_get_array(diag_physics,'swupb'     ,swupb     )
 call mpas_pool_get_array(diag_physics,'swupbc'    ,swupbc    )
 call mpas_pool_get_array(diag_physics,'swupbcln'  ,swupbcln  ) ! Clean Diagnostic
 call mpas_pool_get_array(diag_physics,'swupt'     ,swupt     )
 call mpas_pool_get_array(diag_physics,'swuptc'    ,swuptc    )
 call mpas_pool_get_array(diag_physics,'swuptcln'  ,swuptcln  ) ! Clean Diagnostic
 call mpas_pool_get_array(diag_physics,'swddir'    ,swddir    )
 call mpas_pool_get_array(diag_physics,'swddni'    ,swddni    )
 call mpas_pool_get_array(diag_physics,'swddif'    ,swddif    )
 call mpas_pool_get_array(tend_physics,'rthratensw',rthratensw)
! call mpas_pool_get_array(tend_physics,'rthratenswc',rthratenswc)

 do j = jts,jte
 do i = its,ite
    coszr(i) = coszr_p(i,j)
    gsw(i)    = gsw_p(i,j)
    swcf(i)   = swcf_p(i,j)
    swdnb(i)  = swdnb_p(i,j)
    swdnbc(i) = swdnbc_p(i,j)
    swdnbcln(i) = swdnbcln_p(i,j)
    swdnt(i)  = swdnt_p(i,j)
    swdntc(i) = swdntc_p(i,j)
    swdntcln(i) = swdntcln_p(i,j)
    swupb(i)  = swupb_p(i,j)
    swupbc(i) = swupbc_p(i,j)
    swupbcln(i) = swupbcln_p(i,j)
    swupt(i)  = swupt_p(i,j)
    swuptc(i) = swuptc_p(i,j)
    swuptcln(i) = swuptcln_p(i,j)
    swddir(i) = swddir_p(i,j)
    swddni(i) = swddni_p(i,j)
    swddif(i) = swddif_p(i,j)
 enddo
! enddo

 do k = kts,kte
 do i = its,ite
    rthratensw(k,i) = rthratensw_p(i,k,j)
!    rthratenswc(k,i) = rthratenswc_p(i,k,j)
 enddo
 enddo

 enddo

!call mpas_log_write('--- enter subroutine radiation_sw_to_MPAS:')
!call mpas_log_write(' ')

 end subroutine radiation_sw_to_MPAS

!=================================================================================================================
 subroutine init_radiation_sw(dminfo,configs,mesh,atm_input,diag,diag_physics,state,time_lev)
!=================================================================================================================

!input arguments:
 type(dm_info), intent(in):: dminfo
 type(mpas_pool_type),intent(in):: configs
 type(mpas_pool_type),intent(in),optional:: mesh
 type(mpas_pool_type),intent(in),optional:: diag
 type(mpas_pool_type),intent(in),optional:: diag_physics

 integer,intent(in),optional:: time_lev

!inout arguments:
 type(mpas_pool_type),intent(inout),optional:: atm_input
 type(mpas_pool_type),intent(inout),optional:: state

!local pointers:
 character(len=StrKIND),pointer:: radt_sw_scheme

!-----------------------------------------------------------------------------------------------------------------

 call mpas_pool_get_config(configs,'config_radt_sw_scheme',radt_sw_scheme)

!call to shortwave radiation scheme:
 radiation_sw_select: select case (trim(radt_sw_scheme))
    case ("rrtmg_sw")
       call rrtmg_initsw_forMPAS(dminfo)

    case("cam_sw")
       call camradinit(dminfo,mesh,atm_input,diag,diag_physics,state,time_lev)

    case default
 end select radiation_sw_select

 end subroutine init_radiation_sw

!=================================================================================================================
 subroutine driver_radiation_sw(itimestep,configs,mesh,state,time_lev,diag_physics,atm_input, &
                                sfc_input,tend_physics,xtime_s,its,ite)
!                                sfc_input,tend_physics,config_calc_clean_atm_diag,xtime_s,its,ite)
!=================================================================================================================

!input arguments:
 integer,intent(in):: itimestep
 type(mpas_pool_type),intent(in):: configs
 type(mpas_pool_type),intent(in):: mesh

 integer,intent(in):: its,ite
 integer,intent(in):: time_lev

 real(kind=RKIND),intent(in):: xtime_s

!inout arguments:
 type(mpas_pool_type),intent(inout):: state
 type(mpas_pool_type),intent(inout):: diag_physics
 type(mpas_pool_type),intent(inout):: atm_input
 type(mpas_pool_type),intent(inout):: sfc_input
 type(mpas_pool_type),intent(inout):: tend_physics


!local pointers:
 logical,pointer:: config_o3climatology
 logical,pointer:: config_calc_clean_atm_diag
 character(len=StrKIND),pointer:: radt_sw_scheme
 character(len=StrKIND),pointer:: radt_cld_overlap,radt_cld_dcorrlen

!local variables:
 integer:: o3input
 integer:: cldovrlp,idcor
 real(kind=RKIND):: radt,xtime_m

!-----------------------------------------------------------------------------------------------------------------
!call mpas_log_write(' --- enter subroutine driver_radiation_sw: $i',intArgs=(/itimestep/))

 call mpas_pool_get_config(configs,'config_o3climatology' ,config_o3climatology   )
 call mpas_pool_get_config(configs,'config_radt_sw_scheme',radt_sw_scheme         )
 call mpas_pool_get_config(configs,'config_radt_cld_overlap' ,radt_cld_overlap    )
 call mpas_pool_get_config(configs,'config_radt_cld_dcorrlen',radt_cld_dcorrlen   )
 call mpas_pool_get_config(configs,'config_calc_clean_atm_diag', config_calc_clean_atm_diag)

 call mpas_log_write('=== DEBUG: config_calc_clean_atm_diag check ===')
 if (config_calc_clean_atm_diag) then
   call mpas_log_write('config_calc_clean_atm_diag = TRUE')
 else
   call mpas_log_write('config_calc_clean_atm_diag = FALSE')
 endif
 call mpas_log_write('===============================================')

 xtime_m = xtime_s/60.

!copy MPAS arrays to local arrays:
 call radiation_sw_from_MPAS(configs,mesh,state,time_lev,diag_physics,atm_input,sfc_input,xtime_s,its,ite)
 !call radiation_sw_from_MPAS(configs,mesh,state,time_lev,diag_physics,atm_input,config_calc_clean_atm_diag,sfc_input,xtime_s,its,ite)

! This should be OMP MASTER with barrier afterwards or OMP SINGLE, since declin and solcon
! are global variables in mpas_atmphys_vars.F and race conditions may occur otherwise!
!$OMP SINGLE
!... calculates solar declination:
!call radconst(declin,solcon,julday,degrad,dpd)
 call radconst(declin,solcon,curr_julday,degrad,dpd)
!call mpas_log_write('     ITIMESTEP   = $i', intArgs=(/itimestep/))
!call mpas_log_write('     YEAR        = $i', intArgs=(/year/))
!call mpas_log_write('     JULDAY      = $i', intArgs=(/julday/))
!call mpas_log_write('     GMT         = $r', realArgs=(/gmt/))
!call mpas_log_write('     XTIME_M     = $r', realArgs=(/xtime_m/))
!call mpas_log_write('     CURR_JULDAY = $r', realArgs=(/curr_julday/))
!call mpas_log_write('     SOLCON      = $r', realArgs=(/solcon/))
!call mpas_log_write('     DECLIN      = $r', realArgs=(/declin/))
!$OMP END SINGLE

!... convert the radiation time_step to minutes:
 radt = dt_radtsw/60.
    
!call to shortwave radiation scheme:
 radiation_sw_select: select case (trim(radt_sw_scheme))
    case ("rrtmg_sw")
       o3input = 0
       if(config_o3climatology) o3input = 2

       !initialization of cloud overlap option:
       if(trim(radt_cld_overlap) == "none"              ) cldovrlp = 0
       if(trim(radt_cld_overlap) == "random"            ) cldovrlp = 1
       if(trim(radt_cld_overlap) == "maximum_random"    ) cldovrlp = 2
       if(trim(radt_cld_overlap) == "maximum"           ) cldovrlp = 3
       if(trim(radt_cld_overlap) == "exponential"       ) cldovrlp = 4
       if(trim(radt_cld_overlap) == "exponential_random") cldovrlp = 5

       idcor = 0
       if(trim(radt_cld_overlap)=="exponential" .or. trim(radt_cld_overlap)=="exponential_random") then
          if(trim(radt_cld_dcorrlen) == "constant"        ) idcor = 0
          if(trim(radt_cld_dcorrlen) == "latitude_varying") idcor = 1
       endif

       call mpas_log_write('=== DEBUG: Before calling rrtmg_swrad ===')
       if (config_calc_clean_atm_diag) then
          call mpas_log_write('Passing config_calc_clean_atm_diag = TRUE (logical)')
       else
          call mpas_log_write('Passing config_calc_clean_atm_diag = FALSE (logical)')
       endif
       call mpas_log_write('===========================================')

       call mpas_timer_start('rrtmg_swrad')
       call rrtmg_swrad( &
              p3d        = pres_hyd_p    , p8w      = pres2_hyd_p , pi3d     = pi_p           , &
              t3d        = t_p           , t8w      = t2_p        , dz8w     = dz_p           , &
!             qv3d       = qv_p          , qc3d     = qc_p        , qi3d     = qi_p           , &
!             qs3d       = qs_p          , cldfra3d = cldfrac_p   , tsk      = tsk_p          , &
              qv3d       = qvrad_p       , qc3d     = qcrad_p     , qi3d     = qirad_p        , &
              qs3d       = qsrad_p       , cldfra3d = cldfrac_p   , tsk      = tsk_p          , &
              albedo     = sfc_albedo_p  , xland    = xland_p     , xice     = xice_p         , &
              snow       = snow_p        , coszr    = coszr_p     , xtime    = xtime_m        , &
              gmt        = gmt           , julday   = julday      , radt     = radt           , &
              degrad     = degrad        , declin   = declin      , solcon   = solcon         , &
              xlat       = xlat_p        , xlong    = xlon_p      , icloud   = icloud         , &
              cldovrlp   = cldovrlp      , idcor    = idcor       , o3input  = o3input        , &
              noznlevels = num_oznlevels , pin      = pin_p       , o3clim   = o3clim_p       , &
              gsw        = gsw_p         , swcf     = swcf_p      , rthratensw = rthratensw_p , &
!              rthratenswc = rthratenswc_p, &
              has_reqc   = has_reqc      , has_reqi = has_reqi    , has_reqs   = has_reqs     , &
              re_cloud   = recloud_p     , re_ice   = reice_p     , re_snow    = resnow_p     , &
              aer_opt    = aer_opt       , tauaer3d = tauaer_p    , ssaaer3d   = ssaaer_p     , &
              asyaer3d   = asyaer_p      , swupt    = swupt_p     , swuptc     = swuptc_p     , &
              swuptcln   = swuptcln_p    , swdnt    = swdnt_p  , &
              swdntc     = swdntc_p     , swdntcln   = swdntcln_p    , swupb    = swupb_p  , &
              swupbc     = swupbc_p     , swupbcln   = swupbcln_p    , swdnb    = swdnb_p  , &
              swdnbc     = swdnbc_p     , swdnbcln   = swdnbcln_p    , swupflx  = swupflx_p, &
              swupflxc   = swupflxc_p   , swupflxcln = swupflxcln_p  , swdnflx  = swdnflx_p, &
              swdnflxc   = swdnflxc_p   , swdnflxcln = swdnflxcln_p  , swddir   = swddir_p , &
              swddni     = swddni_p     , swddif     = swddif_p      , &
              config_calc_clean_atm_diag=config_calc_clean_atm_diag, &
              julian     = curr_julday, &
              ids = ids , ide = ide , jds = jds , jde = jde , kds = kds , kde = kde ,           &
              ims = ims , ime = ime , jms = jms , jme = jme , kms = kms , kme = kme ,           &
              its = its , ite = ite , jts = jts , jte = jte , kts = kts , kte = kte             &
                       )
       call mpas_timer_stop('rrtmg_swrad')
       call mpas_log_write('=== DEBUG: Clean diagnostic values after rrtmg_swrad ===')
       call mpas_log_write('swdnbcln_p min/max: $r $r', realArgs=(/minval(swdnbcln_p), maxval(swdnbcln_p)/))
       call mpas_log_write('swuptcln_p min/max: $r $r', realArgs=(/minval(swuptcln_p), maxval(swuptcln_p)/))
       call mpas_log_write('========================================================')

    case ("cam_sw")
       call mpas_timer_start('camrad_sw')
       call camrad( dolw = .false. , dosw = .true. ,                                         &
                p_phy         = pres_hyd_p    , p8w           = pres2_hyd_p   ,              &
                pi_phy        = pi_p          , t_phy         = t_p           ,              &
                z             = zmid_p        , dz8w          = dz_p          ,              &            
                rthratenlw    = rthratenlw_p  , rthratensw    = rthratensw_p  ,              &
                swupt         = swupt_p       , swuptc        = swuptc_p      ,              &
                swdnt         = swdnt_p       , swdntc        = swdntc_p      ,              &
                lwupt         = lwupt_p       , lwuptc        = lwuptc_p      ,              &
                lwdnt         = lwdnt_p       , lwdntc        = lwdntc_p      ,              &
                swupb         = swupb_p       , swupbc        = swupbc_p      ,              &
                swdnb         = swdnb_p       , swdnbc        = swdnbc_p      ,              &
                lwupb         = lwupb_p       , lwupbc        = lwupbc_p      ,              &
                lwdnb         = lwdnb_p       , lwdnbc        = lwdnbc_p      ,              &
                swcf          = swcf_p        , lwcf          = lwcf_p        ,              &
                gsw           = gsw_p         , glw           = glw_p         ,              &
                olr           = olrtoa_p      , cemiss        = cemiss_p      ,              &
                taucldc       = taucldc_p     , taucldi       = taucldi_p     ,              & 
                coszr         = coszr_p       , albedo        = sfc_albedo_p  ,              & 
                emiss         = sfc_emiss_p   , tsk           = tsk_p         ,              & 
                xlat          = xlat_p        , xlong         = xlon_p        ,              &
                rho_phy       = rho_p         , qv3d          = qv_p          ,              & 
                qc3d          = qc_p          , qr3d          = qr_p          ,              &
                qi3d          = qi_p          , qs3d          = qs_p          ,              &
                qg3d          = qg_p          , f_qc          = f_qc          ,              &
                f_qr          = f_qr          , f_qi          = f_qi          ,              &
                f_qs          = f_qs          , f_ice_phy     = f_ice         ,              &
                f_rain_phy    = f_rain        , cldfra        = cldfrac_p     ,              &
                xland         = xland_p       , xice          = xice_p        ,              &
                num_months    = num_months    , levsiz        = num_oznlevels ,              & 
                pin0          = pin_p         , ozmixm        = ozmixm_p      ,              &
                paerlev       = num_aerlevels , naer_c        = num_aerosols  ,              &
                m_psp         = m_psp_p       , m_psn         = m_psn_p       ,              &
                aerosolcp     = aerosolcp_p   , aerosolcn     = aerosolcn_p   ,              &
                m_hybi0       = m_hybi_p      , snow          = snow_p        ,              &
                cam_abs_dim1  = cam_abs_dim1  , cam_abs_dim2  = cam_abs_dim2  ,              &
                gmt           = gmt           , yr            = year          ,              &
                julday        = julday        , julian        = curr_julday   ,              &
                dt            = dt_dyn        , xtime         = xtime_m       ,              &
                declin        = declin        , solcon        = solcon        ,              &
                radt          = radt          , degrad        = degrad        ,              &
                n_cldadv      = 3             , abstot_3d     = abstot_p      ,              &
                absnxt_3d     = absnxt_p      , emstot_3d     = emstot_p      ,              &
                doabsems      = doabsems      ,                                              &
                ids = ids , ide = ide , jds = jds , jde = jde , kds = kds , kde = kde ,      &
                ims = ims , ime = ime , jms = jms , jme = jme , kms = kms , kme = kme ,      &
                its = its , ite = ite , jts = jts , jte = jte , kts = kts , kte = kte        &
                  )
       call mpas_timer_stop('camrad_sw')

    case default
 end select radiation_sw_select

!copy local arrays to MPAS grid:
 call radiation_sw_to_MPAS(diag_physics,tend_physics,its,ite)

!call mpas_log_write('--- end subroutine driver_radiation_sw.')

 end subroutine driver_radiation_sw

!=================================================================================================================
 subroutine radconst(declin,solcon,julian,degrad,dpd)
!=================================================================================================================

!input arguments:
!integer,intent(in):: julian
 real(kind=RKIND),intent(in):: julian
 real(kind=RKIND),intent(in):: degrad,dpd

!output arguments:
 real(kind=RKIND),intent(out):: declin,solcon

!local variables:
 real(kind=RKIND):: obecl,sinob,sxlong,arg,decdeg,djul,rjul,eccfac

!-----------------------------------------------------------------------------------------------------------------

   declin=0.
   solcon=0.

!obecl : obliquity = 23.5 degree.
        
   obecl=23.5*degrad
   sinob=sin(obecl)
        
!calculate longitude of the sun from vernal equinox:        

   if(julian.ge.80.)sxlong=dpd*(julian-80.)
   if(julian.lt.80.)sxlong=dpd*(julian+285.)
   sxlong=sxlong*degrad
   arg=sinob*sin(sxlong)
   declin=asin(arg)
   decdeg=declin/degrad

!solar constant eccentricity factor (paltridge and platt 1976)

   djul=julian*360./365.
   rjul=djul*degrad
   eccfac=1.000110+0.034221*cos(rjul)+0.001280*sin(rjul)+0.000719*  &
          cos(2*rjul)+0.000077*sin(2*rjul)
   solcon=solcon_0*eccfac

 end subroutine radconst

!=================================================================================================================
 end module mpas_atmphys_driver_radiation_sw
!=================================================================================================================
