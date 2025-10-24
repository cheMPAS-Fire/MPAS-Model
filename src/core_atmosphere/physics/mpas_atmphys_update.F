! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!=================================================================================================================
 module mpas_atmphys_update
 use mpas_kind_types
 use mpas_pool_routines

 use mpas_atmphys_driver_convection
 use mpas_atmphys_vars

 implicit none
 private
 public:: physics_update, &
          update_radiation_diagnostics


!Update diagnostics.
!Laura D. Fowler (send comments to laura@ucar.edu).
!2013-05-01.
!
! subroutines in mpas_atmphys_update:
! -----------------------------------
! physics_update              : not used.
! update_radiation_diagnostics: update accumulated radiation diagnostics.
!
! add-ons and modifications to sourcecode:
! ----------------------------------------
! * throughout the sourcecode, replaced all "var_struct" defined arrays by local pointers.
!   Laura D. Fowler (laura@ucar.edu) / 2014-04-22.
! * modified sourcecode to use pools.
!   Laura D. Fowler (laura@ucar.edu) / 2014-05-15.
! * renamed config_conv_deep_scheme to config_convection_scheme.
!   Laura D. Fowler (laura@ucar.edu) / 2014-09-18.


 contains

 
!=================================================================================================================
 subroutine physics_update(domain,dt)
!=================================================================================================================

!input arguments:
 type(domain_type),intent(in):: domain
 real(kind=RKIND),intent(in) :: dt 

!local variables:
 type(block_type),pointer:: block
 integer:: i,j,k

!-----------------------------------------------------------------------------------------------------------------
!call mpas_log_write('')
!call mpas_log_write('--- begin physics_update:')
 block => domain % blocklist
 do while(associated(block))

    !parameterization of convection: update accumulated precipitation.
    !call update_convection_deep(dt,config_bucket_rainc,block%mesh,block%diag_physics)

    block => block % next
 end do
!call mpas_log_write('--- end physics_update:')

 end subroutine physics_update

!=================================================================================================================
 subroutine update_radiation_diagnostics(configs,mesh,diag_physics,its,ite)
!=================================================================================================================

!input arguments:
 type(mpas_pool_type),intent(in):: configs
 type(mpas_pool_type),intent(in):: mesh
 integer,intent(in):: its,ite

!inout arguments:
 type(mpas_pool_type),intent(inout):: diag_physics

!local pointers:
 integer,pointer:: nCellsSolve
 integer,dimension(:),pointer  :: i_acswdnb,i_acswdnbc,i_acswdnt,i_acswdntc, &
                                  i_acswupb,i_acswupbc,i_acswupt,i_acswuptc, &
                                  i_aclwdnb,i_aclwdnbc,i_aclwdnt,i_aclwdntc, &
                                  i_aclwupb,i_aclwupbc,i_aclwupt,i_aclwuptc

 integer,dimension(:),pointer  :: i_acswdnbcln,i_acswdntcln,i_acswupbcln, &
                                  i_acswuptcln,i_aclwdnbcln,i_aclwdntcln, &
                                  i_aclwupbcln,i_aclwuptcln
 real(kind=RKIND),pointer:: bucket_radt
 real(kind=RKIND),dimension(:),pointer:: swdnb,swdnbc,swdnt,swdntc,  &
                                         swupb,swupbc,swupt,swuptc,  &
                                         lwdnb,lwdnbc,lwdnt,lwdntc,  &
                                         lwupb,lwupbc,lwupt,lwuptc
 real(kind=RKIND),dimension(:),pointer:: swdnbcln,swdntcln,swupbcln,swuptcln,&
                                         lwdnbcln,lwdntcln,lwupbcln,lwuptcln
 real(kind=RKIND),dimension(:),pointer:: acswdnb,acswdnbc,acswdnt,acswdntc,  &
                                         acswupb,acswupbc,acswupt,acswuptc,  &
                                         aclwdnb,aclwdnbc,aclwdnt,aclwdntc,  &
                                         aclwupb,aclwupbc,aclwupt,aclwuptc
 real(kind=RKIND),dimension(:),pointer:: acswdnbcln,acswdntcln,acswupbcln,acswuptcln,  &
                                         aclwdnbcln,aclwdntcln,aclwupbcln,aclwuptcln

!local variables and arrays:
 integer:: iCell

!-----------------------------------------------------------------------------------------------------------------

 call mpas_pool_get_config(configs,'config_bucket_radt',bucket_radt)

 call mpas_pool_get_dimension(mesh,'nCellsSolve',nCellsSolve)

 call mpas_pool_get_array(diag_physics,'i_acswdnb' , i_acswdnb )
 call mpas_pool_get_array(diag_physics,'i_acswdnbc', i_acswdnbc)
 call mpas_pool_get_array(diag_physics,'i_acswdnbcln', i_acswdnbcln)
 call mpas_pool_get_array(diag_physics,'i_acswdnt' , i_acswdnt )
 call mpas_pool_get_array(diag_physics,'i_acswdntc', i_acswdntc)
 call mpas_pool_get_array(diag_physics,'i_acswdntcln', i_acswdntcln)
 call mpas_pool_get_array(diag_physics,'i_acswupb' , i_acswupb )
 call mpas_pool_get_array(diag_physics,'i_acswupbc', i_acswupbc)
 call mpas_pool_get_array(diag_physics,'i_acswupbcln', i_acswupbcln)
 call mpas_pool_get_array(diag_physics,'i_acswupt' , i_acswupt )
 call mpas_pool_get_array(diag_physics,'i_acswuptc', i_acswuptc)
 call mpas_pool_get_array(diag_physics,'i_acswuptcln', i_acswuptcln)
 call mpas_pool_get_array(diag_physics,'i_aclwdnb' , i_aclwdnb )
 call mpas_pool_get_array(diag_physics,'i_aclwdnbc', i_aclwdnbc)
 call mpas_pool_get_array(diag_physics,'i_aclwdnbcln', i_aclwdnbcln)
 call mpas_pool_get_array(diag_physics,'i_aclwdnt' , i_aclwdnt )
 call mpas_pool_get_array(diag_physics,'i_aclwdntc', i_aclwdntc)
 call mpas_pool_get_array(diag_physics,'i_aclwdntcln', i_aclwdntcln)
 call mpas_pool_get_array(diag_physics,'i_aclwupb' , i_aclwupb )
 call mpas_pool_get_array(diag_physics,'i_aclwupbc', i_aclwupbc)
 call mpas_pool_get_array(diag_physics,'i_aclwupbcln', i_aclwupbcln)
 call mpas_pool_get_array(diag_physics,'i_aclwupt' , i_aclwupt )
 call mpas_pool_get_array(diag_physics,'i_aclwuptc', i_aclwuptc)
 call mpas_pool_get_array(diag_physics,'i_aclwuptcln', i_aclwuptcln)

 call mpas_pool_get_array(diag_physics,'acswdnb'   , acswdnb   )
 call mpas_pool_get_array(diag_physics,'acswdnbc'  , acswdnbc  )
 call mpas_pool_get_array(diag_physics,'acswdnbcln', acswdnbcln)
 call mpas_pool_get_array(diag_physics,'acswdnt'   , acswdnt   )
 call mpas_pool_get_array(diag_physics,'acswdntc'  , acswdntc  )
 call mpas_pool_get_array(diag_physics,'acswdntcln', acswdntcln)
 call mpas_pool_get_array(diag_physics,'acswupb'   , acswupb   )
 call mpas_pool_get_array(diag_physics,'acswupbc'  , acswupbc  )
 call mpas_pool_get_array(diag_physics,'acswupbcln', acswupbcln)
 call mpas_pool_get_array(diag_physics,'acswupt'   , acswupt   )
 call mpas_pool_get_array(diag_physics,'acswuptc'  , acswuptc  )
 call mpas_pool_get_array(diag_physics,'acswuptcln', acswuptcln)
 call mpas_pool_get_array(diag_physics,'aclwdnb'   , aclwdnb   )
 call mpas_pool_get_array(diag_physics,'aclwdnbc'  , aclwdnbc  )
 call mpas_pool_get_array(diag_physics,'aclwdnbcln', aclwdnbcln)
 call mpas_pool_get_array(diag_physics,'aclwdnt'   , aclwdnt   )
 call mpas_pool_get_array(diag_physics,'aclwdntc'  , aclwdntc  )
 call mpas_pool_get_array(diag_physics,'aclwdntcln', aclwdntcln)
 call mpas_pool_get_array(diag_physics,'aclwupb'   , aclwupb   )
 call mpas_pool_get_array(diag_physics,'aclwupbc'  , aclwupbc  )
 call mpas_pool_get_array(diag_physics,'aclwupbcln', aclwupbcln)
 call mpas_pool_get_array(diag_physics,'aclwupt'   , aclwupt   )
 call mpas_pool_get_array(diag_physics,'aclwuptc'  , aclwuptc  )
 call mpas_pool_get_array(diag_physics,'aclwuptcln', aclwuptcln)

 call mpas_pool_get_array(diag_physics,'swdnb'     , swdnb     )
 call mpas_pool_get_array(diag_physics,'swdnbc'    , swdnbc    )
 call mpas_pool_get_array(diag_physics,'swdnbcln'  , swdnbcln  )
 call mpas_pool_get_array(diag_physics,'swdnt'     , swdnt     )
 call mpas_pool_get_array(diag_physics,'swdntc'    , swdntc    )
 call mpas_pool_get_array(diag_physics,'swdntcln'  , swdntcln  )
 call mpas_pool_get_array(diag_physics,'swupb'     , swupb     )
 call mpas_pool_get_array(diag_physics,'swupbc'    , swupbc    )
 call mpas_pool_get_array(diag_physics,'swupbcln'  , swupbcln  )
 call mpas_pool_get_array(diag_physics,'swupt'     , swupt     )
 call mpas_pool_get_array(diag_physics,'swuptc'    , swuptc    )
 call mpas_pool_get_array(diag_physics,'swuptcln'  , swuptcln  )
 call mpas_pool_get_array(diag_physics,'lwdnb'     , lwdnb     )
 call mpas_pool_get_array(diag_physics,'lwdnbc'    , lwdnbc    )
 call mpas_pool_get_array(diag_physics,'lwdnbcln'  , lwdnbcln  )
 call mpas_pool_get_array(diag_physics,'lwdnt'     , lwdnt     )
 call mpas_pool_get_array(diag_physics,'lwdntc'    , lwdntc    )
 call mpas_pool_get_array(diag_physics,'lwdntcln'  , lwdntcln  )
 call mpas_pool_get_array(diag_physics,'lwupb'     , lwupb     )
 call mpas_pool_get_array(diag_physics,'lwupbc'    , lwupbc    )
 call mpas_pool_get_array(diag_physics,'lwupbcln'  , lwupbcln  )
 call mpas_pool_get_array(diag_physics,'lwupt'     , lwupt     )
 call mpas_pool_get_array(diag_physics,'lwuptc'    , lwuptc    )
 call mpas_pool_get_array(diag_physics,'lwuptcln'  , lwuptcln  )

 do iCell = its, ite
    !short-wave radiation:
    acswdnb(iCell)  = acswdnb (iCell) + swdnb (iCell)*dt_dyn
    acswdnbc(iCell) = acswdnbc(iCell) + swdnbc(iCell)*dt_dyn
    acswdnbcln(iCell) = acswdnbcln(iCell) + swdnbcln(iCell)*dt_dyn
    acswdnt(iCell)  = acswdnt (iCell) + swdnt (iCell)*dt_dyn
    acswdntc(iCell) = acswdntc(iCell) + swdntc(iCell)*dt_dyn
    acswdntcln(iCell) = acswdntcln(iCell) + swdntcln(iCell)*dt_dyn
    acswupb(iCell)  = acswupb (iCell) + swupb (iCell)*dt_dyn
    acswupbc(iCell) = acswupbc(iCell) + swupbc(iCell)*dt_dyn
    acswupbcln(iCell) = acswupbcln(iCell) + swupbcln(iCell)*dt_dyn
    acswupt(iCell)  = acswupt (iCell) + swupt (iCell)*dt_dyn
    acswuptc(iCell) = acswuptc(iCell) + swuptc(iCell)*dt_dyn
    acswuptcln(iCell) = acswuptcln(iCell) + swuptcln(iCell)*dt_dyn
    !long-wave radiation:
    aclwdnb(iCell)  = aclwdnb (iCell) + lwdnb (iCell)*dt_dyn
    aclwdnbc(iCell) = aclwdnbc(iCell) + lwdnbc(iCell)*dt_dyn
    aclwdnbcln(iCell) = aclwdnbcln(iCell) + lwdnbcln(iCell)*dt_dyn
    aclwdnt(iCell)  = aclwdnt (iCell) + lwdnt (iCell)*dt_dyn
    aclwdntc(iCell) = aclwdntc(iCell) + lwdntc(iCell)*dt_dyn
    aclwdntcln(iCell) = aclwdntcln(iCell) + lwdntcln(iCell)*dt_dyn
    aclwupb(iCell)  = aclwupb (iCell) + lwupb (iCell)*dt_dyn
    aclwupbc(iCell) = aclwupbc(iCell) + lwupbc(iCell)*dt_dyn
    aclwupbcln(iCell) = aclwupbcln(iCell) + lwupbcln(iCell)*dt_dyn
    aclwupt(iCell)  = aclwupt (iCell) + lwupt (iCell)*dt_dyn
    aclwuptc(iCell) = aclwuptc(iCell) + lwuptc(iCell)*dt_dyn
    aclwuptcln(iCell) = aclwuptcln(iCell) + lwuptcln(iCell)*dt_dyn
 enddo

 if(l_acradt .and. bucket_radt.gt.0._RKIND) then

    do iCell = its, ite
       !short-wave radiation:
       if(acswdnb(iCell) .gt. bucket_radt) then
          i_acswdnb(iCell) = i_acswdnb(iCell) + 1
          acswdnb(iCell) = acswdnb(iCell) - bucket_radt
       endif   
       if(acswdnbc(iCell) .gt. bucket_radt) then
          i_acswdnbc(iCell) = i_acswdnbc(iCell) + 1
          acswdnbc(iCell) = acswdnbc(iCell) - bucket_radt
       endif   
       if(acswdnbcln(iCell) .gt. bucket_radt) then
          i_acswdnbcln(iCell) = i_acswdnbcln(iCell) + 1
          acswdnbcln(iCell) = acswdnbcln(iCell) - bucket_radt
       endif   
       if(acswdnt(iCell) .gt. bucket_radt) then
          i_acswdnt(iCell) = i_acswdnt(iCell) + 1
          acswdnt(iCell) = acswdnt(iCell) - bucket_radt
       endif   
       if(acswdntc(iCell) .gt. bucket_radt) then
          i_acswdntc(iCell) = i_acswdntc(iCell) + 1
          acswdntc(iCell) = acswdntc(iCell) - bucket_radt
       endif
       if(acswdntcln(iCell) .gt. bucket_radt) then
          i_acswdntcln(iCell) = i_acswdntcln(iCell) + 1
          acswdntcln(iCell) = acswdntcln(iCell) - bucket_radt
       endif
       if(acswupb(iCell) .gt. bucket_radt) then
          i_acswupb(iCell) = i_acswupb(iCell) + 1
          acswupb(iCell) = acswupb(iCell) - bucket_radt
       endif   
       if(acswupbc(iCell) .gt. bucket_radt) then
          i_acswupbc(iCell) = i_acswupbc(iCell) + 1
          acswupbc(iCell) = acswupbc(iCell) - bucket_radt
       endif   
       if(acswupbcln(iCell) .gt. bucket_radt) then
          i_acswupbcln(iCell) = i_acswupbcln(iCell) + 1
          acswupbcln(iCell) = acswupbcln(iCell) - bucket_radt
       endif   
       if(acswupt(iCell) .gt. bucket_radt) then
          i_acswupt(iCell) = i_acswupt(iCell) + 1
          acswupt(iCell) = acswupt(iCell) - bucket_radt
       endif   
       if(acswuptc(iCell) .gt. bucket_radt) then
          i_acswuptc(iCell) = i_acswuptc(iCell) + 1
          acswuptc(iCell) = acswuptc(iCell) - bucket_radt
       endif
       if(acswuptcln(iCell) .gt. bucket_radt) then
          i_acswuptcln(iCell) = i_acswuptcln(iCell) + 1
          acswuptcln(iCell) = acswuptcln(iCell) - bucket_radt
       endif
       !long-wave radiation:
       if(aclwdnb(iCell) .gt. bucket_radt) then
          i_aclwdnb(iCell) = i_aclwdnb(iCell) + 1
          aclwdnb(iCell) = aclwdnb(iCell) - bucket_radt
       endif   
       if(aclwdnbc(iCell) .gt. bucket_radt) then
          i_aclwdnbc(iCell) = i_aclwdnbc(iCell) + 1
          aclwdnbc(iCell) = aclwdnbc(iCell) - bucket_radt
       endif   
       if(aclwdnbcln(iCell) .gt. bucket_radt) then
          i_aclwdnbcln(iCell) = i_aclwdnbcln(iCell) + 1
          aclwdnbcln(iCell) = aclwdnbcln(iCell) - bucket_radt
       endif   
       if(aclwdnt(iCell) .gt. bucket_radt) then
          i_aclwdnt(iCell) = i_aclwdnt(iCell) + 1
          aclwdnt(iCell) = aclwdnt(iCell) - bucket_radt
       endif   
       if(aclwdntc(iCell) .gt. bucket_radt) then
          i_aclwdntc(iCell) = i_aclwdntc(iCell) + 1
          aclwdntc(iCell) = aclwdntc(iCell) - bucket_radt
       endif
       if(aclwdntcln(iCell) .gt. bucket_radt) then
          i_aclwdntcln(iCell) = i_aclwdntcln(iCell) + 1
          aclwdntcln(iCell) = aclwdntcln(iCell) - bucket_radt
       endif
       if(aclwupb(iCell) .gt. bucket_radt) then
          i_aclwupb(iCell) = i_aclwupb(iCell) + 1
          aclwupb(iCell) = aclwupb(iCell) - bucket_radt
       endif   
       if(aclwupbc(iCell) .gt. bucket_radt) then
          i_aclwupbc(iCell) = i_aclwupbc(iCell) + 1
          aclwupbc(iCell) = aclwupbc(iCell) - bucket_radt
       endif   
       if(aclwupbcln(iCell) .gt. bucket_radt) then
          i_aclwupbcln(iCell) = i_aclwupbcln(iCell) + 1
          aclwupbcln(iCell) = aclwupbcln(iCell) - bucket_radt
       endif   
       if(aclwupt(iCell) .gt. bucket_radt) then
          i_aclwupt(iCell) = i_aclwupt(iCell) + 1
          aclwupt(iCell) = aclwupt(iCell) - bucket_radt
       endif   
       if(aclwuptc(iCell) .gt. bucket_radt) then
          i_aclwuptc(iCell) = i_aclwuptc(iCell) + 1
          aclwuptc(iCell) = aclwuptc(iCell) - bucket_radt
       endif
       if(aclwuptcln(iCell) .gt. bucket_radt) then
          i_aclwuptcln(iCell) = i_aclwuptcln(iCell) + 1
          aclwuptcln(iCell) = aclwuptcln(iCell) - bucket_radt
       endif
    enddo

 endif

 end subroutine update_radiation_diagnostics

!=================================================================================================================
 end module mpas_atmphys_update
!=================================================================================================================
