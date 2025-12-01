MODULE FISOC_WM_MOD
  
  USE ESMF
  USE FISOC_WM_Wrapper
  USE FISOC_utils_MOD
  USE FISOC_types_MOD
    
  IMPLICIT NONE
  
  PRIVATE
  
  PUBLIC FISOC_WM_setvm,FISOC_WM_register
  
  CHARACTER(len=ESMF_MAXSTR),SAVE :: WM_impFBname = "WM import fields"
  CHARACTER(len=ESMF_MAXSTR),SAVE :: WM_expFBname = "WM export fields"

CONTAINS
  
  !------------------------------------------------------------------------------
  SUBROUTINE FISOC_WM_setvm(comp, rc)
      
    type(ESMF_GridComp)     :: comp
    integer, intent(out)    :: rc

    CALL ESMF_LogWrite('WM_setvm', ESMF_LOGMSG_INFO)

    CALL ESMF_LogWrite('WM_setvm end', ESMF_LOGMSG_INFO)

    rc=ESMF_SUCCESS
  
  END SUBROUTINE FISOC_WM_setvm

  SUBROUTINE FISOC_WM_register(FISOC_WM, rc)
    
    TYPE(ESMF_GridComp)  :: FISOC_WM
    INTEGER, INTENT(OUT) :: rc
    
    rc = ESMF_FAILURE

    CALL ESMF_GridCompSetEntryPoint(FISOC_WM, ESMF_METHOD_INITIALIZE, &
         userRoutine=FISOC_WM_init_phase1, phase=1, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    CALL ESMF_GridCompSetEntryPoint(FISOC_WM, ESMF_METHOD_INITIALIZE, &
         userRoutine=FISOC_WM_init_phase2, phase=2, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    
    CALL ESMF_GridCompSetEntryPoint(FISOC_WM, ESMF_METHOD_RUN, &
         userRoutine=FISOC_WM_run, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    
    CALL ESMF_GridCompSetEntryPoint(FISOC_WM, ESMF_METHOD_FINALIZE, &
         userRoutine=FISOC_WM_finalise, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    rc = ESMF_SUCCESS
 
  END SUBROUTINE FISOC_WM_register

  !------------------------------------------------------------------------------
  ! Initialisation is implemented in two stages.  The first stage is for independent 
  ! initialisation of the WM and the second stage occurrs after the OM has been 
  ! initialised, to allow inter-component consistency checks or use of the OM 
  ! state to complete WM initalisation.
  SUBROUTINE FISOC_WM_init_phase1(FISOC_WM, WM_ImpSt, WM_ExpSt, FISOC_clock, rc)
    TYPE(ESMF_GridComp)        :: FISOC_WM
    TYPE(ESMF_State)           :: WM_ImpSt, WM_ExpSt
    TYPE(ESMF_Clock)           :: FISOC_clock
    INTEGER, INTENT(OUT)       :: rc

    TYPE(ESMF_config)          :: FISOC_config
    TYPE(ESMF_mesh)            :: WM_mesh
    TYPE(ESMF_grid)            :: WM_grid, OM_grid
    TYPE(ESMF_fieldBundle)     :: WM_ExpFB,OM_ExpFB
    TYPE(ESMF_VM)              :: vm
    TYPE(ESMF_field),ALLOCATABLE :: fieldList(:)
    CHARACTER(len=ESMF_MAXSTR) :: label
    CHARACTER(len=ESMF_MAXSTR),ALLOCATABLE :: WM_DerVarList(:)
    INTEGER                    :: fieldCount
    LOGICAL                    :: WM_UseOMGrid, OM_UseWMGrid

    rc = ESMF_FAILURE

    msg = "WM initialise started"
    CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_INFO, &
       line=__LINE__, file=__FILE__, rc=rc)

    CALL ESMF_GridCompGet(FISOC_WM, config=FISOC_config, vm=vm, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    ! if the WM and OM are on the same grid then the WM export field bundle is 
    ! used as the OM import field bundle
    CALL FISOC_ConfigDerivedAttribute(FISOC_config, WM_UseOMGrid, 'WM_UseOMGrid',rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    CALL FISOC_ConfigDerivedAttribute(FISOC_config, OM_UseWMGrid, 'OM_UseWMGrid',rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    IF (WM_UseOMGrid.OR.OM_UseWMGrid) THEN
      WM_impFBname = "OM export fields"
    END IF

    label = 'FISOC_WM_DerVars:' 
    CALL FISOC_getListFromConfig(FISOC_config, label, WM_DerVarList,rc=rc)
    IF (rc.EQ.ESMF_RC_NOT_FOUND) THEN
      ALLOCATE(WM_DerVarList(0))
    ELSE IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) THEN
      CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    END IF

    CALL FISOC_ConfigDerivedAttribute(FISOC_config, WM_UseOMGrid, 'WM_UseOMGrid',rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    ! create empty field bundle
    WM_ExpFB = ESMF_FieldBundleCreate(name=WM_expFBname, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    
    ! The model-specific initialisation adds the WM vars to the field bundle.
    ! This is followed by adding non-model-specific derived variables. 
# if defined(FISOC_WM_GRID)
    IF (WM_UseOMGrid) THEN
      ! Get grid from field from field list from field bundle from state.
      CALL ESMF_StateGet(WM_ImpSt, WM_impFBname, OM_ExpFB, rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
           line=__LINE__, file=__FILE__)) &
           CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
      CALL ESMF_FieldBundleGet(OM_ExpFB, fieldCount=fieldCount, rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
           line=__LINE__, file=__FILE__)) &
           CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
      ALLOCATE(fieldList(fieldCount))
      CALL ESMF_FieldBundleGet(OM_ExpFB, fieldList=fieldList,rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
           line=__LINE__, file=__FILE__)) &
           CALL ESMF_Finalize(endflag=ESMF_END_ABORT)        
      CALL ESMF_FieldGet(FieldList(1), grid=OM_grid, rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
           line=__LINE__, file=__FILE__)) &
           CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
      
      ! We assign the WM_grid to the OM_grid.  This is not taking a copy, 
      ! but using the same object.
      WM_grid = OM_grid
      
    END IF
    
    CALL FISOC_WM_Wrapper_Init_Phase1(FISOC_config,vm,WM_ExpFB,WM_grid,rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    CALL FISOC_populateFieldBundle(WM_DerVarList,WM_ExpFB,WM_grid,rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    
# elif defined(FISOC_WM_MESH)
    CALL FISOC_WM_Wrapper_Init_Phase1(FISOC_config,vm,WM_ExpFB,WM_mesh,rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    CALL FISOC_populateFieldBundle(WM_DerVarList,WM_ExpFB,WM_mesh,rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

# else
    msg = "ERROR: FISOC does not recognise WM geom type."
    CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_ERROR, &
         line=__LINE__, file=__FILE__, rc=rc)
    CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
# endif

    ! we only add the fields to the import state as a way of letting the coupler get hold of the 
    ! grid.  There must be a better way to do this. 
    ! [edit: there is now, see email from Gerhard Theurich]
    CALL ESMF_StateAdd(WM_ImpSt, (/WM_ExpFB/), rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    CALL ESMF_StateAdd(WM_ExpSt, (/WM_ExpFB/), rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    msg = "WM initialise phase 1 complete"
    CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_INFO, &
       line=__LINE__, file=__FILE__, rc=rc)

    rc = ESMF_SUCCESS

  END SUBROUTINE FISOC_WM_init_phase1


  !------------------------------------------------------------------------------
  SUBROUTINE FISOC_WM_init_phase2(FISOC_WM, WM_ImpSt, WM_ExpSt, FISOC_clock, rc)

    TYPE(ESMF_GridComp)    :: FISOC_WM
    TYPE(ESMF_State)       :: WM_ImpSt, WM_ExpSt
    TYPE(ESMF_Clock)       :: FISOC_clock
    INTEGER, INTENT(OUT)   :: rc

    TYPE(ESMF_config)      :: FISOC_config
    TYPE(ESMF_fieldbundle) :: WM_ImpFB,WM_ExpFB
    TYPE(ESMF_VM)          :: vm

    rc = ESMF_FAILURE

    CALL ESMF_GridCompGet(FISOC_WM, config=FISOC_config, vm=vm, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    CALL ESMF_StateGet(WM_ImpSt, WM_impFBname, WM_ImpFB, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)    

    CALL ESMF_StateGet(WM_ExpSt, WM_expFBname, WM_ExpFB, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)    

    CALL FISOC_WM_Wrapper_Init_Phase2(FISOC_config,vm,WM_ImpFB,WM_ExpFB,rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    msg = "WM initialise phase 2 complete (allows the WM access to the OM initial state) "
    CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_INFO, &
       line=__LINE__, file=__FILE__, rc=rc)

    rc = ESMF_SUCCESS

  END SUBROUTINE FISOC_WM_init_phase2
  

  !------------------------------------------------------------------------------
  SUBROUTINE FISOC_WM_run(FISOC_WM, WM_ImpSt, WM_ExpSt, FISOC_clock, rc)

    TYPE(ESMF_GridComp)    :: FISOC_WM
    TYPE(ESMF_State)       :: WM_ImpSt, WM_ExpSt
    TYPE(ESMF_Clock)       :: FISOC_clock
    INTEGER, INTENT(OUT)   :: rc

    TYPE(ESMF_FileStatus_Flag) :: NC_status
    CHARACTER(len=ESMF_MAXSTR) :: WM_NCfreq
    INTEGER(ESMF_KIND_I8)      :: advanceCount
    INTEGER                    :: advanceCountInt4
    
    TYPE(ESMF_fieldbundle) :: WM_ImpFB,WM_ExpFB
    LOGICAL                :: verbose_coupling,wM_writeNetcdf
    TYPE(ESMF_config)      :: FISOC_config
!    INTEGER                :: localPet
    TYPE(ESMF_VM)          :: vm

    rc = ESMF_FAILURE

    CALL ESMF_GridCompGet(FISOC_WM, config=FISOC_config, vm=vm, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

!    CALL ESMF_VMGet(vm, localPet=localPet, rc=rc)
!    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
!         line=__LINE__, file=__FILE__)) &
!         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    CALL ESMF_ConfigGetAttribute(FISOC_config, verbose_coupling, label='verbose_coupling:', rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    CALL ESMF_ConfigGetAttribute(FISOC_config,  WM_writeNetcdf, label='WM_writeNetcdf:', rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) THEN
       msg="WARNING: WM_writeNetcdf not found in FISOC_config.rc, not writing to NetCDF"
       WM_writeNetcdf = .FALSE.
       CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_WARNING, &
            line=__LINE__, file=__FILE__, rc=rc)
    END IF

    IF (WM_writeNetcdf) THEN      
       CALL ESMF_ConfigGetAttribute(FISOC_config,  WM_NCfreq, label='WM_NCfreq:', rc=rc)
       IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__)) THEN
          msg="WARNING: WM_NCfreq not found in FISOC_config.rc, setting to all"
          WM_NCfreq = "all"
          CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_WARNING, &
               line=__LINE__, file=__FILE__, rc=rc)
       END IF
    END IF

    ! "AdvanceCount" gives the number of (OM) timesteps.  Use it to make NetCDF filename.
    CALL ESMF_ClockGet(FISOC_clock, advanceCount=advanceCount, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)  
    ! extract field bundles from import and export states to send to model-specific wrapper

    CALL ESMF_StateGet(WM_ImpSt, WM_impFBname, WM_ImpFB, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)    

    CALL ESMF_StateGet(WM_ExpSt, WM_expFBname, WM_ExpFB, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)    
    
!!$    CALL FISOC_WM_calcDerivedFields_pre(WM_ExpFB,FISOC_config,FISOC_clock,rc)
!!$    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
!!$         line=__LINE__, file=__FILE__)) &
!!$         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)    
!!$
!!$    CALL FISOC_WM_maskOMfields(FISOC_config,WM_ImpFB,WM_ExpFB,rc=rc)
!!$    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
!!$         line=__LINE__, file=__FILE__)) &
!!$         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    CALL FISOC_WM_Wrapper_Run(FISOC_config,vm,WM_ImpFB=WM_ImpFB, &
         WM_ExpFB=WM_ExpFB,rc_local=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)    

    CALL WM_NetcdfWrapper(FISOC_config,WM_writeNetcdf,WM_NCfreq,advanceCount, &
               WM_ExpFB=WM_ExpFB,WM_ImpFB=WM_ImpFB)   
!!$    CALL FISOC_WM_calcDerivedFields_post(WM_ExpFB,FISOC_config,FISOC_clock,rc)
!!$    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
!!$         line=__LINE__, file=__FILE__)) &
!!$         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)    

    msg = "WM run complete for current timestep"
    CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_INFO, &
       line=__LINE__, file=__FILE__, rc=rc)
    
    rc = ESMF_SUCCESS

  END SUBROUTINE FISOC_WM_run


  !------------------------------------------------------------------------------
  SUBROUTINE FISOC_WM_finalise(FISOC_WM, WM_ImpSt, WM_ExpSt, FISOC_clock, rc)
    TYPE(ESMF_GridComp)  :: FISOC_WM
    TYPE(ESMF_State)     :: WM_ImpSt, WM_ExpSt
    TYPE(ESMF_Clock)     :: FISOC_clock
    INTEGER, INTENT(OUT) :: rc

    INTEGER                      :: localPet
    TYPE(ESMF_VM)                :: vm
    TYPE(ESMF_config)            :: FISOC_config
    TYPE(ESMF_fieldbundle)       :: WM_ImpFB, WM_ExpFB
    INTEGER                      :: FieldCount,ii
    TYPE(ESMF_field),ALLOCATABLE :: ImpFieldList(:)
    TYPE(ESMF_field),ALLOCATABLE :: ExpFieldList(:)
    TYPE(ESMF_mesh)              :: WM_mesh
    LOGICAL                      :: WM_UseOMGrid, OM_UseWMGrid

    CHARACTER(len=ESMF_MAXSTR) :: name
    
    rc = ESMF_FAILURE

    msg = "WM finalise: destroy fields, bundles and mesh"
    CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_INFO, &
       line=__LINE__, file=__FILE__, rc=rc)

    CALL ESMF_GridCompGet(FISOC_WM, config=FISOC_config, vm=vm, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    CALL FISOC_ConfigDerivedAttribute(FISOC_config, WM_UseOMGrid, 'WM_UseOMGrid',rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    CALL FISOC_ConfigDerivedAttribute(FISOC_config, OM_UseWMGrid, 'OM_UseWMGrid',rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    CALL ESMF_VMBarrier(vm, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    
    CALL FISOC_WM_Wrapper_Finalize(FISOC_config,vm,rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)


    ! destroy WM export fields
    CALL ESMF_StateGet(WM_ExpSt, WM_expFBname, WM_ExpFB, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)    
       
    ! ...how many fields?...
    CALL ESMF_FieldBundleGet(WM_ExpFB, fieldCount=FieldCount, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    
    ! ... get list of fields from bundle.
    ALLOCATE(ExpFieldList(FieldCount))
    CALL ESMF_FieldBundleGet(WM_ExpFB, fieldList=ExpFieldList,rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    DO ii = 1,FieldCount
       
       CALL ESMF_FieldGet(field=ExpFieldList(ii),name=name, rc=rc)
       IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__)) &
            CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
       
       CALL ESMF_FieldDestroy(ExpFieldList(ii), rc=rc)
       IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__)) &
            CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    END DO
    
    DEALLOCATE(ExpFieldList)
    
    CALL ESMF_FieldBundleDestroy(WM_ExpFB, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    
    
    IF ( (.NOT.WM_UseOMGrid).AND.(.NOT.OM_UseWMGrid) ) THEN
      ! destroy WM import fields    
      CALL ESMF_StateGet(WM_ImpSt, WM_impFBname, WM_ImpFB, rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
           line=__LINE__, file=__FILE__)) &
           CALL ESMF_Finalize(endflag=ESMF_END_ABORT)    
      
      ! ...how many fields?...
      CALL ESMF_FieldBundleGet(WM_ImpFB, fieldCount=FieldCount, rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
           line=__LINE__, file=__FILE__)) &
           CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
      
      ! ... get list of fields from bundle.
      ALLOCATE(ImpFieldList(FieldCount))
      CALL ESMF_FieldBundleGet(WM_ImpFB, fieldList=ImpFieldList,rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
           line=__LINE__, file=__FILE__)) &
           CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    
!    CALL ESMF_FieldGet(ImpFieldList(1), mesh=WM_mesh, rc=rc)
!    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
!         line=__LINE__, file=__FILE__)) &
!         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    
      DO ii = 1,FieldCount
        CALL ESMF_FieldDestroy(ImpFieldList(ii), rc=rc)
        IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
             line=__LINE__, file=__FILE__)) &
            CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
      END DO
      
!    CALL ESMF_MeshDestroy(WM_mesh,rc=rc)
!    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
!         line=__LINE__, file=__FILE__)) &
!         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    
      DEALLOCATE(ImpFieldList)
    
      CALL ESMF_FieldBundleDestroy(WM_ImpFB, rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
           line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)    
    END IF
      
    rc = ESMF_SUCCESS
    
  END SUBROUTINE FISOC_WM_finalise
  

 !------------------------------------------------------------------------------
  SUBROUTINE WM_NetcdfWrapper(FISOC_config,WM_writeNetcdf,WM_NCfreq,advanceCount,WM_ExpFB,WM_ImpFB)

    TYPE(ESMF_config),INTENT(INOUT)               :: FISOC_config
    LOGICAL,INTENT(IN)                            :: WM_writeNetcdf
    CHARACTER(len=ESMF_MAXSTR),INTENT(IN)         :: WM_NCfreq
    INTEGER(ESMF_KIND_I8),INTENT(IN)              :: advanceCount
    TYPE(ESMF_fieldBundle),INTENT(INOUT),OPTIONAL :: WM_ExpFB, WM_ImpFB

    CHARACTER(len=ESMF_MAXSTR) :: OutputFileName
    INTEGER :: rc

    IF(WM_writeNetcdf) THEN
       
       SELECT CASE (WM_NCfreq)

       ! Writing out at all timesteps: write any field bundles we're given
       CASE("all")          
          IF (PRESENT(WM_ImpFB)) THEN
             msg = "Writing NetCDF output from FISOC on wave grid (WM import)"
             CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_INFO, &
                  line=__LINE__, file=__FILE__, rc=rc)
             WRITE (OutputFileName, "(A14,I0,A3)") "FISOC_WM_imp_t", advanceCount, ".nc"
             CALL FISOC_FB2NC(OutputFileName,WM_ImpFB,FISOC_config)
          END IF
          IF (PRESENT(WM_ExpFB)) THEN
             msg = "Writing NetCDF output from FISOC on wave grid (WM export)"
             CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_INFO, &
                  line=__LINE__, file=__FILE__, rc=rc)
             WRITE (OutputFileName, "(A14,I0,A3)") "FISOC_WM_exp_t", advanceCount, ".nc"
             CALL FISOC_FB2NC(OutputFileName,WM_ExpFB,FISOC_config)
          END IF

       ! Writing out at WM export timesteps: hopefully we're given both field bundles
       CASE("WM")
          IF (PRESENT(WM_ImpFB)) THEN
             IF (PRESENT(WM_ExpFB)) THEN
                msg = "Writing NetCDF output from FISOC on wave grid (WM imp and exp)"
                CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_INFO, &
                     line=__LINE__, file=__FILE__, rc=rc)
                WRITE (OutputFileName, "(A14,I0,A3)") "FISOC_WM_exp_t", advanceCount, ".nc"
                CALL FISOC_FB2NC(OutputFileName,WM_ExpFB,FISOC_config)
                WRITE (OutputFileName, "(A14,I0,A3)") "FISOC_WM_imp_t", advanceCount, ".nc"
                CALL FISOC_FB2NC(OutputFileName,WM_ImpFB,FISOC_config)
          ELSE
                msg = "WM_NCfreq=WM but we have only WM_ImpFB and not WM_ExpFB... :("
                CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_ERROR, &
                     line=__LINE__, file=__FILE__, rc=rc)
                CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
             END IF
          END IF

       CASE DEFAULT
          msg = "WM_NCfreq value not recognised: "//WM_NCfreq
          CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_ERROR, &
               line=__LINE__, file=__FILE__, rc=rc)
          CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

       END SELECT
       
    END IF
    
  END SUBROUTINE WM_NetcdfWrapper

  
END MODULE FISOC_WM_MOD
