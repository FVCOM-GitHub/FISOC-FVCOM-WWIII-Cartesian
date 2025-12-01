MODULE FISOC_coupler_MOD
  
  USE ESMF
  USE FISOC_utils_MOD
    
  IMPLICIT NONE
  
  PRIVATE
  
  PUBLIC FISOC_coupler_register
    
CONTAINS
  
  !------------------------------------------------------------------------------
  SUBROUTINE FISOC_coupler_register(FISOC_coupler, rc)
    
    TYPE(ESMF_CplComp)  :: FISOC_coupler
    INTEGER, INTENT(OUT) :: rc
    
    rc = ESMF_FAILURE

    CALL ESMF_CplCompSetEntryPoint(FISOC_coupler, ESMF_METHOD_INITIALIZE, &
         userRoutine=FISOC_coupler_init_phase1, phase=1, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    CALL ESMF_CplCompSetEntryPoint(FISOC_coupler, ESMF_METHOD_INITIALIZE, &
         userRoutine=FISOC_coupler_init_phase2, phase=2, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    
    CALL ESMF_CplCompSetEntryPoint(FISOC_coupler, ESMF_METHOD_RUN, &
         userRoutine=FISOC_coupler_run_phase1, phase=1, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    
    CALL ESMF_CplCompSetEntryPoint(FISOC_coupler, ESMF_METHOD_RUN, &
         userRoutine=FISOC_coupler_run_phase2, phase=2, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    
    CALL ESMF_CplCompSetEntryPoint(FISOC_coupler, ESMF_METHOD_FINALIZE, &
         userRoutine=FISOC_coupler_finalise, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    rc = ESMF_SUCCESS
 
  END SUBROUTINE FISOC_coupler_register

  !------------------------------------------------------------------------------
  ! Initialisation is implemented in two phases.  The first phase is to regrid 
  ! the WM fields to the OM grid/mesh, and the second phase is to convert the OM 
  ! fields to the WM grid/mesh.
  SUBROUTINE FISOC_coupler_init_phase1(FISOC_coupler, WM_ExpSt, OM_ImpSt, FISOC_clock, rc)

    TYPE(ESMF_CplComp)     :: FISOC_coupler
    TYPE(ESMF_State)       :: WM_ExpSt, OM_ImpSt
    TYPE(ESMF_Clock)       :: FISOC_clock
    INTEGER, INTENT(OUT)   :: rc

    TYPE(ESMF_fieldBundle)        :: WM_ExpFB, OM_ExpFB, OM_ImpFB
    TYPE(ESMF_grid)               :: OM_grid
    TYPE(ESMF_mesh)               :: OM_mesh
    CHARACTER(len=ESMF_MAXSTR),ALLOCATABLE :: fieldNameList(:)
    INTEGER                       :: WM_ExpFieldCount
    TYPE(ESMF_RouteHandle)        :: WM2OM_regridRouteHandle
    TYPE(ESMF_VM)                 :: vm
    LOGICAL                       :: verbose_coupling
    TYPE(ESMF_config)             :: FISOC_config
    TYPE(ESMF_RegridMethod_Flag)  :: Regrid_method
    TYPE(ESMF_ExtrapMethod_Flag)  :: Extrap_method

    rc = ESMF_FAILURE

    ! get some key info from the coupler component and from FISOC config

    CALL ESMF_cplCompGet(FISOC_coupler, vm=vm, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    CALL ESMF_cplCompGet(FISOC_coupler, config=FISOC_config, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    CALL ESMF_ConfigGetAttribute(FISOC_config, verbose_coupling, label='verbose_coupling:', rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    CALL  FISOC_ConfigDerivedAttribute(FISOC_config, Regrid_method, label='WM2OM_regrid:', rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    CALL  FISOC_ConfigDerivedAttribute(FISOC_config, Extrap_method, label='WM2OM_extrap:', rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    ! Create the import bundle for the OM.  This will be populated with 
    ! the WM fields, regridded onto the OM grid or mesh.  It will be added 
    ! to the OM import state.
    OM_ImpFB = ESMF_FieldBundleCreate(name="OM import fields", rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)


    ! get field bundle to be regridded
    CALL ESMF_StateGet(WM_ExpSt, "WM export fields", WM_ExpFB, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    ! get field bundle containing target grid/mesh
    CALL ESMF_StateGet(OM_ImpSt, "OM export fields", OM_ExpFB, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    ! use the first field from the field bundles to make a route handle
    CALL FISOC_makeRHfromFB(WM_ExpFB,OM_ExpFB,        &
         Regrid_method,Extrap_method,verbose_coupling,WM2OM_regridRouteHandle,vm,rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    ! name the route handle to facilitate future extraction from state
    CALL ESMF_RouteHandleSet(WM2OM_regridRouteHandle,name="WM2OM_regridRouteHandle", rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    ! then add the route handle to the WM exp state for future regridding
    ! operations
    CALL ESMF_StateAdd(WM_ExpSt, (/WM2OM_regridRouteHandle/), rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    msg = "Regrid route handle created and added to WM export state "
    CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_INFO, &
         line=__LINE__, file=__FILE__, rc=rc)



    ! We will create OM import fields with names corresponding to the WM 
    ! export field names and use the new routehandle to regrid the WM 
    ! export fields onto the new fields.  Then bundle the OM import fields 
    ! and add them to the OM import state.

    ! get the grid or mesh from the OM exp bundle.  
# if defined(FISOC_OM_GRID)
    CALL FISOC_getGridFromFB(OM_expFB,OM_grid,rc=rc)
# elif defined(FISOC_OM_MESH)
    CALL FISOC_getMeshFromFB(OM_expFB,OM_mesh,rc=rc)
# else
    msg="invalid CPP options for OM geom type"
    CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_WARNING, &
            line=__LINE__, file=__FILE__, rc=rc)          
    CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
# endif
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    ! remove the OM export field bundle from OM export state 
    ! (we only put it there to set up the routehandle).
    CALL ESMF_StateRemove (OM_ImpSt, (/"OM export fields"/), rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    ! get a list of field names from the WM export bundle
    CALL ESMF_FieldBundleGet(WM_ExpFB, fieldCount=WM_ExpFieldCount, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    ALLOCATE(fieldNameList(WM_ExpFieldCount))
    CALL ESMF_FieldBundleGet(WM_expFB, fieldNameList=fieldNameList, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    ! use the OM grid or mesh and the field names from the WM export 
    ! field list to populate the OM import field bundle (fields are 
    ! initially empty or set to zero). 
# if defined(FISOC_OM_GRID)
    CALL FISOC_populateFieldBundle(fieldNameList,OM_impFB,   &
         OM_grid,init_value=REAL(0.0,ESMF_KIND_R8),rc=rc)
# elif defined(FISOC_OM_MESH)
    CALL FISOC_populateFieldBundle(fieldNameList,OM_impFB,    &
         OM_mesh,init_value=REAL(0.0,ESMF_KIND_R8),rc=rc)
# else
    msg="invalid CPP options for OM geom type"
    CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_ERROR, &
            line=__LINE__, file=__FILE__, rc=rc)          
    CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
# endif
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    ! regrid the WM export fields onto the OM grid or mesh to give 
    ! suitable values to the new fields
    CALL FISOC_regridFB(WM_expFB,OM_impFB,WM2OM_regridRouteHandle,rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)


    CALL ESMF_StateAdd(OM_ImpSt, (/OM_ImpFB/), rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    msg = "Regriding complete. Regridded fields stored in OM import state"
    CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_INFO, &
         line=__LINE__, file=__FILE__, rc=rc)

    DEALLOCATE(fieldNameList)

    rc = ESMF_SUCCESS


  END SUBROUTINE FISOC_coupler_init_phase1



  !------------------------------------------------------------------------------
  ! Initialisation is implemented in two phases.  The first phase is to regrid 
  ! the WM fields to the OM grid/mesh, and the second phase is to convert the OM 
  ! fields to the WM grid/mesh.
  SUBROUTINE FISOC_coupler_init_phase2(FISOC_coupler, OM_ExpSt, WM_ImpSt, FISOC_clock, rc)

    TYPE(ESMF_CplComp)     :: FISOC_coupler
    TYPE(ESMF_State)       :: OM_ExpSt, WM_ImpSt
    TYPE(ESMF_Clock)       :: FISOC_clock
    INTEGER, INTENT(OUT)   :: rc

    TYPE(ESMF_fieldBundle)        :: WM_ExpFB, OM_ExpFB, WM_ImpFB
    TYPE(ESMF_grid)               :: WM_grid
    TYPE(ESMF_mesh)               :: WM_mesh
    CHARACTER(len=ESMF_MAXSTR),ALLOCATABLE :: fieldNameList(:)
    INTEGER                       :: OM_ExpFieldCount
    TYPE(ESMF_RouteHandle)        :: OM2WM_regridRouteHandle
    TYPE(ESMF_VM)                 :: vm
    LOGICAL                       :: verbose_coupling
    TYPE(ESMF_config)             :: FISOC_config
    TYPE(ESMF_RegridMethod_Flag)  :: Regrid_method
    TYPE(ESMF_ExtrapMethod_Flag)  :: Extrap_method

    rc = ESMF_FAILURE

    ! get some key info from the coupler component and from FISOC config

    CALL ESMF_cplCompGet(FISOC_coupler, vm=vm, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    CALL ESMF_cplCompGet(FISOC_coupler, config=FISOC_config, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    CALL ESMF_ConfigGetAttribute(FISOC_config, verbose_coupling, label='verbose_coupling:', rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    CALL  FISOC_ConfigDerivedAttribute(FISOC_config, Regrid_method, label='OM2WM_regrid:', rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    CALL  FISOC_ConfigDerivedAttribute(FISOC_config, Extrap_method, label='OM2WM_extrap:', rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    ! Create the import bundle for the WM.  This will be populated with 
    ! the OM fields, regridded onto the WM grid or mesh.  It will be added 
    ! to the WM import state.
    WM_ImpFB = ESMF_FieldBundleCreate(name="WM import fields", rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)


    ! get field bundle to be regridded
    CALL ESMF_StateGet(OM_ExpSt, "OM export fields", OM_ExpFB, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    ! get field bundle containing target grid/mesh
    CALL ESMF_StateGet(WM_ImpSt, "WM export fields", WM_ExpFB, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    ! use the first field from the field bundles to make a route handle
    CALL FISOC_makeRHfromFB(OM_ExpFB,WM_ExpFB,        &
         Regrid_method,Extrap_method,verbose_coupling,OM2WM_regridRouteHandle,vm,rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    ! name the route handle to facilitate future extraction from state
    CALL ESMF_RouteHandleSet(OM2WM_regridRouteHandle,name="OM2WM_regridRouteHandle", rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    ! then add the route handle to the OM exp state for future regridding
    ! operations
    CALL ESMF_StateAdd(OM_ExpSt, (/OM2WM_regridRouteHandle/), rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    msg = "Regrid route handle created and added to OM export state "
    CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_INFO, &
         line=__LINE__, file=__FILE__, rc=rc)



    ! We will create WM import fields with names corresponding to the OM 
    ! export field names and use the new routehandle to regrid the OM 
    ! export fields onto the new fields.  Then bundle the WM import fields 
    ! and add them to the WM import state.

    ! get the grid or mesh from the WM exp bundle.  
# if defined(FISOC_WM_GRID)
    CALL FISOC_getGridFromFB(WM_expFB,WM_grid,rc=rc)
# elif defined(FISOC_WM_MESH)
    CALL FISOC_getMeshFromFB(WM_expFB,WM_mesh,rc=rc)
# else
    msg="invalid CPP options for WM geom type"
    CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_WARNING, &
         line=__LINE__, file=__FILE__, rc=rc)          
    CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
# endif
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    ! remove the WM export field bundle from WM export state 
    ! (we only put it there to set up the routehandle).
    CALL ESMF_StateRemove (WM_ImpSt, (/"WM export fields"/), rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    ! get a list of field names from the OM export bundle
    CALL ESMF_FieldBundleGet(OM_ExpFB, fieldCount=OM_ExpFieldCount, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    ALLOCATE(fieldNameList(OM_ExpFieldCount))
    CALL ESMF_FieldBundleGet(OM_expFB, fieldNameList=fieldNameList, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    ! use the WM grid or mesh and the field names from the OM export 
    ! field list to populate the WM import field bundle (fields are 
    ! initially empty or set to zero). 
# if defined(FISOC_WM_GRID)
    CALL FISOC_populateFieldBundle(fieldNameList,WM_impFB,   &
         WM_grid,init_value=REAL(0.0,ESMF_KIND_R8),rc=rc)
# elif defined(FISOC_WM_MESH)
    CALL FISOC_populateFieldBundle(fieldNameList,WM_impFB,    &
         WM_mesh,init_value=REAL(0.0,ESMF_KIND_R8),rc=rc)
# else
    msg="invalid CPP options for WM geom type"
    CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_ERROR, &
         line=__LINE__, file=__FILE__, rc=rc)          
    CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
# endif
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    ! regrid the OM export fields onto the WM grid or mesh to give 
    ! suitable values to the new fields
    CALL FISOC_regridFB(OM_expFB,WM_impFB,OM2WM_regridRouteHandle,rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)


    CALL ESMF_StateAdd(WM_ImpSt, (/WM_ImpFB/), rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    msg = "Regriding complete. Regridded fields stored in WM import state"
    CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_INFO, &
         line=__LINE__, file=__FILE__, rc=rc)

    DEALLOCATE(fieldNameList)

    rc = ESMF_SUCCESS


  END SUBROUTINE FISOC_coupler_init_phase2

  


  !------------------------------------------------------------------------------
  SUBROUTINE FISOC_coupler_run_phase1(FISOC_coupler, OM_ExpSt, WM_ImpSt, FISOC_clock, rc)
    TYPE(ESMF_CplComp)     :: FISOC_coupler
    TYPE(ESMF_State)       :: OM_ExpSt, WM_ImpSt
    TYPE(ESMF_Clock)       :: FISOC_clock
    INTEGER, INTENT(OUT)   :: rc

    INTEGER,PARAMETER             :: ListLen = 20
    TYPE(ESMF_fieldBundle)        :: WM_ExpFB, OM_ExpFB, WM_ImpFB
    TYPE(ESMF_config)             :: config
    CHARACTER(len=ESMF_MAXSTR)    :: WM_name, OM_name, fieldName, OM_ExpSt_NameList(ListLen), OM2WM_HandleName
    TYPE(ESMF_StateItem_Flag)     :: OM_ExpSt_TypeList(ListLen)
    INTEGER                       :: WM_ImpFieldCount, OM_ExpFieldCount, ii, NumRouteHandleItems, RouteHandleIndex
    TYPE(ESMF_Field),ALLOCATABLE  :: WM_ImpFieldList(:), OM_ExpFieldList(:)
    TYPE(ESMF_RouteHandle)        :: OM2WM_regridRouteHandle
    TYPE(ESMF_TypeKind_Flag)      :: fieldTypeKind
    
    REAL(ESMF_KIND_R8),POINTER    :: optr(:,:),iptr(:)
    INTEGER                       :: nn
    LOGICAL                       :: verbose_coupling, OM_writeNetcdf
    TYPE(ESMF_config)             :: FISOC_config

    rc = ESMF_FAILURE

    CALL ESMF_cplCompGet(FISOC_coupler, config=FISOC_config, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    CALL ESMF_ConfigGetAttribute(FISOC_config, verbose_coupling, label='verbose_coupling:', rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    ! Extract OM2WM regrid routehandle for regridding...
    CALL ESMF_StateGet(OM_ExpSt, "OM2WM_regridRouteHandle", OM2WM_regridRouteHandle, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    ! Extract OM field bundle for regridding...
    CALL ESMF_StateGet(OM_ExpSt, "OM export fields", OM_ExpFB, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    ! ...how many fields?...
    CALL ESMF_FieldBundleGet(OM_ExpFB, fieldCount=OM_ExpFieldCount, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    ! ... get list of fields from bundle.
    ALLOCATE(OM_ExpFieldList(OM_ExpFieldCount))
    CALL ESMF_FieldBundleGet(OM_ExpFB, fieldCount=OM_ExpFieldCount, fieldList=OM_ExpFieldList,rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    IF (verbose_coupling) THEN
       msg = "coupler extracted OM fields from OM export state"
       CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_INFO, &
            line=__LINE__, file=__FILE__, rc=rc)
    END IF
    
    IF (SIZE(OM_ExpFieldList).LT.1) THEN
       msg = "OM field list less than length 1"
       CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_ERROR, &
            line=__LINE__, file=__FILE__, rc=rc)
       CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    END IF

    ! Extract WM field bundle for regridding...
    CALL ESMF_StateGet(WM_ImpSt, "WM import fields", WM_ImpFB, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    ! ...how many fields?...
    CALL ESMF_FieldBundleGet(WM_ImpFB, fieldCount=WM_ImpFieldCount, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    ! ... get list of fields from bundle.
    ALLOCATE(WM_ImpFieldList(WM_ImpFieldCount))
    CALL ESMF_FieldBundleGet(WM_ImpFB, fieldList=WM_ImpFieldList,rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    loop_over_fields: DO ii = 1,OM_ExpFieldCount 

       CALL ESMF_FieldGet(OM_ExpFieldList(ii), name=fieldName, typekind=fieldTypeKind, rc=rc)
       IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__)) &
            CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

       CALL ESMF_FieldRegrid(OM_ExpFieldList(ii),WM_ImpFieldList(ii), &
            routehandle=OM2WM_regridRouteHandle, zeroregion= ESMF_REGION_TOTAL, &
            checkflag=.TRUE.,rc=rc)
       IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__)) &
            CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
       
       IF (verbose_coupling) THEN
          msg = "Regridded field "//fieldName
          CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_INFO, &
               line=__LINE__, file=__FILE__, rc=rc)          
       END IF

    END DO loop_over_fields

    IF (verbose_coupling) THEN
       msg = "Regriding complete. Regridded fields stored in WM import state"
       CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_INFO, &
            line=__LINE__, file=__FILE__, rc=rc)
    END IF

    rc = ESMF_SUCCESS
    
  END SUBROUTINE FISOC_coupler_run_phase1


  !------------------------------------------------------------------------------
  SUBROUTINE FISOC_coupler_run_phase2(FISOC_coupler, WM_ExpSt, OM_ImpSt, FISOC_clock, rc)
    TYPE(ESMF_CplComp)     :: FISOC_coupler
    TYPE(ESMF_State)       :: WM_ExpSt, OM_ImpSt
    TYPE(ESMF_Clock)       :: FISOC_clock
    INTEGER, INTENT(OUT)   :: rc

    INTEGER,PARAMETER             :: ListLen = 20
    TYPE(ESMF_fieldBundle)        :: OM_ExpFB, WM_ExpFB, OM_ImpFB
    TYPE(ESMF_config)             :: config
    CHARACTER(len=ESMF_MAXSTR)    :: OM_name, WM_name, fieldName, WM_ExpSt_NameList(ListLen), WM2OM_HandleName
    TYPE(ESMF_StateItem_Flag)     :: WM_ExpSt_TypeList(ListLen)
    INTEGER                       :: OM_ImpFieldCount, WM_ExpFieldCount, ii, NumRouteHandleItems, RouteHandleIndex
    TYPE(ESMF_Field),ALLOCATABLE  :: OM_ImpFieldList(:), WM_ExpFieldList(:)
    TYPE(ESMF_RouteHandle)        :: WM2OM_regridRouteHandle
    TYPE(ESMF_TypeKind_Flag)      :: fieldTypeKind

    REAL(ESMF_KIND_R8),POINTER    :: optr(:,:),iptr(:)
    INTEGER                       :: nn
    LOGICAL                       :: verbose_coupling
    TYPE(ESMF_config)             :: FISOC_config

    rc = ESMF_FAILURE

    CALL ESMF_cplCompGet(FISOC_coupler, config=FISOC_config, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    CALL ESMF_ConfigGetAttribute(FISOC_config, verbose_coupling, label='verbose_coupling:', rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    ! Extract WM2OM regrid routehandle for regridding...
    CALL ESMF_StateGet(WM_ExpSt, "WM2OM_regridRouteHandle", WM2OM_regridRouteHandle, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    ! Extract WM field bundle for regridding...
    CALL ESMF_StateGet(WM_ExpSt, "WM export fields", WM_ExpFB, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    ! ...how many fields?...
    CALL ESMF_FieldBundleGet(WM_ExpFB, fieldCount=WM_ExpFieldCount, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    ! ... get list of fields from bundle.
    ALLOCATE(WM_ExpFieldList(WM_ExpFieldCount))
    CALL ESMF_FieldBundleGet(WM_ExpFB, fieldCount=WM_ExpFieldCount, fieldList=WM_ExpFieldList,rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    IF (verbose_coupling) THEN
       msg = "coupler extracted WM fields from WM export state"
       CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_INFO, &
            line=__LINE__, file=__FILE__, rc=rc)
    END IF

    IF (SIZE(WM_ExpFieldList).LT.1) THEN
       msg = "WM field list less than length 1"
       CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_ERROR, &
            line=__LINE__, file=__FILE__, rc=rc)
       CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    END IF

    ! Extract OM field bundle for regridding...
    CALL ESMF_StateGet(OM_ImpSt, "OM import fields", OM_ImpFB, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    ! ...how many fields?...
    CALL ESMF_FieldBundleGet(OM_ImpFB, fieldCount=OM_ImpFieldCount, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    ! ... get list of fields from bundle.
    ALLOCATE(OM_ImpFieldList(OM_ImpFieldCount))
    CALL ESMF_FieldBundleGet(OM_ImpFB, fieldList=OM_ImpFieldList,rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    loop_over_fields: DO ii = 1,WM_ExpFieldCount 

       CALL ESMF_FieldGet(WM_ExpFieldList(ii), name=fieldName, typekind=fieldTypeKind, rc=rc)
       IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__)) &
            CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

       CALL ESMF_FieldRegrid(WM_ExpFieldList(ii),OM_ImpFieldList(ii), &
            routehandle=WM2OM_regridRouteHandle, zeroregion= ESMF_REGION_TOTAL, &
            checkflag=.TRUE.,rc=rc)
       IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__)) &
            CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

       IF (verbose_coupling) THEN
          msg = "Regridded field "//fieldName
          CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_INFO, &
               line=__LINE__, file=__FILE__, rc=rc)
       END IF

    END DO loop_over_fields

    IF (verbose_coupling) THEN
       msg = "Regriding complete. Regridded fields stored in OM import state"
       CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_INFO, &
            line=__LINE__, file=__FILE__, rc=rc)
    END IF

    rc = ESMF_SUCCESS

  END SUBROUTINE FISOC_coupler_run_phase2



  !------------------------------------------------------------------------------
  SUBROUTINE FISOC_coupler_finalise(FISOC_coupler, dummy_ExpSt, dummy_ImpSt, FISOC_clock, rc)
    TYPE(ESMF_CplComp)     :: FISOC_coupler
    TYPE(ESMF_State)       :: dummy_ExpSt, dummy_ImpSt
    TYPE(ESMF_Clock)       :: FISOC_clock
    INTEGER, INTENT(OUT)   :: rc

    rc = ESMF_FAILURE

    msg = "FISOC coupler finalise.  Empty routine."
    CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_INFO, &
       line=__LINE__, file=__FILE__, rc=rc)

    rc = ESMF_SUCCESS
    
  END SUBROUTINE FISOC_coupler_finalise


  !------------------------------------------------------------------------------
  TYPE(ESMF_grid) FUNCTION dummyCreateGrid(mesh, rc)

    TYPE(ESMF_mesh),INTENT(IN) :: mesh 
    INTEGER, INTENT(OUT)       :: rc

    TYPE(ESMF_grid)        :: grid 

    rc = ESMF_FAILURE

    dummyCreateGrid = grid

    rc = ESMF_SUCCESS

  END FUNCTION  dummyCreateGrid

END MODULE FISOC_coupler_MOD
