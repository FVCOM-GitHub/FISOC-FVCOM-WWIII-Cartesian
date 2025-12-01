
MODULE  FISOC_parent_MOD
  
  USE ESMF

  USE FISOC_WM_MOD, ONLY : FISOC_WM_setvm,FISOC_WM_register
  USE FISOC_coupler_MOD, ONLY : FISOC_coupler_register
  USE FISOC_OM_MOD, ONLY : FISOC_OM_setvm,FISOC_OM_register
  USE FISOC_utils_MOD
  USE FISOC_types_MOD
  
  IMPLICIT NONE
  
  PRIVATE
  
  PUBLIC FISOC_parent_register
  
  TYPE(ESMF_GridComp), SAVE :: FISOC_WM, FISOC_OM
  TYPE(ESMF_CplComp),  SAVE :: FISOC_coupler
  TYPE(ESMF_State),    SAVE :: WM_ImpSt, WM_ExpSt, OM_ImpSt, OM_ExpSt
  REAL,                SAVE :: OM_time=0.0, WM_time=0.0, AM_time=0.0
  LOGICAL,             SAVE :: WM_UseOMGrid, OM_UseWMGrid

CONTAINS
  
  !------------------------------------------------------------------------------
  SUBROUTINE FISOC_parent_register(FISOC_parent, rc)
    
    TYPE(ESMF_GridComp)  :: FISOC_parent
    INTEGER, INTENT(OUT) :: rc
    
    rc = ESMF_FAILURE

    CALL ESMF_GridCompSetEntryPoint(FISOC_parent, ESMF_METHOD_INITIALIZE, &
         userRoutine=FISOC_init, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    
    CALL ESMF_GridCompSetEntryPoint(FISOC_parent, ESMF_METHOD_RUN, &
         userRoutine=FISOC_run, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    
    CALL ESMF_GridCompSetEntryPoint(FISOC_parent, ESMF_METHOD_FINALIZE, &
         userRoutine=FISOC_finalize, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    
    rc = ESMF_SUCCESS
 
  END SUBROUTINE FISOC_parent_register
  
    
  !------------------------------------------------------------------------------
  SUBROUTINE FISOC_init(FISOC_parent, importState, exportState, FISOC_clock, rc)

    TYPE(ESMF_GridComp)  :: FISOC_parent
    TYPE(ESMF_State)     :: importState, exportState
    TYPE(ESMF_Clock)     :: FISOC_clock
    INTEGER, INTENT(OUT) :: rc

    TYPE(ESMF_State)     :: WM_ImpSt, WM_ExpSt, OM_ImpSt, OM_ExpSt
    TYPE(ESMF_config)    :: config
    INTEGER              :: petCount, localrc, urc, localPet
    INTEGER              :: WM_partitions, OM_partitions
    INTEGER,ALLOCATABLE  :: WM_PetList(:), OM_PetList(:)
    TYPE(ESMF_VM)        :: VM
    LOGICAL              :: verbose_coupling
    character(len=ESMF_MAXSTR) :: parallel_mode
    TYPE(ESMF_Context_Flag) :: OM_ContextFlag, WM_ContextFlag

    rc = ESMF_FAILURE

    CALL ESMF_GridCompGet(FISOC_parent, config=config, localPet=localPet, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    CALL ESMF_ConfigGetAttribute(config, verbose_coupling, label='verbose_coupling:', rc=rc)
    CALL ESMF_ConfigGetAttribute(config, parallel_mode, label='parallel_mode:', rc=rc)

    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    
    IF ((verbose_coupling).AND.(localPet.EQ.0)) THEN
       PRINT*,""
       PRINT*,"************************************************************************************"
       PRINT*,"This is a verbose run, set by the verbose_coupling flag in the FISOC_config.rc file."
       PRINT*,"(this mainly affects printing to screen, but also increases verbosity of the log    "
       PRINT*,"file(s), and dumping of vtk grid files.                                             "
       PRINT*,"************************************************************************************"
       PRINT*,""
    END IF

    IF ((verbose_coupling).AND.(localPet.EQ.0)) THEN
       PRINT*,""
       PRINT*,"******************************************************************************"
       PRINT*,"************    FISOC parent.  Initialise method.       **********************"
       PRINT*,"******************************************************************************"
       PRINT*,""
    END IF

    msg = "Starting FISOC parent initialisation"
    CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_INFO, &
       line=__LINE__, file=__FILE__, rc=rc)

    ! Get config and petCount from the component object (pet is persistent execution thread)
    CALL ESMF_GridCompGet(FISOC_parent, config=config, petCount=petCount, vm=vm, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    
    IF ((verbose_coupling).AND.(localPet.EQ.0)) THEN
       PRINT*,"The VM object contains information about the execution environment of "
       PRINT*,"the Component.  Printing info about parent VM..."
       PRINT*,""
       CALL ESMF_VMPrint(vm, rc=rc)
       IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__)) &
            CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    END IF

    ! Run WM and OM components on a subset of the total partitions if requested.
    CALL ESMF_ConfigGetAttribute(config, WM_partitions, label='WM_partitions:', rc=rc)
    CALL ESMF_ConfigGetAttribute(config, OM_partitions, label='OM_partitions:', rc=rc)
    IF(TRIM(parallel_mode)=='Sequential')THEN
      CALL FISOC_MakePetList(VM, WM_partitions, WM_PetList, WM_ContextFlag)
      CALL FISOC_MakePetList(VM, OM_partitions, OM_PetList, OM_ContextFlag)
    ELSEIF(TRIM(parallel_mode)=='Concurrent')THEN
      CALL FISOC_MakePetList(VM, WM_partitions, WM_PetList, WM_ContextFlag)
      CALL FISOC_MakePetList(VM, OM_partitions, OM_PetList, OM_ContextFlag, WM_partitions)
    ELSE
      msg = "parallel_mode in FISOC_config.rc should be Sequential or Concurrent !"
      CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_INFO, &
           line=__LINE__, file=__FILE__, rc=rc)
      CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    END IF
    IF (localPet.EQ.0) print*,'WM_PetList',WM_PetList
    IF (localPet.EQ.0) print*,'OM_PetList',OM_PetList
    
    ! Create and register child components and routines
    FISOC_WM = ESMF_GridCompCreate(name="Wave Model", config=config, &
         petList = WM_PetList, contextflag=WM_ContextFlag,rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    CALL ESMF_GridCompSetVM(FISOC_WM, userRoutine=FISOC_WM_setvm,          &
           userRc=urc, rc=rc)   
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
       
    FISOC_OM = ESMF_GridCompCreate(name="Ocean Model", config=config, &
         petList = OM_PetList, contextflag=OM_ContextFlag,rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    CALL ESMF_GridCompSetVM(FISOC_OM, userRoutine=FISOC_OM_setvm,          &
           userRc=urc, rc=rc)       
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    FISOC_coupler = ESMF_CplCompCreate(name="FISOC coupler", config=config, &
         contextflag=ESMF_CONTEXT_PARENT_VM,rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    
    CALL ESMF_GridCompSetServices(FISOC_WM, FISOC_WM_register, &
         userRc=urc, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    IF (ESMF_LogFoundError(rcToCheck=urc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
       
    CALL ESMF_GridCompSetServices(FISOC_OM, FISOC_OM_register, &
         userRc=urc, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    IF (ESMF_LogFoundError(rcToCheck=urc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
       
    CALL ESMF_CplCompSetServices(FISOC_coupler, FISOC_coupler_register, &
         userRc=urc, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    IF (ESMF_LogFoundError(rcToCheck=urc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    msg = "FISOC child routines created and registered"
    CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_INFO, &
       line=__LINE__, file=__FILE__, rc=rc)

    ! Check whether regridding is needed... not needed if components are on the same grid.
    CALL FISOC_ConfigDerivedAttribute(config, WM_UseOMGrid, 'WM_UseOMGrid',rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    CALL FISOC_ConfigDerivedAttribute(config, OM_UseWMGrid, 'OM_UseWMGrid',rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    IF (OM_UseWMGrid) THEN
      msg = "OM_UseWMGrid set to true but this option is NYI"
      CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_ERROR, &
           line=__LINE__, file=__FILE__)
           CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    END IF

    IF (OM_UseWMGrid.AND.WM_UseOMGrid) THEN
      msg = "OM_UseWMGrid and WM_UseOMGrid are both set to true."
      CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_ERROR, &
           line=__LINE__, file=__FILE__)
      CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    END IF

    ! If WM and OM are on the same grid, the export states will be imported directly into the 
    ! opposite component.  So we don't need to create separate import states.
    IF ( (.NOT.OM_UseWMGrid) .AND. (.NOT.WM_UseOMGrid) ) THEN
      WM_ImpSt = ESMF_StateCreate(name='WM import state', stateintent=ESMF_STATEINTENT_IMPORT, rc=rc) 
      IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
           line=__LINE__, file=__FILE__)) &
           CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
      
      OM_ImpSt = ESMF_StateCreate(name='OM import state', stateintent=ESMF_STATEINTENT_IMPORT, rc=rc) 
      IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
           line=__LINE__, file=__FILE__)) &
           CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    END IF

    WM_ExpSt = ESMF_StateCreate(name='WM export state', stateintent=ESMF_STATEINTENT_EXPORT, rc=rc) 
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    OM_ExpSt = ESMF_StateCreate(name='OM export state', stateintent=ESMF_STATEINTENT_EXPORT, rc=rc) 
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    msg = "Empty import and export states created"
    CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_INFO, &
       line=__LINE__, file=__FILE__, rc=rc)

    msg = "FISOC parent calling OM init phase 1"
    CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_INFO, &
       line=__LINE__, file=__FILE__, rc=rc)
    IF ( (OM_UseWMGrid) .OR. (WM_UseOMGrid) ) THEN
       CALL ESMF_GridCompInitialize(FISOC_OM, &
            importState=WM_ExpSt, exportState=OM_ExpSt, &
            clock=FISOC_clock, phase=1, rc=rc, userRc=urc)
    ELSE
       CALL ESMF_GridCompInitialize(FISOC_OM, &
            importState=OM_ImpSt, exportState=OM_ExpSt, &
            clock=FISOC_clock, phase=1, rc=rc, userRc=urc)
    END IF
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    IF (ESMF_LogFoundError(rcToCheck=urc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    
!    IF (WM_UseOMGrid) THEN
!      ! Copy OM fields to WM state so that WM can access OM grid.
!      ! The WM will remove the fields after accessing the grid.
!      CALL FISOC_State2StateCopyFB(OM_ExpSt,WM_ImpSt,'OM export fields',rc=rc)
!      IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
!           line=__LINE__, file=__FILE__)) &
!           CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
!    END IF

    msg = "FISOC parent calling WM init phase 1"
    CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_INFO, &
       line=__LINE__, file=__FILE__, rc=rc)
    IF ( (OM_UseWMGrid) .OR. (WM_UseOMGrid) ) THEN
       CALL ESMF_GridCompInitialize(FISOC_WM, &
            importState=OM_ExpSt, exportState=WM_ExpSt, &
            clock=FISOC_clock, phase=1, rc=rc, userRc=urc)
    ELSE 
       CALL ESMF_GridCompInitialize(FISOC_WM, &
            importState=WM_ImpSt, exportState=WM_ExpSt, &
            clock=FISOC_clock, phase=1, rc=rc, userRc=urc)
    END IF
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    IF (ESMF_LogFoundError(rcToCheck=urc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
       
    CALL ESMF_VMBarrier(vm, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)    
    
    IF ( (.NOT.WM_UseOMGrid) .AND. (.NOT.OM_UseWMGrid) ) THEN

       msg = "FISOC parent reconciling states prior to regridding"
       CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_INFO, &
            line=__LINE__, file=__FILE__, rc=rc)
       CALL ESMF_StateReconcile(WM_ExpSt, vm=vm, rc=rc)
       IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__)) &
            CALL ESMF_Finalize(endflag=ESMF_END_ABORT)    
       CALL ESMF_StateReconcile(OM_ImpSt, vm=vm, rc=rc)
       IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__)) &
            CALL ESMF_Finalize(endflag=ESMF_END_ABORT)    
       
       msg = "FISOC parent calling coupler init phase 1"
       CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_INFO, &
            line=__LINE__, file=__FILE__, rc=rc)
       CALL ESMF_cplCompInitialize(FISOC_coupler, &
            importState=WM_ExpSt, exportState=OM_ImpSt, &
            clock=FISOC_clock, phase=1, rc=rc, userRc=urc)
       IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__)) &
            CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
       IF (ESMF_LogFoundError(rcToCheck=urc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__)) &
            CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    END IF
    
    msg = "FISOC parent calling OM init phase 2"
    CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_INFO, &
       line=__LINE__, file=__FILE__, rc=rc)
    IF ( (OM_UseWMGrid) .OR. (WM_UseOMGrid) ) THEN
      CALL ESMF_GridCompInitialize(FISOC_OM, &
           importState=WM_ExpSt, exportState=OM_ExpSt, &
           clock=FISOC_clock, phase=2, rc=rc, userRc=urc)
    ELSE
      CALL ESMF_GridCompInitialize(FISOC_OM, &
           importState=OM_ImpSt, exportState=OM_ExpSt, &
           clock=FISOC_clock, phase=2, rc=rc, userRc=urc)
    END IF
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    IF (ESMF_LogFoundError(rcToCheck=urc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    CALL ESMF_VMBarrier(vm, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)    

    IF ( (.NOT.WM_UseOMGrid) .AND. (.NOT.OM_UseWMGrid) ) THEN

       msg = "FISOC parent reconciling states prior to regridding"
       CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_INFO, &
            line=__LINE__, file=__FILE__, rc=rc)
       CALL ESMF_StateReconcile(OM_ExpSt, rc=rc)
       IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__)) &
            CALL ESMF_Finalize(endflag=ESMF_END_ABORT)    
       CALL ESMF_StateReconcile(WM_ImpSt, rc=rc)
       IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__)) &
            CALL ESMF_Finalize(endflag=ESMF_END_ABORT)    

       msg = "FISOC parent calling coupler init phase 2"
       CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_INFO, &
            line=__LINE__, file=__FILE__, rc=rc)
       CALL ESMF_cplCompInitialize(FISOC_coupler, &
            importState=OM_ExpSt, exportState=WM_ImpSt, &
            clock=FISOC_clock, phase=2, rc=rc, userRc=urc)
       IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__)) &
            CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
       IF (ESMF_LogFoundError(rcToCheck=urc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__)) &
            CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    END IF

    msg = "FISOC parent calling WM init phase 2"
    CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_INFO, &
       line=__LINE__, file=__FILE__, rc=rc)
    IF ( (OM_UseWMGrid) .OR. (WM_UseOMGrid) ) THEN
      CALL ESMF_GridCompInitialize(FISOC_WM, &
           importState=OM_ExpSt, exportState=WM_ExpSt, &
           clock=FISOC_clock, phase=2, rc=rc, userRc=urc)
    ELSE
      CALL ESMF_GridCompInitialize(FISOC_WM, &
           importState=WM_ImpSt, exportState=WM_ExpSt, &
           clock=FISOC_clock, phase=2, rc=rc, userRc=urc)
    END IF
    IF(ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    IF (ESMF_LogFoundError(rcToCheck=urc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    msg = "FISOC initialise completed"
    CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_INFO, &
       line=__LINE__, file=__FILE__, rc=rc)

    ! put child states into the parent import state just so we can make them available during run phases
    IF ( (OM_UseWMGrid) .OR. (WM_UseOMGrid) ) THEN
      CALL ESMF_StateAdd(importstate, (/OM_ExpSt, WM_ExpSt/), rc=rc)
    ELSE
      CALL ESMF_StateAdd(importstate, (/OM_ImpSt, OM_ExpSt, WM_ImpSt, WM_ExpSt/), rc=rc)
    END IF
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    
    rc = ESMF_SUCCESS

    
  END SUBROUTINE FISOC_init
  
  !------------------------------------------------------------------------------
  SUBROUTINE FISOC_run(FISOC_parent, importState, exportState, FISOC_clock, rc)
    TYPE(ESMF_GridComp)  :: FISOC_parent
    TYPE(ESMF_State)     :: importState, exportState
    TYPE(ESMF_Clock)     :: FISOC_clock
    INTEGER, INTENT(OUT) :: rc

    TYPE(ESMF_State)     :: WM_ImpSt, WMExpSt, OM_ImpSt, OM_ExpSt
    TYPE(ESMF_Alarm)     :: alarm_OM, alarm_WM
    TYPE(ESMF_Alarm)     :: alarm_WM_exportAvailable, alarm_OM_output
    LOGICAL              :: verbose_coupling, profiling
    TYPE(ESMF_config)    :: config
    INTEGER              :: urc, localPet
    TYPE(ESMF_VM)        :: vm
    INTEGER(ESMF_KIND_I8):: advanceCount
    CHARACTER(len=ESMF_MAXSTR) :: OM_time_str, WM_time_str
    REAL                 :: preCallTime, postCallTime

    rc = ESMF_FAILURE

    CALL ESMF_GridCompGet(FISOC_parent, config=config, vm=vm, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    CALL ESMF_VMGet(vm, localPet=localPet, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    CALL ESMF_VMBarrier(vm, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)    

    CALL FISOC_ConfigDerivedAttribute(config, profiling, 'profiling', rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    CALL ESMF_ConfigGetAttribute(config, verbose_coupling, label='verbose_coupling:', rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    IF ((verbose_coupling).AND.(localPet.EQ.0)) THEN
       PRINT*,""
       PRINT*,"******************************************************************************"
       PRINT*,"************    FISOC parent.  Run method.              **********************"
       PRINT*,"******************************************************************************"
       PRINT*,""
    END IF

    IF ( (.NOT.WM_UseOMGrid).AND.(.NOT.OM_UseWMGrid) ) THEN
       CALL ESMF_StateGet(importstate, "WM import state", WM_ImpSt, rc=rc)
       IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__)) &
            CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

       CALL ESMF_StateGet(importstate, "OM import state", OM_ImpSt, rc=rc)
       IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__)) &
            CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    END IF

    CALL ESMF_StateGet(importstate, "WM export state", WM_ExpSt, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    CALL ESMF_StateGet(importstate, "OM export state", OM_ExpSt, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    CALL ESMF_ClockGetAlarm(FISOC_clock, "alarm_WM", alarm_WM, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    CALL ESMF_ClockGetAlarm(FISOC_clock, "alarm_OM", alarm_OM, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    msg = "FISOC parent run: got alarms and child states, now start timestepping"
    CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_INFO, &
         line=__LINE__, file=__FILE__, rc=rc)

    CALL ESMF_VMBarrier(vm, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)    

    mainTimeStepping: DO WHILE (.NOT. ESMF_ClockIsStopTime(FISOC_clock, rc=rc))

       IF ((verbose_coupling).AND.(localPet.EQ.0)) THEN
          PRINT*,""
          CALL ESMF_ClockPrint(FISOC_clock, options="advanceCount string isofrac", rc=rc)
          IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
               line=__LINE__, file=__FILE__)) &
               CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
          CALL ESMF_ClockGetAlarm(FISOC_clock, "alarm_OM_output", alarm_OM_output, rc=rc)
          IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
               line=__LINE__, file=__FILE__)) &
               CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
          CALL ESMF_ClockGetAlarm(FISOC_clock, "alarm_WM_exportAvailable", alarm_WM_exportAvailable, rc=rc)
          IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
               line=__LINE__, file=__FILE__)) &
               CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
          PRINT*,"Alarm status: ",ESMF_AlarmIsRinging(alarm_OM),&
               ESMF_AlarmIsRinging(alarm_OM_output),&
               ESMF_AlarmIsRinging(alarm_WM),&
               ESMF_AlarmIsRinging(alarm_WM_exportAvailable)
          PRINT*,""
       END IF

       msg = "FISOC parent starting Ocean Component"
       CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_INFO, &
            line=__LINE__, file=__FILE__, rc=rc)
       IF (profiling) CALL cpu_time(preCallTime)
       CALL ESMF_GridCompRun(FISOC_OM, importState=OM_ImpSt, exportState=OM_ExpSt, &
            clock=FISOC_clock, phase=1, rc=rc, userRc=urc)
       IF (profiling) THEN
          CALL cpu_time(postCallTime)
          OM_time = OM_time + (postCallTime - preCallTime)
       END IF

       msg = "FISOC parent starting Wave Component"
       CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_INFO, &
            line=__LINE__, file=__FILE__, rc=rc)
       IF (profiling) CALL cpu_time(preCallTime)
       CALL ESMF_GridCompRun(FISOC_WM, importState=WM_ImpSt, exportState=WM_ExpSt, &
            clock=FISOC_clock, phase=1, rc=rc, userRc=urc)
       IF (profiling) THEN
          CALL cpu_time(postCallTime)
          WM_time = WM_time + (postCallTime - preCallTime)
       END IF

       msg = "FISOC parent starting cplCompRun "
       CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_INFO, &
            line=__LINE__, file=__FILE__, rc=rc)
       CALL ESMF_cplCompRun(FISOC_coupler, &
            importState=OM_ExpSt, exportState=WM_ImpSt, &
            clock=FISOC_clock, phase=1, rc=rc, userRc=urc)
       IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__)) &
            CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
       IF (ESMF_LogFoundError(rcToCheck=urc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__)) &
            CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

       CALL ESMF_cplCompRun(FISOC_coupler, &
            importState=WM_ExpSt, exportState=OM_ImpSt, &
            clock=FISOC_clock, phase=2, rc=rc, userRc=urc)
       IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__)) &
            CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
       IF (ESMF_LogFoundError(rcToCheck=urc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__)) &
            CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

       CALL ESMF_ClockAdvance(FISOC_clock, rc=rc)   
       IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__)) &
            CALL ESMF_Finalize(endflag=ESMF_END_ABORT)    

       CALL ESMF_ClockGet(FISOC_clock, currTime=FISOC_time, rc=rc)
       IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__)) &
            CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

       CALL ESMF_VMBarrier(vm, rc=rc)
       IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__)) &
            CALL ESMF_Finalize(endflag=ESMF_END_ABORT)    

       IF (profiling) THEN
          WRITE(OM_time_str, *) OM_time
          WRITE(WM_time_str, *) WM_time
          msg = "OM elapsed CPU time is "//OM_time_str
          CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_INFO, &
               line=__LINE__, file=__FILE__, rc=rc)
          msg = "WM elapsed CPU time is "//WM_time_str
          CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_INFO, &
               line=__LINE__, file=__FILE__, rc=rc)
       END IF

       CALL ESMF_LogFlush(rc=rc)
       IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
            line=__LINE__, file=__FILE__)) &
            CALL ESMF_Finalize(endflag=ESMF_END_ABORT)    

    END DO mainTimeStepping

    msg = "FISOC parent run: complete"
    CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_INFO, &
         line=__LINE__, file=__FILE__, rc=rc)

  END SUBROUTINE FISOC_run


  !------------------------------------------------------------------------------
  ! get all states, finalise child components, destroy states.
  SUBROUTINE FISOC_finalize(FISOC_parent, importState, exportState, FISOC_clock, rc)
    TYPE(ESMF_GridComp)  :: FISOC_parent
    TYPE(ESMF_State)     :: importState, exportState
    TYPE(ESMF_Clock)     :: FISOC_clock
    INTEGER, INTENT(OUT) :: rc

    INTEGER              :: urc
    TYPE(ESMF_State)     :: WM_ImpSt, WMExpSt, OM_ImpSt, OM_ExpSt
    TYPE(ESMF_config)    :: config
    LOGICAL              :: verbose_coupling
    INTEGER              :: localPet
    TYPE(ESMF_VM)        :: vm

    rc = ESMF_FAILURE

    CALL ESMF_GridCompGet(FISOC_WM, config=config, vm=vm, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    CALL ESMF_VMGet(vm, localPet=localPet, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    CALL ESMF_ConfigGetAttribute(config, verbose_coupling, label='verbose_coupling:', rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    IF ((verbose_coupling).AND.(localPet.EQ.0)) THEN
       PRINT*,""
       PRINT*,"******************************************************************************"
       PRINT*,"************    FISOC parent.  Finalise method.         **********************"
       PRINT*,"******************************************************************************"
       PRINT*,""
    END IF

    msg = "FISOC parent finalise: getting child states"
    CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_INFO, &
         line=__LINE__, file=__FILE__, rc=rc)

    IF ( (.NOT.WM_UseOMGrid) .AND. (.NOT.OM_UseWMGrid) ) THEN
      CALL ESMF_StateGet(importstate, "WM import state", WM_ImpSt, rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
           line=__LINE__, file=__FILE__)) &
           CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
      
      CALL ESMF_StateGet(importstate, "OM import state", OM_ImpSt, rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
           line=__LINE__, file=__FILE__)) &
           CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    END IF

    CALL ESMF_StateGet(importstate, "WM export state", WM_ExpSt, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    CALL ESMF_StateGet(importstate, "OM export state", OM_ExpSt, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    msg = "FISOC parent finalise: calling child finalise routines"
    CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_INFO, &
         line=__LINE__, file=__FILE__, rc=rc)

    IF ( WM_UseOMGrid .OR. OM_UseWMGrid ) THEN
      CALL ESMF_GridCompFinalize(FISOC_OM, &
           importState=WM_ExpSt, exportState=OM_ExpSt, &
           clock=FISOC_clock, phase=1, rc=rc, userRc=urc)
    ELSE
      CALL ESMF_GridCompFinalize(FISOC_OM, &
           importState=OM_ImpSt, exportState=OM_ExpSt, &
           clock=FISOC_clock, phase=1, rc=rc, userRc=urc)
    END IF
    IF (ESMF_LogFoundError(rcToCheck=urc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)

    IF ( WM_UseOMGrid .OR. OM_UseWMGrid ) THEN
      CALL ESMF_GridCompFinalize(FISOC_WM, &
           importState=OM_ExpSt, exportState=WM_ExpSt, &
           clock=FISOC_clock, phase=1, rc=rc, userRc=urc)
    ELSE
      CALL ESMF_GridCompFinalize(FISOC_WM, &
           importState=WM_ImpSt, exportState=WM_ExpSt, &
           clock=FISOC_clock, phase=1, rc=rc, userRc=urc)
    END IF
    IF (ESMF_LogFoundError(rcToCheck=urc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    
    IF ( (.NOT.WM_UseOMGrid) .AND. (.NOT.OM_UseWMGrid) ) THEN
      CALL ESMF_cplCompFinalize(FISOC_coupler, &
           importState=OM_ExpSt, exportState=WM_ImpSt, &
           clock=FISOC_clock, phase=1, rc=rc, userRc=urc)
      IF (ESMF_LogFoundError(rcToCheck=urc, msg=ESMF_LOGERR_PASSTHRU, &
           line=__LINE__, file=__FILE__)) &
           CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
      IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
           line=__LINE__, file=__FILE__)) &
           CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    END IF
   
    msg = "FISOC parent finalise: destroying states"
    CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_INFO, &
         line=__LINE__, file=__FILE__, rc=rc)

    IF ( (.NOT.WM_UseOMGrid) .AND. (.NOT.OM_UseWMGrid) ) THEN
      CALL ESMF_StateDestroy(WM_ImpSt, rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
           line=__LINE__, file=__FILE__)) &
           CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
      CALL ESMF_StateDestroy(OM_ImpSt, rc=rc)
      IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
           line=__LINE__, file=__FILE__)) &
           CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    END IF
    CALL ESMF_StateDestroy(WM_ExpSt, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    CALL ESMF_StateDestroy(OM_ExpSt, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    CALL ESMF_StateDestroy(exportstate, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    CALL ESMF_StateDestroy(importstate, rc=rc)
    IF (ESMF_LogFoundError(rcToCheck=rc, msg=ESMF_LOGERR_PASSTHRU, &
         line=__LINE__, file=__FILE__)) &
         CALL ESMF_Finalize(endflag=ESMF_END_ABORT)
    
    msg = "FISOC parent finalise: complete"
    CALL ESMF_LogWrite(msg, logmsgFlag=ESMF_LOGMSG_INFO, &
         line=__LINE__, file=__FILE__, rc=rc)

    rc = ESMF_SUCCESS
    
  END SUBROUTINE FISOC_finalize

END MODULE FISOC_parent_MOD
