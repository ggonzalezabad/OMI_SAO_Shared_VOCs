SUBROUTINE read_pcf_file ( pge_error_status )

  USE OMSAO_precision_module
  USE OMSAO_indices_module
  USE OMSAO_errstat_module
  USE OMSAO_metadata_module,       ONLY: &
       pcf_granule_s_time,  pcf_granule_e_time, &
       n_mdata_int, mdata_integer_fields, mdata_integer_values
  USE OMSAO_parameters_module,     ONLY: zerospec_string, str_missval
  USE OMSAO_he5_module,            ONLY: &
       pge_swath_name, process_level, instrument_name, pge_version
  USE OMSAO_omidata_module,        ONLY: l1b_radiance_esdt
  USE OMSAO_variables_module,      ONLY: &
       pge_idx, pge_name, verb_thresh_char, verb_thresh_lev, orbit_number, &
       ecs_version_id, l1b_rad_filename, l1b_irrad_filename, l2_filename,  &
       static_input_fnames, have_amftable, omi_slitfunc_fname,             &
       OMBRO_amf_filename, OMSAO_solcomp_filename, voc_amf_filenames,      &
       refspecs_original, OMSAO_solmonthave_filename,                      &
       OMSAO_refseccor_filename, OMSAO_OMLER_filename,                     &
       OMSAO_refseccor_cld_filename
  USE OMSAO_prefitcol_module,     ONLY: o3_prefit_fname, bro_prefit_fname, &
  																			lqh2o_prefit_fname
  USE OMSAO_radiance_ref_module,  ONLY: l1b_radref_filename
  USE OMSAO_wfamf_module,         ONLY: wfamf_table_lun, climatology_lun,  &
       OMSAO_wfamf_table_filename, OMSAO_climatology_filename

  IMPLICIT NONE

  ! ------------------------------
  ! Name of this module/subroutine
  ! ------------------------------
  CHARACTER (LEN=13), PARAMETER :: modulename = 'read_pcf_file'

  ! ----------------
  ! Output variables
  ! ----------------
  INTEGER (KIND=i4), INTENT (INOUT) :: pge_error_status

  ! ---------------
  ! Local variables
  ! ---------------
  INTEGER   (KIND=i4)                      :: i, j, strlen
  CHARACTER (LEN=maxchlen)                 :: lunstr
  CHARACTER (LEN=PGSd_PC_VALUE_LENGTH_MAX) :: tmpchar

  ! ------------------------
  ! Error handling variables
  ! ------------------------
  INTEGER (KIND=i4) :: errstat, version  

  ! ------------------
  ! External functions
  ! ------------------
  CHARACTER (LEN=maxchlen), EXTERNAL :: int2string
  INTEGER   (KIND=i4),      EXTERNAL :: &
       pgs_pc_getnumberoffiles, pgs_pc_getreference, pgs_pc_getconfigdata, &
       pgs_smf_teststatuslevel

  errstat = pge_errstat_ok

  ! --------------------------------------------------------------------
  ! Read all ConfigData from PCF and set variables associated with them. 
  ! Some, like MoleculeID, might lead to PGE termination if in error.
  ! --------------------------------------------------------------------
  DO i = 1, n_config_luns
     errstat = &
          PGS_PC_GetConfigData ( config_lun_array(i), config_lun_values(i) )
     errstat = PGS_SMF_TestStatusLevel(errstat)
     CALL error_check ( errstat, PGS_SMF_MASK_LEV_S, pge_errstat_warning, OMSAO_W_GETLUN, &
          modulename//f_sep//TRIM(ADJUSTL(config_lun_strings(i))), vb_lev_default, &
          pge_error_status )
     ! -------------------------------------------------------------
     ! Remove any quotes (") that might have made in into the string
     ! -------------------------------------------------------------
     CALL remove_quotes ( config_lun_values(i) )

     SELECT CASE ( config_lun_array(i) )
        ! --------------------------------------------------------------------
        ! Get the verbosity threshold. The setting of this variable determines
        ! how "chatty" in terms of screen output the PGE will be. If the READ
        ! was in error, proceed with DEBUG level
        ! --------------------------------------------------------------------
     CASE ( verbosity_lun )
        IF ( errstat /= PGS_SMF_MASK_LEV_S ) config_lun_values(i) = "1"
        READ (config_lun_values(i), '(I1)') verb_thresh_lev
        ! ----------------------------------------------------------------------
        ! Get the orbit number. This will be checked against L1B MetaData later.
        ! Note that the source for OrbitNumber is the L1B file, NOT the PCF. But
        ! this only matters in case they are different, which will be checked
        ! later. If the ARE different, then the L1B value will take precidence.
        ! ----------------------------------------------------------------------
     CASE ( orbitnumber_lun )
        IF ( errstat /= PGS_SMF_MASK_LEV_S ) config_lun_values(i) = "00000"
        READ (config_lun_values(i), '(I5)') orbit_number

        ! -----------------------------------------------------------------
        ! Get Granule Start and End time. This is provided by the Scheduler
        ! and has to be checked against the L1B MetaData fields for
        ! RangeBeginningTime and RangeEndingTime.
        ! -----------------------------------------------------------------
     CASE ( granule_s_lun )
        IF ( errstat /= PGS_SMF_MASK_LEV_S ) config_lun_values(i) = "T00:00:00.000Z"
        pcf_granule_s_time = TRIM(ADJUSTL(config_lun_values(i)))
     CASE ( granule_e_lun )
        IF ( errstat /= PGS_SMF_MASK_LEV_S ) config_lun_values(i) = "T00:00:00.000Z"
        pcf_granule_e_time = TRIM(ADJUSTL(config_lun_values(i)))

        ! ------------------------------------------------------------------
        ! Get the Process Level
        ! ------------------------------------------------------------------
     CASE ( proclevel_lun )
        IF ( errstat /= PGS_SMF_MASK_LEV_S ) config_lun_values(i) = str_missval
        process_level = TRIM(ADJUSTL(config_lun_values(i)))

        ! ------------------------------------------------------------------
        ! Get the PGE Version
        ! ------------------------------------------------------------------
     CASE ( pge_version_lun )
        IF ( errstat /= PGS_SMF_MASK_LEV_S ) config_lun_values(i) = str_missval
        pge_version = TRIM(ADJUSTL(config_lun_values(i)))

        ! ------------------------------------------------------------------
        ! Get the Instrument Name
        ! ------------------------------------------------------------------
     CASE ( instrument_name_lun )
        IF ( errstat /= PGS_SMF_MASK_LEV_S ) config_lun_values(i) = str_missval
        instrument_name = TRIM(ADJUSTL(config_lun_values(i)))

        ! ------------------------------------------------------------------
        ! Get the HE5 Swath Name; set to MISSING VALUE if it can't be found.
        ! ------------------------------------------------------------------
     CASE ( swathname_lun )
        IF ( errstat /= PGS_SMF_MASK_LEV_S ) config_lun_values(i) = str_missval
        pge_swath_name = TRIM(ADJUSTL(config_lun_values(i)))

        ! ----------------------------------------------------------------------
        ! Get ECS Collection version number and assign it to a Metadata variable
        ! ----------------------------------------------------------------------
     CASE ( versionid_lun )
        READ (config_lun_values(i), '(I1)') ecs_version_id
        getidx: DO j = 1, n_mdata_int
           IF ( mdata_integer_fields(3,j) == "pcf"              .AND. ( &
                INDEX(mdata_integer_fields(1,j), 'VERSIONID') /= 0 .OR. &
                INDEX(mdata_integer_fields(1,j), 'VersionID') /= 0     )  ) THEN
              mdata_integer_values(j) = ecs_version_id
              EXIT getidx
           END IF
        END DO getidx

        ! -------------------------------------------------------------------------
        ! Get the SAO PGE Name string. This string is converted to the PGE
        ! index number (10, 11, 12), which in turn is used to access LUNs and other
        ! elements that are PGE specific. All of this has to be done here, because
        ! we require the PGE index to identify LUNs and other PGE specific items.
        ! >>> We MUST know this, or else we can't execute the PGE <<<
        ! -------------------------------------------------------------------------
     CASE ( pge_molid_lun )
        CALL error_check ( errstat, PGS_SMF_MASK_LEV_S, pge_errstat_fatal, &
             OMSAO_F_GETLUN, modulename//f_sep//TRIM(ADJUSTL(config_lun_strings(i))), &
             vb_lev_default, pge_error_status )
        IF ( pge_error_status >= pge_errstat_error ) RETURN

        ! ---------------------------------------------------------------
        ! Translate molecule string into index (required for LUN look-up)
        ! ---------------------------------------------------------------
        errstat = pge_errstat_ok
        CALL get_pge_ident ( &
             TRIM(ADJUSTL(config_lun_values(i))), pge_name, pge_idx, errstat )
        CALL error_check ( errstat, pge_errstat_ok, pge_errstat_fatal, &
             OMSAO_F_GET_MOLINDEX, modulename, vb_lev_default, pge_error_status )
        IF ( pge_error_status >= pge_errstat_fatal ) RETURN
        CALL error_check ( 0, 1, pge_errstat_ok, OMSAO_S_GET_MOLINDEX, &
             TRIM(ADJUSTL(pge_name)), vb_lev_stmdebug, errstat )  
     END SELECT
  END DO

  errstat = pge_errstat_ok

  ! ---------------------------------------------------------
  ! Static input file with tabulated OMI slit function values
  ! ---------------------------------------------------------
  version = 1
  errstat = PGS_PC_GetReference (omi_slitfunc_lun, version, omi_slitfunc_fname)
  errstat = PGS_SMF_TestStatusLevel(errstat)
  CALL error_check ( errstat, PGS_SMF_MASK_LEV_S, pge_errstat_fatal, OMSAO_F_GETLUN, &
       modulename//f_sep//"OMI_SLITFUNC_LUN", vb_lev_default, pge_error_status )
  IF ( pge_error_status >= pge_errstat_error ) RETURN


  ! ------------------------------------------------------------------
  ! Read names of static input files from PCF. The first entry is the 
  ! algorithm control file (initial fitting values, etc.), the other
  ! entries are reference spectra.
  ! ------------------------------------------------------------------
  static_input_fnames = zerospec_string
  DO i = icf_idx, max_rs_idx
     version = 1
     errstat = PGS_PC_GetReference (pge_static_input_luns(i), version, tmpchar)
     tmpchar = TRIM(ADJUSTL(tmpchar)) ; strlen = LEN(TRIM(ADJUSTL(tmpchar)))

     errstat = PGS_SMF_TestStatusLevel(errstat)
     IF ( errstat /= pgs_smf_mask_lev_s .OR. strlen == 0 ) THEN
        lunstr = int2string ( pge_static_input_luns(i), 1 )
        CALL error_check ( 0, 1, pge_errstat_fatal, OMSAO_F_GETLUN, &
             modulename//f_sep//"PGE_STATIC_INPUT_LUN "//TRIM(ADJUSTL(lunstr)), &
             vb_lev_default, pge_error_status )
        IF ( pge_error_status >= pge_errstat_error ) RETURN
     ELSE IF ( INDEX ( TRIM(ADJUSTL(tmpchar)), zerospec_string ) == 0 ) THEN
        static_input_fnames(i) = TRIM(ADJUSTL(tmpchar))
     END IF
  END DO

  ! ------------------------------------------------------
  ! Save file names with reference spectra to global array
  ! ------------------------------------------------------
  refspecs_original(1:max_rs_idx)%FileName = static_input_fnames(1:max_rs_idx)

  ! ----------------------------------------------------
  ! Read fitting conrol parameters from input file; this
  ! returns L1B_RADIANCE_ESDT, which determines the
  ! ingestion of the L1b radiance file.
  ! ----------------------------------------------------
  errstat = pge_errstat_ok
  CALL read_fitting_control_file ( pge_idx, l1b_radiance_esdt, errstat )
  CALL error_check ( errstat, pge_errstat_ok, pge_errstat_warning, OMSAO_W_SUBROUTINE, &
       modulename//f_sep//"READ_FITTING_CONTROL_FILE.", vb_lev_default, pge_error_status )
  IF ( pge_error_status >= pge_errstat_fatal ) RETURN


  ! -----------------------------
  ! Read Irradiance L1B file name
  ! -----------------------------
  version = 1
  errstat = PGS_PC_GetReference (l1b_irradiance_lun, version, l1b_irrad_filename)
  errstat = PGS_SMF_TestStatusLevel(errstat)
  CALL error_check ( errstat, PGS_SMF_MASK_LEV_S, pge_errstat_fatal, OMSAO_F_GETLUN, &
       modulename//f_sep//"L1B_IRRADIANCE_LUN", vb_lev_default, pge_error_status )
  IF ( pge_error_status >= pge_errstat_error ) RETURN


  ! ---------------------------
  ! Read Radiance L1B file name
  ! ---------------------------
  version = 1
  errstat = PGS_PC_GetReference (l1b_radiance_lun, version, l1b_rad_filename)
  errstat = PGS_SMF_TestStatusLevel(errstat)
  CALL error_check ( errstat, PGS_SMF_MASK_LEV_S, pge_errstat_fatal, OMSAO_F_GETLUN, &
       modulename//f_sep//"L1B_RADIANCE_LUN", vb_lev_default, pge_error_status )
  IF ( pge_error_status >= pge_errstat_error ) RETURN

  l1b_rad_filename = TRIM(ADJUSTL(l1b_rad_filename))


  ! -------------------------------------
  ! Read Radiance Reference L1B file name
  ! -------------------------------------
  version = 1
  errstat = PGS_PC_GetReference (l1b_radianceref_lun, version, l1b_radref_filename)
  errstat = PGS_SMF_TestStatusLevel(errstat)
  CALL error_check ( errstat, PGS_SMF_MASK_LEV_S, pge_errstat_fatal, OMSAO_F_GETLUN, &
       modulename//f_sep//"L1B_RADIANCEREF_LUN", vb_lev_default, pge_error_status )
  IF ( pge_error_status >= pge_errstat_error ) RETURN

  l1b_radref_filename = TRIM(ADJUSTL(l1b_radref_filename))


  ! -------------------------------------------------------------------------
  ! Read name of AMF table file(s). Remember that a missing AMF table is
  ! not a fatal problem, since in that case the slant columns will be written
  ! to the output file.
  ! -------------------------------------------------------------------------
  have_amftable = .TRUE.
  SELECT CASE ( pge_idx )
  CASE ( pge_bro_idx )
     version = 1
     errstat = PGS_PC_GetReference ( OMBRO_amf_lun, version, tmpchar)
     tmpchar = TRIM(ADJUSTL(tmpchar)) ; strlen = LEN(TRIM(ADJUSTL(tmpchar)))
     errstat = PGS_SMF_TestStatusLevel(errstat)
     IF ( (errstat /= pgs_smf_mask_lev_s) .OR. (strlen == 0)  ) THEN
        lunstr = int2string ( OMBRO_amf_lun, 1 )
        CALL error_check ( 0, 1, pge_errstat_warning, OMSAO_W_GETLUN, &
             modulename//f_sep//"PGE_STATIC_INPUT_LUN "//TRIM(ADJUSTL(lunstr)), &
             vb_lev_default, pge_error_status )
        have_amftable = .FALSE.
     ELSE
        OMBRO_amf_filename = TRIM(ADJUSTL(tmpchar))
     END IF
  CASE ( pge_hcho_idx )
     DO i = 1, n_voc_amf_luns
        version = 1
        errstat = PGS_PC_GetReference ( voc_amf_luns(i), version, tmpchar)
        tmpchar = TRIM(ADJUSTL(tmpchar)) ; strlen = LEN(TRIM(ADJUSTL(tmpchar)))
        errstat = PGS_SMF_TestStatusLevel(errstat)
        IF ( (errstat /= pgs_smf_mask_lev_s) .OR. (strlen == 0)  ) THEN
           lunstr = int2string ( voc_amf_luns(i), 1 )
           CALL error_check ( 0, 1, pge_errstat_warning, OMSAO_W_GETLUN, &
                modulename//f_sep//"PGE_STATIC_INPUT_LUN "//TRIM(ADJUSTL(lunstr)), &
                vb_lev_default, pge_error_status )
          have_amftable = .FALSE.
        ELSE
           voc_amf_filenames(i) = TRIM(ADJUSTL(tmpchar))
        END IF
     END DO
  CASE ( pge_gly_idx )
     DO i = 1, n_voc_amf_luns
        version = 1
        errstat = PGS_PC_GetReference ( voc_amf_luns(i), version, tmpchar)
        tmpchar = TRIM(ADJUSTL(tmpchar)) ; strlen = LEN(TRIM(ADJUSTL(tmpchar)))
        errstat = PGS_SMF_TestStatusLevel(errstat)
        IF ( (errstat /= pgs_smf_mask_lev_s) .OR. (strlen == 0)  ) THEN
           lunstr = int2string ( voc_amf_luns(i), 1 )
           CALL error_check ( 0, 1, pge_errstat_warning, OMSAO_W_GETLUN, &
                modulename//f_sep//"PGE_STATIC_INPUT_LUN "//TRIM(ADJUSTL(lunstr)), &
                vb_lev_default, pge_error_status )
          have_amftable = .FALSE.
        ELSE
           voc_amf_filenames(i) = TRIM(ADJUSTL(tmpchar))
        END IF
     END DO
  CASE ( pge_h2o_idx )
     DO i = 1, n_voc_amf_luns
        version = 1
        errstat = PGS_PC_GetReference ( voc_amf_luns(i), version, tmpchar)
        tmpchar = TRIM(ADJUSTL(tmpchar)) ; strlen = LEN(TRIM(ADJUSTL(tmpchar)))
        errstat = PGS_SMF_TestStatusLevel(errstat)
        IF ( (errstat /= pgs_smf_mask_lev_s) .OR. (strlen == 0)  ) THEN
           lunstr = int2string ( voc_amf_luns(i), 1 )
           CALL error_check ( 0, 1, pge_errstat_warning, OMSAO_W_GETLUN, &
                modulename//f_sep//"PGE_STATIC_INPUT_LUN "//TRIM(ADJUSTL(lunstr)), &
                vb_lev_default, pge_error_status )
          have_amftable = .FALSE.
        ELSE
           voc_amf_filenames(i) = TRIM(ADJUSTL(tmpchar))
        END IF
     END DO

  END SELECT

  ! -------------------------------------------------------------------------
  ! gga. The new implementation of the wf amf is intended to be molecule inde
  ! pendent. Therfore I read the file names outside the CASE stament.
  ! -------------------------------------------------------------------------
  ! wavelength dependent AMF table file
  ! -----------------------------------
  version = 1
  errstat = PGS_PC_GetReference ( wfamf_table_lun, version, tmpchar )
  tmpchar = TRIM(ADJUSTL(tmpchar)) ; strlen = LEN(TRIM(ADJUSTL(tmpchar)))
  errstat = PGS_SMF_TestStatusLevel ( errstat )
  IF ( ( errstat /= pgs_smf_mask_lev_s ) .OR. ( strlen == 0 ) ) THEN
     lunstr = int2string ( wfamf_table_lun, 1 )
     CALL error_check ( 0, 1, pge_errstat_warning, OMSAO_W_GETLUN, &
          modulename//f_sep//"PGE_STATIC_INPUT_LUN "//TRIM(ADJUSTL(lunstr)),&
          vb_lev_default, pge_error_status )
  ELSE
     OMSAO_wfamf_table_filename = TRIM ( ADJUSTL (tmpchar))
  END IF

  ! ----------------
  ! Climatology
  ! ----------------
  version = 1
  errstat = PGS_PC_GetReference ( climatology_lun, version, tmpchar )
  tmpchar = TRIM(ADJUSTL(tmpchar)) ; strlen = LEN(TRIM(ADJUSTL(tmpchar)))
  errstat = PGS_SMF_TestStatusLevel ( errstat )
  IF ( ( errstat /= pgs_smf_mask_lev_s ) .OR. ( strlen == 0 ) ) THEN
     lunstr = int2string ( climatology_lun, 1 )
     CALL error_check ( 0, 1, pge_errstat_warning, OMSAO_W_GETLUN, &
          modulename//f_sep//"PGE_STATIC_INPUT_LUN "//TRIM(ADJUSTL(lunstr)),&
          vb_lev_default, pge_error_status )
  ELSE
     OMSAO_climatology_filename = TRIM ( ADJUSTL (tmpchar))
  END IF

  ! -------------------------------------------------------------------------
  ! Read name of file with Solar Spectrum Composite 
  ! (whether we use it or not, since we find that out only after we read the
  !  fitting control file)
  ! -------------------------------------------------------------------------
  version = 1
  errstat = PGS_PC_GetReference ( OMSAO_solcomp_lun, version, tmpchar)
  tmpchar = TRIM(ADJUSTL(tmpchar)) ; strlen = LEN(TRIM(ADJUSTL(tmpchar)))
  errstat = PGS_SMF_TestStatusLevel(errstat)
  IF ( (errstat /= pgs_smf_mask_lev_s) .OR. (strlen == 0)  ) THEN
     lunstr = int2string ( OMSAO_solcomp_lun, 1 )
     CALL error_check ( 0, 1, pge_errstat_warning, OMSAO_W_GETLUN, &
          modulename//f_sep//"PGE_STATIC_INPUT_LUN "//TRIM(ADJUSTL(lunstr)), &
          vb_lev_default, pge_error_status )
  ELSE
     OMSAO_solcomp_filename = TRIM(ADJUSTL(tmpchar))
  END IF

  ! -------------------------------------------------------------------------
  ! Read name of file with Monthly Average Irradiace !gga 
  ! (whether we use it or not, since we find that out only after we read the
  !  fitting control file) !gga
  ! -------------------------------------------------------------------------
  version = 1
  errstat = PGS_PC_GetReference ( OMSAO_solmonthave_lun, version, tmpchar)
  tmpchar = TRIM(ADJUSTL(tmpchar)) ; strlen = LEN(TRIM(ADJUSTL(tmpchar)))
  errstat = PGS_SMF_TestStatusLevel(errstat)
  IF ( (errstat /= pgs_smf_mask_lev_s) .OR. (strlen == 0)  ) THEN
     lunstr = int2string ( OMSAO_solmonthave_lun, 1 )
     CALL error_check ( 0, 1, pge_errstat_warning, OMSAO_W_GETLUN, &
          modulename//f_sep//"PGE_STATIC_INPUT_LUN "//TRIM(ADJUSTL(lunstr)), &
          vb_lev_default, pge_error_status )
  ELSE
     OMSAO_solmonthave_filename = TRIM(ADJUSTL(tmpchar))
  END IF

  ! ------------------------------------------------------------
  ! Read name of file with GEOS-Chem background Reference Sector
  ! concentrations !gga 
  ! ------------------------------------------------------------
  version = 1
  errstat = PGS_PC_GetReference ( OMSAO_refseccor_lun, version, tmpchar)
  tmpchar = TRIM(ADJUSTL(tmpchar)) ; strlen = LEN(TRIM(ADJUSTL(tmpchar)))
  errstat = PGS_SMF_TestStatusLevel(errstat)
  IF ( (errstat /= pgs_smf_mask_lev_s) .OR. (strlen == 0)  ) THEN
     lunstr = int2string ( OMSAO_refseccor_lun, 1 )
     CALL error_check ( 0, 1, pge_errstat_warning, OMSAO_W_GETLUN, &
          modulename//f_sep//"PGE_STATIC_INPUT_LUN "//TRIM(ADJUSTL(lunstr)), &
          vb_lev_default, pge_error_status )
  ELSE
     OMSAO_refseccor_filename = TRIM(ADJUSTL(tmpchar))
  END IF

  ! ---------------------------------------------------------
  ! Read name of file with the radiance reference clouds !gga 
  ! ---------------------------------------------------------
  version = 1
  errstat = PGS_PC_GetReference ( OMSAO_refseccor_cld_lun, version, tmpchar)
  tmpchar = TRIM(ADJUSTL(tmpchar)) ; strlen = LEN(TRIM(ADJUSTL(tmpchar)))
  errstat = PGS_SMF_TestStatusLevel(errstat)
  IF ( (errstat /= pgs_smf_mask_lev_s) .OR. (strlen == 0)  ) THEN
     lunstr = int2string ( OMSAO_refseccor_cld_lun, 1 )
     CALL error_check ( 0, 1, pge_errstat_warning, OMSAO_W_GETLUN, &
          modulename//f_sep//"PGE_STATIC_INPUT_LUN "//TRIM(ADJUSTL(lunstr)), &
          vb_lev_default, pge_error_status )
  ELSE
     OMSAO_refseccor_cld_filename = TRIM(ADJUSTL(tmpchar))
  END IF

  ! ----------------------------------
  ! Read name of OMLER albedo file gga
  ! ----------------------------------
  version = 1
  errstat = PGS_PC_GetReference (OMSAO_OMLER_lun, version, tmpchar)
  tmpchar = TRIM(ADJUSTL(tmpchar)) ; strlen = LEN(TRIM(ADJUSTL(tmpchar)))
  errstat = PGS_SMF_TestStatusLevel(errstat)
  IF ( (errstat /= pgs_smf_mask_lev_s) .OR. (strlen == 0)  ) THEN
     lunstr = int2string ( OMSAO_OMLER_lun, 1 )
     CALL error_check ( 0, 1, pge_errstat_warning, OMSAO_W_GETLUN, &
          modulename//f_sep//"PGE_STATIC_INPUT_LUN "//TRIM(ADJUSTL(lunstr)), &
          vb_lev_default, pge_error_status )
  ELSE
     OMSAO_OMLER_filename = TRIM(ADJUSTL(tmpchar))
  END IF

  ! ---------------------------
  ! Read name of L2 output file
  ! ---------------------------
  version = 1
  errstat = PGS_PC_GetReference (pge_l2_output_lun, version, l2_filename)
  errstat = PGS_SMF_TestStatusLevel(errstat)

  lunstr = int2string ( pge_l2_output_lun, 1 )
  CALL error_check ( errstat, PGS_SMF_MASK_LEV_S, pge_errstat_fatal, OMSAO_F_GETLUN, &
       modulename//f_sep//"PGE_L2_OUTPUT_LUN "//TRIM(ADJUSTL(lunstr)), &
       vb_lev_default, pge_error_status )
  IF ( pge_error_status >= pge_errstat_error ) RETURN


  ! ---------------------------------------------------------
  ! For OMHCHO read HE5 file names with pre-fitted O3 and BrO
  ! ---------------------------------------------------------
  SELECT CASE ( pge_idx )
    
  CASE ( pge_hcho_idx )
     ! -----------
     ! O3 pre-fits
     ! -----------
     version = 1
     errstat = PGS_PC_GetReference (o3_prefit_lun, version, o3_prefit_fname)
     errstat = PGS_SMF_TestStatusLevel(errstat)
     lunstr = int2string ( o3_prefit_lun, 1 )
     CALL error_check ( errstat, PGS_SMF_MASK_LEV_S, pge_errstat_error, OMSAO_E_GETLUN, &
          modulename//f_sep//"O3_PREFIT_LUN "//TRIM(ADJUSTL(lunstr)), &
          vb_lev_default, pge_error_status )
     ! ------------
     ! BrO pre-fits
     ! ------------
     version = 1
     errstat = PGS_PC_GetReference (bro_prefit_lun, version, bro_prefit_fname)
     errstat = PGS_SMF_TestStatusLevel(errstat)
     lunstr = int2string ( bro_prefit_lun, 1 )
     CALL error_check ( errstat, PGS_SMF_MASK_LEV_S, pge_errstat_error, OMSAO_E_GETLUN, &
          modulename//f_sep//"BRO_PREFIT_LUN "//TRIM(ADJUSTL(lunstr)), &
          vb_lev_default, pge_error_status )
  
  CASE ( pge_gly_idx )
     ! ---------------------
     ! Liquid Water pre-fits
     ! ---------------------
     version = 1
     errstat = PGS_PC_GetReference (lqh2o_prefit_lun, version, lqh2o_prefit_fname)
     errstat = PGS_SMF_TestStatusLevel(errstat)
     lunstr = int2string ( lqh2o_prefit_lun, 1 )
     CALL error_check ( errstat, PGS_SMF_MASK_LEV_S, pge_errstat_error, OMSAO_E_GETLUN, &
          modulename//f_sep//"LQH2O_PREFIT_LUN "//TRIM(ADJUSTL(lunstr)), &
          vb_lev_default, pge_error_status )
          
  CASE DEFAULT
  	! Do Nothing
  END SELECT


  RETURN
END SUBROUTINE read_pcf_file
