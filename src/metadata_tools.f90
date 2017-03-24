SUBROUTINE init_metadata ( errstat )

  USE OMSAO_precision_module,   ONLY: i4
  USE OMSAO_indices_module,     ONLY: &
       l1b_irradiance_lun, l1b_radiance_lun, l1b_radianceref_lun, &
       mcf_lun, md_inventory_idx, md_archive_idx,                 &
       pge_hcho_idx, pge_gly_idx, pge_h2o_idx,                    &
       voc_amf_luns, voc_omicld_idx, o3_prefit_lun, bro_prefit_lun
       
  USE OMSAO_metadata_module
  USE OMSAO_errstat_module
  USE OMSAO_parameters_module, ONLY: str_missval, int16_missval, r8_missval
  USE OMSAO_variables_module,  ONLY: l1b_rad_filename, pge_idx, l1br_opf_version
  USE OMSAO_prefitcol_module,  ONLY: yn_o3_prefit, yn_bro_prefit, yn_lqh2o_prefit
  USE OMSAO_he5_module,        ONLY: &
       granule_day, granule_month, granule_year, TAI93At0zOfGranule, l1b_orbitdata
  IMPLICIT NONE

  INCLUDE 'PGS_TD_3.f'

  ! ------------------------------
  ! Name of this module/subroutine
  ! ------------------------------
  CHARACTER (LEN=14), PARAMETER :: modulename = "init_metadata"

  ! ---------------
  ! Output variable
  ! ---------------
  INTEGER (KIND=i4), INTENT (INOUT) :: errstat

  ! ---------------
  ! Local variables
  ! ---------------
  CHARACTER (LEN=3)                         :: mdata_typ, mdata_loc
  CHARACTER (LEN=PGSd_MET_MAX_STRING_SET_L) :: tmp_string, metadata_type
  INTEGER   (KIND=i4)                       :: &
       imd, version, estat, md_stat, locerrstat, mdata_index, jday

  ! ------------------
  ! External functions
  ! ------------------
  INTEGER (KIND=i4), EXTERNAL :: &
       PGS_TD_UTCtoTAI,      PGS_MET_init, &
       PGS_MET_getSetAttr_s, PGS_MET_getSetAttr_i, PGS_MET_getSetAttr_d, &
       PGS_MET_getPCAttr_s,  PGS_MET_getPCAttr_i,  PGS_MET_getPCAttr_d,  &
       day_of_year!, L1B_getGlobalAttrname, 

  ! ------------------------------------------------------------
  ! Variables for the extraction of the L1B radiance OPF version
  ! ------------------------------------------------------------
  INTEGER   (KIND=i4), PARAMETER                     :: n_ipp = 20
  INTEGER   (KIND=i4)                                :: idx
  CHARACTER (LEN=PGSd_MET_NAME_L), DIMENSION (n_ipp) :: l1br_inputp
  CHARACTER (LEN=4)                                  :: opf_string

  ! PGS_MET_getPCAttr_i will fail if LEN is set to anything else
  CHARACTER (LEN=PGSd_MET_NAME_L),           EXTERNAL :: upper_case


  locerrstat = pge_errstat_ok

  ! ----------------------------------------------------------------
  ! Initialize MetaData file (Example taken from PGS Toolkit Primer)
  ! ----------------------------------------------------------------
  mcf_groups = ""
  md_stat = PGS_MET_init ( mcf_lun, mcf_groups )
  CALL error_check ( md_stat, PGS_S_SUCCESS, pge_errstat_fatal, OMSAO_F_METINIT, &
       modulename, vb_lev_default, errstat )
  IF ( errstat >= pge_errstat_fatal ) RETURN

  ! ---------------------------------
  ! Initialize MetaData STRING fields
  ! ---------------------------------
  DO imd = 1, n_mdata_str

     IF ( TRIM(ADJUSTL(mdata_string_fields(3,imd))) == "mcf" .OR. &
          TRIM(ADJUSTL(mdata_string_fields(3,imd))) == "pcf" .OR. &
          TRIM(ADJUSTL(mdata_string_fields(3,imd))) == "l1i" .OR. &
          TRIM(ADJUSTL(mdata_string_fields(3,imd))) == "l1r" .OR. &
          TRIM(ADJUSTL(mdata_string_fields(3,imd))) == "l1R"        ) THEN
        version = 1
        mdata_string_values(imd) = str_missval
      
        ! Short-hand for MetaData category and location
        mdata_typ = TRIM(ADJUSTL(mdata_string_fields(2,imd)))
        mdata_loc = TRIM(ADJUSTL(mdata_string_fields(3,imd)))

        SELECT CASE ( mdata_typ )
        CASE ( "inv" )
           mdata_index = md_inventory_idx
           metadata_type = cmd_str//'.0'
        CASE ( "psa" )
           mdata_index   = md_inventory_idx
           metadata_type = cmd_str//'.0'
           ! Do we have any of these at all?
        CASE ( "arc" )
           mdata_index   = md_archive_idx
           metadata_type = amd_str//'.0'
        CASE DEFAULT
           ! No idea what to do here
        END SELECT

        SELECT CASE ( mdata_loc )
        CASE ( "mcf" )
           md_stat = PGS_MET_getSetAttr_s( mcf_groups(mdata_index), &
                TRIM(ADJUSTL(upper_case(mdata_string_fields(1,imd)))), &
                mdata_string_values(imd) )
           IF ( TRIM(ADJUSTL(mdata_string_fields(1,imd))) == "ShortName" ) &
                mcf_shortname = TRIM(ADJUSTL(mdata_string_values(imd)))
        CASE ( "l1i" )
           md_stat = PGS_MET_getPCAttr_s( &
                l1b_irradiance_lun, version, TRIM(ADJUSTL(metadata_type)), &
                TRIM(ADJUSTL(upper_case(mdata_string_fields(1,imd)))), &
                mdata_string_values(imd) )
        CASE ( "l1r" )
           SELECT CASE ( TRIM(ADJUSTL(mdata_string_fields(2,imd))) )
           CASE ( "arc")
              md_stat = PGS_MET_getPCAttr_s( &
                   l1b_radiance_lun, version, TRIM(ADJUSTL(metadata_type)), &
                   TRIM(ADJUSTL(upper_case(mdata_string_fields(1,imd)))), mdata_string_values(imd) )
           CASE ( "inv")
              md_stat = PGS_MET_getPCAttr_s( &
                   l1b_radiance_lun, version, TRIM(ADJUSTL(metadata_type)), &
                   TRIM(ADJUSTL(upper_case(mdata_string_fields(1,imd)))), &
                   mdata_string_values(imd) )
           END SELECT
           IF ( TRIM(ADJUSTL(mdata_string_fields(1,imd))) == "OrbitData" ) &
                l1b_orbitdata = TRIM(ADJUSTL(mdata_string_values(imd)))
           ! -----------------------------------------------------------
           ! Obtain YEAR, MONTH, and DAY, and convert to Day-of-the-Year
           ! -----------------------------------------------------------
           IF ( TRIM(ADJUSTL(mdata_string_fields(1,imd))) == rbd_str ) THEN
              IF ( TRIM(ADJUSTL(mdata_string_values(imd))) /= str_missval ) THEN
                 READ ( mdata_string_values(imd), '(I4,1X,I2,1X,I2)') &
                      granule_year, granule_month, granule_day
                 jday = day_of_year ( granule_year, granule_month, granule_day )

                 ! ------------------------------------------------------------
                 ! Since we are here, compute Granule Time in seconds; this is
                 ! one of the global file attributes that we need to write out.
                 ! ------------------------------------------------------------
                 tmp_string = ""
                 tmp_string = TRIM(ADJUSTL(mdata_string_values(imd))) //"T00:00:00.000Z"
                 estat = PGS_TD_UTCtoTAI ( TRIM(ADJUSTL(tmp_string)), TAI93At0zOfGranule )
                 IF ( estat /= PGS_S_SUCCESS .AND. estat /= PGSTD_E_NO_LEAP_SECS ) &
                      CALL error_check ( 0, 1, pge_errstat_warning, OMSAO_W_TAI93, &
                      modulename, vb_lev_default, errstat )
                 ! ------------------------------------------------------------
              ELSE 
                 granule_year = 0 ; granule_month = 0 ; granule_day = 0 ; jday = 0
              END IF
           END IF

        CASE ( "l1R" ) ! Radiance reference granule
           SELECT CASE ( TRIM(ADJUSTL(mdata_string_fields(2,imd))) )
           CASE ( "inv")
              md_stat = PGS_MET_getPCAttr_s( &
                   l1b_radianceref_lun, version, TRIM(ADJUSTL(metadata_type)), &
                   TRIM(ADJUSTL(upper_case(mdata_string_fields(1,imd)))), &
                   mdata_string_values(imd) )
           CASE DEFAULT
              ! Nothing to do here
           END SELECT

        CASE DEFAULT
           ! Nothing to do here
        END SELECT

        ! ------------------------------
        ! Error check for initialization
        ! ------------------------------
        CALL error_check ( md_stat, PGS_S_SUCCESS, pge_errstat_warning, OMSAO_W_GETATTR, &
             modulename//f_sep//TRIM(ADJUSTL(mdata_string_fields(1,imd))), &
             vb_lev_default, errstat )

     END IF

  END DO


  ! ----------------------------------------------------
  ! Initialize MetaData STRING fields special for OMHCHO
  ! ----------------------------------------------------
  SELECT CASE ( pge_idx )
  CASE (  pge_hcho_idx )
     DO imd = 1, n_mdata_omhcho
        version = 1
        mdata_omhcho_values(imd) = str_missval
        
        ! Short-hand for MetaData category and location
        mdata_typ = TRIM(ADJUSTL(mdata_omhcho_fields(2,imd)))
        mdata_loc = TRIM(ADJUSTL(mdata_omhcho_fields(3,imd)))

        ! ------------------------------------
        ! We only have Inventory Metadata here
        ! ------------------------------------
        mdata_index = md_inventory_idx

        SELECT CASE ( mdata_loc )
        CASE ( "BRO" )
           IF ( yn_bro_prefit(1) ) THEN
              md_stat = PGS_MET_getPCAttr_s( &
                   bro_prefit_lun, version, cmd_str, &
                   TRIM(ADJUSTL(upper_case(mdata_omhcho_fields(1,imd)))), &
                   mdata_omhcho_values(imd) )
              CALL error_check ( md_stat, PGS_S_SUCCESS, pge_errstat_warning, OMSAO_W_GETATTR, &
                   modulename//f_sep//TRIM(ADJUSTL(mdata_omhcho_fields(1,imd))), &
                   vb_lev_default, errstat )
           END IF
        CASE ( "OOO" )
           IF ( yn_o3_prefit(1) ) THEN
              md_stat = PGS_MET_getPCAttr_s( &
                   o3_prefit_lun, version, cmd_str, &
                   TRIM(ADJUSTL(upper_case(mdata_omhcho_fields(1,imd)))), &
                   mdata_omhcho_values(imd) )
              CALL error_check ( md_stat, PGS_S_SUCCESS, pge_errstat_warning, OMSAO_W_GETATTR, &
                   modulename//f_sep//TRIM(ADJUSTL(mdata_omhcho_fields(1,imd))), &
                   vb_lev_default, errstat )
           END IF
        CASE DEFAULT
           ! No idea what to do here
        END SELECT
     END DO

     DO imd = 1, n_mdata_voc
        version = 1
        mdata_voc_values(imd) = str_missval
        
        ! Short-hand for MetaData category and location
        mdata_typ = TRIM(ADJUSTL(mdata_voc_fields(2,imd)))
        mdata_loc = TRIM(ADJUSTL(mdata_voc_fields(3,imd)))

        ! ------------------------------------
        ! We only have Inventory Metadata here
        ! ------------------------------------
        mdata_index = md_inventory_idx

        SELECT CASE ( mdata_loc )
        CASE ( "CLD" )
           md_stat = PGS_MET_getPCAttr_s( &
                voc_amf_luns(voc_omicld_idx), version, TRIM(ADJUSTL(cmd_str))//'.0', &
                TRIM(ADJUSTL(upper_case(mdata_voc_fields(1,imd)))), &
                mdata_voc_values(imd) )
           CALL error_check ( md_stat, PGS_S_SUCCESS, pge_errstat_warning, OMSAO_W_GETATTR, &
                modulename//f_sep//TRIM(ADJUSTL(mdata_voc_fields(1,imd))), &
                vb_lev_default, errstat )
        CASE DEFAULT
           ! No idea what to do here
        END SELECT
     END DO

  ! ------------------------------------------------------
  ! Initialize MetaData STRING fields special for OMCHOCHO
  ! ------------------------------------------------------
  CASE ( pge_gly_idx)
     DO imd = 1, n_mdata_omchocho
        version = 1
        mdata_voc_values(imd) = str_missval
        
        ! Short-hand for MetaData category and location
        mdata_typ = TRIM(ADJUSTL(mdata_voc_fields(2,imd)))
        mdata_loc = TRIM(ADJUSTL(mdata_voc_fields(3,imd)))

        ! ------------------------------------
        ! We only have Inventory Metadata here
        ! ------------------------------------
        mdata_index = md_inventory_idx

        SELECT CASE ( mdata_loc )
        CASE ( "CLD" )
           md_stat = PGS_MET_getPCAttr_s( &
                voc_amf_luns(voc_omicld_idx), version, TRIM(ADJUSTL(cmd_str))//'.0', &
                TRIM(ADJUSTL(upper_case(mdata_voc_fields(1,imd)))), &
                mdata_voc_values(imd) )
           CALL error_check ( md_stat, PGS_S_SUCCESS, pge_errstat_warning, OMSAO_W_GETATTR, &
                modulename//f_sep//TRIM(ADJUSTL(mdata_voc_fields(1,imd))), &
                vb_lev_default, errstat )
        CASE DEFAULT
           ! No idea what to do here
        END SELECT
     END DO
  ! ----------------------------------------------------------------------------------------------
  ! Initialize MetaData STRING fields special for OMH2O same that for OMHCHO (not completely true)
  ! ----------------------------------------------------------------------------------------------
  CASE ( pge_h2o_idx)
     DO imd = 1, n_mdata_voc
        version = 1
        mdata_voc_values(imd) = str_missval
        
        ! Short-hand for MetaData category and location
        mdata_typ = TRIM(ADJUSTL(mdata_voc_fields(2,imd)))
        mdata_loc = TRIM(ADJUSTL(mdata_voc_fields(3,imd)))

        ! ------------------------------------
        ! We only have Inventory Metadata here
        ! ------------------------------------
        mdata_index = md_inventory_idx

        SELECT CASE ( mdata_loc )
        CASE ( "CLD" )
           md_stat = PGS_MET_getPCAttr_s( &
                voc_amf_luns(voc_omicld_idx), version, TRIM(ADJUSTL(cmd_str))//'.0', &
                TRIM(ADJUSTL(upper_case(mdata_voc_fields(1,imd)))), &
                mdata_voc_values(imd) )
           CALL error_check ( md_stat, PGS_S_SUCCESS, pge_errstat_warning, OMSAO_W_GETATTR, &
                modulename//f_sep//TRIM(ADJUSTL(mdata_voc_fields(1,imd))), &
                vb_lev_default, errstat )
        CASE DEFAULT
           ! No idea what to do here
        END SELECT
     END DO

  END SELECT

  ! ----------------------------------
  ! Initialize MetaData INTEGER fields
  ! ----------------------------------
  DO imd = 1, n_mdata_int
     IF ( TRIM(ADJUSTL(mdata_integer_fields(3,imd))) == "pcf" .OR. &
          TRIM(ADJUSTL(mdata_integer_fields(3,imd))) == "mcf" .OR. &
          TRIM(ADJUSTL(mdata_integer_fields(3,imd))) == "l1r"        ) THEN
  
        version = 1
        mdata_integer_values(imd) = int16_missval
  
        ! Short-hand for MetaData category and location
        mdata_typ = TRIM(ADJUSTL(mdata_integer_fields(2,imd)))
        mdata_loc = TRIM(ADJUSTL(mdata_integer_fields(3,imd)))
  
        SELECT CASE ( mdata_typ )
        CASE ( "inv" )
           mdata_index   = md_inventory_idx
           metadata_type = cmd_str//'.0'
        CASE ( "psa" )
           mdata_index   = md_inventory_idx
           metadata_type = cmd_str//'.0'
           ! Do we have any of these at all?
        CASE ( "arc" )
           mdata_index   = md_archive_idx
           metadata_type = amd_str
        CASE DEFAULT
           ! No idea what to do here
        END SELECT
  
        SELECT CASE ( mdata_loc )
        CASE ( "mcf" )
           md_stat = PGS_MET_getSetAttr_i( mcf_groups(mdata_index), &
                TRIM(ADJUSTL(upper_case(mdata_integer_fields(1,imd)))), &
                mdata_integer_values(imd) )
        CASE ( "pcf" )
           md_stat = PGS_MET_getSetAttr_i( mcf_groups(mdata_index), &
                TRIM(ADJUSTL(upper_case(mdata_integer_fields(1,imd)))), &
                mdata_integer_values(imd) )
           IF ( TRIM(ADJUSTL(mdata_integer_fields(1,imd))) == vid_str ) &
                mcf_versionid = mdata_integer_values(imd)
        CASE ( "l1r" )
           md_stat = PGS_MET_getPCAttr_i( &
                l1b_radiance_lun, version, TRIM(ADJUSTL(metadata_type)), &
                TRIM(ADJUSTL(upper_case(mdata_integer_fields(1,imd)))), &
                mdata_integer_values(imd) )
        CASE ( "l1R" )
           md_stat = PGS_MET_getPCAttr_i( &
                l1b_radiance_lun, version, TRIM(ADJUSTL(metadata_type)), &
                TRIM(ADJUSTL(upper_case(mdata_integer_fields(1,imd)))), &
                mdata_integer_values(imd) )
        CASE DEFAULT
           ! No idea what to do here
        END SELECT
  
        ! ------------------------------
        ! Error check for initialization
        ! ------------------------------
        CALL error_check ( &
             md_stat, PGS_S_SUCCESS, pge_errstat_error, OMSAO_E_GETATTR, &
             modulename//f_sep//TRIM(ADJUSTL(upper_case(mdata_integer_fields(1,imd)))), &
             vb_lev_default, errstat )

     END IF
  END DO


  ! ----------------------------------
  ! Initialize MetaData DOUBLE fields
  ! ----------------------------------
  DO imd = 1, n_mdata_dbl
     IF ( TRIM(ADJUSTL(mdata_double_fields(3,imd))) == "mcf" .OR. &
          TRIM(ADJUSTL(mdata_double_fields(3,imd))) == "l1r"        ) THEN
  
        version = 1
        mdata_double_values(imd) = r8_missval
  
        ! Short-hand for MetaData category and location
        mdata_typ = TRIM(ADJUSTL(mdata_double_fields(2,imd)))
        mdata_loc = TRIM(ADJUSTL(mdata_double_fields(3,imd)))
  
        SELECT CASE ( mdata_typ )
        CASE ( "inv" )
           mdata_index   = md_inventory_idx
           metadata_type = cmd_str//'.0'
        CASE ( "psa" )
           mdata_index   = md_inventory_idx
           metadata_type = cmd_str//'.0'
           ! Do we have any of these at all?
        CASE ( "arc" )
           mdata_index   = md_archive_idx
           metadata_type = amd_str//'.0'
        CASE DEFAULT
           ! No idea what to do here
        END SELECT
  
        SELECT CASE ( mdata_loc )
        CASE ( "mcf" )
           md_stat = PGS_MET_getSetAttr_d( mcf_groups(mdata_index), &
                TRIM(ADJUSTL(upper_case(mdata_double_fields(1,imd)))), &
                mdata_double_values(imd) )
        CASE ( "l1r" )
           md_stat = PGS_MET_getPCAttr_d( &
                l1b_radiance_lun, version, TRIM(ADjUSTL(metadata_type)), &
                TRIM(ADJUSTL(upper_case(mdata_double_fields(1,imd)))), &
                mdata_double_values(imd) )
        CASE ( "l1R" )
           md_stat = PGS_MET_getPCAttr_d( &
                l1b_radiance_lun, version, TRIM(ADjUSTL(metadata_type)), &
                TRIM(ADJUSTL(upper_case(mdata_double_fields(1,imd)))), &
                mdata_double_values(imd) )
        CASE DEFAULT
           ! No idea what to do here
        END SELECT
  
        ! ------------------------------
        ! Error check for initialization
        ! ------------------------------
        CALL error_check ( &
             md_stat, PGS_S_SUCCESS, pge_errstat_error, OMSAO_E_GETATTR, &
             modulename//f_sep//TRIM(ADJUSTL(upper_case(mdata_double_fields(1,imd)))), &
             vb_lev_default, errstat )
  
     END IF
  END DO
  
  ! -----------------------------
  ! Read L1B radiance OPF version
  ! -----------------------------
  ! (somehow this fails when called as subroutine)
  ! ----------------------------------------------
  md_stat = PGS_MET_getPCAttr_s ( &
       l1b_radiance_lun, version, TRIM(ADJUSTL(cmd_str//'.0')), 'INPUTPOINTER', l1br_inputp )
  l1br_opf_version = -1 ; version = 1
  get_opf: DO imd = 1, n_ipp
     tmp_string = l1br_inputp(imd)
     idx        = INDEX ( tmp_string, 'OML1BOPF' )
     IF ( idx > 0 ) THEN
        opf_string = tmp_string(idx+10:idx+14)
        READ ( opf_string, '(I4)') l1br_opf_version
     END IF
     IF ( l1br_opf_version > 0 ) EXIT get_opf
  END DO get_opf


  RETURN
END SUBROUTINE init_metadata

SUBROUTINE check_metadata_consistency ( errstat )

  ! ================================================================
  ! At this point we compare common Metadata values from PCF, MCF,
  ! and L1B files. Any inconstencies lead to error exits of the PGE.
  ! ================================================================

  USE OMSAO_precision_module, ONLY: i4
  USE OMSAO_variables_module, ONLY: orbit_number
  USE OMSAO_metadata_module
  USE OMSAO_errstat_module

  IMPLICIT NONE

  ! ------------------------------
  ! Name of this module/subroutine
  ! ------------------------------
  CHARACTER (LEN=26), PARAMETER :: modulename = "check_metadata_consistency"

  ! ---------------
  ! Output variable
  ! ---------------
  INTEGER (KIND=i4), INTENT (INOUT) :: errstat


  ! ---------------
  ! Local variables
  ! ---------------
  INTEGER (KIND=i4) :: i, k, locerrstat

  locerrstat = pge_errstat_ok

  ! -----------------------------------------------------------
  ! Checking for Orbit Number consistency. Sources: PCF and L1B.
  ! -------------------------------------------------------------------
  ! If the values from PCF and L1B are different, L1B takes precedence,
  ! and ORBIT_NUMBER, which at this point carries the value from the
  ! PCF, must be set to what was found in the L1B.
  ! -------------------------------------------------------------------
  getidx_1: DO i = 1, n_mdata_int
     IF ( mdata_integer_fields(3,i) == "l1r"              .AND. &
          mdata_integer_fields(2,i) == "inv"              .AND. &
          INDEX(mdata_integer_fields(1,i), orbn_str) /= 0 .AND. &
          orbit_number /= mdata_integer_values(i)           ) THEN
        orbit_number = mdata_integer_values(i)
        CALL error_check ( &
             1, 0, pge_errstat_warning, OMSAO_W_MDMISMATCH, &
             modulename//f_sep//orbn_str, vb_lev_default, errstat )
        EXIT getidx_1
     END IF
  END DO getidx_1

  ! -----------------------------------------------------------
  ! Checking for String consistency. Sources: MCF and L1B
  ! --------------------------------------------------------------
  ! (first find STRING string in L1B and MCF Metadata arrays,
  !  but exclude 'ShortName', which is the name of the PGE and
  !  thus different between L1B and L2)
  ! --------------------------------------------------------------
  DO k = 1, n_mdata_str
     IF ( mdata_string_fields(3,k) == "mcf" .AND. &
          mdata_string_fields(2,k) == "inv" .AND. &
          mdata_string_fields(1,k) /= "ShortName"  ) THEN
        getidx_2: DO i = 1, n_mdata_str
           IF ( mdata_string_fields(3,i) == "l1r" .AND. &
                mdata_string_fields(2,i) == "inv" .AND. &
                TRIM(ADJUSTL(mdata_string_fields(1,i))) == TRIM(ADJUSTL(mdata_string_fields(1,k))) .AND. &
                TRIM(ADJUSTL(mdata_string_values(k))) /= &
                TRIM(ADjUSTL(mdata_string_values(i))) ) THEN
              CALL error_check ( &
                   1, 0, pge_errstat_warning, OMSAO_W_MDMISMATCH, &
                   modulename//f_sep//TRIM(ADJUSTL(mdata_string_fields(1,k))), &
                   vb_lev_default, errstat )
           END IF
        END DO getidx_2
     END IF
  END DO

  ! -------------------------------------------------------------------------
  ! Checking for consistency of Granule Start/End times. Sources: L1B and PCF
  ! -------------------------------------------------------------------------
  DO k = 1, n_mdata_str
     IF ( mdata_string_fields(3,k) == "l1r" .AND. &
          mdata_string_fields(2,k) == "inv" ) THEN

        IF ( TRIM(ADJUSTL(mdata_string_fields(1,k))) == "RangeBeginningTime" .AND. &
             TRIM(ADJUSTL(mdata_string_values(k))) /= TRIM(ADJUSTL(pcf_granule_s_time)) ) THEN 
           CALL error_check ( &
                1, 0, pge_errstat_warning, OMSAO_W_MDMISMATCH, &
                modulename//f_sep//TRIM(ADJUSTL(mdata_string_fields(1,k))), &
                vb_lev_default, errstat )
        END IF

        IF ( TRIM(ADJUSTL(mdata_string_fields(1,k))) == "RangeEndingTime" .AND. &
             TRIM(ADJUSTL(mdata_string_values(k))) /= TRIM(ADJUSTL(pcf_granule_e_time)) ) THEN
           CALL error_check ( &
                1, 0, pge_errstat_warning, OMSAO_W_MDMISMATCH, &
                modulename//f_sep//TRIM(ADJUSTL(mdata_string_fields(1,k))), &
                vb_lev_default, errstat )
        END IF
     END IF
  END DO


  RETURN
END SUBROUTINE check_metadata_consistency


SUBROUTINE set_l2_metadata ( errstat )

  USE OMSAO_precision_module,  ONLY: i4
  USE OMSAO_indices_module,    ONLY: md_inventory_idx, md_archive_idx, &
       pge_l2_output_lun
  USE OMSAO_he5_module,        ONLY: n_lun_inp_max, n_lun_inp, lun_input
  USE OMSAO_metadata_module
  USE OMSAO_errstat_module
  USE OMSAO_parameters_module, ONLY: maxchlen
  USE OMSAO_variables_module,  ONLY: pge_idx, l2_filename

  IMPLICIT NONE

  ! ------------------------------
  ! Name of this module/subroutine
  ! ------------------------------
  CHARACTER (LEN=15), PARAMETER :: modulename = "set_l2_metadata"

  ! ---------------
  ! Output variable
  ! ---------------
  INTEGER (KIND=i4), INTENT (INOUT) :: errstat

  ! ---------------
  ! Local variables
  ! ---------------
  CHARACTER (LEN=PGS_SMF_MAX_MSG_SIZE)                                :: md_msg
  CHARACTER (LEN=PGSd_MET_MAX_STRING_SET_L), DIMENSION(n_lun_inp_max) :: input_pointer
  CHARACTER (LEN=maxchlen)                                            :: lgid
  INTEGER   (KIND=i4) :: ip, md_stat, version, locerrstat, ics, he5_met_id


  ! ------------------
  ! External functions
  ! ------------------
  INTEGER (KIND=i4),               EXTERNAL :: &
       OMI_localGranuleID, PGS_PC_GetReference, PGS_MET_setAttr_s,    &
       PGS_MET_setAttr_d, PGS_MET_setAttr_i, PGS_MET_setMultiAttr_s, PGS_MET_sfstart, &
       PGS_MET_SFend, PGS_MET_Remove, PGS_MET_write
  CHARACTER (LEN=maxchlen),        EXTERNAL :: compose_localgranuleid
  CHARACTER (LEN=PGSd_MET_NAME_L), EXTERNAL :: upper_case
  CHARACTER (LEN=PGSd_MET_MAX_STRING_SET_L), EXTERNAL :: int2string

  locerrstat = pge_errstat_ok

  ! ----------------
  ! Get InputPointer
  ! ----------------
  input_pointer(1:n_lun_inp_max) = PGSd_MET_STR_END
  DO ip = 1, n_lun_inp
     version = 1
     md_stat = PGS_PC_GetReference( lun_input(ip), version, md_msg )
     CALL error_check ( &
          md_stat, PGS_S_SUCCESS, pge_errstat_error, OMSAO_E_GETREF, &
          modulename//f_sep//"lun_input_pointer", vb_lev_default, errstat )
     ! ------------------------------------------------------------------
     ! Extract file name only from relative/absolute file path:
     ! find first occurrance of '/' from the end of the string backwards;
     ! anything forward of it is the file name.
     ! ------------------------------------------------------------------
     ics = INDEX ( md_msg, '/', BACK = .TRUE. ) + 1
     CALL slice_string (LEN_TRIM(ADJUSTL(md_msg)), TRIM(ADJUSTL(md_msg)), ics, input_pointer(ip) )
     !input_pointer(ip) = TRIM(ADJUSTL( md_msg(ics:) ))
  END DO

  ! --------------------------
  ! Compose the LocalGranuleID
  ! --------------------------
  lgid = compose_localgranuleid ( locerrstat )
  IF ( locerrstat /= pge_errstat_ok ) errstat = MAX ( errstat, locerrstat )

  ! -------------------------------------
  ! Now set all the L2 INVENTORY Metadata
  ! -------------------------------------

  ! --------------
  ! LocalGranuleID
  ! --------------
  md_stat = PGS_MET_setAttr_s ( &
       mcf_groups(md_inventory_idx), TRIM(ADJUSTL(lgid_str)), TRIM(ADJUSTL(lgid)) )
  CALL error_check ( &
       md_stat, PGS_S_SUCCESS, pge_errstat_error, OMSAO_E_MDL2INV, &
       modulename//f_sep//TRIM(ADJUSTL(lgid_str)), vb_lev_default, errstat )

  ! ------------------------------------------------------------
  ! Set various non-PGE generated Inventory and Archive MetaData
  ! ------------------------------------------------------------
  CALL set_str_int_dbl_metadata ( errstat )

  ! -------------
  ! Input Pointer
  ! -------------
  md_stat = PGS_MET_setMultiAttr_s ( &
       mcf_groups(md_inventory_idx), "InputPointer", n_lun_inp, &
       input_pointer(1:n_lun_inp) )
  CALL error_check ( &
       md_stat, PGS_S_SUCCESS, pge_errstat_error, OMSAO_E_MDL2INV, &
       modulename//f_sep//"InputPointer", vb_lev_default, errstat )

  ! --------------------------------------
  ! Set PGE Inventory and Archive MetaData
  ! --------------------------------------
  CALL set_pge_metadata ( errstat )

  ! -----------------------------------------
  ! Copy PSAs from L1B file to L2 output file
  ! -----------------------------------------
  CALL copy_psa_metadata ( errstat )

  ! ****************************************************************************
  ! Plenty of things are pending, including setting QA flags, Archived Metadata,
  ! and the copying of Metadata objects from L1B to L2 (according to K. Yang).
  ! ****************************************************************************

  ! ---------------------------------
  ! Write Metatdata to L2 output file
  ! ---------------------------------

  ! * Obtain L2 output filename
  !   (replicate what's being done in the READ_PCF_FILE subroutine)
  version = 1
  md_stat = PGS_PC_GetReference (pge_l2_output_lun, version, l2_filename)
  IF ( md_stat /= PGS_S_SUCCESS .OR. LEN(TRIM(ADJUSTL(l2_filename))) == 0 ) THEN
     CALL error_check ( &
          1, 0, pge_errstat_error, OMSAO_E_GETLUN, &
          modulename//f_sep//"L2 Output LUN", vb_lev_default, errstat )
  ENDIF

  ! ----------------------------------
  ! Initiate file for writing Metadata
  ! ----------------------------------
  md_stat = PGS_MET_sfstart ( TRIM(ADJUSTL(l2_filename)), HDF5_ACC_RDWR, he5_met_id )
  CALL error_check ( &
       md_stat, PGS_S_SUCCESS, pge_errstat_error, OMSAO_E_MDINIT, &
       modulename, vb_lev_default, errstat )

  ! -------------------
  ! Write CoreMetatData
  ! -------------------
  md_stat = PGS_MET_write ( mcf_groups(md_inventory_idx), cmd_str, he5_met_id )
  CALL error_check ( &
       md_stat, PGS_S_SUCCESS, pge_errstat_warning, OMSAO_W_CMDWRT, &
       modulename, vb_lev_default, errstat )

  ! -----------------------
  ! Write ArchivedMetatData
  ! -----------------------
  md_stat = PGS_MET_write ( mcf_groups(md_archive_idx), amd_str, he5_met_id )
  CALL error_check ( &
       md_stat, PGS_S_SUCCESS, pge_errstat_error, OMSAO_W_AMDWRT, &
       modulename, vb_lev_default, errstat )

  ! ---------------------------------------------------
  ! Close the L2 Metadata file and free Metadata memory
  ! ---------------------------------------------------
  md_stat = PGS_MET_SFend ( he5_met_id )
  md_stat = PGS_MET_Remove ( )

  RETURN
END SUBROUTINE set_l2_metadata


SUBROUTINE copy_psa_metadata ( errstat )

  USE OMSAO_precision_module,  ONLY: i4
  USE OMSAO_indices_module, ONLY: md_inventory_idx, &
       l1b_radiance_lun
  USE OMSAO_metadata_module
  USE OMSAO_errstat_module

  IMPLICIT NONE

  ! ---------------
  ! Output variable
  ! ---------------
  INTEGER (KIND=i4), INTENT (INOUT) :: errstat

  ! ---------------
  ! Local variables
  ! ---------------
  CHARACTER (LEN=5)                          :: cntr_str
  CHARACTER (LEN=PGSd_MET_NAME_L)            :: tmp_field
  CHARACTER (LEN=PGSd_MET_MAX_STRING_SET_L)  :: tmp_value

  INTEGER   (KIND=i4), PARAMETER                                 :: nstr_max = 100
  CHARACTER (LEN=PGSd_MET_MAX_STRING_SET_L), DIMENSION(nstr_max) :: str_array

  INTEGER   (KIND=i4) :: &
       imd, md_l1b_stat, md_stat, version, locerrstat, l1b_cntr, l2_cntr, nstr

  ! --------------------------------------
  ! Character strings for Attribute fields
  ! --------------------------------------
  CHARACTER (LEN=24), PARAMETER :: aan_str = "ADDITIONALATTRIBUTENAME."
  CHARACTER (LEN=15), PARAMETER :: pva_str = "PARAMETERVALUE."

  ! ------------------
  ! External functions
  ! ------------------
  CHARACTER (LEN=PGSd_MET_NAME_L),           EXTERNAL :: upper_case
  CHARACTER (LEN=PGSd_MET_MAX_STRING_SET_L), EXTERNAL :: int2string
  INTEGER (KIND=i4),                         EXTERNAL :: &
       PGS_MET_setAttr_s, PGS_MET_getPCAttr_s, PGS_MET_setMultiAttr_s


  locerrstat = pge_errstat_ok

  ! ----------------------------
  ! Copy PSA MetaData to L2 file
  ! -----------------------------------------------------------------------
  ! We read ADDITIONALATTRIBUTENAME.l1b_cntr from the L1B file, and write
  ! this as ADDITIONALATTRIBUTENAME.l2_cntr to the L2 file. We do this util
  ! l2_cntr equals the total number of PSAs to be copied, or until we run
  ! out of Attributes to read from the L1B file. The latter condition is
  ! checked via the md_l1b status variable
  ! -----------------------------------------------------------------------
  l1b_cntr = 1; l2_cntr  = 1; md_l1b_stat = pgs_s_success
  cp_psa: DO WHILE ( l2_cntr <= n_mdata_psa )
     IF ( l2_cntr > n_mdata_psa  ) EXIT cp_psa

     ! -------------------------------
     ! Blank counter and field strings
     ! -------------------------------
     cntr_str = ""; tmp_field = ""

     ! ---------------------------------------
     ! Compose attribute field name to be read
     ! ---------------------------------------
     cntr_str  = int2string ( l1b_cntr, 1 )
     tmp_field = aan_str//TRIM(ADJUSTL(cntr_str))

     ! ------------------------------------------------
     ! Read attribute from L1B file; exit loop if error
     ! ------------------------------------------------
     version = 1
     md_l1b_stat = PGS_MET_getPCAttr_s ( &
          l1b_radiance_lun, version, cmd_str//".0", &
          TRIM(ADJUSTL(upper_case(tmp_field))), tmp_value )

     IF ( md_l1b_stat /= pgs_s_success ) EXIT cp_psa

     write_psa: DO imd = 1, n_mdata_psa
        IF ( mdata_psa_fields(3,imd) == "l1r" .AND. &
             mdata_psa_fields(2,imd) == "psa" .AND. &
             TRIM(ADJUSTL(mdata_psa_fields(1,imd))) == TRIM(ADJUSTL(tmp_value)) ) THEN

           tmp_field = ""
           tmp_field = pva_str//TRIM(ADJUSTL(cntr_str))

           ! -----------------
           ! Blank value array
           ! -----------------
           str_array(:) = ""

           ! ---------------------
           ! Read attribute values
           ! ---------------------
           version = 1
           md_stat = PGS_MET_getPCAttr_s ( &
                l1b_radiance_lun, version, cmd_str//".0", &
                TRIM(ADJUSTL(upper_case(tmp_field))), str_array )
           IF ( md_stat /= pgs_s_success ) EXIT write_psa
           
           ! -----------------------------------------------------
           ! Find number of entries in str_array (non-Zero length)
           ! -----------------------------------------------------
           get_nstr: DO nstr = 1, nstr_max
              IF ( LEN(TRIM(ADJUSTL(str_array(nstr)))) <= 0 ) EXIT get_nstr
           END DO get_nstr

           ! -------------------------------------------------------------
           ! NSTR-1 is the last entry of non-Zero length. Adjust and check
           ! whether we have anything left to write. Exit loop if not.
           ! -------------------------------------------------------------
           nstr = nstr - 1;  IF ( nstr <= 0 ) EXIT write_psa

           ! -------------------
           ! Copy to L2 met file
           ! -------------------

           ! -------------------------------------
           ! Set L2 attribute name with L2 counter
           ! -------------------------------------
           cntr_str = "";  cntr_str  = int2string ( l2_cntr, 1 )
           tmp_field = ""; tmp_field = aan_str//TRIM(ADJUSTL(cntr_str))

           md_stat = PGS_MET_setAttr_s ( mcf_groups(md_inventory_idx), &
                TRIM(ADJUSTL(tmp_field)), TRIM(ADJUSTL(tmp_value)) )

           tmp_field = ""; tmp_field = pva_str//TRIM(ADJUSTL(cntr_str))

           md_stat = PGS_MET_setMultiAttr_s ( mcf_groups(md_inventory_idx), &
                TRIM(ADJUSTL(tmp_field)), nstr, str_array(1:nstr) )

           l2_cntr = l2_cntr + 1
           EXIT write_psa
        END IF
     END DO write_psa

     l1b_cntr = l1b_cntr + 1
  END DO cp_psa
  
  errstat = MAX ( locerrstat, errstat )

  RETURN
END SUBROUTINE copy_psa_metadata

SUBROUTINE set_str_int_dbl_metadata ( errstat )

  USE OMSAO_precision_module, ONLY: i4
  USE OMSAO_indices_module, ONLY: md_inventory_idx, md_archive_idx
  USE OMSAO_metadata_module
  USE OMSAO_errstat_module
  USE OMSAO_variables_module, ONLY: pge_idx
  IMPLICIT NONE

  ! ------------------------------
  ! Name of this module/subroutine
  ! ------------------------------
  CHARACTER (LEN=24), PARAMETER :: modulename = "set_str_int_dbl_metadata"

  ! ---------------
  ! Output variable
  ! ---------------
  INTEGER (KIND=i4), INTENT (INOUT) :: errstat

  ! ---------------
  ! Local variables
  ! ---------------
  CHARACTER (LEN=3)                         :: mdata_typ, mdata_loc
  INTEGER   (KIND=i4) :: imd, imd2, md_stat, locerrstat, mdata_index
  LOGICAL             :: write_field


  ! ------------------
  ! External functions
  ! ------------------
  CHARACTER (LEN=PGSd_MET_NAME_L), EXTERNAL :: upper_case
  INTEGER (KIND=i4),               EXTERNAL :: &
       PGS_MET_setAttr_s, PGS_MET_setAttr_d, PGS_MET_setAttr_i

  locerrstat = pge_errstat_ok

  ! ------------------------------------------------------------
  ! NOTE that ArchiveMetadata have to become global attributes
  ! and are written in HE5_WRITE_GLOBAL_ATTRIBUTES, i.e., BEFORE
  ! the HE5 file is closed.
  ! ------------------------------------------------------------

  ! -----------------------------
  ! Write STRING MetaData to file
  ! -----------------------------
  DO imd = 1, n_mdata_str
     IF ( TRIM(ADJUSTL(mdata_string_fields(3,imd))) == "l1r" .AND. &
          TRIM(ADJUSTL(mdata_string_fields(2,imd))) /= "psa"       ) THEN

        ! Short-hand for MetaData category and location
        mdata_typ = TRIM(ADJUSTL(mdata_string_fields(2,imd)))
        mdata_loc = TRIM(ADJUSTL(mdata_string_fields(3,imd)))

        write_field = .TRUE.

        SELECT CASE ( mdata_typ )
        CASE ( "inv" )
           mdata_index = md_inventory_idx
        CASE ( "arc" )
           write_field = .FALSE.
        CASE DEFAULT
           write_field = .FALSE.
        END SELECT

        ! -------------------------------------------------------------
        ! Check whether the current L1B inventory field has a match in
        ! the MCF fields or in the PGE. If so, skip writing this value.
        ! -------------------------------------------------------------
        fieldmatch: DO imd2 = 1, n_mdata_str
           IF ( ( &
                TRIM(ADJUSTL(mdata_string_fields(3,imd2))) == "mcf" .OR. &
                TRIM(ADJUSTL(mdata_string_fields(3,imd2))) == "pge"     ) .AND. &
                INDEX ( TRIM(ADJUSTL(mdata_string_fields(1,imd2))),       &
                TRIM(ADJUSTL(mdata_string_fields(1,imd))) ) /= 0 ) THEN
              write_field = .FALSE.
              EXIT fieldmatch
           END IF
        END DO fieldmatch

        IF ( write_field ) THEN
           IF ( mdata_index == md_archive_idx ) THEN
              !md_stat = HE5_EHwrglatt ( &
              !     pge_swath_file_id, TRIM(ADJUSTL(mdata_string_fields(1,imd))),  &
              !     HE5T_NATIVE_CHAR, LEN_TRIM(ADJUSTL(mdata_string_values(imd))), &
              !     TRIM(ADJUSTL(mdata_string_values(imd))) )
           ELSE
              md_stat = PGS_MET_setAttr_s ( mcf_groups(mdata_index), &
                   TRIM(ADJUSTL(mdata_string_fields(1,imd))), TRIM(ADJUSTL(mdata_string_values(imd))) )
           END IF
           CALL error_check ( &
                md_stat, PGS_S_SUCCESS, pge_errstat_warning, OMSAO_W_MDL2INV, &
                modulename//f_sep//TRIM(ADJUSTL(mdata_string_fields(1,imd))), &
                vb_lev_default, errstat )
        END IF

     END IF
  END DO

  ! ----------------------------------
  ! Write L1B INTEGER MetaData to file
  ! ----------------------------------
  DO imd = 1, n_mdata_int
     IF ( TRIM(ADJUSTL(mdata_integer_fields(3,imd))) == "l1r" ) THEN

        ! Short-hand for MetaData category and location
        mdata_typ = TRIM(ADJUSTL(mdata_integer_fields(2,imd)))
        mdata_loc = TRIM(ADJUSTL(mdata_integer_fields(3,imd)))

        write_field = .TRUE.

        SELECT CASE ( mdata_typ )
        CASE ( "inv" )
           mdata_index = md_inventory_idx
        CASE ( "arc" )
           write_field = .FALSE.
        CASE DEFAULT
           write_field = .FALSE.
        END SELECT

        IF ( write_field ) THEN
           md_stat = PGS_MET_setAttr_i ( mcf_groups(mdata_index), &
                TRIM(ADJUSTL(mdata_integer_fields(1,imd))), mdata_integer_values(imd) )
           CALL error_check ( &
                md_stat, PGS_S_SUCCESS, pge_errstat_warning, OMSAO_W_MDL2INV, &
                modulename//f_sep//TRIM(ADJUSTL(mdata_integer_fields(1,imd))), &
                vb_lev_default, errstat )
        END IF

     END IF
  END DO

  ! -------------------
  ! Set DOUBLE MetaData
  ! -------------------
  DO imd = 1, n_mdata_dbl
     IF ( TRIM(ADJUSTL(mdata_double_fields(3,imd))) == "l1r" ) THEN
        ! Short-hand for MetaData category and location
        mdata_typ = TRIM(ADJUSTL(mdata_double_fields(2,imd)))
        mdata_loc = TRIM(ADJUSTL(mdata_double_fields(3,imd)))
        write_field = .TRUE.

        SELECT CASE ( mdata_typ )
        CASE ( "inv" )
           mdata_index = md_inventory_idx
        CASE ( "arc" )
           write_field = .FALSE.
        CASE DEFAULT
           write_field = .FALSE.
        END SELECT

        IF ( write_field ) THEN
           md_stat = PGS_MET_setAttr_d ( mcf_groups(mdata_index), &
                TRIM(ADJUSTL(mdata_double_fields(1,imd))), mdata_double_values(imd) )
           CALL error_check ( &
                md_stat, PGS_S_SUCCESS, pge_errstat_warning, OMSAO_W_MDL2INV, &
                modulename//f_sep//TRIM(ADJUSTL(mdata_double_fields(1,imd))), &
                vb_lev_default, errstat )
        END IF
     END IF
  END DO

  RETURN
END SUBROUTINE set_str_int_dbl_metadata


SUBROUTINE set_pge_metadata ( errstat )

  USE OMSAO_indices_module, ONLY: md_inventory_idx
  USE OMSAO_metadata_module
  USE OMSAO_errstat_module
  USE OMSAO_variables_module, ONLY: pge_idx

  IMPLICIT NONE

  ! ------------------------------
  ! Name of this module/subroutine
  ! ------------------------------
  CHARACTER (LEN=16), PARAMETER :: modulename = "set_pge_metadata"

  ! ---------------
  ! Output variable
  ! ---------------
  INTEGER (KIND=i4), INTENT (INOUT) :: errstat

  ! ---------------
  ! Local variables
  ! ---------------
  CHARACTER (LEN=5)                         :: ipstring
  CHARACTER (LEN=PGSd_MET_MAX_STRING_SET_L) :: parname
  INTEGER   (KIND=i4)                       :: imd, md_stat, locerrstat

  ! ------------------
  ! External functions
  ! ------------------
  CHARACTER (LEN=PGSd_MET_MAX_STRING_SET_L), EXTERNAL :: int2string
  INTEGER (KIND=i4),                         EXTERNAL :: &
       PGS_MET_setAttr_s, PGS_MET_setAttr_i


  locerrstat = pge_errstat_ok

  ! ---------------------------------------------------------
  ! THE PGE output parameter, and associated inventory fields
  ! ---------------------------------------------------------
  DO imd = 1, n_l2_output_paras

     ipstring = int2string ( imd, 1 )

     parname = "ParameterName."//TRIM(ADJUSTL(ipstring)) 
     md_stat = PGS_MET_setAttr_s ( mcf_groups(md_inventory_idx), &
          TRIM(ADJUSTL(parname)), TRIM(ADJUSTL(pge_parameter_names(imd,pge_idx))) )
     CALL error_check ( &
          md_stat, PGS_S_SUCCESS, pge_errstat_warning, OMSAO_W_MDL2INV, &
          modulename//f_sep//TRIM(ADJUSTL(parname)), vb_lev_default, errstat )

     parname = "AutomaticQualityFlag."//TRIM(ADJUSTL(ipstring))
     md_stat = PGS_MET_setAttr_s ( mcf_groups(md_inventory_idx), &
          TRIM(ADJUSTL(parname)), TRIM(ADJUSTL(AutomaticQualityFlag)) )
     CALL error_check ( &
          md_stat, PGS_S_SUCCESS, pge_errstat_warning, OMSAO_W_MDL2INV, &
          modulename//f_sep//TRIM(ADJUSTL(parname)), vb_lev_default, errstat )

     parname = "AutomaticQualityFlagExplanation."//TRIM(ADJUSTL(ipstring))
     md_stat = PGS_MET_setAttr_s ( mcf_groups(md_inventory_idx), &
          TRIM(ADJUSTL(parname)), TRIM(ADJUSTL(AutomaticQualityFlagExplanation)) )     
     CALL error_check ( &
          md_stat, PGS_S_SUCCESS, pge_errstat_warning, OMSAO_W_MDL2INV, &
          modulename//f_sep//TRIM(ADJUSTL(parname)), vb_lev_default, errstat )

     parname = "QAPercentMissingData."//TRIM(ADJUSTL(ipstring))
     md_stat = PGS_MET_setAttr_i ( mcf_groups(md_inventory_idx), &
          TRIM(ADJUSTL(parname)), QAPercentMissingData )
     CALL error_check ( &
          md_stat, PGS_S_SUCCESS, pge_errstat_warning, OMSAO_W_MDL2INV, &
          modulename//f_sep//TRIM(ADJUSTL(parname)), vb_lev_default, errstat )

     parname = "QAPercentOutofBoundsData."//TRIM(ADJUSTL(ipstring))
     md_stat = PGS_MET_setAttr_i ( mcf_groups(md_inventory_idx), &
          TRIM(ADJUSTL(parname)), QAPercentOutofBoundsData )
     CALL error_check ( &
          md_stat, PGS_S_SUCCESS, pge_errstat_warning, OMSAO_W_MDL2INV, &
          modulename//f_sep//TRIM(ADJUSTL(parname)), vb_lev_default, errstat )

  END DO

  RETURN
END SUBROUTINE set_pge_metadata


CHARACTER (LEN=*) FUNCTION compose_localgranuleid ( errstat ) RESULT ( local_granule_id )

  ! ====================================
  ! Function to compose local granule id
  ! ====================================

  USE OMSAO_precision_module,  ONLY: i4
  USE OMSAO_parameters_module, ONLY: maxchlen
  USE OMSAO_metadata_module
  USE OMSAO_errstat_module
  USE OMSAO_variables_module,  ONLY: orbit_number

  IMPLICIT NONE

  ! ---------------------------------------
  ! Name of this module/subroutine/function
  ! ---------------------------------------
  CHARACTER (LEN=22), PARAMETER :: modulename = "compose_localgranuleid"

  ! ---------------
  ! Output variable
  ! ---------------
  INTEGER (KIND=i4), INTENT (INOUT) :: errstat

  ! ---------------
  ! Local variables
  ! ---------------
  INTEGER   (KIND=i4)        :: &
       year, month, day, julday, hour, minute, second, clen, &
       j1, j2, imd, locerrstat
  CHARACTER (LEN=maxchlen)   :: &
       year_str, month_str, day_str, hour_str, minute_str, second_str, &
       orbnum_str, mcf_vid_str, proc_dt_str, rbdt_str, ipsn_str, tmpstr1, tmpstr2

  ! ------------------
  ! External functions
  ! ------------------
  CHARACTER (LEN=maxchlen), EXTERNAL :: int2string

  ! --------------------------------------
  ! Compose Processing Date and Time (UTC)
  ! --------------------------------------

  locerrstat = pge_errstat_ok

  ! * get UTC date
  CALL utc_julian_date_and_time ( year, month, day, julday, hour, minute, second )

  ! * convert date numbers to STRING; do only those required for LocalGranuleID
  year_str   = int2string ( year,   4 )  ;  month_str  = int2string ( month,  2 )
  day_str    = int2string ( day,    2 )  ;  hour_str   = int2string ( hour,   2 )
  minute_str = int2string ( minute, 2 )  ;  second_str = int2string ( second, 2)

  ! * concatinate date strings 
  proc_dt_str = ""
  proc_dt_str = TRIM(ADJUSTL(year_str)) // &
       "m" // TRIM(ADJUSTL(month_str)) // TRIM(ADJUSTL(day_str)) // &
       "t" // TRIM(ADJUSTL(hour_str)) // TRIM(ADJUSTL(minute_str)) // &
       TRIM(ADJUSTL(second_str))

  ! ---------------------------------------------------------------------
  ! Construct "Range Beginning Date&Time", which is used in construction
  ! of LocalGranuleID. First find indices for Date and Time in Metadata
  ! field array, then concatinate. 
  ! ---------------------------------------------------------------------
  rbdt_str = ""  ;  j1 = 0  ;  j2 = 0

  get_idx1: DO imd = 1, n_mdata_str
     IF ( mdata_string_fields(3,imd) == "l1r" .AND. &
          mdata_string_fields(2,imd) == "inv"         ) THEN
        IF ( INDEX ( TRIM(ADJUSTL(mdata_string_fields(1,imd))), rbd_str ) /= 0 ) j1 = imd
        IF ( INDEX ( TRIM(ADJUSTL(mdata_string_fields(1,imd))), rbt_str ) /= 0 ) j2 = imd
        IF ( MIN(j1,j2) > 0 ) EXIT get_idx1
     END IF
  END DO get_idx1

  IF ( MIN(j1,j2) > 0 ) THEN     
     tmpstr1 = TRIM(ADJUSTL(mdata_string_values(j1)))
     tmpstr2 = TRIM(ADJUSTL(mdata_string_values(j2)))
     rbdt_str = tmpstr1(1:4) // "m" // tmpstr1(6:7) // tmpstr1(9:10) // &
          "t" // tmpstr2(1:2) // tmpstr2(4:5) // tmpstr2(7:8)
  ELSE
     CALL error_check ( &
          1, 0, pge_errstat_error, OMSAO_E_LGIDCOMP,  modulename, vb_lev_default, errstat )
  END IF

  ! ---------------------------------------------------------------------
  ! Construct "InstrumentPlatformName", which is used in construction
  ! of LocalGranuleID. It follows the same principle as the construction
  ! of "Range Beginning Date&Time".
  ! ---------------------------------------------------------------------
  ipsn_str = ""  ;  j1 = 0  ;  j2 = 0
  get_idx2: DO imd = 1, n_mdata_str
     IF ( mdata_string_fields(3,imd) == "l1r" .AND. &
          mdata_string_fields(2,imd) == "inv"         ) THEN
        IF ( INDEX ( TRIM(ADJUSTL(mdata_string_fields(1,imd))), aisn_str ) /= 0 ) j1 = imd
        IF ( INDEX ( TRIM(ADJUSTL(mdata_string_fields(1,imd))), apsn_str ) /= 0 ) j2 = imd
        IF ( MIN(j1,j2) > 0 ) EXIT get_idx2
     END IF
  END DO get_idx2
  IF ( MIN(j1,j2) > 0 ) THEN
     ipsn_str = TRIM(ADJUSTL(mdata_string_values(j1))) // "-" // &
          TRIM(ADJUSTL(mdata_string_values(j2)))
  ELSE
     CALL error_check ( &
          1, 0, pge_errstat_error, OMSAO_E_LGIDCOMP,  modulename, vb_lev_default, errstat )
  END IF

  ! ---------------------------------------------------
  ! Still missing are PCF OrbitNumber and MCF VersionID
  ! ------------------------------------------------------------------
  ! Note that, at this point, ORBIT_NUMBER is consistent with what was
  ! found in the L1B file. This is taken care of by the call to
  ! CHECK_METADATA_CONSITENCY earlier in the flow of the PGE.
  ! ------------------------------------------------------------------
  orbnum_str  = int2string ( orbit_number,  5 )
  mcf_vid_str = int2string ( mcf_versionid, 3 )

  ! -------------------------------------------------------------------
  ! Now bring it all together for LocalGranuleID, which is of the form
  !
  !  OMI-Aura_L2-OMBRO_2004m0131t0732-o01696_v001-2004m0425t124736.he5
  ! -------------------------------------------------------------------
  ! (remember that RBDT_STR contains seconds, which we don't want in
  !  the string for LocalGranuleID)
  ! -------------------------------------------------------------------  
  clen = LEN(TRIM(ADJUSTL(rbdt_str))) ;  rbdt_str = TRIM(ADJUSTL(rbdt_str))
  local_granule_id = &
       TRIM(ADJUSTL(ipsn_str)) // "_L2-" // TRIM(ADJUSTL(mcf_shortname)) // &
       "_" // rbdt_str(1:clen-2) // "-o" // TRIM(ADJUSTL(orbnum_str)) // &
       "_v" // TRIM(ADJUSTL(mcf_vid_str)) // "-" // TRIM(ADJUSTL(proc_dt_str)) // ".he5"

  RETURN
END FUNCTION compose_localgranuleid


SUBROUTINE set_automatic_quality_flag ( PercentGoodOutputSamples )

  USE OMSAO_precision_module
  USE OMSAO_metadata_module, ONLY: AutomaticQualityFlag, AutomaticQualityFlagExplanation
  USE OMSAO_omidata_module,  ONLY: qa_percent_passed, qa_percent_suspect

  IMPLICIT NONE

  ! ---------------
  ! Input Variables
  ! ---------------
  REAL (KIND=r4), INTENT (IN) :: PercentGoodOutputSamples

  ! ---------------
  ! Local Variables
  ! ---------------
  INTEGER   (KIND=i4) :: ipgos
  CHARACTER (LEN=3)   :: pstr, sstr

  ! ------------------
  ! External functions
  ! ------------------
  CHARACTER (LEN=2), EXTERNAL :: int2string


  pstr = int2string ( qa_percent_passed,  1 )
  sstr = int2string ( qa_percent_suspect, 1 )

  AutomaticQualityFlagExplanation = &
        "Flag set to Passed if NrofGoodOutputSamples >= "//TRIM(ADJUSTL(pstr))//          &
        "% of NrofGoodInputSamples, to Suspect if NrofGoodOutputSamples >= "//            &
        TRIM(ADJUSTL(sstr))//"% but < "//TRIM(ADJUSTL(pstr))//"% of NrofGoodInputSamples, "// &
        "and to Failed if NrofGoodOutputSamples < "//TRIM(ADJUSTL(sstr))//                &
        "% of NrofGoodInputSamples."

  ipgos = ANINT ( PercentGoodOutputSamples, KIND=i4 )

  IF ( ipgos >= qa_percent_passed ) THEN
     AutomaticQualityFlag = "Passed"
  ELSE IF ( ipgos < qa_percent_suspect ) THEN
     AutomaticQualityFlag = "Failed"
  ELSE
     AutomaticQualityFlag = "Suspect"
  END IF

  RETURN
END SUBROUTINE set_automatic_quality_flag

SUBROUTINE get_l1br_opf_version ( l1br_opf_version )

  USE OMSAO_precision_module, ONLY: i4
  USE OMSAO_indices_module, ONLY: l1b_radiance_lun
  USE OMSAO_metadata_module
  IMPLICIT NONE

  ! ---------------
  ! Output variable
  ! ---------------
  INTEGER (KIND=i4), INTENT (OUT) :: l1br_opf_version

  ! ------------------------------------------------------------
  ! Variables for the extraction of the L1B radiance OPF version
  ! ------------------------------------------------------------
  INTEGER   (KIND=i4), PARAMETER  :: n_ipp = 20
  INTEGER   (KIND=i4)             :: idx, md_stat, imd, version

  CHARACTER (LEN=PGSd_MET_NAME_L), DIMENSION (n_ipp) :: l1br_inputp
  CHARACTER (LEN=PGSd_MET_NAME_L)                    :: tmp_string
  CHARACTER (LEN=4)                                  :: opf_string

  INTEGER (KIND=i4), EXTERNAL :: PGS_MET_getPCAttr_s


  ! -----------------------------
  ! Read L1B radiance OPF version
  ! -----------------------------
  md_stat = PGS_MET_getPCAttr_s ( &
       l1b_radiance_lun, version, TRIM(ADJUSTL(cmd_str//'.0')), 'INPUTPOINTER', l1br_inputp )

  l1br_opf_version = -1 ; version = 1
  get_opf: DO imd = 1, n_ipp
     tmp_string = l1br_inputp(imd)

     idx        = INDEX ( tmp_string, 'OML1BOPF' )
     IF ( idx > 0 ) THEN
        opf_string = tmp_string(idx+10:idx+14)
        READ ( opf_string, '(I4)') l1br_opf_version
     END IF
     IF ( l1br_opf_version > 0 ) EXIT get_opf
  END DO get_opf

  RETURN
END SUBROUTINE get_l1br_opf_version
