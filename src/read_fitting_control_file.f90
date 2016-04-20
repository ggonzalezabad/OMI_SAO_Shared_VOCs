SUBROUTINE read_fitting_control_file ( pge_idx, l1b_radiance_esdt, pge_error_status )

  ! ***********************************************************
  !
  !   Read fitting control parameters from input control file
  !
  ! ***********************************************************

  USE OMSAO_precision_module
  USE OMSAO_indices_module, ONLY: &
       max_rs_idx, calfit_strings, max_calfit_idx, radfit_strings, mxs_idx, n_max_fitpars, &
       hwe_idx, refspec_strings,    icf_idx, pge_static_input_luns,                        &
       genline_str, socline_str, racline_str, rrsline_str, procline_str,                   &
       rafline_str, rspline_str, iofline_str, molline_str, eoi3str, us1_idx, us2_idx,      &
       pge_oclo_idx,  pge_hcho_idx, pge_no2_idx, pge_gly_idx, pge_no2d_idx,                &
       solcal_idx, radcal_idx, radref_idx, radfit_idx, wavwindow_str, fitresconst_str,     &
       destriping_str, scpline_str, nrmline_str, comline_str, o3amf_str, maxgoodcol_str,   &
       comm_idx, procmode_diag, procmode_prod, solmonthave_str, wfmod_amf_str,             &
       newshift_str, refseccor_str, scattweight_str
  USE OMSAO_parameters_module,   ONLY: maxchlen, forever, n_fit_winwav, max_mol_fit
  USE OMSAO_variables_module,    ONLY: &
       fitcol_idx, n_mol_fit, max_itnum_sol, max_itnum_rad, yn_smooth, yn_doas,            &
       fitvar_sol_init, fitvar_rad_init,  fitvar_rad_saved, szamax,                        &
       zatmos, lo_sunbnd, up_sunbnd, lo_radbnd, up_radbnd, n_refspec,                      &
       yn_use_labslitfunc, radwavcal_freq, tol, epsrel,  epsabs,  epsx, pm_one, phase,     &
       fit_winwav_lim, fit_winexc_lim, pixnum_lim, radfit_latrange,                        &
       static_input_fnames, fitvar_rad_str, verb_thresh_lev, winwav_min, winwav_max,       &
       have_undersampling, l1b_rad_filename,                                               &
       n_fitres_loop, fitres_range, l1b_channel, yn_solar_i0,                              &
       yn_solar_comp, solar_comp_typ, solar_comp_orb, yn_spectrum_norm, yn_common_iter,    &
       common_fitpos, common_fitvar, common_latrange, yn_o3amf_cor, yn_diagnostic_run,     &
       max_good_col, yn_solmonthave, yn_newshift, yn_refseccor, yn_sw
  USE OMSAO_prefitcol_module, ONLY: yn_bro_prefit, yn_o3_prefit, yn_lqh2o_prefit
  USE OMSAO_omidata_module, ONLY: vis_lower_wvl, omi_ccdpix_selection, nxtrack_max
  USE OMSAO_radiance_ref_module, ONLY: &
       yn_radiance_reference, radref_latrange, yn_remove_target, target_npol
  USE OMSAO_destriping_module, ONLY: &
       ctr_pol_base, ctr_pol_scal, ctr_pol_patt, ctr_nloop, ctrdst_latrange, ctr_nblocks, &
       ctr_fitfunc_calls, ctr_maxcol, yn_remove_ctrbias, ctr_bias_pol, yn_run_destriping
  USE OMSAO_casestring_module, ONLY: lower_case
  USE OMSAO_errstat_module
  USE OMSAO_wfamf_module, ONLY: &
       yn_amf_wfmod, amf_wfmod_idx, amf_alb_lnd, amf_alb_sno, amf_wvl, amf_wvl2, &
       amf_alb_cld, amf_max_sza

  IMPLICIT NONE

  ! ---------------
  ! Input variables
  ! ---------------
  INTEGER (KIND=i4), INTENT (IN) :: pge_idx

  ! ---------------
  ! Output variable
  ! ---------------
  CHARACTER (LEN=maxchlen), INTENT (OUT)   :: l1b_radiance_esdt
  INTEGER   (KIND=i4),      INTENT (INOUT) :: pge_error_status

  ! ---------------
  ! Local variables
  ! ---------------
  INTEGER   (KIND=i4)      :: i, j, k, file_read_stat, sidx, ridx
  INTEGER   (KIND=i4)      :: fit_ctrl_unit
  CHARACTER (LEN=maxchlen) :: tmpchar
  CHARACTER (LEN=3)        :: idxchar
  LOGICAL                  :: yn_eoi
  REAL      (KIND=r8)      :: vartmp, lotmp, uptmp

  ! ==============================
  ! Name of this module/subroutine
  ! ==============================
  CHARACTER (LEN=25), PARAMETER :: modulename = 'read_fitting_control_file'

  ! ========================
  ! Error handling variables
  ! ========================
  INTEGER (KIND=i4) :: errstat, version

  ! =================================
  ! External OMI and Toolkit routines
  ! =================================
  INTEGER (KIND=i4) :: pgs_smf_teststatuslevel, pgs_io_gen_openf, pgs_io_gen_closef

  errstat = pge_errstat_ok

  ! ---------------------------
  ! Initialize output variables
  ! ---------------------------
  ctr_pol_base    = 0      ; ctr_pol_scal = 0 ; ctr_pol_patt      = 0
  ctr_nloop       = 0      ; ctr_nloop    = 0 ; ctr_fitfunc_calls = 0
  ctrdst_latrange = 0.0_r4 ; ctr_maxcol   = 0.0_r8
  ctr_bias_pol    = 0      ; yn_remove_ctrbias = .FALSE.

  ! -------------------------
  ! Open fitting control file
  ! -------------------------
  version = 1
  errstat = PGS_IO_GEN_OPENF ( &
       pge_static_input_luns(icf_idx), PGSd_IO_Gen_RSeqFrm, 0, &
       fit_ctrl_unit, version )
  errstat = PGS_SMF_TESTSTATUSLEVEL(errstat)
  CALL error_check ( &
       errstat, PGS_SMF_MASK_LEV_S, pge_errstat_fatal, OMSAO_F_OPEN_FITCTRL_FILE, &
       modulename//f_sep//TRIM(ADJUSTL(static_input_fnames(icf_idx))), &
       vb_lev_default, pge_error_status )
  IF ( pge_error_status >= pge_errstat_error ) RETURN

  ! -----------------------------------------------
  ! Position cursor to read molecule name(s) to fit
  ! -----------------------------------------------
  REWIND ( fit_ctrl_unit )
  CALL skip_to_filemark ( fit_ctrl_unit, molline_str, tmpchar, file_read_stat )
  CALL error_check ( &
       file_read_stat, file_read_ok, pge_errstat_fatal, OMSAO_F_READ_FITCTRL_FILE, &
       modulename//f_sep//molline_str, vb_lev_default, pge_error_status )
  IF ( pge_error_status >= pge_errstat_error ) RETURN

  READ (fit_ctrl_unit, '(A)') tmpchar
  CALL get_mols_for_fitting ( tmpchar, n_mol_fit, fitcol_idx, errstat )
  CALL error_check ( &
       errstat, pge_errstat_ok, pge_errstat_fatal, OMSAO_F_GET_MOLFITNAME, &
       modulename, vb_lev_default, pge_error_status )
  IF ( pge_error_status >= pge_errstat_error ) RETURN

  ! --------------------------------------------------------------
  ! Position cursor to read processing mode. Set YN_DIAGNOSTIC_RUN
  ! to .TRUE. if "diagnostic" is selected as processing mode.
  ! --------------------------------------------------------------
  REWIND ( fit_ctrl_unit )
  CALL skip_to_filemark ( fit_ctrl_unit, procline_str, tmpchar, file_read_stat )
  CALL error_check ( &
       file_read_stat, file_read_ok, pge_errstat_fatal, OMSAO_F_READ_FITCTRL_FILE, &
       modulename//f_sep//procline_str, vb_lev_default, pge_error_status )
  IF ( pge_error_status >= pge_errstat_error ) RETURN

  READ (fit_ctrl_unit, '(A)') tmpchar
  tmpchar = lower_case ( TRIM(ADJUSTL(tmpchar)) )
  IF ( TRIM(ADJUSTL(tmpchar)) == procmode_diag ) THEN
     yn_diagnostic_run = .TRUE.
  ELSE
     yn_diagnostic_run = .FALSE.
  END IF

  ! ------------------------------------------------
  ! Position cursor to read general input parameters
  ! ------------------------------------------------
  REWIND (fit_ctrl_unit)
  CALL skip_to_filemark ( fit_ctrl_unit, genline_str, tmpchar, file_read_stat )
  CALL error_check ( &
       file_read_stat, file_read_ok, pge_errstat_fatal, OMSAO_F_READ_FITCTRL_FILE, &
       modulename//f_sep//genline_str, vb_lev_default, pge_error_status )
  IF ( pge_error_status >= pge_errstat_error ) RETURN

  READ (fit_ctrl_unit, *) yn_smooth
  READ (fit_ctrl_unit, *) yn_doas
  READ (fit_ctrl_unit, *) pixnum_lim(1:2)
  READ (fit_ctrl_unit, *) pixnum_lim(3:4)
  READ (fit_ctrl_unit, *) radfit_latrange(1:2)
  READ (fit_ctrl_unit, *) tol
  READ (fit_ctrl_unit, *) epsrel
  READ (fit_ctrl_unit, *) epsabs
  READ (fit_ctrl_unit, *) epsx

  ! ------------------------------------------------
  ! Check for consistency of pixel limits to process
  ! ------------------------------------------------
  IF ( ANY ( pixnum_lim(3:4) < 1 ) ) pixnum_lim(3:4) = (/ 1, nxtrack_max /)
  IF ( pixnum_lim(1) > pixnum_lim(2) ) pixnum_lim(1) = pixnum_lim(2)
  IF ( pixnum_lim(3) > pixnum_lim(4) ) pixnum_lim(3) = pixnum_lim(4)

  ! -------------------------------------------------
  ! Position cursor to read solar composite selection
  ! -------------------------------------------------
  REWIND (fit_ctrl_unit)
  CALL skip_to_filemark ( fit_ctrl_unit, scpline_str, tmpchar, file_read_stat )
  CALL error_check ( &
       file_read_stat, file_read_ok, pge_errstat_fatal, OMSAO_F_READ_FITCTRL_FILE, &
       modulename//f_sep//scpline_str, vb_lev_default, pge_error_status )
  IF ( pge_error_status >= pge_errstat_error ) RETURN

  READ (fit_ctrl_unit, *) yn_solar_comp
  IF ( yn_solar_comp ) READ (fit_ctrl_unit,*) solar_comp_typ

  ! -----------------------------------------------------
  ! Position cursor to read solar monthly average section
  ! -----------------------------------------------------
  REWIND (fit_ctrl_unit)
  CALL skip_to_filemark ( fit_ctrl_unit, solmonthave_str, tmpchar, file_read_stat )
  CALL error_check ( &
       file_read_stat, file_read_ok, pge_errstat_fatal, OMSAO_F_READ_FITCTRL_FILE, &
       modulename//f_sep//scpline_str, vb_lev_default, pge_error_status )
  IF ( pge_error_status >= pge_errstat_error ) RETURN

  READ (fit_ctrl_unit, *) yn_solmonthave

  ! -------------------------------------------------------
  ! Position cursor to read spectum normalization selection
  ! -------------------------------------------------------
  REWIND (fit_ctrl_unit)
  CALL skip_to_filemark ( fit_ctrl_unit, nrmline_str, tmpchar, file_read_stat )
  CALL error_check ( &
       file_read_stat, file_read_ok, pge_errstat_fatal, OMSAO_F_READ_FITCTRL_FILE, &
       modulename//f_sep//nrmline_str, vb_lev_default, pge_error_status )
  IF ( pge_error_status >= pge_errstat_error ) RETURN

  READ (fit_ctrl_unit, *) yn_spectrum_norm


  ! ---------------------------------------------------------------------
  ! Position cursor to read common mode iteration. .TRUE. will lead to a
  ! second run through the orbit, using the common mode spectrum that has
  ! been created during the first pass.
  ! ---------------------------------------------------------------------
  REWIND (fit_ctrl_unit)
  CALL skip_to_filemark ( fit_ctrl_unit, comline_str, tmpchar, file_read_stat )
  CALL error_check ( &
       file_read_stat, file_read_ok, pge_errstat_fatal, OMSAO_F_READ_FITCTRL_FILE, &
       modulename//f_sep//comline_str, vb_lev_default, pge_error_status )
  IF ( pge_error_status >= pge_errstat_error ) RETURN

  READ (fit_ctrl_unit, *) yn_common_iter
  IF ( yn_common_iter ) READ (fit_ctrl_unit, *) common_latrange


  ! ---------------------------------------------------------------------
  ! Position cursor to read radiance reference settings.
  ! ---------------------------------------------------------------------
  REWIND (fit_ctrl_unit)
  CALL skip_to_filemark ( fit_ctrl_unit, rrsline_str, tmpchar, file_read_stat )
  CALL error_check ( &
       file_read_stat, file_read_ok, pge_errstat_fatal, OMSAO_F_READ_FITCTRL_FILE, &
       modulename//f_sep//rrsline_str, vb_lev_default, pge_error_status )
  IF ( pge_error_status >= pge_errstat_error ) RETURN

  READ (fit_ctrl_unit, *) yn_radiance_reference
  READ (fit_ctrl_unit, *) yn_remove_target, target_npol
  READ (fit_ctrl_unit, *) radref_latrange(1:2)

  ! ----------------------------------------------------------
  ! Position cursor to read solar calibration input parameters
  ! ----------------------------------------------------------
  REWIND (fit_ctrl_unit)
  CALL skip_to_filemark ( fit_ctrl_unit, socline_str, tmpchar, file_read_stat )
  CALL error_check ( &
       file_read_stat, file_read_ok, pge_errstat_fatal, OMSAO_F_READ_FITCTRL_FILE, &
       modulename//f_sep//socline_str, vb_lev_default, pge_error_status )
  IF ( pge_error_status >= pge_errstat_error ) RETURN

  READ (fit_ctrl_unit, *) max_itnum_sol

  yn_use_labslitfunc = .TRUE.
  solpars: DO i = 1, max_calfit_idx

     READ (fit_ctrl_unit, *) idxchar, vartmp, lotmp, uptmp
     ! ---------------------------------------------------------
     ! Check for consitency of bounds and adjust where necessary
     ! ---------------------------------------------------------
     IF ( lotmp > vartmp .OR. uptmp < vartmp ) THEN
        lotmp = vartmp ; uptmp = vartmp
     END IF
     IF ( lotmp == uptmp .AND. lotmp /= vartmp ) THEN
        uptmp = vartmp ; lotmp = vartmp
     END IF

     IF ( idxchar == eoi3str ) EXIT solpars

     CALL string2index ( calfit_strings, max_calfit_idx, idxchar, sidx )
     IF ( sidx > 0 ) THEN
        fitvar_sol_init(sidx) = vartmp
        lo_sunbnd(sidx) = lotmp ; up_sunbnd(sidx) = uptmp
        ! ----------------------------------------------------
        ! Check whether we will be doing slit function fitting
        ! ----------------------------------------------------
        IF ( sidx == hwe_idx .AND. fitvar_sol_init(sidx) > 0.0_r8 ) yn_use_labslitfunc = .FALSE.
     END IF
  END DO solpars

  ! -------------------------------------------------------------
  ! Position cursor to read radiance calibration input parameters
  ! -------------------------------------------------------------
  REWIND (fit_ctrl_unit)
  CALL skip_to_filemark ( fit_ctrl_unit, racline_str, tmpchar, file_read_stat )
  CALL error_check ( &
       file_read_stat, file_read_ok, pge_errstat_fatal, OMSAO_F_READ_FITCTRL_FILE, &
       modulename//f_sep//racline_str, vb_lev_default, pge_error_status )
  IF ( pge_error_status >= pge_errstat_error ) RETURN

  ! --------------------------------------------------------------------
  ! Special logicals for HCHO - O3 and BrO prefits. Two-dimensional:
  ! (1) Use O3/BrO prefitted columns?
  ! (2) Vary those prefitted columns within their fitting uncertainties?
  !
  ! Special logicals for CHOCHO - Liquid Water prefits. Two-dimensional:
  ! (1) Use Liquid Water prefitted columns?
  ! (2) Vary those prefitted columns within their fitting uncertainties?
  ! --------------------------------------------------------------------
  SELECT CASE ( pge_idx )
  CASE ( pge_hcho_idx )
     READ (fit_ctrl_unit, *) yn_o3_prefit (1:2)
     READ (fit_ctrl_unit, *) yn_bro_prefit(1:2)
  CASE ( pge_gly_idx )
     READ (fit_ctrl_unit, *) yn_lqh2o_prefit (1:2)
  CASE DEFAULT
  END SELECT

  READ (fit_ctrl_unit, *) yn_solar_i0
  READ (fit_ctrl_unit, *) max_itnum_rad
  READ (fit_ctrl_unit, *) radwavcal_freq
  READ (fit_ctrl_unit, *) szamax
  READ (fit_ctrl_unit, *) zatmos
  READ (fit_ctrl_unit, *) phase
  fitvar_rad_init = 0.0_r8
  radpars: DO i = 1, max_calfit_idx
     READ (fit_ctrl_unit, *) idxchar, vartmp, lotmp, uptmp
     ! ---------------------------------------------------------
     ! Check for consitency of bounds and adjust where necessary
     ! ---------------------------------------------------------
     IF ( lotmp > vartmp .OR. uptmp < vartmp ) THEN
        lotmp = vartmp ; uptmp = vartmp
     END IF
     IF ( lotmp == uptmp .AND. lotmp /= vartmp ) THEN
        uptmp = vartmp ; lotmp = vartmp
     END IF

     IF ( idxchar == eoi3str ) EXIT radpars
     CALL string2index ( calfit_strings, max_calfit_idx, idxchar, sidx )
     IF ( sidx > 0 ) THEN
        fitvar_rad_init(sidx) = vartmp
        fitvar_rad_str (sidx) = TRIM(ADJUSTL(idxchar))
        lo_radbnd(sidx) = lotmp ; up_radbnd(sidx) = uptmp
     END IF
  END DO radpars
  fitvar_rad_saved = fitvar_rad_init

  ! ---------------------------------------------------------------------
  ! Check the latitude for the radiance reference, a.k.a. wavelength
  ! calibration spectrum. 
  ! If larger than 90 deg, set to 0.0, i.e., the Equator.
  ! ---------------------------------------------------------------------
  WHERE ( ABS(radref_latrange) > 90.0_r4 )
     radref_latrange = 0.0_r4
  END WHERE

  ! ---------------------------------------------------------
  ! Position cursor to read WFmodified AMF logical
  ! ---------------------------------------------------------
  REWIND (fit_ctrl_unit)
  CALL skip_to_filemark ( fit_ctrl_unit, wfmod_amf_str, tmpchar, file_read_stat )
  CALL error_check ( &
       file_read_stat, file_read_ok, pge_errstat_fatal, OMSAO_F_READ_FITCTRL_FILE, &
       modulename//f_sep//wfmod_amf_str, vb_lev_default, pge_error_status )
  IF ( pge_error_status >= pge_errstat_error ) RETURN
  READ (fit_ctrl_unit, *) yn_amf_wfmod, amf_wfmod_idx
  READ (fit_ctrl_unit, *) amf_alb_lnd, amf_alb_sno, amf_alb_cld
  IF ( .NOT. yn_amf_wfmod ) READ (fit_ctrl_unit, *) amf_wvl, amf_wvl2, amf_max_sza

  ! ---------------------------------------------------------
  ! Position cursor to read O3 AMF correction logical
  ! ---------------------------------------------------------
  REWIND (fit_ctrl_unit)
  CALL skip_to_filemark ( fit_ctrl_unit, o3amf_str, tmpchar, file_read_stat )
  CALL error_check ( &
       file_read_stat, file_read_ok, pge_errstat_fatal, OMSAO_F_READ_FITCTRL_FILE, &
       modulename//f_sep//o3amf_str, vb_lev_default, pge_error_status )
  IF ( pge_error_status >= pge_errstat_error ) RETURN
  READ (fit_ctrl_unit, *) yn_o3amf_cor



  ! ---------------------------------------------------------
  ! Position cursor to read radiance fitting input parameters
  ! ---------------------------------------------------------
  REWIND (fit_ctrl_unit)
  CALL skip_to_filemark ( fit_ctrl_unit, rafline_str, tmpchar, file_read_stat )
  CALL error_check ( &
       file_read_stat, file_read_ok, pge_errstat_fatal, OMSAO_F_READ_FITCTRL_FILE, &
       modulename//f_sep//rafline_str, vb_lev_default, pge_error_status )
  IF ( pge_error_status >= pge_errstat_error ) RETURN


  ! ------------------------------------------------------------------------
  ! By default we set the undersampling spectrum to FALSE. Only if we select
  ! it to be included in the fitting does it become trues. This way we save
  ! computation time for cases where we don't include the undersampling.
  ! ------------------------------------------------------------------------
  have_undersampling = .FALSE.
  getpars: DO j = 1, max_rs_idx
     READ (UNIT=fit_ctrl_unit, FMT='(A)', IOSTAT=errstat) tmpchar
     CALL error_check ( &
          errstat, file_read_ok, pge_errstat_fatal, OMSAO_F_READ_FITCTRL_FILE, &
          modulename//f_sep//"radiance fitting parameters", vb_lev_default, pge_error_status )
     IF ( pge_error_status >= pge_errstat_error ) RETURN

     CALL check_for_endofinput ( TRIM(ADJUSTL(tmpchar)), yn_eoi )
     IF ( yn_eoi ) EXIT getpars
     CALL string2index ( refspec_strings, max_rs_idx, tmpchar, ridx )

     ! ------------------------------------------------------------------------
     ! The loop goes over MXS_IDX+1 to catch any "eoi" that may have been added
     ! after a three-element block initialization (O3 for polynomial-dependent
     ! cross-section correction).
     ! ------------------------------------------------------------------------
     getblock: DO k = 1, mxs_idx+1

        READ (fit_ctrl_unit, *) idxchar, vartmp, lotmp, uptmp
        IF ( TRIM(ADJUSTL(idxchar)) == 'eoi' ) EXIT getblock

        ! ---------------------------------------------------------
        ! Check for consitency of bounds and adjust where necessary
        ! ---------------------------------------------------------
        IF ( lotmp > vartmp .OR. uptmp < vartmp ) THEN
           lotmp = vartmp ; uptmp = vartmp
        END IF
        IF ( lotmp == uptmp .AND. lotmp /= vartmp ) THEN
           uptmp = vartmp ; lotmp = vartmp
        END IF

        CALL string2index ( radfit_strings, mxs_idx, idxchar, sidx )

        IF ( sidx > 0 .AND. ridx > 0 ) THEN
           i = max_calfit_idx + (ridx-1)*mxs_idx + sidx
           fitvar_rad_init(i) = vartmp
           fitvar_rad_str (i) = TRIM(ADJUSTL(tmpchar))
           lo_radbnd (i) = lotmp ; up_radbnd (i) = uptmp
           IF ( (ridx == us1_idx .OR. ridx == us2_idx) .AND. &
                ANY ( (/ vartmp,lotmp,uptmp /) /= 0.0_r8 ) ) have_undersampling(ridx) = .TRUE.

           ! --------------------------------------------------
           ! Check for on-line Common Mode spectrum computation
           ! --------------------------------------------------
           IF ( ridx == comm_idx .AND. yn_common_iter ) THEN
              common_fitvar(1:3) = (/ vartmp, lotmp, uptmp /)
              common_fitpos      = sidx
           END IF
        END IF
     END DO getblock
  END DO getpars

  ! ------------------------------------------------------------------
  ! Formerly wavelength array specifications, now superceeded by the
  ! CCD pixel specifications. Included here for backward compatibility
  ! and easy of human reading (how else would we know the extend of
  ! the CCD pixel slice in terms of wavelengths?)
  ! ------------------------------------------------------------------
  REWIND (fit_ctrl_unit)
  CALL skip_to_filemark ( fit_ctrl_unit, wavwindow_str, tmpchar, file_read_stat )
  CALL error_check ( &
       file_read_stat, file_read_ok, pge_errstat_fatal, OMSAO_F_READ_FITCTRL_FILE, &
       modulename//f_sep//wavwindow_str, vb_lev_default, pge_error_status )
  IF ( pge_error_status >= pge_errstat_error ) RETURN
  READ (fit_ctrl_unit, *) l1b_channel
  READ (fit_ctrl_unit, *) fit_winwav_lim(1:n_fit_winwav), fit_winexc_lim(1:2)

  ! ------------------------------------------------------------------
  ! Acceptable window for the fitting residual in multiples of its
  ! standard deviation, and the number of iterations  we will perform 
  ! to get the residual within these bounds. Negative values stand for
  ! "no constraints" - this is mainly for the radiance wavelength
  ! calibration, which consists of the fitting of the solar reference to a
  ! radiance spectrum and hence produces a large residual anyway.  Smaller
  ! windows and more iterations lead to larger execution times.
  ! ------------------------------------------------------------------
  REWIND (fit_ctrl_unit)
  CALL skip_to_filemark ( fit_ctrl_unit, fitresconst_str, tmpchar, file_read_stat )
  CALL error_check ( &
       file_read_stat, file_read_ok, pge_errstat_fatal, OMSAO_F_READ_FITCTRL_FILE, &
       modulename//f_sep//fitresconst_str, vb_lev_default, pge_error_status )
  IF ( pge_error_status >= pge_errstat_error ) RETURN
  READ (fit_ctrl_unit, *) fitres_range(solcal_idx), n_fitres_loop(solcal_idx)
  READ (fit_ctrl_unit, *) fitres_range(radcal_idx), n_fitres_loop(radcal_idx)
  READ (fit_ctrl_unit, *) fitres_range(radref_idx), n_fitres_loop(radref_idx)
  READ (fit_ctrl_unit, *) fitres_range(radfit_idx), n_fitres_loop(radfit_idx)


  ! ---------------------------------------------------------------------------
  ! Position cursor to read maximum good column amount
  ! ---------------------------------------------------------------------------
  REWIND (fit_ctrl_unit)
  CALL skip_to_filemark ( fit_ctrl_unit, maxgoodcol_str, tmpchar, file_read_stat )
  CALL error_check ( &
       file_read_stat, file_read_ok, pge_errstat_fatal, OMSAO_F_READ_FITCTRL_FILE, &
       modulename//f_sep//maxgoodcol_str, vb_lev_default, pge_error_status )
  IF ( pge_error_status >= pge_errstat_error ) RETURN
  READ (fit_ctrl_unit, *) max_good_col

  ! ---------------------------------------------------------------------------
  ! Position cursor to read destriping parameters:
  ! * Order of cross track polynomials (Baseline, Scaling; XTR pattern)
  ! * Number of swath lines (or latitude limits) to be averaged for XTR pattern
  ! * Maximum number of calls to fitting function
  ! * Number of iterations for destriping
  ! * Absolute maximum column value (+/- range) to include in averaging
  ! ---------------------------------------------------------------------------
  REWIND (fit_ctrl_unit)
  CALL skip_to_filemark ( fit_ctrl_unit, destriping_str, tmpchar, file_read_stat )
  CALL error_check ( &
       file_read_stat, file_read_ok, pge_errstat_fatal, OMSAO_F_READ_FITCTRL_FILE, &
       modulename//f_sep//destriping_str, vb_lev_default, pge_error_status )
  IF ( pge_error_status >= pge_errstat_error ) RETURN
  READ (fit_ctrl_unit, *) yn_run_destriping
  READ (fit_ctrl_unit, *) yn_remove_ctrbias, ctr_bias_pol
  READ (fit_ctrl_unit, *) ctr_pol_base, ctr_pol_scal, ctr_pol_patt
  READ (fit_ctrl_unit, *) ctr_nblocks, ctrdst_latrange
  READ (fit_ctrl_unit, *) ctr_fitfunc_calls
  READ (fit_ctrl_unit, *) ctr_nloop

  ! --------------------------------------------------------
  ! Position cursor to read new shift and squeeze option gga
  ! --------------------------------------------------------
  REWIND (fit_ctrl_unit)
  CALL skip_to_filemark ( fit_ctrl_unit, newshift_str, tmpchar, file_read_stat )
  CALL error_check ( &
       file_read_stat, file_read_ok, pge_errstat_fatal, OMSAO_F_READ_FITCTRL_FILE, &
       modulename//f_sep//destriping_str, vb_lev_default, pge_error_status )
  IF ( pge_error_status >= pge_errstat_error ) RETURN
  READ (fit_ctrl_unit, *) yn_newshift

  ! -------------------------------------------------------------------
  ! Position cursor to read logical for Reference Sector Correction gga
  ! -------------------------------------------------------------------
  REWIND (fit_ctrl_unit)
  CALL skip_to_filemark ( fit_ctrl_unit, refseccor_str, tmpchar, file_read_stat )
  CALL error_check ( &
       file_read_stat, file_read_ok, pge_errstat_fatal, OMSAO_F_READ_FITCTRL_FILE, &
       modulename//f_sep//destriping_str, vb_lev_default, pge_error_status )
  IF ( pge_error_status >= pge_errstat_error ) RETURN
  READ (fit_ctrl_unit, *) yn_refseccor

  ! ------------------------------------------------------------------
  ! Position cursor to read logical for Scattering Weights and ... gga
  ! ------------------------------------------------------------------
  REWIND (fit_ctrl_unit)
  CALL skip_to_filemark ( fit_ctrl_unit, scattweight_str, tmpchar, file_read_stat )
  CALL error_check ( &
       file_read_stat, file_read_ok, pge_errstat_fatal, OMSAO_F_READ_FITCTRL_FILE, &
       modulename//f_sep//destriping_str, vb_lev_default, pge_error_status )
  IF ( pge_error_status >= pge_errstat_error ) RETURN
  READ (fit_ctrl_unit, *) yn_sw

  ! -------------------------------------------------------------------
  ! Unless we come up with a reason against it, the maximum good column
  ! also applies to the destriping procedure.
  ! -------------------------------------------------------------------
  ctr_maxcol = max_good_col

  ! -------------------------------------------------------------------------
  ! Determine minimum and maximum wavelength in selected read/fitting windows
  ! -------------------------------------------------------------------------
  winwav_min = MINVAL((/ fit_winwav_lim(1:n_fit_winwav) /))
  winwav_max = MAXVAL((/ fit_winwav_lim(1:n_fit_winwav) /))

  errstat = pge_errstat_ok

  ! ----------------------------------------------------------------------
  ! At this point we check whether we have to read from UV2 or VIS file.
  ! This is set throught the 3-letter string preceeding the start and end
  ! wavelengths of the fitting window. Of course, the PCF must contain the
  ! corresponding file under the LUN associated with the radiance file.
  ! ----------------------------------------------------------------------
  SELECT CASE ( l1b_channel )
  CASE ( 'UV1' )
     l1b_radiance_esdt  = 'OML1BRUG'
  CASE ( 'UV2' )
     l1b_radiance_esdt  = 'OML1BRUG'
  CASE ( 'VIS' )
     l1b_radiance_esdt  = 'OML1BRVG'
  END SELECT

  CALL find_radiance_fitting_variables ( errstat )

  IF ( yn_doas ) THEN
     pm_one     = -1.0_r8
  ELSE
     pm_one     =  1.0_r8
  END IF

  ! -----------------------------------------------
  ! Close fitting control file, report SUCCESS read
  ! -----------------------------------------------
  errstat = PGS_IO_GEN_CLOSEF ( fit_ctrl_unit )
  errstat = PGS_SMF_TESTSTATUSLEVEL(errstat)
  CALL error_check ( &
       errstat, pgs_smf_mask_lev_s, pge_errstat_warning, OMSAO_W_CLOSE_FITCTRL_FILE, &
       modulename//f_sep//TRIM(ADJUSTL(static_input_fnames(icf_idx))),               &
       vb_lev_default, pge_error_status )

  CALL error_check ( &
       0, 1, pge_errstat_ok, OMSAO_S_READ_FITCTRL_FILE,                &
       modulename//f_sep//TRIM(ADJUSTL(static_input_fnames(icf_idx))), &
       vb_lev_default, pge_error_status )

  RETURN
END SUBROUTINE read_fitting_control_file

SUBROUTINE get_mols_for_fitting ( tmpchar, n_mol_fit, fitcol_idx, errstat )

  USE OMSAO_precision_module,  ONLY: i4
  USE OMSAO_indices_module,    ONLY: refspec_strings, max_rs_idx
  USE OMSAO_parameters_module, ONLY: max_mol_fit
  USE OMSAO_errstat_module,    ONLY: pge_errstat_error, pge_errstat_ok

  IMPLICIT NONE

  ! ===============
  ! Input variables
  ! ===============
  CHARACTER (LEN=*), INTENT (IN) :: tmpchar

  ! ================
  ! Output variables
  ! ================
  INTEGER (KIND=i4),                          INTENT (OUT) :: errstat, n_mol_fit
  INTEGER (KIND=i4), DIMENSION (max_mol_fit), INTENT (OUT) :: fitcol_idx

  ! ===============
  ! Local variables
  ! ===============
  INTEGER   (KIND=i4)          :: i, ncl, sstart, sidx
  LOGICAL                      :: yn_eoc
  CHARACTER (LEN=LEN(tmpchar)) :: tmpsub, local_tmpchar

  errstat = pge_errstat_ok

  ! ----------------------------
  ! Initialize output quantities
  ! ----------------------------
  n_mol_fit = 0  ;  fitcol_idx = 0

  ! -----------------------------
  ! Copy tmpchar to local_tmpchar
  ! -----------------------------
  local_tmpchar = tmpchar

  ! -----------------------------------------------
  ! Get names and indices for main molecules to fit
  ! -----------------------------------------------
  sstart = 0 ; ncl = 0 ; yn_eoc = .FALSE.
  getmolnames: DO i = 1, max_mol_fit
     ! ---------------------------------------------------------
     ! Extract index string, find index, then extract file name.
     ! ---------------------------------------------------------
     CALL get_substring ( local_tmpchar, sstart, tmpsub, ncl, yn_eoc )

     IF ( ncl > 0 ) THEN
        CALL string2index ( refspec_strings, max_rs_idx, tmpsub(1:ncl), sidx )
        IF ( sidx > 0 ) THEN
           n_mol_fit = n_mol_fit + 1
           fitcol_idx(n_mol_fit) = sidx
        END IF
     END IF
     IF ( yn_eoc ) EXIT getmolnames
  END DO getmolnames
  IF ( n_mol_fit == 0 .OR. ALL(fitcol_idx == 0) ) errstat = pge_errstat_error

  RETURN
END SUBROUTINE get_mols_for_fitting


SUBROUTINE find_radiance_fitting_variables ( errstat )

  USE OMSAO_precision_module
  USE OMSAO_errstat_module
  USE OMSAO_indices_module,      ONLY:                     &
       max_rs_idx, max_calfit_idx, mns_idx, mxs_idx,       &
       calfit_titles,  radfit_titles,  refspec_titles,     &
       calfit_strings, radfit_strings, refspec_strings,    &
       o3_t1_idx, o3_t2_idx, o3_t3_idx, bro_idx, comm_idx, &
       hwe_idx, asy_idx, ad1_idx, ad2_idx
  USE OMSAO_variables_module,    ONLY: &
       n_fitvar_rad, all_radfit_idx, mask_fitvar_rad, fitvar_rad_init,         &
       lo_radbnd, up_radbnd, n_fincol_idx, n_mol_fit, fitcol_idx, fincol_idx,  &
       yn_common_iter, common_fitpos, ozone_idx, ozone_log
  USE OMSAO_prefitcol_module, ONLY:                                            &
       o3_prefit_fitidx, bro_prefit_fitidx, yn_bro_prefit, yn_o3_prefit,       &
       n_prefit_vars
  USE OMSAO_omidata_module,      ONLY: &
       correlation_names, correlation_names_concat, nclenfit

  IMPLICIT NONE

  ! ---------------------------------------------------------------------
  ! Correlation Matrix Output: Only the correlations
  ! with the first main fitting element - FINCOL_IDX(1,1) - are written
  ! to file. This is only relevant for OMSAO3, where the final column can
  ! be composed of more than one fitting index.
  ! ---------------------------------------------------------------------

  ! ----------------
  ! Output Variables
  ! ----------------
  INTEGER (KIND=i4), INTENT (INOUT) :: errstat

  ! ---------------
  ! Local Variables
  ! ---------------
  INTEGER (KIND=i4) :: i, j, k, idx, locerrstat

  CHARACTER (LEN=maxchlen), EXTERNAL :: int2string

  ! CCM Fit Lineshape for every spectrum
  LOGICAL :: yn_fit_itf = .TRUE.

  locerrstat = pge_errstat_ok

  ! -------------------------------------------------------------
  ! Find the indices of those variables that are actually varied
  ! during the fitting, and save those in MASK_FITVAR_RAD. Save
  ! the number of varied parameters in N_FITVAR_RAD.
  !
  ! In addition, we need to determine the fitting indices that
  ! will make up the final fitted column of the molecule(s) in
  ! question. For this we have to jump through a double loop:
  ! Since we are compressing the fitting parameter array to 
  ! include only the varied parameters, the final covariance 
  ! matrix, which is crucial for determining the uncertainties,
  ! only knows the compressed indices. Therefore we have to have
  ! an "index of an index" type array, that remembers the index
  ! position of indices of the final molecule(s), AS THEY APPEAR
  ! IN THE COMPRESSED FITTING PARAMETER LIST.
  !
  ! For the latter task it is easier to split the loops into
  ! calibration parameters (unrelated to the final fitted column)
  ! and reference spectra parameters.
  !
  ! Later addition: ALL_RADFIT_IDX contains the indices of the
  ! fitting parameters without, in the case of the trace gases,
  ! sub-indices, but with MAX_CALFIT_IDX added for uniqueness.
  ! Remembering these values here saves us some headaches when
  ! we collect the columns of all fitting parameters and want to
  ! match those with reference spectra normalization factors.
  ! -------------------------------------------------------------
  n_fitvar_rad = 0 ; mask_fitvar_rad = 0; correlation_names = ""; nclenfit = 0
  all_radfit_idx = 0

  ! --------------------------------
  ! First the calibration parameters
  ! --------------------------------
  DO i = 1, max_calfit_idx

     IF ( i == hwe_idx .OR. i == asy_idx ) THEN

        IF ( .NOT. yn_fit_itf ) CYCLE

     END IF
     
     IF ( lo_radbnd(i) < up_radbnd(i) ) THEN
        n_fitvar_rad = n_fitvar_rad + 1
        mask_fitvar_rad(n_fitvar_rad) = i
        all_radfit_idx (n_fitvar_rad) = i
        correlation_names(n_fitvar_rad) = &
             TRIM(ADJUSTL(calfit_strings(i)))//": "//TRIM(ADJUSTL(calfit_titles(i)))
     END IF

  END DO

  ! ------------------------------------
  ! Now the reference spectra parameters
  ! ------------------------------------
  n_fincol_idx  = 0
  n_prefit_vars = 0
  DO i = 1, max_rs_idx
     idx = max_calfit_idx + (i-1)*mxs_idx
     DO j = mns_idx, mxs_idx        
        ! ----------------------------------------------
        ! Assign only entries that are varied in the fit
        ! ----------------------------------------------
        IF ( lo_radbnd(idx+j) < up_radbnd(idx+j) ) THEN
           n_fitvar_rad = n_fitvar_rad + 1
           mask_fitvar_rad(n_fitvar_rad) = idx+j
           all_radfit_idx (n_fitvar_rad) = max_calfit_idx + i
           correlation_names(n_fitvar_rad) = &
                TRIM(ADJUSTL(refspec_strings(i)))//' - '//TRIM(ADJUSTL(radfit_strings(j)))//&
                ": "//TRIM(ADJUSTL(refspec_titles(i)))//" "//TRIM(ADJUSTL(radfit_titles(j)))

           ! --------------------------------------------------------------------------
           ! And here the loop over the final column molecules. We have to match the
           ! FITCOL_IDX with the current molecule index, and then remember the position
           ! of the fitting index in the MASK_FITVAR_RAD array. The second index
           ! remembers the reference spectrum that is associated with this molecule, so
           ! that we can easily access its normalization factor.
           ! --------------------------------------------------------------------------
           getfincol: DO k = 1, n_mol_fit
              IF ( fitcol_idx(k) == i ) THEN
                 n_fincol_idx = n_fincol_idx + 1
                 fincol_idx (1,n_fincol_idx) = n_fitvar_rad
                 fincol_idx (2,n_fincol_idx) = i
                 EXIT getfincol
              END IF
           END DO getfincol
           ! ----------------------------------------------------------------
           ! To select which TOMS profile is to be used in AMF calculation we
           ! need an idea of O3 column
           ! ----------------------------------------------------------------
           IF ( i .EQ. o3_t1_idx .OR. i .EQ. o3_t2_idx .OR. i .EQ. o3_t1_idx &
                .AND. ozone_log )  THEN
              ozone_idx = n_fitvar_rad
              ozone_log = .FALSE.
           ENDIF

        END IF

        ! -------------------------------------------------------------------------
        ! Check for any indices that correspond to pre-fitted columns
        ! -------------------------------------------------------------------------
        SELECT CASE ( i )
        CASE ( bro_idx )
           IF ( yn_bro_prefit(1) .AND. &
                (fitvar_rad_init(idx+j) /= 0.0_r8 .OR. &
                (lo_radbnd(idx+j)       < up_radbnd(idx+j))) ) bro_prefit_fitidx   = idx+j
        CASE ( o3_t1_idx:o3_t3_idx )
           IF ( yn_o3_prefit(1) .AND. &
                (fitvar_rad_init(idx+j) /= 0.0_r8 .OR. &
                (lo_radbnd(idx+j)       < up_radbnd(idx+j))) ) o3_prefit_fitidx(i) = idx+j
        CASE DEFAULT
        END SELECT

     END DO

     ! ---------------------------------------------------------------
     ! Finally, check for Common Mode Iteration. If selected, we have
     ! to add the fitting parameter information here but keep it fixed
     ! at Zero for the first round. This is a bit of a kludge, and it
     ! wastes time because we are carrying along a Zero value fitting
     ! parameter, but since the initialization of the output file and
     ! also the set-up of temporary arrays require all dimensions, it
     ! is really the easiest way to do things.
     ! ---------------------------------------------------------------
     IF ( i == comm_idx .AND. yn_common_iter ) THEN

        ! ------------------------------------------------------------------
        ! The Common Mode fitting position has been saved from the
        ! fitting control file
        ! ------------------------------------------------------------------
        j = common_fitpos
        idx = max_calfit_idx + (i-1)*mxs_idx+j

        !n_fitvar_rad = n_fitvar_rad + 1
        !mask_fitvar_rad(n_fitvar_rad) = idx
        !all_radfit_idx (n_fitvar_rad) = max_calfit_idx + i
        !correlation_names(n_fitvar_rad) = &
        !     TRIM(ADJUSTL(refspec_strings(i)))//' - '//TRIM(ADJUSTL(radfit_strings(j)))//&
        !     ": "//TRIM(ADJUSTL(refspec_titles(i)))//" "//TRIM(ADJUSTL(radfit_titles(j)))

        lo_radbnd      (idx) = 0.0_r8
        up_radbnd      (idx) = 0.0_r8
        fitvar_rad_init(idx) = 0.0_r8
     END IF

  END DO

  ! -------------------------------------------------
  ! Find out which indices hold information for ozone
  ! -------------------------------------------------
  IF (ozone_log)       ozone_log = .FALSE.
  IF (.NOT. ozone_log) ozone_log = .TRUE.

  ! ------------------------------------------------------------------------------------
  ! Concatinate the names of fitting elements that will appear in the correlation matrix
  ! ------------------------------------------------------------------------------------
  correlation_names_concat = ''
  DO i = 1, n_fitvar_rad
     correlation_names_concat = &
          TRIM(ADJUSTL(correlation_names_concat))//TRIM(ADJUSTL(correlation_names(i)))
     IF ( i < n_fitvar_rad )    &
          correlation_names_concat = TRIM(ADJUSTL(correlation_names_concat))//','
  END DO
  nclenfit = LEN_TRIM(ADJUSTL(correlation_names_concat))

  errstat = MAX ( errstat, locerrstat )

  RETURN
END SUBROUTINE find_radiance_fitting_variables
