SUBROUTINE xtrack_radiance_wvl_calibration (             &
     yn_radiance_reference, yn_solar_comp,               &
     first_pix, last_pix, n_max_rspec, n_comm_wvl_out, errstat )

  USE OMSAO_precision_module
  USE OMSAO_indices_module,    ONLY: &
       wvl_idx, spc_idx, sig_idx, max_calfit_idx, max_rs_idx, hwe_idx, asy_idx,  &
       shi_idx, squ_idx, oclo_idx, o2o2_idx, bro_idx, solar_idx, pge_bro_idx,    &
       pge_oclo_idx, ccd_idx, radcal_idx
  USE OMSAO_parameters_module, ONLY: maxchlen, downweight, normweight
  USE OMSAO_variables_module,  ONLY:  &
       verb_thresh_lev, hw1e, e_asym, n_rad_wvl, curr_rad_spec, rad_spec_wavcal, &
       sol_wav_avg, database, fitvar_cal, fitvar_cal_saved,                      &
       fitvar_rad_init, pge_idx, rad_wght_wavcal,         &
       n_fitres_loop, fitres_range, yn_diagnostic_run
  USE OMSAO_slitfunction_module, ONLY: saved_shift, saved_squeeze
  USE OMSAO_omidata_module ! , exept_this_one => n_comm_wvl
  USE OMSAO_radiance_ref_module, ONLY: radref_latrange
  USE OMSAO_errstat_module
  USE EZspline_obj
  USE EZspline

  IMPLICIT NONE

  ! ---------------
  ! Input variables
  ! ---------------
  INTEGER (KIND=i4), INTENT (IN) :: first_pix, last_pix, n_max_rspec
  LOGICAL,           INTENT (IN) :: yn_radiance_reference, yn_solar_comp

  ! ---------------
  ! Output variable
  ! ---------------
  INTEGER (KIND=i4), INTENT (OUT) :: n_comm_wvl_out

  ! -----------------
  ! Modified variable
  ! -----------------
  INTEGER (KIND=i4), INTENT (INOUT) :: errstat

  ! ---------------
  ! Local variables
  ! ---------------
  INTEGER   (KIND=i2)      :: radcal_itnum
  INTEGER   (KIND=i4)      :: locerrstat, ipix, radcal_exval, i, imax, n_ref_wvl !, nxtloc, xtr_add
  REAL      (KIND=r8)      :: chisquav, rad_spec_avg, rms
  LOGICAL                  :: yn_skip_pix, yn_bad_pixel, yn_full_range
  CHARACTER (LEN=maxchlen) :: addmsg
  INTEGER (KIND=i4), DIMENSION (4)             :: select_idx
  INTEGER (KIND=i4), DIMENSION (2)             :: exclud_idx
  REAL    (KIND=r8), ALLOCATABLE, DIMENSION(:) :: fitres, fitspec
  REAL    (KIND=r8), DIMENSION (n_max_rspec)   :: ref_wvl, ref_spc, ref_wgt, rad_wvl

  ! ------------------------------
  ! Name of this module/subroutine
  ! ------------------------------
  CHARACTER (LEN=31), PARAMETER :: modulename = 'xtrack_radiance_wvl_calibration'


  locerrstat = pge_errstat_ok

  fitvar_cal_saved(1:max_calfit_idx) = fitvar_rad_init(1:max_calfit_idx)

  ! -------------------------------------------------
  ! Set the number of wavelengths for the common mode
  ! -------------------------------------------------
  n_comm_wvl_out = MAXVAL ( omi_nwav_radref(first_pix:last_pix) )
  IF ( MAXVAL(omi_nwav_rad(first_pix:last_pix,0)) > n_comm_wvl_out ) &
    n_comm_wvl_out = MAXVAL(omi_nwav_rad(first_pix:last_pix,0))

  ! --------------------------------
  ! Loop over cross-track positions. 
  ! --------------------------------
  XTrackWavCal: DO ipix = first_pix, last_pix

     locerrstat = pge_errstat_ok

     curr_xtrack_pixnum = ipix

     ! ---------------------------------------------------------------------
     ! If we already determined that this cross track pixel position carries
     ! an error, we don't even have to start processing.
     ! ---------------------------------------------------------------------
     IF ( omi_cross_track_skippix(ipix) ) CYCLE


     ! ---------------------------------------------------------------------------
     ! For each cross-track position we have to initialize the saved Shift&Squeeze
     ! ---------------------------------------------------------------------------
     saved_shift = -1.0e+30_r8 ; saved_squeeze = -1.0e+30_r8

     ! ----------------------------------------------------
     ! Assign number of radiance and irradiance wavelengths
     ! ----------------------------------------------------
     n_omi_irradwvl = omi_nwav_irrad(ipix  )
     n_omi_radwvl   = omi_nwav_rad  (ipix,0)

     ! -----------------------------------------------------------------
     ! tpk: Should the following be "> n_fitvar_rad"??? No, because that
     !      value is set only inside OMI_ADJUST_RADIANCE_DATA!!!
     ! -----------------------------------------------------------------
     IF ( n_omi_irradwvl <= 0 .OR. n_omi_radwvl <= 0 ) CYCLE

     ! ---------------------------------------------------------------
     ! Restore solar fitting variables for across-track reference in
     ! Earthshine fitting. Use the Radiance References if appropriate.
     ! ---------------------------------------------------------------
     sol_wav_avg = omi_sol_wav_avg(ipix)
     hw1e        = omi_solcal_pars(hwe_idx,ipix)
     e_asym      = omi_solcal_pars(asy_idx,ipix)

     ! -----------------------------------------------------
     ! Assign (hopefully predetermined) "reference" weights.
     ! -----------------------------------------------------
     IF ( .NOT. yn_solar_comp ) THEN
        n_omi_irradwvl            = omi_nwav_irrad(ipix)
        ref_wgt(1:n_omi_irradwvl) = omi_irradiance_wght(1:n_omi_irradwvl,ipix)

       ! -----------------------------------------------------
       ! Catch the possibility that N_OMI_RADWVL > N_OMI_IRRADWVL
       ! -----------------------------------------------------
       IF ( n_omi_radwvl > n_omi_irradwvl ) THEN
          i = n_omi_radwvl - n_omi_irradwvl
          ref_wgt(n_omi_irradwvl+1:n_omi_irradwvl+i) = downweight
          n_omi_irradwvl = n_omi_radwvl
       END IF
     ELSE
        n_omi_irradwvl          = n_omi_radwvl
        ref_wgt(1:n_omi_radwvl) = normweight
     END IF

     ! ---------------------------------------------------------------
     ! If a Radiance Reference is being used, then it must be calibrated
     ! rather than the swath line that has been read.
     ! ---------------------------------------------------------------
     IF ( yn_radiance_reference ) THEN
        omi_radiance_wavl(1:n_omi_radwvl,ipix,0) = omi_radref_wavl(1:n_omi_radwvl,ipix)
        omi_radiance_spec(1:n_omi_radwvl,ipix,0) = omi_radref_spec(1:n_omi_radwvl,ipix)
        omi_radiance_qflg(1:n_omi_radwvl,ipix,0) = omi_radref_qflg(1:n_omi_radwvl,ipix)
     END IF

     ! ---------------------------------------------------------------------------
     ! Set up generic fitting arrays. Remember that OMI_RADIANCE_XXX arrays are
     ! 3-dim with the last dimension being the scan line numbers. For the radiance
     ! wavelength calibration we only have one scan line at index "0".
     ! ---------------------------------------------------------------------------
     select_idx(1:4) = omi_ccdpix_selection(ipix,1:4)
     exclud_idx(1:2) = omi_ccdpix_exclusion(ipix,1:2)


     CALL omi_adjust_radiance_data ( &           ! Set up generic fitting arrays
          select_idx(1:4), exclud_idx(1:2),            &
          n_omi_radwvl,                                &
          omi_radiance_wavl  (1:n_omi_radwvl,ipix,0),  &
          omi_radiance_spec  (1:n_omi_radwvl,ipix,0),  &
          omi_radiance_qflg  (1:n_omi_radwvl,ipix,0),  &
          omi_radiance_ccdpix(1:n_omi_radwvl,ipix,0),  &
          n_omi_irradwvl, ref_wgt(1:n_omi_irradwvl),   &
          n_rad_wvl, curr_rad_spec(wvl_idx:ccd_idx,1:n_omi_radwvl), rad_spec_avg, &
          yn_skip_pix )

     ! ------------------------------------------------------------------------------------
     IF ( yn_skip_pix .OR. locerrstat >= pge_errstat_error ) THEN
        errstat = MAX ( errstat, locerrstat )
        omi_cross_track_skippix (ipix) = .TRUE.
        addmsg = ''
        WRITE (addmsg, '(A,I2)') 'SKIPPING cross track pixel #', ipix
        CALL error_check ( 0, 1, pge_errstat_warning, OMSAO_W_SKIPPIX, &
             modulename//f_sep//TRIM(ADJUSTL(addmsg)), vb_lev_default, &
             locerrstat )
        CYCLE
     END IF

     ! -----------------------------------------------------
     ! Assign the solar average wavelength - the wavelength
     ! calibration will not converge without it!
     ! -----------------------------------------------------
     sol_wav_avg = &
          SUM ( curr_rad_spec(wvl_idx,1:n_omi_radwvl) ) / REAL(n_omi_radwvl,KIND=r8)
     yn_bad_pixel = .FALSE.

     ALLOCATE( fitspec(1:n_rad_wvl))
     ALLOCATE( fitres(1:n_rad_wvl))
     CALL radiance_wavcal ( &                       ! Radiance wavelength calibration
          ipix, n_fitres_loop(radcal_idx), fitres_range(radcal_idx),       &
          n_rad_wvl, curr_rad_spec(wvl_idx:ccd_idx,1:n_rad_wvl),           &
          radcal_exval, radcal_itnum, chisquav, yn_bad_pixel, rms,         &
          fitspec, fitres, locerrstat )

     IF ( yn_bad_pixel .OR. locerrstat >= pge_errstat_error ) THEN
        errstat = MAX ( errstat, locerrstat )
        omi_cross_track_skippix (ipix) = .TRUE.
        addmsg = ''
        WRITE (addmsg, '(A,I2)') 'SKIPPING cross track pixel #', ipix
        CALL error_check ( 0, 1, pge_errstat_warning, OMSAO_W_SKIPPIX, &
             modulename//f_sep//TRIM(ADJUSTL(addmsg)), vb_lev_default, &
             locerrstat )
        CYCLE
     END IF
     ! ------------------------------------------------------------------------------------
           
     addmsg = ''
     WRITE (addmsg, '(A,I2,4(A,1PE10.3),2(A,I5))') 'RADIANCE Wavcal    #', ipix, &
          ': hw 1/e = ', hw1e, '; e_asy = ', e_asym, '; shift = ', &
          fitvar_cal(shi_idx), '; squeeze = ', fitvar_cal(squ_idx), &
         '; exit val = ', radcal_exval, '; iter num = ', radcal_itnum
     CALL error_check ( &
          0, 1, pge_errstat_ok, OMSAO_S_PROGRESS, TRIM(ADJUSTL(addmsg)), &
          vb_lev_omidebug, locerrstat )
     IF ( verb_thresh_lev >= vb_lev_screen ) WRITE (*, '(A)') TRIM(ADJUSTL(addmsg))

     ! ------------------------------------------
     ! Write out radiance calibration diagnostics
     ! ------------------------------------------
     CALL he5_radcal_write( &
          ipix, n_rad_wvl, radcal_exval, radcal_itnum, fitvar_cal(shi_idx), &
          hw1e, e_asym, rms, curr_rad_spec(wvl_idx, 1:n_rad_wvl),           &
          curr_rad_spec(spc_idx, 1:n_rad_wvl),                              &
          fitspec(1:n_rad_wvl), fitres(1:n_rad_wvl), radref_latrange, errstat)            

     DEALLOCATE(fitspec)
     DEALLOCATE(fitres)

     ! ---------------------------------
     ! Save crucial variables for output
     ! ---------------------------------
     omi_radcal_pars (1:max_calfit_idx,ipix) = fitvar_cal(1:max_calfit_idx)
     omi_radcal_xflag(ipix)                  = INT (radcal_exval, KIND=i2)
     omi_radcal_itnum(ipix)                  = INT (radcal_itnum, KIND=i2)
     omi_radcal_chisq(ipix)                  = chisquav
     ! -----------------------------------------------------------------------

     IF ( .NOT. (yn_radiance_reference) ) THEN
        n_ref_wvl = n_omi_irradwvl
        ref_wvl(1:n_ref_wvl) = omi_irradiance_wavl(1:n_ref_wvl,ipix)
        ref_spc(1:n_ref_wvl) = omi_irradiance_spec(1:n_ref_wvl,ipix)
        ref_wgt(1:n_ref_wvl) = omi_irradiance_wght(1:n_ref_wvl,ipix)
     ELSE
        n_ref_wvl = n_rad_wvl
        ref_wvl(1:n_ref_wvl) = curr_rad_spec(wvl_idx,1:n_rad_wvl)
        ref_spc(1:n_ref_wvl) = curr_rad_spec(spc_idx,1:n_rad_wvl)
        ref_wgt(1:n_ref_wvl) = curr_rad_spec(sig_idx,1:n_rad_wvl)

        omi_nwav_radref(ipix)             = n_ref_wvl
        omi_radref_wavl(1:n_ref_wvl,ipix) = curr_rad_spec(wvl_idx,1:n_rad_wvl)
        omi_radref_spec(1:n_ref_wvl,ipix) = curr_rad_spec(spc_idx,1:n_rad_wvl)
        omi_radref_wght(1:n_ref_wvl,ipix) = curr_rad_spec(sig_idx,1:n_rad_wvl)
     END IF

     ! ----------------------------------------------------
     ! Spline reference spectra to current wavelength grid.
     ! ----------------------------------------------------
     rad_wvl(1:n_rad_wvl) = curr_rad_spec(wvl_idx,1:n_rad_wvl)
     Call prepare_databases ( &
          ipix, n_ref_wvl, ref_wvl(1:n_ref_wvl), ref_spc(1:n_ref_wvl), &
          n_rad_wvl, rad_wvl(1:n_rad_wvl), n_max_rspec, locerrstat )
     ! --------------------------------------------------------------------------------

     IF ( locerrstat >= pge_errstat_error ) EXIT XTrackWavCal

     ! ---------------------------------------------------------
     ! Save DATABASE in OMI_DATABASE for radiance fitting loops.
     ! ---------------------------------------------------------
     omi_database (1:max_rs_idx,1:n_rad_wvl,ipix) = database (1:max_rs_idx,1:n_rad_wvl)
     n_omi_database_wvl(ipix)                     = n_rad_wvl
     omi_database_wvl(1:n_rad_wvl, ipix)          = curr_rad_spec(wvl_idx,1:n_rad_wvl)
     
     ! ----------------------------------------------------------------------
     ! Update the radiance reference with the wavelength calibrated values.
     ! ----------------------------------------------------------------------
     IF ( yn_radiance_reference ) THEN
        omi_radref_wavl(1:n_rad_wvl,ipix) = curr_rad_spec(wvl_idx,1:n_rad_wvl)
        omi_radref_spec(1:n_rad_wvl,ipix) = curr_rad_spec(spc_idx,1:n_rad_wvl)
        omi_radref_wght(n_rad_wvl+1:nwavel_max,ipix) = downweight

        ! --------------------------------------------------------
        ! Update the solar spectrum entry in OMI_DATABASE. First
        ! re-sample the solar reference spectrum to the OMI grid
        ! then assign to data base.
        !
        ! We need to keep the irradiance spectrum because we still
        ! have to fit the radiance reference, and we can't really
        ! do that against itself. In a later module the irradiance
        ! is replaced by the radiance reference.
        ! --------------------------------------------------------

        ! ------------------------------------------------------------------
        ! Prevent failure of interpolation by finding the maximum wavelength
        ! of the irradiance wavelength array.
        ! ------------------------------------------------------------------
        imax = MAXVAL ( MAXLOC ( omi_irradiance_wavl(1:n_omi_irradwvl,ipix) ) )
        !imin = MINVAL ( MINLOC ( omi_irradiance_wavl(1:imax,          ipix) ) )

        CALL interpolation ( &
             modulename, &
             imax, omi_irradiance_wavl(1:imax,ipix),                     &
             omi_irradiance_spec(1:imax,ipix),                           &
             n_rad_wvl, omi_database_wvl(1:n_rad_wvl,ipix),              &
             omi_database(solar_idx,1:n_rad_wvl,ipix),                   &
             'endpoints', 0.0_r8, yn_full_range, locerrstat )

        IF ( locerrstat >= pge_errstat_error ) THEN
           errstat = MAX ( errstat, locerrstat )
           omi_cross_track_skippix (ipix) = .TRUE.
           addmsg = ''
           WRITE (addmsg, '(A,I2)') 'SKIPPING cross track pixel #', ipix
           CALL error_check ( 0, 1, pge_errstat_warning, OMSAO_W_SKIPPIX, &
                modulename//f_sep//TRIM(ADJUSTL(addmsg)), vb_lev_default, &
                locerrstat )
           CYCLE
        END IF

     END IF

  END DO XTrackWavCal

  ! CCM Write splined/convolved databases if necessary
  IF( yn_diagnostic_run ) THEN
     ! omi_database maybe omi_database_wvl?
     CALL he5_write_omi_database(omi_database(1:max_rs_idx,1:n_rad_wvl,1:nxtrack_max), &
          omi_database_wvl(1:n_rad_wvl, 1:nxtrack_max), &
          max_rs_idx, n_rad_wvl, nxtrack_max, errstat)
     
  ENDIF

  errstat = MAX ( errstat, locerrstat )

  RETURN
END SUBROUTINE xtrack_radiance_wvl_calibration


SUBROUTINE xtrack_radiance_fitting_loop (                             &
     n_max_rspec, first_pix, last_pix, pge_idx, iloop,                &
     ctr_maxcol, n_fitvar_rad, allfit_cols, allfit_errs, corr_matrix, &
     target_var, errstat, fitspc_out, fitspc_out_dim0                 )

  USE OMSAO_precision_module
  USE OMSAO_indices_module,    ONLY: &
       wvl_idx, spc_idx, sig_idx, o3_t1_idx, o3_t3_idx, hwe_idx, asy_idx, shi_idx, squ_idx, &
       pge_o3_idx, pge_hcho_idx, n_max_fitpars, solar_idx, ccd_idx, radfit_idx, bro_idx,    &
       pge_gly_idx
  USE OMSAO_parameters_module, ONLY: &
       i2_missval, i4_missval, r4_missval, r8_missval, maxchlen, elsunc_less_is_noise
  USE OMSAO_variables_module,  ONLY:  &
       database, curr_sol_spec, n_rad_wvl, curr_rad_spec, sol_wav_avg, hw1e, e_asym,     &
       verb_thresh_lev, fitvar_rad_saved, fitvar_rad_init, n_database_wvl, &
       fitvar_rad, rad_wght_wavcal, n_fitres_loop, fitres_range, refspecs_original,      &
       yn_solar_comp, max_itnum_rad, szamax, n_fincol_idx, ozone_idx, ozone_log
  USE OMSAO_radiance_ref_module, ONLY: yn_radiance_reference, yn_reference_fit
  USE OMSAO_slitfunction_module, ONLY: saved_shift, saved_squeeze
  USE OMSAO_prefitcol_module, ONLY: &
       yn_o3_prefit,  o3_prefit_col,  o3_prefit_dcol,  &
       yn_bro_prefit, bro_prefit_col, bro_prefit_dcol, &
       yn_lqh2o_prefit, lqh2o_prefit_col, lqh2o_prefit_dcol
  USE OMSAO_omidata_module  ! nxtrack_max, ...
  USE OMSAO_errstat_module
     
  IMPLICIT NONE

  ! ---------------
  ! Input Variables
  ! ---------------
  REAL    (KIND=r8), INTENT (IN) :: ctr_maxcol
  INTEGER (KIND=i4), INTENT (IN) :: &
       pge_idx, iloop, first_pix, last_pix, n_max_rspec, n_fitvar_rad, &
       fitspc_out_dim0

  ! -----------------
  ! Modified variable
  ! -----------------
  INTEGER (KIND=i4), INTENT (INOUT) :: errstat
  REAL    (KIND=r8), INTENT (OUT  ), DIMENSION (n_fitvar_rad,first_pix:last_pix) :: &
       allfit_cols, allfit_errs, corr_matrix

  ! ---------------------------------------------------------
  ! Optional output variable (fitted variable for target gas)
  ! ---------------------------------------------------------
  REAL (KIND=r8), DIMENSION(n_fincol_idx,first_pix:last_pix), INTENT (OUT) :: target_var

  ! CCM Output fit spectra
  !REAL (KIND=r8), DIMENSION(n_comm_wvl,nxtrack_max,4), INTENT (OUT) :: fitspc_out	
  REAL (KIND=r8), DIMENSION(fitspc_out_dim0,nxtrack_max,4), INTENT (OUT) :: fitspc_out

  ! ---------------
  ! Local variables
  ! ---------------
  INTEGER (KIND=i4) :: locerrstat, ipix, radfit_exval, radfit_itnum
  REAL    (KIND=r8) :: fitcol, rms, dfitcol, chisquav, rad_spec_avg
  REAL    (KIND=r8) :: brofit_col, brofit_dcol
  REAL    (KIND=r8) :: lqh2ofit_col, lqh2ofit_dcol
  REAL    (KIND=r8), DIMENSION (o3_t1_idx:o3_t3_idx) :: o3fit_cols, o3fit_dcols
  LOGICAL                                     :: yn_skip_pix, yn_cycle_this_pix
  LOGICAL                                     :: yn_bad_pixel
  INTEGER (KIND=i4), DIMENSION (4)            :: select_idx
  INTEGER (KIND=i4), DIMENSION (2)            :: exclud_idx
  INTEGER (KIND=i4)                           :: n_solar_pts
  REAL    (KIND=r8), DIMENSION (n_max_rspec)  :: solar_wvl

  ! CCM Array for holding fitted spectra
  REAL    (KIND=r8), DIMENSION (fitspc_out_dim0) :: fitspc
  REAL    (KIND=i4)                           :: id

  CHARACTER (LEN=28), PARAMETER :: modulename = 'xtrack_radiance_fitting_loop'

  locerrstat = pge_errstat_ok

  !!!fitvar_rad_saved = fitvar_rad_init

  XTrackPix: DO ipix = first_pix, last_pix

     curr_xtrack_pixnum = ipix

     ! ---------------------------------------------------------------------
     ! If we already determined that this cross track pixel position carries
     ! an error, we don't even have to start processing.
     ! ---------------------------------------------------------------------
     IF ( omi_cross_track_skippix(ipix) .OR. szamax < omi_szenith(ipix,iloop) ) CYCLE
    
     locerrstat = pge_errstat_ok

     n_database_wvl = n_omi_database_wvl(ipix)
     n_omi_radwvl   = omi_nwav_rad      (ipix,iloop)

     ! ---------------------------------------------------------------------------
     ! For each cross-track position we have to initialize the saved Shift&Squeeze
     ! ---------------------------------------------------------------------------
     saved_shift = -1.0e+30_r8 ; saved_squeeze = -1.0e+30_r8

     ! ----------------------------------------------------------------------------
     ! Assign the solar wavelengths. Those should be current in the DATABASE array
     ! and can be taken from there no matter which case - YN_SOLAR_COMP and/or
     ! YN_RADIANCE_REFRENCE we are processing.
     ! ----------------------------------------------------------------------------
     n_solar_pts              = n_omi_database_wvl(ipix)
     if (n_solar_pts < 1) cycle  ! JED fix

     solar_wvl(1:n_solar_pts) = omi_database_wvl  (1:n_solar_pts, ipix)
     n_omi_irradwvl           = n_solar_pts

     CALL check_wavelength_overlap ( &
          n_fitvar_rad,                                                &
          n_solar_pts,          solar_wvl (1:n_solar_pts),             &
          n_omi_radwvl, omi_radiance_wavl (1:n_omi_radwvl,ipix,iloop), &
          yn_cycle_this_pix )

     IF ( yn_cycle_this_pix                .OR. &
          (n_database_wvl <= 0) .OR. (n_omi_radwvl <= 0) ) CYCLE
          !(n_database_wvl <= n_fitvar_rad) .OR. (n_omi_radwvl <= n_fitvar_rad) ) CYCLE

     ! ----------------------------------------------
     ! Restore DATABASE from OMI_DATABASE (see above)
     ! ----------------------------------------------
     database (1:max_rs_idx,1:n_database_wvl) = omi_database (1:max_rs_idx,1:n_database_wvl,ipix)
                 
     ! ---------------------------------------------------------------------------------
     ! Restore solar fitting variables for across-track reference in Earthshine fitting.
     ! Note that, for the YN_SOLAR_COMP case, some variables have been assigned already
     ! in the XTRACK_RADIANCE_WAVCAL loop.
     ! ---------------------------------------------------------------------------------
     sol_wav_avg                             = omi_sol_wav_avg(ipix)
     hw1e                                    = omi_solcal_pars(hwe_idx,ipix)
     e_asym                                  = omi_solcal_pars(asy_idx,ipix)
     curr_sol_spec(wvl_idx,1:n_database_wvl) = omi_database_wvl(1:n_database_wvl,ipix)
     curr_sol_spec(spc_idx,1:n_database_wvl) = omi_database    (solar_idx,1:n_database_wvl,ipix)
     ! --------------------------------------------------------------------------------

     omi_xtrackpix_no = ipix

     ! -------------------------------------------------------------------------
     select_idx(1:4) = omi_ccdpix_selection(ipix,1:4)
     exclud_idx(1:2) = omi_ccdpix_exclusion(ipix,1:2)
     CALL omi_adjust_radiance_data ( &           ! Set up generic fitting arrays
          select_idx(1:4), exclud_idx(1:2),                        &
          n_omi_radwvl,                                            &
          omi_radiance_wavl  (1:n_omi_radwvl,ipix,iloop),          &
          omi_radiance_spec  (1:n_omi_radwvl,ipix,iloop),          &
          omi_radiance_qflg  (1:n_omi_radwvl,ipix,iloop),          &
          omi_radiance_ccdpix(1:n_omi_radwvl,ipix,iloop),          &
          n_omi_radwvl, omi_radref_wght(1:n_omi_radwvl,ipix),      &
          n_rad_wvl, curr_rad_spec(wvl_idx:ccd_idx,1:n_omi_radwvl),&
          rad_spec_avg, yn_skip_pix )

     SELECT CASE ( pge_idx )
		 
     CASE (pge_hcho_idx)
        o3fit_cols (o3_t1_idx:o3_t3_idx) = o3_prefit_col (o3_t1_idx:o3_t3_idx,ipix,iloop)
        o3fit_dcols(o3_t1_idx:o3_t3_idx) = o3_prefit_dcol(o3_t1_idx:o3_t3_idx,ipix,iloop)
        brofit_col                       = bro_prefit_col (ipix,iloop)
        brofit_dcol                      = bro_prefit_dcol(ipix,iloop)
     
     CASE ( pge_gly_idx )
        lqh2ofit_col                     = lqh2o_prefit_col (ipix,iloop)
        lqh2ofit_dcol                    = lqh2o_prefit_dcol(ipix,iloop)
     CASE DEFAULT
     		! Nothing
     END SELECT

     ! --------------------
     ! The radiance fitting
     ! --------------------
     fitcol       = r8_missval
     dfitcol      = r8_missval
     radfit_exval = INT(i2_missval, KIND=i4)
     radfit_itnum = INT(i2_missval, KIND=i4)
     rms          = r8_missval

     yn_reference_fit = .FALSE.
     IF ( MAXVAL(curr_rad_spec(spc_idx,1:n_rad_wvl)) > 0.0_r8 .AND.     &
          n_rad_wvl > n_fitvar_rad .AND. (.NOT. yn_skip_pix)          ) THEN

        yn_bad_pixel = .FALSE.
        CALL radiance_fit ( &
             pge_idx, ipix, n_fitres_loop(radfit_idx), fitres_range(radfit_idx),   &
             yn_reference_fit,                                                     &
             n_rad_wvl, curr_rad_spec(wvl_idx:ccd_idx,1:n_rad_wvl),                &
             fitcol, rms, dfitcol, radfit_exval, radfit_itnum, chisquav,           &
             o3fit_cols, o3fit_dcols, brofit_col, brofit_dcol,                     &
             lqh2ofit_col, lqh2ofit_dcol,                                          &
             target_var(1:n_fincol_idx,ipix),                                      &
             allfit_cols(1:n_fitvar_rad,ipix), allfit_errs(1:n_fitvar_rad,ipix),   &
             corr_matrix(1:n_fitvar_rad,ipix), yn_bad_pixel, fitspc(1:n_rad_wvl) )

        IF ( yn_bad_pixel ) CYCLE
     END IF
    
     ! -----------------------------------
     ! Assign pixel values to final arrays
     ! -----------------------------------
     omi_fitconv_flag (ipix,iloop) = INT (radfit_exval, KIND=i2)
     omi_itnum_flag   (ipix,iloop) = INT (radfit_itnum, KIND=i2)
     omi_radfit_chisq (ipix,iloop) = chisquav
     omi_fit_rms      (ipix,iloop) = rms
     omi_column_amount(ipix,iloop) = fitcol
     omi_column_uncert(ipix,iloop) = dfitcol
     IF (ozone_log) THEN
        omi_ozone_amount(ipix,omi_iline)  = allfit_cols(ozone_idx,ipix)
     ENDIF

     ! CCM assign fit residual
     fitspc_out(1:n_rad_wvl,ipix,1) = fitspc(1:n_rad_wvl)
     fitspc_out(1:n_rad_wvl,ipix,2) = curr_rad_spec(spc_idx,1:n_rad_wvl)
     fitspc_out(1:n_rad_wvl,ipix,3) = curr_rad_spec(wvl_idx,1:n_rad_wvl)
     fitspc_out(1:n_rad_wvl,ipix,4) = curr_rad_spec(sig_idx,1:n_rad_wvl)

     IF ( pge_idx == pge_o3_idx ) THEN
        omi_o3_amount(o3_t1_idx:o3_t3_idx,ipix,iloop) = o3fit_cols (o3_t1_idx:o3_t3_idx)
        omi_o3_uncert(o3_t1_idx:o3_t3_idx,ipix,iloop) = o3fit_dcols(o3_t1_idx:o3_t3_idx)
     END IF

  END DO XTrackPix

  errstat = MAX ( errstat, locerrstat )

  RETURN
END SUBROUTINE xtrack_radiance_fitting_loop
