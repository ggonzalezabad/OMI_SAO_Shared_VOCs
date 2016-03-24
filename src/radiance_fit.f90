SUBROUTINE radiance_fit ( &
     pge_idx, ipix, n_fitres_loop, fitres_range, yn_reference_fit, &
     n_rad_wvl, curr_rad_spec,                                     &
     fitcol, rms, dfitcol, radfit_exval, radfit_itnum, chisquav,   &
     o3fit_cols, o3fit_dcols, brofit_cols, brofit_dcols,           &
     lqh2ofit_cols, lqh2ofit_dcols,                                &
     target_var, allfit, allerr, corrmat, yn_bad_pixel, fitspc_out )

  ! ***************************************************************
  !
  !   Perform solar wavelength calibration and slit width fitting
  !
  ! ***************************************************************

  USE OMSAO_precision_module
  USE OMSAO_indices_module,      ONLY: &
       solar_idx, n_max_fitpars, wvl_idx, spc_idx, sig_idx, ccd_idx, &
       max_calfit_idx, o3_t1_idx, o3_t2_idx, o3_t3_idx, bro_idx,     &
       pge_o3_idx, pge_hcho_idx, pge_gly_idx
  USE OMSAO_parameters_module,   ONLY: &
       r8_missval, i4_missval, i2_missval, downweight, elsunc_less_is_noise
  USE OMSAO_variables_module,    ONLY:                                    &
       n_mol_fit, fitcol_idx, n_fincol_idx, fincol_idx, pm_one, database, &
       n_database_wvl, yn_doas, yn_smooth, rad_wav_avg, fitvar_rad,       &
       fitvar_rad_init, fitvar_rad_saved,  n_fitvar_rad,                  &
       lo_radbnd, up_radbnd, lobnd, upbnd, fitweights, currspec, fitwavs, &
       fit_winwav_idx, mask_fitvar_rad, max_itnum_rad, refspecs_original, &
       xtrack_fitres_limit, all_radfit_idx, yn_o3amf_cor, &
       n_rad_wvl_max
  USE OMSAO_prefitcol_module, ONLY:                                       &
       n_prefit_vars, yn_o3_prefit, yn_bro_prefit, bro_prefit_var,        &
       o3_prefit_var, o3_prefit_fitidx, bro_prefit_fitidx,                &
       yn_lqh2o_prefit,lqh2o_prefit_var,lqh2o_prefit_fitidx
  USE OMSAO_omidata_module

  IMPLICIT NONE

  ! *******************************************************************
  ! CAREFUL: Assumes that radiance and solar wavelength arrays have the
  ! same number of points. That must not be the case if we read in a
  ! general EL1 file. Examine and adjust! (tpk, note to himself)
  ! *******************************************************************

  ! ---------------
  ! Input variables
  ! ---------------
  INTEGER (KIND=i4), INTENT (IN) :: pge_idx, ipix, n_rad_wvl, n_fitres_loop, fitres_range
  REAL    (KIND=r8), INTENT (IN) :: brofit_cols, brofit_dcols
  REAL    (KIND=r8), INTENT (IN) :: lqh2ofit_cols, lqh2ofit_dcols
  LOGICAL,           INTENT (IN) :: yn_reference_fit

  ! -----------------------------
  ! (Possibly) Modified Variables
  ! -----------------------------
  REAL (KIND=r8), DIMENSION (ccd_idx, n_rad_wvl),  INTENT (INOUT) :: curr_rad_spec
  REAL (KIND=r8), DIMENSION (o3_t1_idx:o3_t3_idx), INTENT (INOUT) :: o3fit_cols, o3fit_dcols

  ! ------------------
  ! Modified variables
  ! ------------------
  LOGICAL,                                     INTENT (OUT) :: yn_bad_pixel
  INTEGER (KIND=i4),                           INTENT (OUT) :: radfit_exval, radfit_itnum
  REAL    (KIND=r8),                           INTENT (OUT) :: fitcol, rms, dfitcol, chisquav
  REAL    (KIND=r8), DIMENSION (n_fitvar_rad), INTENT (OUT) :: allfit, allerr, corrmat
  REAL    (KIND=r8), DIMENSION (n_fincol_idx), INTENT (OUT) :: target_var

  ! CCM Also return fitted spectrum
  ! REAL    (KIND=r8), DIMENSION (n_rad_wvl),    INTENT (OUT) :: fitspc_out 
  REAL    (KIND=r8), DIMENSION (n_rad_wvl_max), INTENT (OUT) :: fitspc_out 

  ! ---------------
  ! Local variables
  ! ---------------
  INTEGER (KIND=i4) :: i, j, idx, j1, j2, k1, k2, l, ll_rad, lu_rad, index
  INTEGER (KIND=i4) :: n_fitwav_rad, locitnum, n_nozero_wgt
  REAL    (KIND=r8)                                         :: asum, ssum
  REAL    (KIND=r8)                                         :: mean, mdev, sdev, loclim, normfac, mfac
  REAL    (KIND=r8), DIMENSION (n_rad_wvl)                  :: fitres, fitspec, tmp
  REAL    (KIND=r8), DIMENSION (n_max_fitpars)              :: fitvar
  REAL    (KIND=r8), DIMENSION (n_fitvar_rad, n_fitvar_rad) :: covar_matrix
  REAL    (KIND=r8), DIMENSION (:,:), ALLOCATABLE           :: covar

  REAL    (KIND=r8) :: fitcol_saved

  CHARACTER(LEN=2)  :: n_fitvar_rad_str

  EXTERNAL specfit_func, specfit_func_o3exp
  SAVE fitcol_saved

  yn_bad_pixel = .FALSE.

  radfit_exval = INT(i2_missval, KIND=i4)
  radfit_itnum = INT(i2_missval, KIND=i4)
  chisquav     = r8_missval

  IF ( pge_idx == pge_o3_idx ) THEN
     o3fit_cols = r8_missval  ;  o3fit_dcols = r8_missval
  END IF

  ! ============================================================
  ! Assign LL_RAD, LU_RAD, and SIG for each earthshine radiance
  ! ============================================================
  ll_rad = fit_winwav_idx(2)  ;  lu_rad = fit_winwav_idx(3)

  fitwavs   (1:n_rad_wvl) = curr_rad_spec(wvl_idx,1:n_rad_wvl)
  currspec  (1:n_rad_wvl) = curr_rad_spec(spc_idx,1:n_rad_wvl)
  fitweights(1:n_rad_wvl) = curr_rad_spec(sig_idx,1:n_rad_wvl)

  ! ---------------------------------------------------------------
  ! High pass filtering for DOAS. First, take log (rad/irrad), then
  ! filter by subtracting a cubic, then re-add the log (irradiance).
  ! This way we are fitting the log (rad) with the proper filtering
  ! having been done. The spectrum subroutine will start by
  ! subtracting the irradiance spectrum, where the proper shifting
  ! and squeezing can take place. In this high-pass filtering, we
  ! are ignoring the small extra wavelength calibration change, as
  ! for Ring effect, above.
  ! ---------------------------------------------------------------
  IF ( yn_doas ) THEN
     currspec(1:n_rad_wvl) = LOG ( currspec(1:n_rad_wvl) / database(solar_idx,1:n_rad_wvl) )
     CALL subtract_cubic_meas (fitwavs(1:n_rad_wvl), n_rad_wvl, currspec(1:n_rad_wvl), ll_rad, lu_rad)
     currspec(1:n_rad_wvl) = currspec(1:n_rad_wvl) + LOG ( database(solar_idx,1:n_rad_wvl) )
  END IF

  ! --------------------------------------------------------------------
  ! Apply smoothing (1/16,1/4,3/8,1/4,1/16); 2/98 uhe/ife recommendation
  ! --------------------------------------------------------------------
  IF ( yn_smooth ) THEN
     tmp(1:n_rad_wvl) = currspec(1:n_rad_wvl)
     currspec (3:n_rad_wvl-2) = 0.375_r8 * tmp (3:n_rad_wvl-2) +  &
          0.25_r8   * (tmp (4:n_rad_wvl-1) + tmp (2:n_rad_wvl-3)) +  &
          0.0625_r8 * (tmp (5:n_rad_wvl) + tmp (1:n_rad_wvl-4))
  END IF

  ! ---------------------------------------
  ! Compute average of radiance wavelengths
  ! ---------------------------------------
  asum = SUM ( fitwavs(1:n_rad_wvl) * ( fitweights(1:n_rad_wvl)*fitweights(1:n_rad_wvl) ) )
  ssum = SUM (             1.0_r8   * ( fitweights(1:n_rad_wvl)*fitweights(1:n_rad_wvl) ) )
  rad_wav_avg = asum / ssum
  
  radfit_exval = 0
     
  ! --------------------------------------------------------------------
  ! Initialize the fitting variables with the initial guess. Rather than
  ! subjecting ourselves to the vagarities of a wrong convergence, we
  ! bite the computationally more expensive bullet of starting from
  ! scratch each and every time.
  ! --------------------------------------------------------------------

  ! -----------------------------------------------------------
  ! Initialize the fitting variables. FITVAR_RAD_SAVED has been
  ! set to the initial values in the calling routine outside the
  ! pixel loop. Here we use FITVAR_RAD_SAVED, which will be 
  ! updated with current values from the previous fit if that 
  ! fit has gone well.
  ! -----------------------------------------------------------
  fitvar_rad(1:n_max_fitpars) = fitvar_rad_saved(1:n_max_fitpars)

  ! ---------------------------------------------------
  ! For HCHO we check whether we will be doing pre-fits
  ! ---------------------------------------------------------------------------
  ! Note that the conversion to slant columns (BrO only) and the multiplication
  ! with the cross section normalization has already been done at the time when
  ! the pre-fitted columns are read in.
  ! ---------------------------------------------------------------------------
  SELECT CASE ( pge_idx )
  CASE( pge_hcho_idx )
     ! ---
     ! BrO
     ! ---
     IF ( yn_bro_prefit(1) ) THEN
        bro_prefit_var = 0.0_r8
        IF ( brofit_cols > r8_missval ) THEN
           fitvar_rad(bro_prefit_fitidx) = brofit_cols
           IF ( yn_bro_prefit(2) ) THEN
              lo_radbnd(bro_prefit_fitidx) = brofit_cols - 1.0_r8 * brofit_dcols
              up_radbnd(bro_prefit_fitidx) = brofit_cols + 1.0_r8 * brofit_dcols
           ELSE
              lo_radbnd(bro_prefit_fitidx) = brofit_cols
              up_radbnd(bro_prefit_fitidx) = brofit_cols
              bro_prefit_var               = brofit_cols
           END IF
        END IF
     END IF
     ! ---
     ! O3
     ! ---
     IF ( yn_o3_prefit(1) ) THEN
        o3_prefit_var(o3_t1_idx:o3_t3_idx) = 0.0_r8
        DO j = o3_t1_idx, o3_t3_idx
           IF ( o3fit_cols(j) > r8_missval ) THEN
              fitvar_rad(o3_prefit_fitidx(j)) = o3fit_cols(j)
              IF ( yn_o3_prefit(2) ) THEN
                 lo_radbnd(o3_prefit_fitidx(j)) = o3fit_cols(j) - 2.0_r8 * o3fit_dcols(j)
                 up_radbnd(o3_prefit_fitidx(j)) = o3fit_cols(j) + 2.0_r8 * o3fit_dcols(j)
              ELSE
                 lo_radbnd(o3_prefit_fitidx(j)) = o3fit_cols(j)
                 up_radbnd(o3_prefit_fitidx(j)) = o3fit_cols(j)
                 o3_prefit_var(j)               = o3fit_cols(j)
              END IF
           END IF
        END DO
     END IF
  
  CASE( pge_gly_idx )
  
  	! ------------
    ! Liquid Water
    ! ------------
    IF ( yn_lqh2o_prefit(1) ) THEN
       lqh2o_prefit_var = 0.0_r8
       IF ( lqh2ofit_cols > r8_missval ) THEN
          fitvar_rad(lqh2o_prefit_fitidx) = lqh2ofit_cols
          IF ( yn_lqh2o_prefit(2) ) THEN
             lo_radbnd(lqh2o_prefit_fitidx) = lqh2ofit_cols - 1.0_r8 * lqh2ofit_dcols
             up_radbnd(lqh2o_prefit_fitidx) = lqh2ofit_cols + 1.0_r8 * lqh2ofit_dcols
          ELSE
             lo_radbnd(lqh2o_prefit_fitidx) = lqh2ofit_cols
             up_radbnd(lqh2o_prefit_fitidx) = lqh2ofit_cols
             lqh2o_prefit_var               = lqh2ofit_cols
          END IF
       END IF
    END IF
    
  CASE DEFAULT
  	! Do nothing
  END SELECT

  ! -----------------------------------------------------------------
  ! Create a condensed array of fitting variables that only contains
  ! the ones that are varied. This considerably reduces the execution
  ! time of the fitting routine.
  ! -----------------------------------------------------------------
  fitvar = 0.0_r8 ; lobnd = 0.0_r8 ; upbnd = 0.0_r8
  DO i = 1, n_fitvar_rad
     idx       = mask_fitvar_rad(i)
     fitvar(i) = fitvar_rad(idx)
     lobnd (i) = lo_radbnd(idx)
     upbnd (i) = up_radbnd(idx)
  END DO

  ! --------------------------------------------------------------------
  ! Check whether we enough spectral points to carry out the fitting. If
  ! not, call it a bad pixel and return.
  ! --------------------------------------------------------------------
  IF ( (n_fitvar_rad-n_prefit_vars) >= n_rad_wvl ) THEN
     yn_bad_pixel = .TRUE.  ;  RETURN
  END IF

  ! -------------------------------------
  ! Allocate memory for COVARIANCE MATRIX
  ! -------------------------------------
  ALLOCATE ( covar(1:n_fitvar_rad,1:n_fitvar_rad) )
  covar = r8_missval
  ! ------------------------------------------------------------------------------
  ! For certain PGEs we use a fitting function that corrects the O3 cross sections
  ! with a quadratic polynomial in wavelength, to account for the fact that O3
  ! absorption can vary greatly in magnitude over the fitting window (HH bands).
  ! ------------------------------------------------------------------------------
  IF ( yn_o3amf_cor ) THEN
     CALL specfit (                                                    &
          n_fitvar_rad, fitvar(1:n_fitvar_rad), n_rad_wvl,             &
          lobnd(1:n_fitvar_rad), upbnd(1:n_fitvar_rad), max_itnum_rad, &
          covar(1:n_fitvar_rad, 1:n_fitvar_rad), fitspec(1:n_rad_wvl), &
          fitres(1:n_rad_wvl), radfit_exval, locitnum, specfit_func_o3exp )
  ELSE
     CALL specfit (                                                    &
          n_fitvar_rad, fitvar(1:n_fitvar_rad), n_rad_wvl,             &
          lobnd(1:n_fitvar_rad), upbnd(1:n_fitvar_rad), max_itnum_rad, &
          covar(1:n_fitvar_rad, 1:n_fitvar_rad), fitspec(1:n_rad_wvl), &
          fitres(1:n_rad_wvl), radfit_exval, locitnum, specfit_func      )
  END IF

  covar_matrix(1:n_fitvar_rad,1:n_fitvar_rad) = covar(1:n_fitvar_rad,1:n_fitvar_rad)

  ! ------------------------------------------
  ! Assign iteration number from the first fit
  ! ------------------------------------------
  radfit_itnum = locitnum
   
  n_nozero_wgt = MAX ( INT ( ANINT ( SUM(fitweights(1:n_rad_wvl)) ) ), 1 )
  mean         = SUM  ( fitres(1:n_rad_wvl) )                 / REAL(n_nozero_wgt, KIND=r8)
  sdev         = SQRT ( SUM ( (fitres(1:n_rad_wvl)-mean)**2 ) / REAL(n_nozero_wgt-1, KIND=r8) )
  mdev         = SUM  ( ABS(fitres(1:n_rad_wvl)-mean) )       / REAL(n_nozero_wgt, KIND=r8)
  loclim       = mean + REAL(fitres_range, KIND=r8)*sdev

  ! ----------------------
  ! Fitting RMS and CHI**2
  ! ----------------------
  IF ( n_nozero_wgt > 0 ) THEN
     rms     = SQRT ( SUM ( fitres(1:n_rad_wvl)**2 ) / REAL(n_nozero_wgt, KIND=r8) )
     ! ---------------------------------------------
     ! This gives the same CHI**2 as the NR routines
     ! ---------------------------------------------
     chisquav = SUM  ( fitres(1:n_rad_wvl)**2 )
  ELSE
     rms      = r8_missval
     chisquav = r8_missval
  END IF

  IF ( yn_reference_fit ) xtrack_fitres_limit(ipix) = loclim

  IF ( ( n_fitres_loop                    >  0             ) .AND. &
       ( loclim                           >  0.0_r8        ) .AND. &
       ( MAXVAL(ABS(fitres(1:n_rad_wvl))) >= loclim        ) .AND. &
       ( n_nozero_wgt                     >  n_fitvar_rad  )          ) THEN

     fitloop: DO j = 1, n_fitres_loop

        WHERE ( ABS(fitres(1:n_rad_wvl)) >= loclim )
           fitweights(1:n_rad_wvl) = downweight
        END WHERE
        
        IF ( yn_o3amf_cor ) THEN 
           CALL specfit (                                                     &
                n_fitvar_rad, fitvar(1:n_fitvar_rad), n_rad_wvl,              &
                lobnd(1:n_fitvar_rad), upbnd(1:n_fitvar_rad), max_itnum_rad,  &
                covar(1:n_fitvar_rad, 1:n_fitvar_rad), fitspec(1:n_rad_wvl),  &
                fitres(1:n_rad_wvl), radfit_exval, locitnum, specfit_func_o3exp )
        ELSE
           CALL specfit (                                                     &
                n_fitvar_rad, fitvar(1:n_fitvar_rad), n_rad_wvl,              &
                lobnd(1:n_fitvar_rad), upbnd(1:n_fitvar_rad), max_itnum_rad,  &
                covar(1:n_fitvar_rad, 1:n_fitvar_rad), fitspec(1:n_rad_wvl),  &
                fitres(1:n_rad_wvl), radfit_exval, locitnum, specfit_func       )
        END IF
        covar_matrix(1:n_fitvar_rad,1:n_fitvar_rad) = covar(1:n_fitvar_rad,1:n_fitvar_rad)

        ! -----------------------------
        ! Add any subsequent iterations
        ! -----------------------------
        radfit_itnum = radfit_itnum + locitnum

        n_nozero_wgt = MAX ( INT ( ANINT ( SUM(fitweights(1:n_rad_wvl)) ) ), 1 )
        mean         = SUM  ( fitres(1:n_rad_wvl) )                 / REAL(n_nozero_wgt, KIND=r8)
        sdev         = SQRT ( SUM ( (fitres(1:n_rad_wvl)-mean)**2 ) / REAL(n_nozero_wgt-1, KIND=r8) )
        loclim       = mean + REAL(fitres_range, KIND=r8)*sdev

        ! ----------------------
        ! Fitting RMS and CHI**2
        ! ----------------------
        IF ( n_nozero_wgt > 0 ) THEN
           rms     = SQRT ( SUM ( fitres(1:n_rad_wvl)**2 ) / REAL(n_nozero_wgt, KIND=r8) )
           ! ---------------------------------------------
           ! This gives the same CHI**2 as the NR routines
           ! ---------------------------------------------
           chisquav = SUM  ( fitres(1:n_rad_wvl)**2 )
        ELSE
           rms      = r8_missval
           chisquav = r8_missval
        END IF

        ! ---------------------------------------------------------------
        ! Exit iteration loop if fitting residual is within contstraints.
        ! Save the (spike-adjusted) radiance reference fitting weights .
        ! ---------------------------------------------------------------
        IF ( MAXVAL(ABS(fitres(1:n_rad_wvl))) < loclim ) EXIT fitloop

     END DO fitloop
  END IF

  ! ----------------------------------------
  ! De-allocate memory for COVARIANCE MATRIX
  ! ----------------------------------------
  IF ( ALLOCATED (covar) ) DEALLOCATE (covar)

  ! ---------------------------------------------------------------------
  ! Save correlation information from covariance matrix !gga to real corr
  ! elation. Only correlation of the main variable (retrieved molecule) w
  ! ith the other variables is kept
  ! ---------------------------------------------------------------------
  corrmat = r8_missval
  index   = fincol_idx(1,1)
  DO i = 1, n_fitvar_rad
     IF (covar_matrix (i,i) .EQ. 0.0 .OR. covar_matrix(index,index) .EQ. 0.0) CYCLE
     corrmat(i) =       covar_matrix ( index, i               ) / &
                  SQRT( covar_matrix ( index, index ) *           &
                        covar_matrix ( i              , i )   )
  END DO

  ! --------------------------------------------------------------------
  ! Save fitting weights for possible use through radiance reference fit
  ! --------------------------------------------------------------------
  curr_rad_spec(sig_idx,1:n_rad_wvl) = fitweights(1:n_rad_wvl)

  ! CCM save fitted spectrum
  fitspc_out(1:n_rad_wvl) = fitspec(1:n_rad_wvl)

  ! --------------------------------------------------------------------
  ! At this point we have to make sure we only report the fitted columns
  ! if EXVAL, the exit variable from the fitting routine, is >=0. All
  ! negative values mean trouble and are likely to have produced strange
  ! uncertainties. In any case, the values for the fitted 
  ! columns should not be trusted.
  ! --------------------------------------------------------------------
  SELECT CASE ( radfit_exval )
  CASE ( 0: )     

     ! ---------------------------
     ! Update common mode spectrum
     ! ---------------------------
     IF ( .NOT. yn_reference_fit ) THEN
        !fitres(1:n_rad_wvl) = fitres(1:n_rad_wvl) !/ fitweights(1:n_rad_wvl)
        CALL compute_common_mode ( &
             yn_reference_fit, ipix, n_rad_wvl, fitwavs(1:n_rad_wvl), &
             fitres(1:n_rad_wvl), .FALSE. )
     END IF

     ! =====================================================================
     ! Compute the actual number of radiance wavelengths used in the fitting
     ! =====================================================================
     n_fitwav_rad = &
          INT ( SUM ( 1.0_r8 * (fitweights(1:n_rad_wvl)*fitweights(1:n_rad_wvl)) ) )

     ! --------------------------------------------------------------------------
     ! Assign total column. We have done the preliminary work with FITCOL_IDX and
     ! PM_ONE, so we don't have to "IF DOAS" here.
     !
     ! Reminder about the (admittedly confusing) variable names:
     !
     !   FITCOL:       Fitted column
     !   DFITCOL:      Uncertainty of fitted column
     !   N_FINCOL_IDX: Number of fitting variables that make up the final FITCOL
     !                 (e.g., O3 at more than one temperature).
     !   FINCOL_IDX:   Array of dimension (2,N_FINCOL_IDX*MXS_IDX), where MXS_IDX
     !                 is the maximum fitting sub-index (e.g., AD1_IDX, LBE_IDX).
     !                 FINCOL_IDX(1,*) carries the relative indices of the varied
     !                 final column variables in the array that was passed to
     !                 the fitting routine. FINCOL_IDX(2,*) contains the index
     !                 for the associated reference spectrum; this we need for
     !                 access to the normalization factor. See subroutine
     !                 READ_CONTROL_FILE for assignment of these indices.
     !
     !  TARGET_VAR     Saved fitting parameter values for possible use to remove
     !                 target gas from radiance reference
     ! --------------------------------------------------------------------------
     fitcol = 0.0_r8  ;  dfitcol = 0.0_r8 ; target_var = 1.0_r8
     DO i = 1, n_fincol_idx
        ! --------------------------------------------------
        ! First add the contribution of the diagonal element
        ! --------------------------------------------------
        j1 = fincol_idx(1, i) ; k1 = fincol_idx(2,i)

        fitcol  = fitcol  + pm_one * fitvar(j1) / refspecs_original(k1)%NormFactor
        dfitcol = dfitcol + covar_matrix(j1,j1) / refspecs_original(k1)%NormFactor**2

        IF ( pge_idx == pge_o3_idx ) THEN
           o3fit_cols (k1) = fitvar(j1)   / refspecs_original(k1)%NormFactor
           ! -------------------------------------------------------------------------
           ! Just as for the combined column uncertainty, we can't compute O3FIT_DCOLS
           ! if any of the major elements are missing.
           ! -------------------------------------------------------------------------
           IF ( (n_fitwav_rad - n_fitvar_rad + n_prefit_vars) <= 0 .OR. &
                (rms == r8_missval) .OR. (n_rad_wvl == 0)                 ) THEN
              o3fit_dcols(k1) = r8_missval
           ELSE
              o3fit_dcols(k1) = rms * &
                   SQRT ( covar_matrix(j1,j1) / refspecs_original(k1)%NormFactor**2 *  &
                   REAL ( n_fitwav_rad, KIND=r8 ) /                                    &
                   REAL ( n_fitwav_rad-n_fitvar_rad+n_prefit_vars, KIND=r8 ) )
           END IF
        END IF

        ! -------------------------------------------------------------------------
        ! Then add the contributions from off-diagonal elements (correlations)
        ! -------------------------------------------------------------------------
        DO l = i+1, n_fincol_idx
           j2 = fincol_idx(1,l) ; k2 = fincol_idx(2,l)
           dfitcol = dfitcol + 2.0_r8 * covar_matrix(j1,j2) / &
                (refspecs_original(k1)%NormFactor*refspecs_original(k2)%NormFactor)
        END DO

        ! ------------------------------------------
        ! Save fitting parameter value of target gas
        ! ------------------------------------------
        target_var(i) = fitvar(j1)

     END DO

     ! ------------------------------------------------------------
     ! Fitting columns and uncertainties for all fitting parameters
     ! ------------------------------------------------------------
     DO i = 1, n_fitvar_rad
        allfit(i) = fitvar(i)
        allerr(i) = covar_matrix(i,i)
        j = all_radfit_idx(i)
        IF ( j > max_calfit_idx ) THEN
           j = j - max_calfit_idx
           normfac = refspecs_original(j)%NormFactor
           IF ( normfac == 0.0_r8 ) normfac = 1.0_r8
           allfit(i) = allfit(i) / normfac
           allerr(i) = allerr(i) / normfac**2
        END IF
     END DO

     ! ---------------------------------------------------------------------
     ! Computing rationale for DFITCOL:
     !
     ! We are in the CASE branch where iteration has converged to something,
     ! but we don't know whether that something actually makes sense. It is
     ! well possible that some quantities are missing because of numerical
     ! overflows, to name just one. By setting DFITCOL to "missing" we 
     ! acknowledge that we have no idea what the fitting uncertainties are
     ! for the current pixel. Setting DFITCOL to a large value (+1E+30 for
     ! example) would be an alternative, but that is just another way of
     ! expressing our ingnorance about what really is going on.
     ! ---------------------------------------------------------------------
     IF ( ( (n_fitwav_rad-n_fitvar_rad+n_prefit_vars) <= 0 ) .OR. &
          ( rms == r8_missval                              ) .OR. &
          ( n_rad_wvl == 0                                 )        ) THEN
        dfitcol                = r8_missval
        allerr(1:n_fitvar_rad) = r8_missval
     ELSE
        ! --------------------------------------------------------
        ! Adjust fitting errors according to formulas given in 
        ! Numerical Recipes, $15.5: Multiply by SQRT[chi**2/(n-2)]
        ! --------------------------------------------------------
        mfac = SQRT ( chisquav / REAL(n_fitwav_rad-n_fitvar_rad+n_prefit_vars, KIND=r8) )
        dfitcol                = SQRT ( dfitcol )                * mfac
        allerr(1:n_fitvar_rad) = SQRT ( allerr(1:n_fitvar_rad) ) * mfac
     END IF

     ! -------------------------------------------------------------------
     ! The following assignment makes sense only because FITVAR_RAD is
     ! updated with FITVAR (using the proper mask) in SPECTRUM_EARTHSHINE.
     !
     ! HOWEVER: The OMI data are rather noisy, and fitting uncertainties
     !          are generally large. Hence this assignment is rather
     !          dangerous. By going back to the initial guess, we may
     !          lose some speed, but we gain predictability.
     ! -------------------------------------------------------------------
     !IF ( (radfit_exval >= INT(elsunc_less_is_noise, KIND=i4)) .AND. &
     !     (fitcol+1.0_r8*dfitcol >= 0.0_r8)                    .AND. &
     !     ( .NOT. yn_reference_fit )                                  )  THEN
     !   fitvar_rad_saved(1:n_max_fitpars) = fitvar_rad(1:n_max_fitpars)
     !ELSE
        fitvar_rad_saved(1:n_max_fitpars) = fitvar_rad_init(1:n_max_fitpars)
     !END IF
     ! ======================================================
     ! Anything but EXVAL > 0 means that trouble has occurred, 
     ! and the fit most likely is screwed.
     ! ======================================================
  CASE ( :-1 )
     fitcol = r8_missval ; dfitcol = r8_missval
  CASE DEFAULT
     ! -------------------------------------------------------------------
     ! We should never reach here, because the above CASE statements cover
     ! all possible values of EXVAL. But better safe than sorry.
     fitcol = r8_missval ; dfitcol = r8_missval
     ! -------------------------------------------------------------------
  END SELECT

  fitcol_saved = fitcol

  RETURN
END SUBROUTINE radiance_fit

