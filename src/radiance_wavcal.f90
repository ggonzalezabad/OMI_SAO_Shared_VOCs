SUBROUTINE radiance_wavcal (                              &
     ipix, n_fitres_loop, fitres_range, n_rad_wvl,        &
     curr_rad_spec, radcal_exval, radcal_itnum, chisquav, &
     yn_bad_pixel, rms, fitspec, fitres, errstat )

  USE OMSAO_precision_module,     ONLY: i2, i4, r8
  USE OMSAO_parameters_module,    ONLY: &
       i2_missval, i4_missval, r8_missval, downweight, elsunc_less_is_noise
  USE OMSAO_indices_module,       ONLY: &
       max_calfit_idx, shi_idx, squ_idx, wvl_idx, spc_idx, sig_idx, ccd_idx, &
       hwe_idx, asy_idx
  USE OMSAO_variables_module,   ONLY: &
       fitwavs, fitweights, currspec,                                 &
       fitvar_cal, fitvar_rad_init, fitvar_sol_init, fitvar_cal_saved,&
       mask_fitvar_cal, lo_sunbnd, up_sunbnd, n_fitvar_cal,           &
       fitvar_rad, n_fitvar_rad, lo_radbnd, up_radbnd, lobnd, upbnd,  &
       max_itnum_sol, hw1e, e_asym, rad_wav_avg, yn_newshift, sol_wav_avg
  USE OMSAO_errstat_module

  IMPLICIT NONE

  ! ----------------
  ! Input parameters
  ! ----------------
  INTEGER (KIND=i4), INTENT (IN) :: ipix, n_fitres_loop, n_rad_wvl, fitres_range

  ! -------------------
  ! Modified parameters
  ! -------------------
  LOGICAL,                                          INTENT (OUT)   :: yn_bad_pixel
  INTEGER (KIND=i2),                                INTENT (OUT)   :: radcal_itnum
  INTEGER (KIND=i4),                                INTENT (OUT)   :: radcal_exval
  REAL    (KIND=r8),                                INTENT (OUT)   :: chisquav, rms
  REAL    (KIND=r8), DIMENSION (n_rad_wvl),         INTENT (OUT)   :: fitres, fitspec
  INTEGER (KIND=i4),                                INTENT (INOUT) :: errstat
  REAL    (KIND=r8), DIMENSION (ccd_idx,n_rad_wvl), INTENT (INOUT) :: curr_rad_spec
  
  ! ---------------
  ! Local variables
  ! ---------------
  INTEGER (KIND=i4)  :: i, locerrstat, locitnum, n_nozero_wgt
  REAL    (KIND=r8)  :: mean, mdev, sdev, loclim
  REAL    (KIND=r8), DIMENSION (max_calfit_idx)    :: fitvar
  REAL    (KIND=r8), DIMENSION (:,:), ALLOCATABLE  :: covar

  EXTERNAL specfit_func_sol


  yn_bad_pixel = .FALSE.


  ! Select and wavelength calibrate radiance spectrum

  locerrstat = pge_errstat_ok

  radcal_exval = i4_missval
  radcal_itnum = i2_missval
  chisquav     = r8_missval

  ! ----------
  ! ELSUNC fit
  ! ----------
  fitwavs   (1:n_rad_wvl) = curr_rad_spec(wvl_idx, 1:n_rad_wvl)
  currspec  (1:n_rad_wvl) = curr_rad_spec(spc_idx, 1:n_rad_wvl)
  fitweights(1:n_rad_wvl) = curr_rad_spec(sig_idx, 1:n_rad_wvl)

  ! ------------------------------------------
  ! Update wavelegths for common mode spectrum
  ! (the .TRUE. in the call below selects the
  !  "wavelength update only" branch)
  ! ------------------------------------------
  CALL compute_common_mode ( &
       .TRUE., ipix, n_rad_wvl, fitwavs(1:n_rad_wvl), currspec(1:n_rad_wvl), .FALSE. )        

  ! -------------------------------------------------------------
  ! Initialize the fitting variables. FITVAR_CAL_SAVED has been
  ! set to the initial values in the calling routine. outside the
  ! pixel loop. Here we use FITVAR_CAL_SAVED, which will be 
  ! updated with current values from the previous fit if that fit
  ! has gone well.
  ! -------------------------------------------------------------
  !fitvar_cal(1:max_calfit_idx) = fitvar_cal_saved(1:max_calfit_idx)
  !fitvar_cal(1:max_calfit_idx) = fitvar_rad_init(1:max_calfit_idx)
  fitvar_cal(1:max_calfit_idx) = fitvar_sol_init(1:max_calfit_idx)

  ! -------------------------------------------------------------------------
  ! Keep the slit function variables from solar fit fixed. Remember to reduce
  ! the number of solar fitting variables if previously varied.
  ! -------------------------------------------------------------------------
  fitvar = 0.0_r8 ; lobnd = 0.0_r8 ; upbnd = 0.0_r8 ; n_fitvar_cal = 0
  DO i = 1, max_calfit_idx
     SELECT CASE ( i )
     CASE ( hwe_idx )
        fitvar_cal(hwe_idx) = hw1e
        lo_radbnd (hwe_idx) = hw1e
        up_radbnd (hwe_idx) = hw1e
     CASE ( asy_idx )
        fitvar_cal(asy_idx) = e_asym
        lo_radbnd (asy_idx) = e_asym
        up_radbnd (asy_idx) = e_asym
     CASE DEFAULT
        IF (lo_radbnd(i) < up_radbnd(i) ) THEN
           n_fitvar_cal  = n_fitvar_cal + 1
           mask_fitvar_cal(n_fitvar_cal) = i
           fitvar(n_fitvar_cal) = fitvar_cal(i)
           lobnd (n_fitvar_cal) = lo_radbnd(i)
           upbnd (n_fitvar_cal) = up_radbnd(i)
        END IF
     END SELECT
  END DO

  ! --------------------------------------------------------------------
  ! Check whether we enough spectral points to carry out the fitting. If
  ! not, call it a bad pixel and return.
  ! --------------------------------------------------------------------
  IF ( n_fitvar_cal >= n_rad_wvl ) THEN
     yn_bad_pixel = .TRUE.  ;  RETURN
  END IF

  ! -------------------------------------
  ! Allocate memory for COVARIANCE MATRIX
  ! -------------------------------------
  ALLOCATE ( covar(1:n_fitvar_cal,1:n_fitvar_cal) )

  !WRITE (*,'(A9,20(F12.6:))') '>>>>>>>>>', fitvar(1:n_fitvar_cal)
  CALL specfit (                                                    &
       n_fitvar_cal, fitvar(1:n_fitvar_cal), n_rad_wvl,             &
       lobnd(1:n_fitvar_cal), upbnd(1:n_fitvar_cal), max_itnum_sol, &
       covar(1:n_fitvar_cal,1:n_fitvar_cal), fitspec(1:n_rad_wvl),  &
       fitres(1:n_rad_wvl), radcal_exval, locitnum, specfit_func_sol )
  !WRITE (*,'(I6,I3,20(F12.6:))') radcal_exval, locitnum, fitvar(1:n_fitvar_cal)

  radcal_itnum = INT ( locitnum, KIND=i2 )

  ! ---------------------------------------------------------------------
  ! Attempt to standardize the re-iteration with spectral points excluded
  ! that have fitting residuals larger than a pre-set window. Needs more
  ! thinking before it can replace a simple window determined empirically
  ! from fitting lots of spectra.
  ! ---------------------------------------------------------------------
  n_nozero_wgt = INT ( ANINT ( SUM(fitweights(1:n_rad_wvl)) ) )
  mean         = SUM  ( fitres(1:n_rad_wvl) )                 / REAL(n_nozero_wgt,   KIND=r8)
  sdev         = SQRT ( SUM ( (fitres(1:n_rad_wvl)-mean)**2 ) / REAL(n_nozero_wgt-1, KIND=r8) )
  mdev         = SUM  ( ABS(fitres(1:n_rad_wvl)-mean) )       / REAL(n_nozero_wgt,   KIND=r8)
  loclim       = REAL (fitres_range, KIND=r8)*sdev

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

  IF ( ( n_fitres_loop                    >  0             ) .AND. &
       ( loclim                           >  0.0_r8        ) .AND. &
       ( MAXVAL(ABS(fitres(1:n_rad_wvl))) >= loclim        ) .AND. &
       ( n_nozero_wgt                     >  n_fitvar_cal  )  ) THEN

     fitloop: DO i = 1, n_fitres_loop
        WHERE ( ABS(fitres(1:n_rad_wvl)) > loclim )
           fitweights(1:n_rad_wvl) = downweight
        END WHERE

        !WRITE (*,'(A9,20(F12.6:))') '>>>>>>>>>', fitvar(1:n_fitvar_cal)
        CALL specfit (                                                    &
             n_fitvar_cal, fitvar(1:n_fitvar_cal), n_rad_wvl,             &
             lobnd(1:n_fitvar_cal), upbnd(1:n_fitvar_cal), max_itnum_sol, &
             covar(1:n_fitvar_cal,1:n_fitvar_cal), fitspec(1:n_rad_wvl),  &
             fitres(1:n_rad_wvl), radcal_exval, locitnum, specfit_func_sol )
        !WRITE (*,'(I6,I3,20(F12.6:))') radcal_exval, locitnum, fitvar(1:n_fitvar_cal)

        ! ----------------------
        ! Fitting RMS and CHI**2
        ! ----------------------
        n_nozero_wgt = INT ( ANINT ( SUM(fitweights(1:n_rad_wvl)) ) )
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
        
        radcal_itnum = radcal_itnum + INT ( locitnum, KIND=i2 )

        IF ( MAXVAL(ABS(fitres(1:n_rad_wvl))) <= loclim ) EXIT fitloop
     END DO fitloop
  END IF

  ! ----------------------------------------
  ! De-allocate memory for COVARIANCE MATRIX
  ! ----------------------------------------
  IF ( ALLOCATED (covar) ) DEALLOCATE (covar)

  ! ------------------------------------------------------------------
  ! The following assignment makes sense only because FITVAR_CAL is
  ! updated with FITVAR (using the proper mask) in SPECTRUM_SOLAR.
  ! ------------------------------------------------------------------
  IF ( radcal_exval >= INT(elsunc_less_is_noise, KIND=i4) ) THEN
     fitvar_cal_saved(1:max_calfit_idx) = fitvar_cal(1:max_calfit_idx)
  ELSE
     fitvar_cal_saved(1:max_calfit_idx) = fitvar_rad_init(1:max_calfit_idx)
  END IF


  ! -----------------------------------------------------------
  ! Reality check: set SQUEEZE to 1.0 to avoid division by Zero
  ! -----------------------------------------------------------
  IF ( fitvar_cal(squ_idx) == -1.0_r8 ) fitvar_cal(squ_idx) = 0.0_r8

  ! ---------------------
  ! Perform Shift&Squueze
  ! ---------------------
  ! gga to include Xiong comments
  IF (yn_newshift .EQV. .true.) THEN
     curr_rad_spec(wvl_idx,1:n_rad_wvl) = ( &
          curr_rad_spec(wvl_idx,1:n_rad_wvl) - fitvar_cal(shi_idx) + sol_wav_avg * fitvar_cal(squ_idx)) / &
          (1.0_r8 + fitvar_cal(squ_idx))
  ELSE
     curr_rad_spec(wvl_idx,1:n_rad_wvl) = ( &
          curr_rad_spec(wvl_idx,1:n_rad_wvl) - fitvar_cal(shi_idx) ) / (1.0_r8 + fitvar_cal(squ_idx))
  END IF

  ! -----------------------------------------------------------------
  ! We haven't implemented any error checks in this subroutine, hence
  ! it doesn't make sense yet to update the error variable. We'll do
  ! it anyway to make ourselves feel better.
  ! -----------------------------------------------------------------
  IF ( locerrstat /= pge_errstat_ok ) errstat = MAX ( errstat, locerrstat )

  RETURN
END SUBROUTINE radiance_wavcal
