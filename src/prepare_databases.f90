SUBROUTINE prepare_databases ( &
     xtrack_pix, n_sol_wvl, sol_wvl, sol_spc, n_rad_wvl, curr_rad_wvl, n_max_rspec, errstat )

  ! ===========================================
  !
  ! Called from XTRACK_RADIANCE_WVL_CALIBRATION
  !
  ! ===========================================


  USE OMSAO_precision_module
  USE OMSAO_variables_module, ONLY: hw1e, e_asym, phase, have_undersampling, database
  USE OMSAO_errstat_module
  IMPLICIT NONE

  ! ---------------
  ! Input variables
  ! ---------------
  INTEGER (KIND=i4),                        INTENT (IN) :: xtrack_pix, n_rad_wvl, n_sol_wvl, n_max_rspec
  REAL    (KIND=r8), DIMENSION (n_rad_wvl), INTENT (IN) :: curr_rad_wvl
  REAL    (KIND=r8), DIMENSION (n_sol_wvl), INTENT (IN) :: sol_wvl, sol_spc

  ! ---------------
  ! Output variable
  ! ---------------
  INTEGER (KIND=i4), INTENT (INOUT) :: errstat

  ! ---------------
  ! Local variables
  ! ---------------
  INTEGER (KIND=i4) :: locerrstat

  locerrstat = pge_errstat_ok

  ! ---------------------------------------------------------
  ! Spline external reference spectra to common radiance grid
  ! ---------------------------------------------------------
  CALL dataspline ( xtrack_pix, n_rad_wvl, curr_rad_wvl(1:n_rad_wvl), n_max_rspec, locerrstat )
  errstat = MAX ( errstat, locerrstat )
  IF ( errstat >= pge_errstat_error ) RETURN

  ! -----------------------------------
  ! Calculate the undersampled spectrum
  ! -----------------------------------
  IF ( ANY (have_undersampling) ) &
       CALL undersample (                                                               &
       xtrack_pix, n_rad_wvl, curr_rad_wvl(1:n_rad_wvl), hw1e, e_asym, phase, locerrstat )
  errstat = MAX ( errstat, locerrstat )
  IF ( errstat >= pge_errstat_error ) RETURN
  !IF ( ANY (have_undersampling) ) &
  !     CALL undersample_new (                                                             &
  !     xtrack_pix, n_rad_wvl, curr_rad_wvl(1:n_rad_wvl), n_sol_wvl, sol_wvl(1:n_sol_wvl), &
  !     hw1e, e_asym, locerrstat )
  !errstat = MAX ( errstat, locerrstat )
  !IF ( errstat >= pge_errstat_error ) RETURN

  ! -------------------------------------------------------------------------------------
  ! Calculate the splined fitting database. This will be saved into a nXtrack-dimensional
  ! array in a calling routine higher up (inside the fitting loop). 
  ! -------------------------------------------------------------------------------------
  CALL prepare_solar_refspec ( &
       n_sol_wvl, sol_wvl(1:n_sol_wvl), sol_spc(1:n_sol_wvl), n_rad_wvl, &
       curr_rad_wvl(1:n_rad_wvl), locerrstat )
  errstat = MAX ( errstat, locerrstat )
  IF ( errstat >= pge_errstat_error ) RETURN

  RETURN
END SUBROUTINE prepare_databases


SUBROUTINE prepare_solar_refspec ( &
     n_solpts, sol_wvl, sol_spc, n_radpts, curr_rad_wvl, errstat )

  ! ***********************************************************
  !
  !   Spline the solar measured solar spectrum to the radiance
  !   wavelength grid
  !
  !   Called from PREPARE_DATABASE
  !
  ! ***********************************************************

  USE OMSAO_precision_module
  USE OMSAO_indices_module, ONLY: max_rs_idx, solar_idx, ring_idx
  USE OMSAO_variables_module,  ONLY: yn_doas, fit_winwav_idx, &
       yn_smooth, n_refspec, database
  USE OMSAO_errstat_module

  IMPLICIT NONE

  ! ---------------
  ! Input variables
  ! ---------------
  INTEGER (KIND=i4),                       INTENT (IN)    :: n_radpts, n_solpts
  REAL    (KIND=r8), DIMENSION (n_solpts), INTENT (IN)    :: sol_wvl, sol_spc
  REAL    (KIND=r8), DIMENSION (n_radpts), INTENT (IN)    :: curr_rad_wvl
  INTEGER (KIND=i4),                       INTENT (INOUT) :: errstat

  ! ---------------
  ! Local variables
  ! ---------------
  LOGICAL                                 :: yn_full_range
  INTEGER (KIND=i4)                       :: j, ll_rad, lu_rad, locerrstat, n_sol_tmp
  REAL    (KIND=r8), DIMENSION (n_radpts) :: spline_sun
  REAL    (KIND=r8), DIMENSION (n_solpts) :: tmp_sol_spec, tmp_sol_wvl

  ! ------------------------------
  ! Name of this subroutine/module
  ! ------------------------------
  CHARACTER (LEN=21), PARAMETER :: modulename = 'prepare_solar_refspec'


  locerrstat = pge_errstat_ok

  ! ---------------------------------------------
  ! Spline irradiance spectrum onto radiance grid
  ! ---------------------------------------------
  ! There are two issues here:
  ! * The radiance and irradiance grids are likely to be off-set, i.e., at one
  !   or both end points we will run into an extrapolation.
  ! * Any bad pixel in the solar spectrum must be avoided in the interpolation.

  ! ----------------------------------------------------------------------
  ! First create an array with solar spectrum and wavelengths that is free
  ! of any bad pixels. "Bad" pixels are set to -1.0_r8
  ! ----------------------------------------------------------------------
  n_sol_tmp = 0
  DO j = 1, n_solpts
     IF ( sol_spc(j) > 0.0_r8 ) THEN
        n_sol_tmp = n_sol_tmp + 1
        tmp_sol_wvl (n_sol_tmp) = sol_wvl(j)
        tmp_sol_spec(n_sol_tmp) = sol_spc(j)
        ! ------------------------------------------------
        ! Check whether wavelengths are strictly ascending
        ! ------------------------------------------------
        IF ( n_sol_tmp > 1 ) THEN
           IF ( tmp_sol_wvl (n_sol_tmp) <= tmp_sol_wvl (n_sol_tmp-1) ) n_sol_tmp = n_sol_tmp - 1
        END IF
     END IF
  END DO

  ! ----------------------------------------------------------------
  ! Interpolate (spectral overlap is taken care of in INTERPOLATION)
  ! ----------------------------------------------------------------
  CALL interpolation ( &
       modulename,                                                                 &
       n_sol_tmp, tmp_sol_wvl(1:n_sol_tmp), tmp_sol_spec(1:n_sol_tmp),             &
       n_radpts, curr_rad_wvl(1:n_radpts), spline_sun(1:n_radpts),                 &
       'endpoints', 0.0_r8, yn_full_range, locerrstat )
  CALL error_check ( &
       locerrstat, pge_errstat_ok, pge_errstat_error, OMSAO_E_INTERPOL, &
       modulename, vb_lev_default, errstat )
  IF ( errstat >= pge_errstat_error ) RETURN
  IF ( .NOT. yn_full_range )   CALL error_check ( 0, 1, pge_errstat_warning, &
       OMSAO_W_INTERPOL_RANGE, modulename, vb_lev_develop, errstat )

  ! --------------------------------------------------------------------------
  ! Finally, we have to check whether we encountered any bad solar pixels. The
  ! adopted strategy is to exclude any such pixels from the fit, plus a window
  ! of +/- 0.1 nm around them.
  ! --------------------------------------------------------------------------
  DO j = 1, n_solpts
     IF ( sol_spc(j) < 0.0_r8 ) THEN
        WHERE ( ABS(curr_rad_wvl(1:n_radpts)-sol_wvl(j)) <= 0.1_r8)
           spline_sun(1:n_radpts) = -1.0_r8
        END WHERE
     END IF
  END DO

  ! --------------------------------------------------------
  ! Now we are ready to assign the splined solar spectrum to
  ! its final array - DATABASE(solar_idx,*).
  ! --------------------------------------------------------
  database(solar_idx,1:n_radpts) = spline_sun(1:n_radpts)

  ! =================================================================
  ! Note that the UNDERSAMPLING spectrum has already been assigned to
  ! DATABASE(us1/2_idx,*) in the UNDERSPEC routine.
  ! =================================================================

  ! ---------------------------------------------------------------------
  ! Set up Ring spectrum for high pass filtering of divided spectrum
  ! (afterward, add solar spectrum back in order to divide by solar
  ! spectrum with altered wavelength calibration in subroutine spectrum).
  ! ---------------------------------------------------------------------
  IF ( yn_doas ) database(ring_idx, 1:n_radpts) = &
       database(ring_idx, 1:n_radpts) / spline_sun(1:n_radpts)

  ! ---------------------------------------------------------
  ! For the DOAS case, high-pass filter the reference spectra
  ! ---------------------------------------------------------
  IF ( yn_doas ) THEN
     ll_rad = fit_winwav_idx(2) ; lu_rad = fit_winwav_idx(3)
     CALL subtract_cubic (curr_rad_wvl, n_radpts, ll_rad, lu_rad)
     database(ring_idx, 1:n_radpts) = &
          database(ring_idx, 1:n_radpts) * spline_sun(1:n_radpts)
  END IF

  IF ( yn_smooth ) database(1:max_rs_idx, 3:n_radpts-2) = &
       0.375_r8  * database(1:max_rs_idx, 3:n_radpts-2)  +  &
       0.25_r8   * (database(1:max_rs_idx, 4:n_radpts-1) + &
       database(1:max_rs_idx, 2:n_radpts-3)) + &
       0.0625_r8 * (database(1:max_rs_idx, 5:n_radpts)   + &
       database(1:max_rs_idx, 1:n_radpts-4))

  errstat = MAX ( errstat, locerrstat )

  RETURN
END SUBROUTINE prepare_solar_refspec
