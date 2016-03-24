SUBROUTINE spectrum_solar ( &
     npoints, nfitvar, smooth, sol_wav_avg, locwvl, fit, fitvar, doas )

  USE OMSAO_precision_module
  USE OMSAO_indices_module, ONLY: &
       wvl_idx, spc_idx,      &
       max_rs_idx, solar_idx, &
       bl0_idx, bl1_idx, bl2_idx, bl3_idx, bl4_idx, bl5_idx,&
       sc0_idx, sc1_idx, sc2_idx, sc3_idx, sc4_idx, sc5_idx,&
       sin_idx, hwe_idx, asy_idx, shi_idx, squ_idx
  USE OMSAO_parameters_module, ONLY: max_spec_pts
  USE OMSAO_variables_module,  ONLY: &
       refspecs_original, phase, fitwavs, solar_spec_convolved, yn_use_labslitfunc, &
       fitvar_cal, mask_fitvar_cal, yn_spectrum_norm, yn_newshift
  USE OMSAO_omidata_module,  ONLY: curr_xtrack_pixnum
  USE OMSAO_slitfunction_module
  USE OMSAO_errstat_module
  USE OMSAO_solcomp_module
  USE EZspline_obj
  USE EZspline

  IMPLICIT NONE


  INTEGER (KIND=i4),                      INTENT (IN)    :: npoints, nfitvar
  LOGICAL,                                INTENT (INOUT) :: smooth, doas
  REAL    (KIND=r8),                      INTENT (IN)    :: sol_wav_avg
  REAL    (KIND=r8), DIMENSION (nfitvar), INTENT (INOUT) :: fitvar
  REAL    (KIND=r8), DIMENSION (npoints), INTENT (INOUT) :: locwvl, fit

  ! ---------------
  ! Local variables
  ! ---------------
  LOGICAL                                                 :: yn_full_range
  REAL    (KIND=r8), DIMENSION (npoints)                  :: del, sunspec_ss
  INTEGER (KIND=i4)                                       :: i, idx, npts, errstat
  ! Shorthands for solar reference spectrum
  REAL    (KIND=r8), DIMENSION (refspecs_original(solar_idx)%nPoints) :: &
       solar_pos, solar_spec


  ! ------------------------------
  ! Name of this subroutine/module
  ! ------------------------------
  CHARACTER (LEN=14), PARAMETER :: modulename = 'spectrum_solar'


  errstat = pge_errstat_ok

  ! -----------------------------------------------------------------------------------
  ! First, we have to undo the compression of the FITVAR_SOL array. This compression
  ! is performed in the SOLAR_FIT and RADIANCE_WAVCAL subroutines and accelerates the
  ! fitting process, because ELSUNC has to handle less indices. But here we require the
  ! original layout, otherwise the index assingment is screwed.
  ! -----------------------------------------------------------------------------------
  DO i = 1, nfitvar
     idx = mask_fitvar_cal(i)
     fitvar_cal(idx) = fitvar(i)
  END DO

  npts               = refspecs_original(solar_idx)%nPoints
  solar_pos (1:npts) = refspecs_original(solar_idx)%RefSpecWavs(1:npts)
  solar_spec(1:npts) = refspecs_original(solar_idx)%RefSpecData(1:npts)
  IF ( .NOT. yn_spectrum_norm ) &
       solar_spec(1:npts) = solar_spec(1:npts) * refspecs_original(solar_idx)%NormFactor

  ! =========================================================================
  !     Spectrum Calculation for Solar and Radiance Wavelength Calibration
  ! =========================================================================
  
  !     Calculate the spectrum:
  !     First do the shift and squeeze. Shift by FITVAR(SHI_IDX), squeeze by
  !     1 + FITVAR(SQU_IDX); do in absolute sense, to make it easy to back-convert
  !     OMI data.
  !     Now, after Xiong recommendation if yn_newfit equal true then (gga):
  !     Lambda = Lambda * (1 + squeeze) + shift - sol_wav_avg * squeeze


  IF (yn_newshift .EQ. .true.) THEN ! gga
     solar_pos(1:npts) = solar_pos(1:npts) * (1.0_r8 + fitvar_cal(squ_idx)) + fitvar_cal(shi_idx) - sol_wav_avg * fitvar_cal(squ_idx)
  ELSE ! gga
     solar_pos(1:npts) = solar_pos(1:npts) * (1.0_r8 + fitvar_cal(squ_idx)) + fitvar_cal(shi_idx)
  END IF

  ! ----------------------------------------------
  ! Convolve only if we don't do a solar composite
  ! ----------------------------------------------
  IF ( yn_use_labslitfunc ) THEN
     ! ------------------------------------------------------------------------
     ! Only if either SHIFT or SQUEEZE have changed from the last iteration do
     ! we need to reconvolve the solar spectrum (its convolved value is saved
     ! in SOLAR_SPEC_CONVOLVED through MODULE association.
     !
     ! The choice of OMI lab slit function vs. Gaussian is made in the fitting
     ! control file: If the initial value of FITVAR(hwe_idx) is 0.0 then we are
     ! using the lab measurements, otherwise the Gaussian.
     ! ------------------------------------------------------------------------
     IF ( fitvar_cal(squ_idx) /= saved_squeeze .OR. &
          fitvar_cal(shi_idx) /= saved_shift ) THEN
        saved_squeeze = fitvar_cal(squ_idx)
        saved_shift   = fitvar_cal(shi_idx)
        solar_spec_convolved = 0.0_r8
        CALL omi_slitfunc_convolve (                                  &
             curr_xtrack_pixnum, npts, solar_pos(1:npts),             &
             solar_spec(1:npts), solar_spec_convolved(1:npts), errstat )
        CALL error_check ( &
             errstat, pge_errstat_ok, pge_errstat_error, OMSAO_E_INTERPOL, &
             modulename//f_sep//'Convolution', vb_lev_default, errstat )
        IF ( errstat >= pge_errstat_error ) RETURN
     END IF
  ELSE
     CALL asymmetric_gaussian_sf (                                           &
          npts, fitvar_cal(hwe_idx), fitvar_cal(asy_idx),                    &
          solar_pos(1:npts), solar_spec(1:npts), solar_spec_convolved(1:npts) )
  END IF
  

  ! =============================================
  ! Broadening and re-sampling of solar spectrum:
  ! =============================================
  ! Case for wavelength fitting of irradiance and radiance
  ! Broaden the solar reference by the hw1e value
  ! ------------------------------------------------------

  ! ------------------------------------------------------
  ! Re-sample the solar reference spectrum to the OMI grid
  ! ------------------------------------------------------
  CALL interpolation ( &
       modulename, npts, solar_pos(1:npts), solar_spec_convolved(1:npts), &
       npoints, locwvl(1:npoints), sunspec_ss(1:npoints), 'endpoints', 0.0_r8, &
       yn_full_range, errstat )
  CALL error_check ( &
       errstat, pge_errstat_ok, pge_errstat_error, OMSAO_E_INTERPOL, &
       modulename//f_sep//'Resampling to Solar Grid', vb_lev_default, errstat )
  IF ( errstat >= pge_errstat_error ) RETURN

  ! --------------------------------------------------------------------
  ! Add up the contributions, with solar intensity as FITVAR (SIN_IDX),
  ! to include possible linear and Beer's law forms.  Do these as 
  ! linear-Beer's-linear. In order to do DOAS we only need to be careful
  ! to include just linear contributions, since I already high-pass 
  ! filtered them.
  ! --------------------------------------------------------------------

  ! -----------
  !  Doing BOAS
  ! -----------
  fit(1:npoints) = fitvar_cal(sin_idx) * sunspec_ss(1:npoints)

  ! ----------------
  ! Add the scaling.
  ! ----------------
  del(1:npoints) = locwvl(1:npoints) - sol_wav_avg
  fit(1:npoints) = fit(1:npoints) * ( &
       fitvar_cal(sc0_idx)                                               + &
       fitvar_cal(sc1_idx) * del(1:npoints)                              + &
       fitvar_cal(sc2_idx) * del(1:npoints)*del(1:npoints)               + &
       fitvar_cal(sc3_idx) * del(1:npoints)*del(1:npoints)*del(1:npoints)+ &
       fitvar_cal(sc4_idx) * del(1:npoints)**4.0                         + &
       fitvar_cal(sc5_idx) * del(1:npoints)**5.0                          )
  
  ! ------------------------
  ! Add baseline parameters.
  ! ------------------------
  fit(1:npoints) = fit(1:npoints) + &
       fitvar_cal(bl0_idx)                                               + &
       fitvar_cal(bl1_idx) * del(1:npoints)                              + &
       fitvar_cal(bl2_idx) * del(1:npoints)*del(1:npoints)               + &
       fitvar_cal(bl3_idx) * del(1:npoints)*del(1:npoints)*del(1:npoints)+ &
       fitvar_cal(bl4_idx) * del(1:npoints)**4.0                         + &
       fitvar_cal(bl5_idx) * del(1:npoints)**5.0                          


  RETURN
END SUBROUTINE spectrum_solar

SUBROUTINE spectrum_earthshine ( &
     npts, n_fitvar, smooth, rad_wav_avg, locwvl, fit, fitvar, database, doas )

  USE OMSAO_precision_module
  USE OMSAO_indices_module, ONLY: &
       max_rs_idx, max_calfit_idx, solar_idx, ring_idx, ad1_idx, &
       lbe_idx, ad2_idx, mxs_idx, wvl_idx, spc_idx,                   &
       bl0_idx, bl1_idx, bl2_idx, bl3_idx, bl4_idx, bl5_idx, sc0_idx, &
       sc1_idx, sc2_idx, sc3_idx, sc4_idx, sc5_idx, sin_idx, hwe_idx, &
       asy_idx, shi_idx, squ_idx, bro_idx, o3_t1_idx, o3_t3_idx, io_idx
  USE OMSAO_parameters_module, ONLY: max_spec_pts, downweight
  USE OMSAO_variables_module,  ONLY: &
       n_database_wvl, curr_sol_spec, fitvar_rad, mask_fitvar_rad, fitweights, &
       yn_solar_comp, yn_spectrum_norm, yn_newshift
  USE OMSAO_prefitcol_module,  ONLY:                             &
       o3_prefit_fitidx,    yn_o3_prefit,    o3_prefit_var,      &
       bro_prefit_fitidx,   yn_bro_prefit,   bro_prefit_var,     &
       lqh2o_prefit_fitidx, yn_lqh2o_prefit, lqh2o_prefit_var
  USE OMSAO_omidata_module,      ONLY: curr_xtrack_pixnum, omi_solcal_pars, nwavel_max
  USE OMSAO_slitfunction_module, ONLY: saved_shift, saved_squeeze
  USE OMSAO_radiance_ref_module, ONLY: yn_radiance_reference, yn_reference_fit
  USE OMSAO_errstat_module
  USE OMSAO_solcomp_module
  USE EZspline_obj
  USE EZspline

  IMPLICIT NONE


  ! ===============
  ! Input variables
  ! ===============
  LOGICAL,                                                INTENT (IN) :: smooth, doas
  INTEGER (KIND=i4),                                      INTENT (IN) :: npts, n_fitvar
  REAL    (KIND=r8),                                      INTENT (IN) :: rad_wav_avg
  REAL    (KIND=r8), DIMENSION (n_fitvar),                INTENT (IN) :: fitvar
  REAL    (KIND=r8), DIMENSION (npts),                    INTENT (IN) :: locwvl
  REAL    (KIND=r8), DIMENSION (max_rs_idx,nwavel_max),   INTENT (IN) :: database

  ! ================
  ! Output variables
  ! ================
  REAL (KIND=r8), DIMENSION (npts), INTENT (OUT) :: fit

  ! ===============
  ! Local variables
  ! ===============
  REAL    (KIND=r8), PARAMETER                  :: expmax = REAL(MAXEXPONENT(1.0_r4), KIND=r8)
  REAL    (KIND=r8), PARAMETER                  :: expmin = REAL(MINEXPONENT(1.0_r4), KIND=r8)
  LOGICAL                                       :: yn_full_range, yn_solsynth
  INTEGER (KIND=i4)                             :: i, j, idx, errstat, j1, j2, n_sunpos
  REAL    (KIND=r8)                             :: shift, squeeze, soco_shi
  REAL    (KIND=r8), DIMENSION (npts)           :: del, sunspec_ss, tmpexp, sumexp
  REAL    (KIND=r8), DIMENSION (npts)           :: locwvl_shift
  REAL    (KIND=r8), DIMENSION (max_spec_pts)   :: sunpos_ss, sunspec_loc, sunspec_save

  ! ------------------------------
  ! Name of this subroutine/module
  ! ------------------------------
  CHARACTER (LEN=19), PARAMETER :: modulename = 'spectrum_earthshine'

  SAVE sunspec_save

  !     Calculate the spectrum:
  !     First do the shift and squeeze. Shift by FITVAR(SHI_IDX), squeeze by
  !     1 + FITVAR(SQU_IDX); do in absolute sense, to make it easy to back-convert
  !     OMI data.


  errstat = pge_errstat_ok

  ! ----------------------------------------------------------------------------
  ! Here is a logical to determine whether we need to compute a "sythetic"
  ! solar spectrum from the solar composite. The cases for YES are
  !
  ! (1) We are using the solar composite and are NOT doing a radiance reference
  ! (2) We are using the solar composite and ARE doing a radiance reference, and
  !     this happens to be the radiance reference fit.
  ! ----------------------------------------------------------------------------
  ! ---------------------------------------------------------------------
  ! The solar composite spectrum may have an additional shift, which was
  ! determined during the solar wavelength calibration. This needs to be
  ! taken into account when computing the spectra. But careful: It should
  ! NOT be added to the local wavelength array, since that is related to
  ! the radiance only. The Solar Composite shift must be subtracted from
  ! the wavelength array, hence the negative sign.
  ! ---------------------------------------------------------------------
  IF ( ( yn_solar_comp .AND. (.NOT. yn_radiance_reference) ) .OR. &
       ( yn_solar_comp .AND. (yn_radiance_reference .AND. yn_reference_fit) ) ) THEN
     yn_solsynth = .TRUE.
     soco_shi = -omi_solcal_pars(shi_idx,curr_xtrack_pixnum)
  ELSE
     yn_solsynth = .FALSE.
     soco_shi = 0.0_r8
  END IF

  ! -----------------------------------------------------------------------------------
  ! First, we have to undo the compression of the FITVAR_RAD array. This compression
  ! is performed in the RADIANCE_FIT subroutine and accelerates the fitting process,
  ! because ELSUNC has to handle less paraemters. But here we require the original
  ! layout, otherwise the index assingment is screwed.
  ! -----------------------------------------------------------------------------------
  DO i = 1, n_fitvar
     idx = mask_fitvar_rad(i)
     fitvar_rad(idx) = fitvar(i)
  END DO
  shift   = fitvar_rad(shi_idx)
  squeeze = fitvar_rad(squ_idx)

	! -------------------------------------
	! Dealing with any pre-fitted variables 
	! -------------------------------------
	! (1) OMCHOCHO
	! ------------
	IF ( yn_lqh2o_prefit(1) .AND. (.NOT. yn_lqh2o_prefit(2)) .AND. lqh2o_prefit_var /= 0.0_r8 ) &
       fitvar_rad(lqh2o_prefit_fitidx) = lqh2o_prefit_var
  ! ----------
  ! (2) OMHCHO
  !-----------
  IF ( yn_bro_prefit(1) .AND. (.NOT. yn_bro_prefit(2)) .AND. bro_prefit_var /= 0.0_r8 ) &
       fitvar_rad(bro_prefit_fitidx) = bro_prefit_var
  IF ( yn_o3_prefit(1)  .AND. (.NOT. yn_o3_prefit(2))  ) THEN
     DO j = o3_t1_idx, o3_t3_idx
        IF ( o3_prefit_var(j) /= 0.0_r8 ) fitvar_rad(o3_prefit_fitidx(j)) = o3_prefit_var(j)
     END DO
  END IF

  ! -----------------------------------------------------------------------------------------
  ! Assign current solar spectrum to local arrays. This depends on whether we are using
  ! actual measured solar spectra or solar composites. Since there is no point to interpolate
  ! already interpolated spectra, we use the original solar composites here as base for the
  ! interpolation to the final radiance wavelengths.
  ! -----------------------------------------------------------------------------------------
  n_sunpos                = n_database_wvl
  sunpos_ss  (1:n_sunpos) = curr_sol_spec(wvl_idx,1:n_sunpos)
  sunspec_loc(1:n_sunpos) = curr_sol_spec(spc_idx,1:n_sunpos)

  ! ----------------------------------------------
  ! Sort local arrays - important to pass EZspline
  ! ----------------------------------------------
  CALL array_sort_r8 ( n_sunpos, sunpos_ss(1:n_sunpos), sunspec_loc(1:n_sunpos) )

  ! ---------------------------------------------------------------------
  ! Apply Shift&Squeeze
  ! Changed to include Xiong comments (gga) if yn_newshift equal .true. :
  ! Lambda = Lambda * (1 + squeeze) + shift - sol_wav_avg * squeeze
  ! ---------------------------------------------------------------------
  j1 = -1; j2 = -1
  IF ( squeeze == 0.0_r8 .AND. yn_solsynth ) THEN
     locwvl_shift(1:npts) = locwvl(1:npts) - shift
     CALL array_locate_r8 ( npts, locwvl(1:npts), locwvl_shift(   1), 'GE', j1 )
     CALL array_locate_r8 ( npts, locwvl(1:npts), locwvl_shift(npts), 'LE', j2 )
  ELSE IF (yn_newshift) THEN !gga
     sunpos_ss(1:n_sunpos) = sunpos_ss(1:n_sunpos) * (1.0_r8 + squeeze) +       &
                             shift - rad_wav_avg * squeeze
     CALL array_locate_r8 ( npts, locwvl(1:npts), sunpos_ss(       1), 'GE', j1 )
     CALL array_locate_r8 ( npts, locwvl(1:npts), sunpos_ss(n_sunpos), 'LE', j2 ) !gga
  ELSE
     sunpos_ss(1:n_sunpos) = sunpos_ss(1:n_sunpos) * (1.0_r8 + squeeze) + shift
     CALL array_locate_r8 ( npts, locwvl(1:npts), sunpos_ss(       1), 'GE', j1 )
     CALL array_locate_r8 ( npts, locwvl(1:npts), sunpos_ss(n_sunpos), 'LE', j2 )
  END IF

  ! ---------------------------------------------------------------------
  ! Re-sample the solar reference spectrum to the current radiance grid
  ! ---------------------------------------------------------------------
  ! 
  ! The endpoints may be problematic due to no-strict ascendence. If that
  ! happens, exclude end-points.
  ! ---------------------------------------------------------------------

  IF ( j1 <= 0 .OR. j2 <= 0 ) THEN
     CALL error_check ( &
       0, 1, pge_errstat_warning, OMSAO_W_INTERPOL_RANGE, &
       modulename//f_sep//'Resampling to Radiance Grid -- no solar spectrum!!!', &
       vb_lev_default, errstat )
  ELSE

     IF ( squeeze /= saved_squeeze .OR. shift /= saved_shift ) THEN

        IF ( squeeze == 0.0_r8 .AND. yn_solsynth ) THEN
           CALL soco_compute ( &
                yn_spectrum_norm, curr_xtrack_pixnum, npts, &
                locwvl_shift(1:npts)+soco_shi, sunspec_ss(1:npts) )
        ELSE
           CALL interpolation (                                                 &
                modulename,                                                     &
                n_sunpos, sunpos_ss(1:n_sunpos), sunspec_loc(1:n_sunpos),       &
                npts, locwvl(1:npts), sunspec_ss(1:npts), 'endpoints', 0.0_r8,  &
                yn_full_range, errstat                                            )
           CALL error_check ( &
                errstat, pge_errstat_ok, pge_errstat_error, OMSAO_E_INTERPOL,      &
                modulename//f_sep//'Resampling to Radiance Grid -- interpolation', &
                vb_lev_default, errstat )
        END IF
        sunspec_save(1:npts) = sunspec_ss(1:npts)
        saved_shift          = shift
        saved_squeeze        = squeeze
     ELSE
        sunspec_ss(1:npts) = sunspec_save(1:npts)
     END IF

  END IF

  ! Add up the contributions, with solar intensity as FITVAR_RAD(sin_idx), trace
  ! species beginning at FITVAR_RAD(SQU_IDX+1), to include possible linear and
  ! Beer's law forms.  Do these as linear-Beer's-linear. In order to
  ! do DOAS I only need to be careful to include just linear
  ! contributions, since I already high-pass filtered them.

  IF ( j1 > 1    )  fitweights(1:j1-1)    = downweight
  IF ( j2 < npts )  fitweights(j2+1:npts) = downweight

  fit = 0.0_r8
  
  ! ==================================================================
  ! For BOAS or any wavelength calibration, we have the following line
  ! ==================================================================


  fit(j1:j2) = fitvar_rad(sin_idx) * sunspec_ss(j1:j2)

  !     DOAS here - the spectrum to be fitted needs to be re-defined:
  IF ( doas ) THEN

     i = max_calfit_idx + (ring_idx-1)*mxs_idx + ad1_idx

     fit(j1:j2) = &
          ! For DOAS, FITVAR_RAD(SIN_IDX) should == 1., and not be varied
          fitvar_rad(sin_idx) * LOG ( sunspec_ss(j1:j2) ) + &
          ! Ring adjustment
          fitvar_rad(i) * (database(ring_idx,j1:j2) / sunspec_ss (j1:j2))

     DO j = 1, max_rs_idx
        IF ( j /= solar_idx .AND. j /= ring_idx ) THEN
           i = max_calfit_idx + (j-1)*mxs_idx + ad1_idx
           fit(j1:j2) = fit(j1:j2) + fitvar_rad(i) * database(j,j1:j2)
        END IF
     END DO

  ELSE
     ! -----------------------------
     ! Initial add-on contributions.
     ! -----------------------------
     DO j = 1, max_rs_idx
        IF ( j /= solar_idx ) THEN
           i = max_calfit_idx + (j-1)*mxs_idx + ad1_idx
           fit(j1:j2) = fit(j1:j2) + fitvar_rad(i) * database(j,j1:j2)
        END IF
     END DO
     ! -----------------------------
     ! Beer's law contributions.
     ! -----------------------------
     ! ---------------------------------------------------------------
     ! We sum over all contributions and take the EXP only at the end.
     ! This should shave a few seconds off the execution time.
     ! ---------------------------------------------------------------
     sumexp(j1:j2) = 0.0_r8
     DO j = 1, max_rs_idx
        IF ( j /= solar_idx ) THEN
           i = max_calfit_idx + (j-1)*mxs_idx + lbe_idx
           tmpexp = 0.0_r8  ;  tmpexp(j1:j2) = fitvar_rad(i)*database(j,j1:j2)
           WHERE ( tmpexp(j1:j2) >= expmax )
              tmpexp(j1:j2) = expmax
           ENDWHERE
           WHERE ( tmpexp(j1:j2) <= expmin )
              tmpexp(j1:j2) = expmin
           ENDWHERE
           sumexp(j1:j2) = sumexp(j1:j2) - tmpexp(j1:j2)

        END IF
     END DO
     WHERE ( sumexp(j1:j2) >= expmax )
        sumexp(j1:j2) = expmax
     ENDWHERE
     WHERE ( sumexp(j1:j2) <= expmin )
        sumexp(j1:j2) = expmin
     ENDWHERE
     fit(j1:j2) = fit(j1:j2) * EXP(sumexp(j1:j2))

     ! Final add-on contributions.
     DO j = 1, max_rs_idx
        IF ( j /= solar_idx ) THEN !.AND. database(j,j1:j2) /= 0.0_r8 ) THEN
           i = max_calfit_idx + (j-1)*mxs_idx + ad2_idx
           fit(j1:j2) = fit(j1:j2) + fitvar_rad(i) * database(j,j1:j2)
        END IF
     END DO

  END IF


  ! Add the scaling.
  del(j1:j2) = locwvl(j1:j2) - rad_wav_avg
  fit(j1:j2) = fit(j1:j2) * ( &
       fitvar_rad(sc0_idx)                                       + &
       fitvar_rad(sc1_idx) * del(j1:j2)                          + &
       fitvar_rad(sc2_idx) * del(j1:j2)*del(j1:j2)               + &
       fitvar_rad(sc3_idx) * del(j1:j2)*del(j1:j2)*del(j1:j2)    + &
       fitvar_rad(sc4_idx) * del(j1:j2)**4.0                     + &
       fitvar_rad(sc5_idx) * del(j1:j2)**5.0                     )
  ! Add baseline parameters.
  fit(j1:j2) = fit(j1:j2)                                        + &
       fitvar_rad(bl0_idx)                                       + &
       fitvar_rad(bl1_idx) * del(j1:j2)                          + &
       fitvar_rad(bl2_idx) * del(j1:j2)*del(j1:j2)               + &
       fitvar_rad(bl3_idx) * del(j1:j2)*del(j1:j2)*del(j1:j2)    + &
       fitvar_rad(bl4_idx) * del(j1:j2)**4.0                     + &
       fitvar_rad(bl5_idx) * del(j1:j2)**5.0                     

  ! ----------------------------------------------------------------
  ! Final sanity check: If the various multiplications and additions
  ! have lead to NaN values, we set those to ZERO. This is somewhat
  ! experimental, and if we come up with a better way of doing this,
  ! then the logic below should be changed accordingly.
  ! ----------------------------------------------------------------
  WHERE ( .NOT. ( fit(j1:j2) > -HUGE(1.0_r8) .AND. fit(j1:j2) < HUGE(1.0_r8) ) )
     fit(j1:j2) = 0.0_r8
  END WHERE

  RETURN
END SUBROUTINE spectrum_earthshine

SUBROUTINE spectrum_earthshine_o3exp ( &
     npts, n_fitvar, smooth, rad_wav_avg, locwvl, fit, fitvar, database, doas )

  USE OMSAO_precision_module
  USE OMSAO_indices_module, ONLY: &
       max_rs_idx, max_calfit_idx, solar_idx, ring_idx, ad1_idx, &
       lbe_idx, ad2_idx, mxs_idx, wvl_idx, spc_idx,                   &
       bl0_idx, bl1_idx, bl2_idx, bl3_idx, bl4_idx, bl5_idx, sc0_idx, &
       sc1_idx, sc2_idx, sc3_idx, sc4_idx, sc5_idx, sin_idx, hwe_idx, &
       asy_idx, shi_idx, squ_idx, bro_idx, o3_t1_idx, o3_t2_idx, o3_t3_idx
  USE OMSAO_parameters_module, ONLY: max_spec_pts, downweight
  USE OMSAO_variables_module,  ONLY: &
       n_database_wvl, curr_sol_spec, fitvar_rad, mask_fitvar_rad, fitweights, &
       yn_solar_comp, yn_spectrum_norm, yn_newshift
  USE OMSAO_prefitcol_module,  ONLY:                                           &
       bro_prefit_fitidx, o3_prefit_fitidx, yn_bro_prefit, bro_prefit_var,     &
       yn_o3_prefit, o3_prefit_var
  USE OMSAO_omidata_module,      ONLY: curr_xtrack_pixnum, omi_solcal_pars
  USE OMSAO_slitfunction_module, ONLY: saved_shift, saved_squeeze
  USE OMSAO_radiance_ref_module, ONLY: yn_radiance_reference, yn_reference_fit
  USE OMSAO_errstat_module
  USE OMSAO_solcomp_module
  USE EZspline_obj
  USE EZspline

  IMPLICIT NONE


  ! ===============
  ! Input variables
  ! ===============
  LOGICAL,                                                INTENT (IN) :: smooth, doas
  INTEGER (KIND=i4),                                      INTENT (IN) :: npts, n_fitvar
  REAL    (KIND=r8),                                      INTENT (IN) :: rad_wav_avg
  REAL    (KIND=r8), DIMENSION (n_fitvar),                INTENT (IN) :: fitvar
  REAL    (KIND=r8), DIMENSION (npts),                    INTENT (IN) :: locwvl
  REAL    (KIND=r8), DIMENSION (max_rs_idx,max_spec_pts), INTENT (IN) :: database

  ! ================
  ! Output variables
  ! ================
  REAL (KIND=r8), DIMENSION (npts), INTENT (OUT) :: fit

  ! ===============
  ! Local variables
  ! ===============
  REAL    (KIND=r8), PARAMETER                  :: expmax = REAL(MAXEXPONENT(1.0_r4), KIND=r8)
  REAL    (KIND=r8), PARAMETER                  :: expmin = REAL(MINEXPONENT(1.0_r4), KIND=r8)
  LOGICAL                                       :: yn_full_range, yn_solsynth
  INTEGER (KIND=i4)                             :: i, j, idx, errstat, j1, j2, n_sunpos, k1, k2
  REAL    (KIND=r8)                             :: shift, squeeze, soco_shi
  REAL    (KIND=r8), DIMENSION (npts)           :: del, sunspec_ss, tmpexp, sumexp
  REAL    (KIND=r8), DIMENSION (npts)           :: locwvl_shift
  REAL    (KIND=r8), DIMENSION (max_spec_pts)   :: sunpos_ss, sunspec_loc, sunspec_save

  ! ------------------------------
  ! Name of this subroutine/module
  ! ------------------------------
  CHARACTER (LEN=19), PARAMETER :: modulename = 'spectrum_earthshine'

  SAVE sunspec_save

  !  Calculate the spectrum:
  !  First do the shift and squeeze. Shift by FITVAR(SHI_IDX), squeeze by
  !  1 + FITVAR(SQU_IDX); do in absolute sense, to make it easy to back-convert
  !  OMI data.


  errstat = pge_errstat_ok

  ! ----------------------------------------------------------------------------
  ! Here is a logical to determine whether we need to compute a "sythetic"
  ! solar spectrum from the solar composite. The cases for YES are
  !
  ! (1) We are using the solar composite and are NOT doing a radiance reference
  ! (2) We are using the solar composite and ARE doing a radiance reference, and
  !     this happens to be the radiance reference fit.
  ! ----------------------------------------------------------------------------
  ! ---------------------------------------------------------------------
  ! The solar composite spectrum may have an additional shift, which was
  ! determined during the solar wavelength calibration. This needs to be
  ! taken into account when computing the spectra. But careful: It should
  ! NOT be added to the local wavelength array, since that is related to
  ! the radiance only. The Solar Composite shift must be subtracted from
  ! the wavelength array, hence the negative sign.
  ! ---------------------------------------------------------------------
  IF ( ( yn_solar_comp .AND. (.NOT. yn_radiance_reference) ) .OR. &
       ( yn_solar_comp .AND. (yn_radiance_reference .AND. yn_reference_fit) ) ) THEN
     yn_solsynth = .TRUE.
     soco_shi = -omi_solcal_pars(shi_idx,curr_xtrack_pixnum)
  ELSE
     yn_solsynth = .FALSE.
     soco_shi = 0.0_r8
  END IF

  ! -----------------------------------------------------------------------------------
  ! First, we have to undo the compression of the FITVAR_RAD array. This compression
  ! is performed in the RADIANCE_FIT subroutine and accelerates the fitting process,
  ! because ELSUNC has to handle less paraemters. But here we require the original
  ! layout, otherwise the index assingment is screwed.
  ! -----------------------------------------------------------------------------------
  DO i = 1, n_fitvar
     idx = mask_fitvar_rad(i)
     fitvar_rad(idx) = fitvar(i)
  END DO
  shift   = fitvar_rad(shi_idx)
  squeeze = fitvar_rad(squ_idx)

  IF ( yn_bro_prefit(1) .AND. (.NOT. yn_bro_prefit(2)) .AND. bro_prefit_var /= 0.0_r8 ) &
       fitvar_rad(bro_prefit_fitidx) = bro_prefit_var
  IF ( yn_o3_prefit(1)  .AND. (.NOT. yn_o3_prefit(2))  ) THEN
     DO j = o3_t1_idx, o3_t3_idx
        IF ( o3_prefit_var(j) /= 0.0_r8 ) fitvar_rad(o3_prefit_fitidx(j)) = o3_prefit_var(j)
     END DO
  END IF

  ! ---------------------------------------------------------------------------------
  ! Assign current solar spectrum to local arrays. This depends on whether we are
  ! using actual measured solar spectra or solar composites. Since there is no point
  ! to interpolate already interpolated spectra, we use the original solar composites
  ! here as base for the interpolation to the final radiance wavelengths.
  ! ---------------------------------------------------------------------------------
  n_sunpos                = n_database_wvl
  sunpos_ss  (1:n_sunpos) = curr_sol_spec(wvl_idx,1:n_sunpos)
  sunspec_loc(1:n_sunpos) = curr_sol_spec(spc_idx,1:n_sunpos)


  ! ----------------------------------------------
  ! Sort local arrays - important to pass EZspline
  ! ----------------------------------------------
  CALL array_sort_r8 ( n_sunpos, sunpos_ss(1:n_sunpos), sunspec_loc(1:n_sunpos) )

  ! ---------------------------------------------------------------------
  ! Apply Shift&Squeeze
  ! Changed to include Xiong comments (gga) if yn_newshift equal .true. :
  ! Lambda = Lambda * (1 + squeeze) + shift - sol_wav_avg * squeeze
  ! ---------------------------------------------------------------------
  j1 = -1; j2 = -1
  IF ( squeeze == 0.0_r8 .AND. yn_solsynth ) THEN
     locwvl_shift(1:npts) = locwvl(1:npts) - shift
     CALL array_locate_r8 ( npts, locwvl(1:npts), locwvl_shift(   1), 'GE', j1 )
     CALL array_locate_r8 ( npts, locwvl(1:npts), locwvl_shift(npts), 'LE', j2 )
  ELSE IF (yn_newshift .EQ. .true.) THEN !gga
     sunpos_ss(1:n_sunpos) = sunpos_ss(1:n_sunpos) * (1.0_r8 + squeeze) +       &
                             shift - rad_wav_avg * squeeze
     CALL array_locate_r8 ( npts, locwvl(1:npts), sunpos_ss(       1), 'GE', j1 )
     CALL array_locate_r8 ( npts, locwvl(1:npts), sunpos_ss(n_sunpos), 'LE', j2 ) !gga
  ELSE
     sunpos_ss(1:n_sunpos) = sunpos_ss(1:n_sunpos) * (1.0_r8 + squeeze) + shift
     CALL array_locate_r8 ( npts, locwvl(1:npts), sunpos_ss(       1), 'GE', j1 )
     CALL array_locate_r8 ( npts, locwvl(1:npts), sunpos_ss(n_sunpos), 'LE', j2 )
  END IF


  ! ---------------------------------------------------------------------
  ! Re-sample the solar reference spectrum to the current radiance grid
  ! ---------------------------------------------------------------------
  ! 
  ! The endpoints may be problematic due to no-strict ascendence. If that
  ! happens, exclude end-points.
  ! ---------------------------------------------------------------------

  IF ( j1 <= 0 .OR. j2 <= 0 ) THEN
     CALL error_check ( &
       0, 1, pge_errstat_warning, OMSAO_W_INTERPOL_RANGE, &
       modulename//f_sep//'Resampling to Radiance Grid -- no solar spectrum!!!', &
       vb_lev_default, errstat )
  ELSE

     IF ( squeeze /= saved_squeeze .OR. shift /= saved_shift ) THEN

        IF ( squeeze == 0.0_r8 .AND. yn_solsynth ) THEN
           CALL soco_compute ( &
                yn_spectrum_norm, curr_xtrack_pixnum, npts, &
                locwvl_shift(1:npts)+soco_shi, sunspec_ss(1:npts) )
        ELSE
           CALL interpolation (                                                         &
                modulename,                                                             &
                n_sunpos, sunpos_ss(1:n_sunpos), sunspec_loc(1:n_sunpos),               &
                npts, locwvl(1:npts), sunspec_ss(1:npts), 'endpoints', 0.0_r8, &
                yn_full_range, errstat                                                   )
           CALL error_check ( &
                errstat, pge_errstat_ok, pge_errstat_error, OMSAO_E_INTERPOL, &
                modulename//f_sep//'Resampling to Radiance Grid -- interpolation', &
                vb_lev_default, errstat )
        END IF

        sunspec_save(1:npts) = sunspec_ss(1:npts)
        saved_shift   = shift
        saved_squeeze = squeeze
     ELSE
        sunspec_ss(1:npts) = sunspec_save(1:npts)
     END IF

  END IF

  !     Add up the contributions, with solar intensity as FITVAR_RAD(sin_idx), trace
  !     species beginning at FITVAR_RAD(SQU_IDX+1), to include possible linear and
  !     Beer's law forms.  Do these as linear-Beer's-linear. In order to
  !     do DOAS I only need to be careful to include just linear
  !     contributions, since I already high-pass filtered them.

  IF ( j1 > 1    )  fitweights(1:j1-1)       = downweight
  IF ( j2 < npts )  fitweights(j2+1:npts) = downweight

  fit = 0.0_r8

  ! --------------------------------------------------------
  ! Compute abcissae for exponential x-section modification:
  ! Values between -1 and +1 on the fitting wavelength grid.
  ! --------------------------------------------------------
  del(j1:j2) = (locwvl(j1:j2) - locwvl(j1))/(locwvl(j2)-locwvl(j1)) - 0.5_r8

  ! ==================================================================
  ! For BOAS or any wavelength calibration, we have the following line
  ! ==================================================================


  fit(j1:j2) = fitvar_rad(sin_idx) * sunspec_ss(j1:j2)

  !     DOAS here - the spectrum to be fitted needs to be re-defined:
  IF ( doas ) THEN

     i = max_calfit_idx + (ring_idx-1)*mxs_idx + ad1_idx

     fit(j1:j2) = &
          ! For DOAS, FITVAR_RAD(SIN_IDX) should == 1., and not be varied
          fitvar_rad(sin_idx) * LOG ( sunspec_ss(j1:j2) ) + &
          ! Ring adjustment
          fitvar_rad(i) * (database(ring_idx,j1:j2) / sunspec_ss (j1:j2))

     DO j = 1, max_rs_idx
        IF ( j /= solar_idx .AND. j /= ring_idx  .AND. &
             j /= o3_t1_idx .AND. j /= o3_t2_idx .AND. j /= o3_t3_idx ) THEN
           i = max_calfit_idx + (j-1)*mxs_idx + ad1_idx
           fit(j1:j2) = fit(j1:j2) + fitvar_rad(i) * database(j,j1:j2)
        END IF
     END DO

  ELSE
     ! -----------------------------
     ! Initial add-on contributions.
     ! -----------------------------
     DO j = 1, max_rs_idx
        IF ( j /= solar_idx .AND. &
             j /= o3_t1_idx .AND. j /= o3_t2_idx .AND. j /= o3_t3_idx ) THEN
           i = max_calfit_idx + (j-1)*mxs_idx + ad1_idx
           fit(j1:j2) = fit(j1:j2) + fitvar_rad(i) * database(j,j1:j2)
        END IF
     END DO
     ! -----------------------------
     ! Beer's law contributions.
     ! -----------------------------
     ! ---------------------------------------------------------------
     ! We sum over all contributions and take the EXP only at the end.
     ! This should shave a few seconds off the execution time.
     ! ---------------------------------------------------------------
     sumexp(j1:j2) = 0.0_r8
     DO j = 1, max_rs_idx
        IF ( j /= solar_idx ) THEN
           tmpexp = 0.0_r8
           i = max_calfit_idx + (j-1)*mxs_idx + lbe_idx
           IF ( j == o3_t1_idx .OR. j == o3_t2_idx .OR. j == o3_t3_idx ) THEN
              k1 = max_calfit_idx + (j-1)*mxs_idx + ad1_idx
              k2 = max_calfit_idx + (j-1)*mxs_idx + ad2_idx
              tmpexp(j1:j2) = fitvar_rad(i)*database(j,j1:j2) *  &
                   (1.0_r8 + fitvar_rad(k1)*del(j1:j2) + &
                   fitvar_rad(k2)*del(j1:j2)*del(j1:j2))                   
           ELSE
              tmpexp(j1:j2) = fitvar_rad(i)*database(j,j1:j2)
           END IF

           WHERE ( tmpexp(j1:j2) >= expmax )
              tmpexp(j1:j2) = expmax
           ENDWHERE
           WHERE ( tmpexp(j1:j2) <= expmin )
              tmpexp(j1:j2) = expmin
           ENDWHERE
           sumexp(j1:j2) = sumexp(j1:j2) - tmpexp(j1:j2)
        END IF
     END DO
     WHERE ( sumexp(j1:j2) >= expmax )
        sumexp(j1:j2) = expmax
     ENDWHERE
     WHERE ( sumexp(j1:j2) <= expmin )
        sumexp(j1:j2) = expmin
     ENDWHERE
     fit(j1:j2) = fit(j1:j2) * EXP(sumexp(j1:j2))

     ! Final add-on contributions.
     DO j = 1, max_rs_idx
        IF ( j /= solar_idx .AND. &
             j /= o3_t1_idx .AND. j /= o3_t2_idx .AND. j /= o3_t3_idx ) THEN
           i = max_calfit_idx + (j-1)*mxs_idx + ad2_idx
           fit(j1:j2) = fit(j1:j2) + fitvar_rad(i) * database(j,j1:j2)
        END IF
     END DO

  END IF


  ! ----------------------------------------
  ! Compute abcissae for closure polynomials
  ! ----------------------------------------
  del(j1:j2) = locwvl(j1:j2) - rad_wav_avg

!   ! Add the scaling.
  del(j1:j2) = locwvl(j1:j2) - rad_wav_avg
  fit(j1:j2) = fit(j1:j2) * ( &
       fitvar_rad(sc0_idx)                                       + &
       fitvar_rad(sc1_idx) * del(j1:j2)                          + &
       fitvar_rad(sc2_idx) * del(j1:j2)*del(j1:j2)               + &
       fitvar_rad(sc3_idx) * del(j1:j2)*del(j1:j2)*del(j1:j2)    + &
       fitvar_rad(sc4_idx) * del(j1:j2)**4.0                     + &
       fitvar_rad(sc5_idx) * del(j1:j2)**5.0                     )
  ! Add baseline parameters.
  fit(j1:j2) = fit(j1:j2)                                        + &
       fitvar_rad(bl0_idx)                                       + &
       fitvar_rad(bl1_idx) * del(j1:j2)                          + &
       fitvar_rad(bl2_idx) * del(j1:j2)*del(j1:j2)               + &
       fitvar_rad(bl3_idx) * del(j1:j2)*del(j1:j2)*del(j1:j2)    + &
       fitvar_rad(bl4_idx) * del(j1:j2)**4.0                     + &
       fitvar_rad(bl5_idx) * del(j1:j2)**5.0                     

  ! ----------------------------------------------------------------
  ! Final sanity check: If the various multiplications and additions
  ! have lead to NaN values, we set those to ZERO. This is somewhat
  ! experimental, and if we come up with a better way of doing this,
  ! then the logic below should be changed accordingly.
  ! ----------------------------------------------------------------
  WHERE ( .NOT. ( fit(j1:j2) > -HUGE(1.0_r8) .AND. fit(j1:j2) < HUGE(1.0_r8) ) )
     fit(j1:j2) = 0.0_r8
  END WHERE

  RETURN
END SUBROUTINE spectrum_earthshine_o3exp
