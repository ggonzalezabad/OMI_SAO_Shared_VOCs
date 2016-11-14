MODULE OMSAO_radiance_ref_module

  USE OMSAO_precision_module
  USE OMSAO_parameters_module,   ONLY: maxchlen
  USE OMSAO_variables_module, ONLY: n_rad_wvl_max

  IMPLICIT NONE

  ! -------------------------------------------------------
  ! Variables connected with  a radiance reference spectrum
  ! -------------------------------------------------------
  LOGICAL                          :: yn_radiance_reference, yn_reference_fit, yn_remove_target
  INTEGER (KIND=i4)                :: target_npol
  INTEGER (KIND=i4), DIMENSION (2) :: radiance_reference_lnums
  REAL    (KIND=r4), DIMENSION (2) :: radref_latrange

  CHARACTER (LEN=maxchlen) :: l1b_radref_filename

CONTAINS


  SUBROUTINE omi_get_radiance_reference ( l1bfile, radwcal_lines, errstat )

    USE OMSAO_parameters_module, ONLY: &
         r4_missval, r8_missval, downweight, normweight, deg2rad, rad2deg
    USE OMSAO_indices_module,    ONLY: &
       qflg_mis_idx, qflg_bad_idx, qflg_err_idx, qflg_tra_idx, qflg_rts_idx, qflg_sat_idx
    USE OMSAO_variables_module,  ONLY: pixnum_lim, fit_winwav_lim, fit_winexc_lim
    USE OMSAO_omidata_module,    ONLY: &
         omi_ccdpix_selection, omi_radiance_qflg, omi_radiance_spec, omi_radiance_wavl, &
         omi_szenith, omi_vzenith, omi_nwav_radref, omi_radref_spec, omi_radref_wavl,   &
         omi_radref_qflg, omi_radref_sza, omi_radref_vza, omi_radref_wght,              &
         omi_ccdpix_exclusion, nlines_max, n_comm_wvl, omi_sol_wav_avg, &
         omi_nwav_rad
    USE OMSAO_errstat_module
    USE L1B_Reader_class

    IMPLICIT NONE

    ! ------------------------------
    ! Name of this module/subroutine
    ! ------------------------------
    CHARACTER (LEN=28), PARAMETER :: modulename = 'omi_get_radiance_reference' ! JED fixed

    ! ---------------
    ! Input variables
    ! ---------------
    CHARACTER (LEN=*), INTENT (IN) :: l1bfile

    ! -----------------------------
    ! Output and Modified variables
    ! -----------------------------
    INTEGER (KIND=i4), DIMENSION (2), INTENT (OUT)   :: radwcal_lines
    INTEGER (KIND=i4),                INTENT (INOUT) :: errstat

    ! ---------------
    ! Local variables
    ! ---------------
    INTEGER (KIND=i2), PARAMETER :: nbits = 16
    LOGICAL                      :: yn_have_scanline, yn_dummy
    LOGICAL, DIMENSION (2)       :: yn_have_limits
    CHARACTER (LEN=maxchlen)     :: l1bswath

    INTEGER (KIND=i4) :: ntrr, nxrr, nwrr
    INTEGER (KIND=i4) :: fpix, lpix, midpt_line
    INTEGER (KIND=i4) :: nloop, j1, iline, ix, iloop, imin, imax, icnt
    REAL    (KIND=r4) :: lat_midpt
    REAL    (KIND=r8) :: specsum

!!$    INTEGER (KIND=i4), DIMENSION (0:ntrr-1,2)     :: xtrange
!!$    INTEGER (KIND=i1), DIMENSION (0:ntrr-1)       :: binfac
!!$    LOGICAL,           DIMENSION (0:ntrr-1)       :: ynzoom
!!$
!!$    REAL    (KIND=r4), DIMENSION (nxrr)           :: szacount
!!$    REAL    (KIND=r4), DIMENSION (nxrr,0:ntrr-1)  :: latr4
!!$    REAL    (KIND=r8), DIMENSION (nxrr, nwrr)     :: radref_spec, radref_wavl
!!$    REAL    (KIND=r8), DIMENSION (nwrr)           :: radref_wavl_ix
!!$    REAL    (KIND=r8), DIMENSION (nxrr, nwrr)     :: allcount, dumcount
!!$    REAL    (KIND=r8), DIMENSION (nwrr      )     :: cntr8
!!$    INTEGER (KIND=i2), DIMENSION (nwrr,0:nbits-1) :: qflg_bit
!!$    INTEGER (KIND=i2), DIMENSION (nwrr)           :: qflg_mask

    INTEGER (KIND=i4), ALLOCATABLE, DIMENSION (:,:) :: xtrange
    INTEGER (KIND=i1), ALLOCATABLE, DIMENSION (:)   :: binfac
    LOGICAL,           ALLOCATABLE, DIMENSION (:)   :: ynzoom

    REAL    (KIND=r4), ALLOCATABLE, DIMENSION (:)   :: szacount
    REAL    (KIND=r4), ALLOCATABLE, DIMENSION (:,:) :: latr4
    REAL    (KIND=r8), ALLOCATABLE, DIMENSION (:,:) :: radref_spec, radref_wavl
    REAL    (KIND=r8), ALLOCATABLE, DIMENSION (:)   :: radref_wavl_ix
    REAL    (KIND=r8), ALLOCATABLE, DIMENSION (:,:) :: allcount, dumcount
    REAL    (KIND=r8), ALLOCATABLE, DIMENSION (:)   :: cntr8
    INTEGER (KIND=i2), ALLOCATABLE, DIMENSION (:,:) :: qflg_bit
    INTEGER (KIND=i2), ALLOCATABLE, DIMENSION (:)   :: qflg_mask


    ! ------------------------------
    ! Initialize some some variables
    ! ------------------------------
    radiance_reference_lnums = -1  ! This will be written to file, hence needs a value
    lat_midpt = SUM ( radref_latrange ) / 2.0_r4

    ! -------------------------
    ! Retrieve swath dimensions
    ! -------------------------
    CALL omi_read_radiance_paras ( &
         l1bfile, ntrr, nxrr, nwrr, l1bswath, errstat )

    ! ----------------------------------------------------------
    ! Allocate variables after reading radiance swath dimensions
    ! ----------------------------------------------------------
    ALLOCATE (xtrange(0:ntrr-1,2))      ; ALLOCATE (binfac(0:ntrr-1))
    ALLOCATE (ynzoom(0:ntrr-1))         ; ALLOCATE (szacount(nxrr))
    ALLOCATE (latr4(nxrr,0:ntrr-1))
    ALLOCATE (radref_spec(nxrr,nwrr))   ; ALLOCATE (radref_wavl(nxrr,nwrr))
    ALLOCATE (radref_wavl_ix(nwrr))     ; ALLOCATE (allcount(nxrr,nwrr))
    ALLOCATE (dumcount(nxrr,nwrr))      ; ALLOCATE (cntr8(nwrr))
    ALLOCATE (qflg_bit(nwrr,0:nbits-1)) ; ALLOCATE (qflg_mask(nwrr))

    CALL omi_read_binning_factor ( &
         l1bfile, l1bswath, ntrr, binfac(0:ntrr-1), ynzoom(0:ntrr-1), errstat )

    CALL omi_set_xtrpix_range ( &
         ntrr, nxrr, pixnum_lim(3:4), binfac(0:ntrr-1), &
         xtrange(0:ntrr-1,1:2), fpix, lpix, errstat    )

    CALL read_latitude (l1bfile, l1bswath, ntrr, nxrr, latr4(1:nxrr,0:ntrr-1) )

    ! ----------------------------------------------------------------------
    ! Locate the swath line numbers corresponding the center of the latitude
    ! range to average into radiance reference spectrum.
    ! ----------------------------------------------------------------------
    CALL find_swathline_by_latitude ( &
         nxrr, 0, ntrr-1, latr4(1:nxrr,0:ntrr-1), lat_midpt, &
         xtrange(0:ntrr-1,1:2), midpt_line, yn_have_scanline )

    ! --------------------------------------------------------------------
    ! If lower and upper bounds of the radiance reference block to average
    ! are identical, then we keep the midpoint line number as the only
    ! reference. Else locate the corresponding swath line numbers.
    ! --------------------------------------------------------------------
    IF ( radref_latrange(1) == radref_latrange(2) ) THEN
       radiance_reference_lnums(1:2) = midpt_line
       yn_have_limits(1:2)           = .TRUE.
    ELSE
       CALL find_swathline_by_latitude ( &
            nxrr, 0, midpt_line, latr4(1:nxrr,0:midpt_line), radref_latrange(1), &
            xtrange(0:midpt_line,1:2), radiance_reference_lnums(1), yn_have_limits(1)   )
       CALL find_swathline_by_latitude ( &
            nxrr, midpt_line, ntrr-1, latr4(1:nxrr,midpt_line:ntrr-1), radref_latrange(2), &
            xtrange(midpt_line:ntrr-1,1:2), radiance_reference_lnums(2), yn_have_limits(2) )
    END IF

    print*, midpt_line, yn_have_scanline, lat_midpt, radref_latrange
    print*, radiance_reference_lnums(1), radiance_reference_lnums(2)
    stop

    ! -----------------------------------------------------
    ! If we don't find a working scan line, we have to fold
    ! -----------------------------------------------------
    IF ( ( .NOT. yn_have_scanline )               .OR. &
         ( ANY ( .NOT. yn_have_limits(1:2) ) )    .OR. &
         ( midpt_line < 0 )                       .OR. &
         ( ANY ( radiance_reference_lnums < 0 ) )        ) THEN
       CALL error_check ( 1, 0, pge_errstat_fatal, OMSAO_E_READ_L1B_FILE, &
            modulename//f_sep//"Failed to find working radiance spectrum.", vb_lev_default, errstat )
       RETURN
    END IF

    ! ------------------------------------------------------------------
    ! A kludge for now, or maybe not. Assign the radiance reference
    ! swath line numbers to the radiance wavelength calibration swath
    ! line numbers.
    ! ------------------------------------------------------------------
    radwcal_lines(1:2) = radiance_reference_lnums(1:2)

    ! ---------------------------------------------------------------------------
    ! Fudge the selection of the complete spectrum; this is required because the
    ! routine that does the radiance line reading expects it to be available from
    ! after the Irradiance reading (which we are skipping). We will be setting up
    ! the proper numbers below.
    ! ---------------------------------------------------------------------------
    omi_ccdpix_selection(1:nxrr,1) = 1 ; omi_ccdpix_selection(1:nxrr,4) = nwrr

    ! --------------------------------------------------------------------
    ! Now we can average the spectra and the wavelength arrays. Loop over
    ! the block of swath lines in multiples of NLINES_MAX (100 by default)
    ! --------------------------------------------------------------------
    allcount    = 0.0_r8  ;  dumcount    = 0.0_r8  ;  szacount = 0.0_r4
    radref_wavl = 0.0_r8  ;  radref_spec = 0.0_r8
    omi_radref_sza = 0.0_r4 ; omi_radref_vza = 0.0_r4

    DO iline = radiance_reference_lnums(1), radiance_reference_lnums(2), nlines_max

       ! --------------------------------------------------------
       ! Check if loop ends before n_times_loop max is exhausted, 
       ! or if we are outside the FIRST_LINE -> LAST_LINE range.
       ! --------------------------------------------------------
       nloop = MIN( nlines_max, radiance_reference_lnums(2)-radiance_reference_lnums(1)+1 )

       IF ( (iline+nloop) > ntrr  )  nloop = ntrr - iline

       ! ------------------------------
       ! Get NTIMES_LOOP radiance lines
       ! ------------------------------
       CALL omi_read_radiance_lines (              &
            l1bfile, iline, nxrr, nloop, nwrr, errstat )

       ! Global used to set the dimension of fitspc
       n_rad_wvl_max = MAXVAL(omi_nwav_rad(:,0))

       DO iloop = 0, nloop-1

          ! ------------------------------------------------------
          ! Skip this cross-track position if there isn't any data
          ! ------------------------------------------------------
          fpix = xtrange(iline+iloop,1)
          lpix = xtrange(iline+iloop,2)

          DO ix = 1, nxrr

             IF ( (ix < fpix) .OR. (ix > lpix) ) CYCLE

             ! ----------------------------
             ! Find the pixel quality flags
             ! ----------------------------
             ! -------------------------------------------------------------------
             ! CAREFUL: Only 15 flags/positions (0:14) can be returned or else the
             !          conversion will result in a numeric overflow.
             ! -------------------------------------------------------------------
             CALL convert_2bytes_to_16bits ( nbits-1, nwrr, &
                  omi_radiance_qflg(1:nwrr,ix,iloop), qflg_bit(1:nwrr,0:nbits-2) )
             ! --------------------------------------------------------------------
             ! Add contributions from various quality flags. Any CCD pixel that has
             ! a cumulative quality flag > 0 will be excluded form the averaging.
             ! --------------------------------------------------------------------
             qflg_mask(1:nwrr) = 0_i2
             qflg_mask(1:nwrr) = &
                  qflg_bit(1:nwrr,qflg_mis_idx) + &   ! Missing pixel
                  qflg_bit(1:nwrr,qflg_bad_idx) + &   ! Bad pixel
                  qflg_bit(1:nwrr,qflg_err_idx) !+ &   ! Processing error
                  !qflg_bit(1:nwrr,qflg_tra_idx) + &   ! Transient pixel
                  !qflg_bit(1:nwrr,qflg_rts_idx) + &   ! RTS pixel
                  !qflg_bit(1:nwrr,qflg_sat_idx)       ! Saturation

             cntr8(1:nwrr) = 1.0_r8
             WHERE ( qflg_mask(1:nwrr) > 0_i2 )
                cntr8(1:nwrr) = 0.0_r8
             END WHERE

             ! ------------------------------------
             ! Only proceed if we have a good value
             ! ------------------------------------
             IF ( ANY ( cntr8(1:nwrr) > 0.0_r8 ) ) THEN

                omi_radiance_spec(1:nwrr,ix,iloop) = &
                     omi_radiance_spec(1:nwrr,ix,iloop)*cntr8(1:nwrr) * cntr8(1:nwrr)

                specsum = SUM ( omi_radiance_spec(1:nwrr,ix,iloop) ) / SUM ( cntr8(1:nwrr) )
                IF ( specsum == 0.0_r8 ) specsum = 1.0_r8

                specsum = 1.0_r8

                radref_spec(ix,1:nwrr) = &
                     radref_spec(ix,1:nwrr) + omi_radiance_spec(1:nwrr,ix,iloop)/specsum
                radref_wavl(ix,1:nwrr) = &
                     radref_wavl(ix,1:nwrr) + omi_radiance_wavl(1:nwrr,ix,iloop)
                allcount(ix,1:nwrr) = allcount(ix,1:nwrr) + cntr8(1:nwrr)
                dumcount(ix,1:nwrr) = dumcount(ix,1:nwrr) + 1.0_r8

                IF ( omi_szenith(ix,iloop) /= r4_missval .AND. &
                     omi_vzenith(ix,iloop) /= r4_missval         ) THEN
                   omi_radref_sza(ix) = omi_radref_sza(ix) + omi_szenith(ix,iloop)
                   omi_radref_vza(ix) = omi_radref_vza(ix) + omi_vzenith(ix,iloop)
                   szacount      (ix) = szacount  (ix) + 1.0_r4
                END IF

             END IF

          END DO

       END DO

    END DO

    ! -----------------------------------------------------------
    ! Now for the actual averaging and assignment of final arrays
    ! -----------------------------------------------------------
    n_comm_wvl = 0
    DO ix = 1, nxrr

       ! -----------------------------------
       ! Average the wavelengths and spectra
       ! -----------------------------------
       WHERE ( allcount(ix,1:nwrr) /= 0.0_r8 )
          radref_spec(ix,1:nwrr) = radref_spec(ix,1:nwrr) / allcount(ix,1:nwrr)
       END WHERE
       WHERE ( dumcount(ix,1:nwrr) /= 0.0_r8 )
          radref_wavl(ix,1:nwrr) = radref_wavl(ix,1:nwrr) / dumcount(ix,1:nwrr)
       END WHERE

       ! -------------------------------------------
       ! Average the Solar and Viewing Zenith Angles
       ! -------------------------------------------
       IF ( szacount(ix) > 0.0_r4 ) THEN
          omi_radref_sza(ix) = omi_radref_sza(ix) / szacount(ix)
          omi_radref_vza(ix) = omi_radref_vza(ix) / szacount(ix)
       ELSE
          omi_radref_sza(ix) = r4_missval
          omi_radref_vza(ix) = r4_missval             
       END IF

       ! -------------------------------------------------------------------------------
       ! Determine the CCD pixel numbers based on the selected wavelength fitting window
       ! -------------------------------------------------------------------------------

       radref_wavl_ix = radref_wavl (ix, 1:nwrr)
       DO j1 = 1, 3, 2
          CALL array_locate_r8 ( &
               nwrr, radref_wavl_ix, fit_winwav_lim(j1  ), 'LE', &
               omi_ccdpix_selection(ix,j1  ) )
          CALL array_locate_r8 ( &
               nwrr, radref_wavl_ix, fit_winwav_lim(j1+1), 'GE', &
               omi_ccdpix_selection(ix,j1+1) )
       END DO

       imin = omi_ccdpix_selection(ix,1)
       imax = omi_ccdpix_selection(ix,4)

       icnt = imax - imin + 1
       omi_nwav_radref(       ix) = icnt
       omi_radref_spec(1:icnt,ix) = radref_spec(ix,imin:imax)
       omi_radref_wavl(1:icnt,ix) = radref_wavl_ix(imin:imax)
       omi_radref_qflg(1:icnt,ix) = 0_i2
       omi_radref_wght(1:icnt,ix) = normweight

       ! -----------------------------------------------------------------
       ! Re-assign the average solar wavelength variable, sinfe from here
       ! on we are concerned with radiances.
       ! -----------------------------------------------------------------
       omi_sol_wav_avg(ix) =  SUM( omi_radref_wavl(1:icnt,ix) ) / REAL(icnt, KIND=r8)

       ! ------------------------------------------------------------------
       ! Set weights and quality flags to "bad" for missing spectral points
       ! ------------------------------------------------------------------
       allcount(ix,1:icnt) = allcount(ix,imin:imax)
       WHERE ( allcount(ix,1:icnt) == 0.0_r8 )
          omi_radref_qflg(1:icnt,ix) = 7_i2
          omi_radref_wght(1:icnt,ix) = downweight
       END WHERE
          
       ! ------------------------------------------------------------------------------
       ! If any window is excluded, find the corresponding indices. This has to be done
       ! after the array assignements above because we need to know which indices to
       ! exclude from the final arrays, not the complete ones read from the HE4 file.
       ! ------------------------------------------------------------------------------
       omi_ccdpix_exclusion(ix,1:2) = -1
       IF ( MINVAL(fit_winexc_lim(1:2)) > 0.0_r8 ) THEN
          CALL array_locate_r8 ( &
               nwrr, radref_wavl_ix, fit_winexc_lim(1), 'GE', &
               omi_ccdpix_exclusion(ix,1) )
          CALL array_locate_r8 ( &
               nwrr, radref_wavl_ix, fit_winexc_lim(2), 'LE', &
               omi_ccdpix_exclusion(ix,2) )
       END IF

       ! ----------------------------------------
       ! Update the maximum number of wavelengths
       ! ----------------------------------------
       n_comm_wvl = MAX ( n_comm_wvl, icnt )

    END DO

    ! --------------------
    ! Deallocate variables
    ! --------------------
    DEALLOCATE (xtrange)        ; DEALLOCATE (binfac)
    DEALLOCATE (ynzoom)         ; DEALLOCATE (szacount)
    DEALLOCATE (latr4)
    DEALLOCATE (radref_spec)    ; DEALLOCATE (radref_wavl)
    DEALLOCATE (radref_wavl_ix) ; DEALLOCATE (allcount)
    DEALLOCATE (dumcount)       ; DEALLOCATE (cntr8)
    DEALLOCATE (qflg_bit)       ; DEALLOCATE (qflg_mask)

    RETURN
  END SUBROUTINE omi_get_radiance_reference



  SUBROUTINE xtrack_radiance_reference_loop (     &
       yn_radiance_reference, yn_remove_target, nx, nw, fpix, lpix, pge_idx, errstat )

    USE OMSAO_indices_module,    ONLY: &
         wvl_idx, spc_idx, sig_idx, o3_t1_idx, o3_t3_idx, hwe_idx, asy_idx, shi_idx, squ_idx, &
         pge_bro_idx, pge_o3_idx, pge_hcho_idx, n_max_fitpars, solar_idx, ccd_idx, radref_idx,&
         bro_idx
    USE OMSAO_parameters_module, ONLY:  &
         i2_missval, i4_missval, r4_missval, r8_missval, downweight, normweight
    USE OMSAO_variables_module,  ONLY:  &
         database, curr_sol_spec, n_rad_wvl, curr_rad_spec, sol_wav_avg,                  &
         hw1e, e_asym, n_fitvar_rad, verb_thresh_lev, fitvar_rad_saved, fitvar_rad_init,  &
         n_database_wvl, fitvar_rad, n_fincol_idx, fincol_idx,                            &
         n_fitres_loop, fitres_range, refspecs_original, xtrack_fitres_limit, &
         yn_solar_comp, n_rad_wvl_max

    USE OMSAO_prefitcol_module, ONLY:                                                     &
         yn_o3_prefit, o3_prefit_col, o3_prefit_dcol,                                     &
         yn_bro_prefit, bro_prefit_col, bro_prefit_dcol,                                  &
         yn_lqh2o_prefit, lqh2o_prefit_col, lqh2o_prefit_dcol
    USE OMSAO_slitfunction_module, ONLY: saved_shift, saved_squeeze
    USE OMSAO_omidata_module, ONLY: &
         omi_nwav_irrad, omi_irradiance_wght, omi_nwav_rad, n_omi_database_wvl, &
         omi_cross_track_skippix, curr_xtrack_pixnum, n_omi_radwvl, max_rs_idx, &
         omi_database, omi_database_wvl, omi_sol_wav_avg, omi_solcal_pars,      &
         omi_radiance_wavl, omi_radref_wavl, omi_radiance_spec, omi_radref_spec,&
         omi_radiance_qflg, omi_radref_qflg, omi_radiance_spec,                 &
         omi_ccdpix_selection, omi_radiance_ccdpix, omi_ccdpix_exclusion,       &
         omi_xtrackpix_no, omi_radref_wght, omi_radref_pars, max_calfit_idx,    &
         omi_radref_xflag, omi_radref_itnum, omi_radref_chisq, omi_radref_col,  &
         omi_radref_rms, omi_radref_dcol, omi_radref_xtrcol
    USE OMSAO_errstat_module

    IMPLICIT NONE

    ! ---------------
    ! Input Variables
    ! ---------------
    INTEGER (KIND=i4), INTENT (IN) :: pge_idx, nx, nw, fpix, lpix
    LOGICAL,           INTENT (IN) :: yn_radiance_reference, yn_remove_target

    ! -----------------
    ! Modified variable
    ! -----------------
    INTEGER (KIND=i4), INTENT (INOUT) :: errstat

    ! ---------------
    ! Local variables
    ! ---------------
    INTEGER (KIND=i4) :: locerrstat, ipix, jpix, radfit_exval, radfit_itnum, i
    REAL    (KIND=r8) :: fitcol, rms, dfitcol, chisquav, rad_spec_avg
    REAL    (KIND=r8), DIMENSION (o3_t1_idx:o3_t3_idx) :: o3fit_cols, o3fit_dcols
    REAL    (KIND=r8), DIMENSION (n_fitvar_rad)        :: corr_matrix_tmp, allfit_cols_tmp, allfit_errs_tmp
    LOGICAL                  :: yn_skip_pix
    CHARACTER (LEN=maxchlen) :: addmsg
    LOGICAL                                          :: yn_bad_pixel
    INTEGER (KIND=i4), DIMENSION (4)                 :: select_idx
    INTEGER (KIND=i4), DIMENSION (2)                 :: exclud_idx
    INTEGER (KIND=i4)                                :: n_solar_pts
    REAL    (KIND=r8), DIMENSION (1:nw)              :: solar_wgt
    REAL    (KIND=r8), DIMENSION (n_fincol_idx,1:nx) :: target_var 

    CHARACTER (LEN=30), PARAMETER :: modulename = 'xtrack_radiance_reference_loop'

    ! CCM fitted spectrum now returned from radiance_fit.f90
    !REAL    (KIND=r8), DIMENSION (n_rad_wvl)         :: fitspctmp
    ! The above will not work since size varies with the loop index --JED
    REAL    (KIND=r8), DIMENSION (n_rad_wvl_max) :: fitspctmp

    ! -------------------------
    ! Initialize some variables
    ! -------------------------
    locerrstat          = pge_errstat_ok
    xtrack_fitres_limit = 0.0_r8
    target_var          = r8_missval
    fitvar_rad_saved    = fitvar_rad_init


    ! ---------------------------------------------------
    ! Note that this initialization will overwrite valid
    ! results on any second call to this subroutine. This
    ! happens, for example, when YN_RADIANCE_REFERENCE
    ! and YN_REMOVE_TARGET are selected simultaneously.
    ! In that case, however, we write the results to file
    ! before the second call.
    ! ---------------------------------------------------
    omi_radref_pars  (1:max_calfit_idx,1:nx) = r8_missval
    omi_radref_xflag (1:nx)                  = i2_missval
    omi_radref_itnum (1:nx)                  = i2_missval
    omi_radref_chisq (1:nx)                  = r8_missval
    omi_radref_col   (1:nx)                  = r8_missval
    omi_radref_dcol  (1:nx)                  = r8_missval
    omi_radref_rms   (1:nx)                  = r8_missval
    omi_radref_xtrcol(1:nx)                  = r8_missval


    XTrackPix: DO jpix = fpix, lpix


       locerrstat = pge_errstat_ok

       ! -----------------------------------------------------------
       ! The current cross-track pixel number is required further on
       ! when the slit function is computed: The CCD position based
       ! hyper-parameterization requires the knowledge of the row#.
       ! -----------------------------------------------------------
       ipix = jpix

       curr_xtrack_pixnum = ipix

       ! ---------------------------------------------------------------------
       ! If we already determined that this cross track pixel position carries
       ! an error, we don't even have to start processing.
       ! ---------------------------------------------------------------------
       IF ( omi_cross_track_skippix(ipix) ) CYCLE

       n_database_wvl = n_omi_database_wvl(ipix)
       n_omi_radwvl   = omi_nwav_rad      (ipix,0)

       ! ---------------------------------------------------------------------------
       ! For each cross-track position we have to initialize the saved Shift&Squeeze
       ! ---------------------------------------------------------------------------
       saved_shift = -1.0e+30_r8 ; saved_squeeze = -1.0e+30_r8

       ! ------------------------------------------------------------------
       ! Assign number of irradiance wavelengths and the fitting weights
       ! from the solar wavelength calibration. Why? gga
       ! ------------------------------------------------------------------
       n_solar_pts              = omi_nwav_irrad(ipix)
       !       IF (yn_solar_comp) THEN !gga
       !          solar_wgt(1:n_solar_pts) = 1.0_r8
       !       ELSE
       !          solar_wgt(1:n_solar_pts) = omi_irradiance_wght(1:n_solar_pts,ipix)
       !       END IF !gga
       solar_wgt(1:n_solar_pts) = normweight

       ! -----------------------------------------------------
       ! Catch the possibility that N_OMI_RADWVL > N_SOLAR_PTS
       ! -----------------------------------------------------
       IF ( n_omi_radwvl > n_solar_pts ) THEN
          i = n_omi_radwvl - n_solar_pts
          solar_wgt(n_solar_pts+1:n_solar_pts+i) = downweight
          n_solar_pts = n_omi_radwvl
       END IF

       ! tpk: Should this be "> n_fitvar_rad"?
       IF ( n_database_wvl > 0 .AND. n_omi_radwvl > 0 ) THEN

          ! ----------------------------------------------
          ! Restore DATABASE from OMI_DATABASE (see above)
          ! ----------------------------------------------
          database (1:max_rs_idx,1:n_database_wvl) = omi_database (1:max_rs_idx,1:n_database_wvl,ipix)

          ! -----------------------------------------------------------------------
          ! Restore solar fitting variables for across-track reference in Earthshine fitting
          ! --------------------------------------------------------------------------------
          sol_wav_avg                             = omi_sol_wav_avg(ipix)
          hw1e                                    = omi_solcal_pars(hwe_idx,ipix)
          e_asym                                  = omi_solcal_pars(asy_idx,ipix)
          curr_sol_spec(wvl_idx,1:n_database_wvl) = omi_database_wvl(1:n_database_wvl,ipix)
          curr_sol_spec(spc_idx,1:n_database_wvl) = omi_database    (solar_idx,1:n_database_wvl,ipix)
          ! --------------------------------------------------------------------------------


          ! ---------------------------------------------------------------
          ! If a Radiance Reference is being used, then it must be calibrated
          ! rather than the swath line that has been read.
          ! ---------------------------------------------------------------
          IF ( yn_radiance_reference ) THEN
             omi_radiance_wavl(1:n_omi_radwvl,ipix,0) = omi_radref_wavl(1:n_omi_radwvl,ipix)
             omi_radiance_spec(1:n_omi_radwvl,ipix,0) = omi_radref_spec(1:n_omi_radwvl,ipix)
             omi_radiance_qflg(1:n_omi_radwvl,ipix,0) = omi_radref_qflg(1:n_omi_radwvl,ipix)
          END IF

          omi_xtrackpix_no = ipix

          ! -------------------------------------------------------------------------
          select_idx(1:4) = omi_ccdpix_selection(ipix,1:4)
          exclud_idx(1:2) = omi_ccdpix_exclusion(ipix,1:2)
          CALL omi_adjust_radiance_data ( &           ! Set up generic fitting arrays
               select_idx(1:4), exclud_idx(1:2),            &
               n_omi_radwvl,                                &
               omi_radiance_wavl   (1:n_omi_radwvl,ipix,0), &
               omi_radiance_spec   (1:n_omi_radwvl,ipix,0), &
               omi_radiance_qflg   (1:n_omi_radwvl,ipix,0), &
               omi_radiance_ccdpix (1:n_omi_radwvl,ipix,0), &
               n_solar_pts, solar_wgt(1:n_solar_pts),       &
               n_rad_wvl, curr_rad_spec(wvl_idx:ccd_idx,1:n_omi_radwvl), rad_spec_avg, &
               yn_skip_pix )

          ! -------------------------------------------------------------------------

          ! --------------------------------------------------------------------
          ! Update the weights for the Reference/Wavelength Calibration Radiance
          ! --------------------------------------------------------------------
          omi_radref_wght(1:n_omi_radwvl,ipix) = curr_rad_spec(sig_idx,1:n_omi_radwvl)


          IF ( pge_idx == pge_hcho_idx .AND. &
               ( .NOT. yn_radiance_reference )  ) THEN
             o3fit_cols (o3_t1_idx:o3_t3_idx) = o3_prefit_col (o3_t1_idx:o3_t3_idx,ipix,0)
             o3fit_dcols(o3_t1_idx:o3_t3_idx) = o3_prefit_dcol(o3_t1_idx:o3_t3_idx,ipix,0)
          ELSE
             o3fit_cols (o3_t1_idx:o3_t3_idx) = 0.0_r8
             o3fit_dcols(o3_t1_idx:o3_t3_idx) = 0.0_r8
          END IF

          ! --------------------
          ! The radiance fitting
          ! --------------------
          fitcol       = r8_missval
          dfitcol      = r8_missval
          radfit_exval = INT(i2_missval, KIND=i4)
          radfit_itnum = INT(i2_missval, KIND=i4)
          rms          = r8_missval

          addmsg = ''
          IF ( MAXVAL(curr_rad_spec(spc_idx,1:n_rad_wvl)) > 0.0_r8 .AND.     &
               n_rad_wvl > n_fitvar_rad .AND. (.NOT. yn_skip_pix)              ) THEN

             yn_bad_pixel     = .FALSE.

             CALL radiance_fit ( &
                  pge_idx, ipix, n_fitres_loop(radref_idx), fitres_range(radref_idx),       &
                  yn_reference_fit,                                                         &
                  n_rad_wvl, curr_rad_spec(wvl_idx:ccd_idx,1:n_rad_wvl),                    &
                  fitcol, rms, dfitcol, radfit_exval, radfit_itnum, chisquav,               &
                  o3fit_cols, o3fit_dcols, bro_prefit_col(ipix,0), bro_prefit_dcol(ipix,0), &
                  lqh2o_prefit_col(ipix,0), lqh2o_prefit_dcol(ipix,0),                      &
                  target_var(1:n_fincol_idx,ipix),                                          &
                  allfit_cols_tmp(1:n_fitvar_rad), allfit_errs_tmp(1:n_fitvar_rad),         &
                  corr_matrix_tmp(1:n_fitvar_rad), yn_bad_pixel, fitspctmp )

             yn_reference_fit = .FALSE.

             IF ( yn_bad_pixel ) CYCLE

             WRITE (addmsg, '(A,I2,4(A,1PE10.3),2(A,I5))') 'RADIANCE Reference #', ipix, &
                  ': hw 1/e = ', hw1e, '; e_asy = ', e_asym, '; shift = ', &
                  fitvar_rad(shi_idx), '; squeeze = ', fitvar_rad(squ_idx),&
                  '; exit val = ', radfit_exval, '; iter num = ', radfit_itnum
          ELSE
             WRITE (addmsg, '(A,I2,A)') 'RADIANCE Reference #', ipix, ': Skipped!'
          END IF

          ! ------------------
          ! Report on progress
          ! ------------------
          CALL error_check ( &
               0, 1, pge_errstat_ok, OMSAO_S_PROGRESS, TRIM(ADJUSTL(addmsg)), vb_lev_omidebug, errstat )
          IF ( verb_thresh_lev >= vb_lev_screen ) WRITE (*, '(A)') TRIM(ADJUSTL(addmsg))

          ! -----------------------------------
          ! Assign pixel values to final arrays
          ! -----------------------------------
          omi_radref_pars (1:max_calfit_idx,ipix) = fitvar_rad(1:max_calfit_idx)
          omi_radref_xflag(ipix)                  = INT (radfit_exval, KIND=i2)
          omi_radref_itnum(ipix)                  = INT (radfit_itnum, KIND=i2)
          omi_radref_chisq(ipix)                  = chisquav
          omi_radref_col  (ipix)                  = fitcol
          omi_radref_dcol (ipix)                  = dfitcol
          omi_radref_rms  (ipix)                  = rms

          ! -------------------------------------------------------------------------
          ! Remember weights for the reference radiance, to be used as starting point
          ! in the regular radiance fitting
          ! -------------------------------------------------------------------------
          omi_radref_wght(1:n_rad_wvl,ipix) = curr_rad_spec(sig_idx,1:n_rad_wvl)

          ! -----------------------------------------------
          ! Update the solar spectrum entry in OMI_DATABASE
          ! -----------------------------------------------
          IF ( yn_radiance_reference ) &
               omi_database (solar_idx,1:n_rad_wvl,ipix) = omi_radref_spec(1:n_rad_wvl,ipix)

       END IF

    END DO XTrackPix

    ! -----------------------------------------
    ! Remove target gas from radiance reference
    ! -----------------------------------------

    IF ( yn_radiance_reference .AND. yn_remove_target ) THEN
       ! ----------------------------------------------------------------
       ! Removing the target gas from the radiance reference will alter
       ! OMI_RADREF_SPEC (1:NWVL,FPIX:LPIX). This is being passed to the
       ! subroutine via MODULE use rather than through the argument list.
       ! ----------------------------------------------------------------
       CALL remove_target_from_radiance (                                  &
            nw, fpix, lpix, n_fincol_idx, fincol_idx(1:2,1:n_fincol_idx),  &
            target_npol, target_var(1:n_fincol_idx,fpix:lpix), omi_radref_xtrcol(fpix:lpix) )

    END IF

    errstat = MAX ( errstat, locerrstat )

    RETURN
  END SUBROUTINE xtrack_radiance_reference_loop


  SUBROUTINE remove_target_from_radiance (       &
       nw, ipix, jpix, n_fincol_idx, fincol_idx, &
       target_npol, target_var, target_fit         )

    USE OMSAO_indices_module,   ONLY: solar_idx
    USE OMSAO_parameters_module,ONLY: downweight, r8_missval
    USE OMSAO_omidata_module,   ONLY: omi_radref_spec, n_omi_database_wvl, omi_database, omi_nwav_radref
    USE OMSAO_variables_module, ONLY: refspecs_original
    USE OMSAO_median_module,    ONLY: median

    IMPLICIT NONE


    ! ---------------
    ! Input Variables
    ! ---------------
    INTEGER (KIND=i4),                                     INTENT (IN) :: &
         nw, ipix, jpix, n_fincol_idx, target_npol
    INTEGER (KIND=i4), DIMENSION (2,n_fincol_idx),         INTENT (IN) :: fincol_idx

    ! ------------------
    ! Modified Variables
    ! ------------------
    REAL (KIND=r8), DIMENSION (n_fincol_idx,ipix:jpix), INTENT (INOUT) :: target_var
    REAL (KIND=r8), DIMENSION              (ipix:jpix), INTENT (OUT)   :: target_fit

    ! ------------------------------
    ! Local Variables and Parameters
    ! ------------------------------
    INTEGER (KIND=i4)                 :: i, j, k, l, nwvl
    REAL    (KIND=r8)                 :: yfloc
    REAL    (KIND=r8), DIMENSION (nw) :: tmpexp

    CHARACTER (LEN=27), PARAMETER :: modulename = 'remove_target_from_radiance'

    ! ----------------
    ! DPOLFt variables
    ! ----------------
    INTEGER (KIND=i4)            :: ndeg, ierr, nx, npol, nfit
    REAL    (KIND=r8)            :: eps
    REAL    (KIND=r8), DIMENSION (jpix-ipix+1) :: x, y, yf, w
    REAL    (KIND=r8), DIMENSION (3*((jpix-ipix+1)+target_npol+1)) :: a

    target_fit = 0.0_r8

    npol = target_npol
    nx   = jpix-ipix+1
    DO j = 1, n_fincol_idx
       k = fincol_idx(2,j)
     
       ! ----------------------------------------------------------------------
       ! If we can/have to, fit a cross-track polynomial to the fitted columns,
       ! we do this individually for each FINCOL_IDX and remove the smoothed 
       ! column loading rather than the originally fitted one. In any case, YF
       ! will contain the column values to be removed. Hence the outer loop 
       ! over N_FINCOL_IDX rather than cross-track position.
       ! ----------------------------------------------------------------------
       nfit = 0
       DO i = ipix, jpix
          IF ( target_var(j,i) > r8_missval ) THEN
             nfit    = nfit + 1
             y(nfit) = target_var(j,i)
             w(nfit) = 1.0_r8
          END IF
       END DO

       IF ( nfit /= (jpix-ipix+1) .AND. nfit > 0 ) THEN
          WHERE ( target_var(j,ipix:jpix) <= r8_missval )
             target_var(j,ipix:jpix) = SUM(y(1:nfit))/REAL(nfit,KIND=r8)
          ENDWHERE
       END IF

       ! ----------------------------------------------------
       ! We either fit a polynomial or simply use the Median.
       ! The distinction is made depending on
       !
       ! (a) the order of the cross-track polynomial, and
       ! (b) the number of cross-track points we can fit.
       ! ----------------------------------------------------
       IF ( npol > 0 .AND. jpix-ipix+1 > npol ) THEN

          !IF ( npol >=0 .AND. nfit > npol ) THEN
          !eps  = -0.1_r8  ! Chose the best-fitting order
           
          eps =  0.0_r8  ! Fit the complete NPOL polynomial
          ndeg = npol
          x(1:nx) = (/ ( REAL(i-nx/2, KIND=r8), i = 1, nx ) /) / REAL(nx/2, KIND=r8)
          y(1:nx) = target_var(j,ipix:jpix)
          WHERE ( y(1:nx) > r8_missval )
             w(1:nx) = 1.0_r8
          ELSEWHERE
             w(1:nx) = downweight
          END WHERE
          CALL dpolft (&
               nx, x(1:nx), y(1:nx), w(1:nx), npol, ndeg, eps, yf(1:nx), ierr, a )

          !CALL dpolft (&
          !     nfit, x(1:fit), y(1:nfit), w(1:nfit), npol, ndeg, eps, yf(1:nfit), ierr, a )

       ELSE
          ! -----------------------------------------------------------------
          ! The Median is a better choice than the Mean, since the former
          ! is less sensitive to outliers. The Mean may be skewed towards
          ! abnormally high values at the edges of the swath.
          ! -----------------------------------------------------------------
          nfit = nx
          yf(1:nx) = median(nx, target_var(j,ipix:jpix))
       END IF

       DO i = ipix, jpix

          l = i - ipix + 1

          nwvl = omi_nwav_radref(i)

          IF ( yf(l) > r8_missval ) THEN
             yfloc = yf(l)
          ELSE
             yfloc = 0.0_r8
          END IF

          tmpexp(1:nwvl) = yfloc * omi_database(k,1:nwvl,i)
          WHERE ( tmpexp >= MAXEXPONENT(1.0_r8) )
             tmpexp = MAXEXPONENT(1.0_r8) - 1.0_r8
          ENDWHERE
          WHERE ( tmpexp <= MINEXPONENT(1.0_r8) )
             tmpexp = MINEXPONENT(1.0_r8) + 1.0_r8
          ENDWHERE

          omi_radref_spec(1:nwvl,i) = omi_radref_spec(1:nwvl,i) * &
               EXP(+tmpexp(1:nwvl))

          target_fit(i) = target_fit(i) + yfloc/refspecs_original(k)%NormFactor

       END DO
    END DO

    RETURN
  END SUBROUTINE remove_target_from_radiance


END MODULE OMSAO_radiance_ref_module
