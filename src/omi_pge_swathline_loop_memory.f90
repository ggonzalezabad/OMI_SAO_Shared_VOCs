SUBROUTINE omi_pge_swathline_loop_memory (                               &
     pge_idx, nt, nx, nccd, n_max_rspec, yn_process,                     &
     xtrange, yn_radiance_reference, yn_remove_target, ntargpol,         &
     yn_commit, mem_column_amount, mem_column_uncertainty, mem_rms,      &
     mem_fit_flag, mem_xtrflg, mem_latitude, mem_longitude, mem_sza,     &
     mem_vza, mem_height, errstat)


  USE OMSAO_precision_module,  ONLY: i4, r8
  USE OMSAO_parameters_module, ONLY: i2_missval, r8_missval, main_qa_missing, &
       maxchlen, r4_missval
  USE OMSAO_indices_module,    ONLY: pge_hcho_idx, n_max_fitpars, solar_idx
  USE OMSAO_variables_module,  ONLY:  &
       n_fitvar_rad, fitvar_rad_init, fitvar_rad_saved, l1b_rad_filename, &
       verb_thresh_lev, n_fincol_idx, fincol_idx, n_rad_wvl, n_rad_wvl_max
  USE OMSAO_omidata_module,    ONLY:  &
       nlines_max, nUTCdim, omi_scanline_no, omi_blockline_no,                  &
       omi_itnum_flag, omi_fitconv_flag, omi_column_amount,                     &
       omi_column_uncert, omi_time_utc, omi_time, omi_latitude, omi_fit_rms,    &
       omi_radiance_errstat, omi_nwav_radref, omi_radref_spec, omi_radref_wavl, &
       omi_szenith, omi_vzenith, omi_longitude, omi_xtrflg, omi_height
  USE OMSAO_destriping_module, ONLY: ctr_maxcol
  USE OMSAO_prefitcol_module ! brings in nxtrack_max, nlines_max, ...
  USE OMSAO_errstat_module
  USE OMSAO_radiance_ref_module, ONLY: remove_target_from_radiance

  IMPLICIT NONE


  ! ---------------
  ! Input variables
  ! ---------------
  INTEGER (KIND=i4), INTENT (IN) :: pge_idx, nx, nt, nccd, ntargpol, n_max_rspec
  INTEGER (KIND=i4), DIMENSION (0:nt-1,1:2),  INTENT (IN) :: xtrange
  LOGICAL,           DIMENSION (0:nt-1),      INTENT (IN) :: yn_process
  LOGICAL,           INTENT (IN) :: yn_commit, yn_radiance_reference, yn_remove_target

  ! ------------------
  ! Modified variables
  ! ------------------
  INTEGER (KIND=i4), INTENT (INOUT) :: errstat

  ! ---------------
  ! Local variables
  ! ---------------
  INTEGER   (KIND=i4)      :: iline, iloop, nblock, fpix, lpix, ipix, estat, locerrstat
  CHARACTER (LEN=maxchlen) :: addmsg

  ! CCM for looping
  INTEGER (KIND=i4) :: i,j,k

  ! ---------------------------------------------------------------
  ! Variables to remove target gas from radiance reference spectrum
  ! ---------------------------------------------------------------
  REAL (KIND=r8), DIMENSION (n_fincol_idx,1:nx) :: target_var, targsum, targcnt
  REAL (KIND=r8), DIMENSION (1:nx)              :: target_fit, target_col

  ! ---------------------------------------------------------------------------------
  ! CCM Array to hold (1) Fitted Spec (2) Observed Spec (3) Spec Pos (4) Weight flags
  ! ---------------------------------------------------------------------------------
  REAL (KIND=r8), DIMENSION (n_rad_wvl_max,nxtrack_max,4) :: fitspc_tmp
!!$  REAL (KIND=r8), DIMENSION (n_rad_wvl_max,nxtrack_max,4,0:nt-1) :: omi_fitspc

  ! -------------------------------------
  ! Correlations with main output product 
  ! -------------------------------------
  REAL (KIND=r8), DIMENSION (n_fitvar_rad,nx,0:nlines_max-1) :: &
       all_fitted_columns, all_fitted_errors, correlation_columns

  ! -----------------------------
  ! Keeping the results in memory 
  ! -----------------------------
  REAL    (KIND=r8), DIMENSION (nx,0:nt-1), INTENT(INOUT) :: &
       mem_column_amount, mem_column_uncertainty, mem_rms
  REAL    (KIND=r4), DIMENSION (nx,0:nt-1), INTENT(INOUT) :: &
       mem_latitude, mem_longitude, mem_sza, mem_vza, mem_height
  INTEGER (KIND=i2), DIMENSION (nx,0:nt-1), INTENT(INOUT) :: mem_fit_flag, mem_xtrflg

  locerrstat = pge_errstat_ok

  ! --------------------------------
  ! Initialize fitting output arrays
  ! --------------------------------
  all_fitted_columns  = r8_missval
  all_fitted_errors   = r8_missval
  correlation_columns = r8_missval

  ! ------------------------
  ! Initialize memory arrays
  ! ------------------------
  mem_column_amount      = r8_missval
  mem_column_uncertainty = r8_missval
  mem_rms                = r8_missval
  mem_latitude           = r4_missval
  mem_longitude          = r4_missval
  mem_sza                = r4_missval
  mem_vza                = r4_missval
  mem_fit_flag           = i2_missval

  IF ( yn_radiance_reference .AND. yn_remove_target ) THEN
     target_var = 0.0_r8
     targsum    = 0.0_r8
     targcnt    = 0.0_r8
     target_fit = 0.0_r8
     target_col = 0.0_r8
  END IF

  ! ---------------------------------------------------------------------
  ! Loop over all scan lines, in multiples of NLINES_MAX (100 by default)
  ! ---------------------------------------------------------------------
  ScanLines: DO iline = 0, nt-1, nlines_max

     ! ---------------------------------------------------------
     ! Check if loop ends before n_times_loop max is exhausted.
     ! Not a serious problem but it saves a few bytes of memory.
     ! ---------------------------------------------------------
     nblock = nlines_max
     IF ( (iline+nblock) > nt ) nblock = nt - iline
     ! ------------------------------
     ! Get NBLOCK radiance lines
     ! ------------------------------
     CALL omi_read_radiance_lines (                   &
          l1b_rad_filename, iline, nx, nblock, nccd, locerrstat )
     ! -----------------------------------------------------------------------------------

     ! ------------------------------------------
     ! Initialize output fields with MissingValue
     ! ------------------------------------------
     omi_itnum_flag   (1:nx,     0:nblock-1) = i2_missval
     omi_fitconv_flag (1:nx,     0:nblock-1) = i2_missval
     omi_column_amount(1:nx,     0:nblock-1) = r8_missval
     omi_column_uncert(1:nx,     0:nblock-1) = r8_missval
     omi_fit_rms      (1:nx,     0:nblock-1) = r8_missval
     omi_time_utc     (1:nUTCdim,0:nblock-1) = i2_missval

     ! -----------------------------------------
     ! Skip if we don't have anything to process
     ! -----------------------------------------
     IF ( .NOT. ( ANY ( yn_process(iline:iline+nblock-1) ) ) ) CYCLE

     ! --------------------------------
     ! Read pre-fitted molecule columns
     ! --------------------------------
     IF ( ( .NOT. yn_radiance_reference ) .AND. ( pge_idx == pge_hcho_idx ) .AND. &
          ANY( (/yn_o3_prefit(1), yn_bro_prefit(1)/) ) ) THEN
        CALL read_prefit_columns ( pge_idx, nx, nblock, iline, locerrstat )
        errstat = MAX ( errstat, locerrstat )
        IF ( errstat >= pge_errstat_error ) RETURN
     END IF

     ! -------------------------------------
     ! Re-initialize saved fitting variables
     ! -------------------------------------
     fitvar_rad_saved(1:n_max_fitpars ) = fitvar_rad_init(1:n_max_fitpars)

     ! -----------------------------------------------
     ! Loops over all scan lines in current data block
     ! -----------------------------------------------
     ScanLineBlock: DO iloop = 0, nblock-1

        ! --------------------------------------------------------------------
        ! Further down, in deeper layers of the algorithm, we require both the
        ! current line in the data block and the absolute swath line number.
        ! Both values are initialized here.
        ! --------------------------------------------------------------------
        omi_blockline_no = iloop
        omi_scanline_no  = iline+iloop

        IF ( omi_scanline_no > nt-1 ) EXIT ScanLines

        ! ----------------------------------------------------------
        ! Skip this line if it isn't in the list of those to process
        ! ----------------------------------------------------------
        IF ( .NOT. yn_process(omi_scanline_no) ) CYCLE

        ! ------------------
        ! Report on progress
        ! ------------------
        addmsg = ''
        WRITE (addmsg,'(A,I5)') 'Working on scan line', omi_scanline_no
        estat = OMI_SMF_setmsg ( OMSAO_S_PROGRESS, TRIM(ADJUSTL(addmsg)), " ", vb_lev_omidebug )
        !IF ( verb_thresh_lev >= vb_lev_screen ) WRITE (*, '(A)') TRIM(ADJUSTL(addmsg))

        IF ( omi_radiance_errstat(iloop) /= pge_errstat_error ) THEN

           fpix = xtrange(omi_scanline_no,1)
           lpix = xtrange(omi_scanline_no,2)

           CALL xtrack_radiance_fitting_loop (                         &
                n_max_rspec, fpix, lpix, pge_idx, iloop,               &
                ctr_maxcol, n_fitvar_rad,                              &
                all_fitted_columns (1:n_fitvar_rad,fpix:lpix,iloop),   & !gga (1:nx to fpix:lpix)
                all_fitted_errors  (1:n_fitvar_rad,fpix:lpix,iloop),   & !gga (1:nx to fpix:lpix)
                correlation_columns(1:n_fitvar_rad,fpix:lpix,iloop),   & !gga (1:nx to fpix:lpix)
                target_var(1:n_fincol_idx,fpix:lpix), locerrstat, fitspc_tmp, n_rad_wvl_max)
           ipix = (fpix+lpix)/2
           addmsg = ''
           WRITE (addmsg,'(I5, I3,3(1PE15.5),I5)') omi_scanline_no, ipix, &
                omi_column_amount(ipix, iloop), omi_column_uncert(ipix, iloop), &
                omi_fit_rms   (ipix, iloop), MAX(-1,omi_itnum_flag(ipix, iloop))
           estat = OMI_SMF_setmsg ( OMSAO_S_PROGRESS, TRIM(addmsg), " ", vb_lev_omidebug )
           IF ( verb_thresh_lev >= vb_lev_screen ) WRITE (*, '(A)') TRIM(addmsg)

           ! CCM Add omi_fitspc - Assignment problem - do an inefficient loop for now
!!$           DO i=1,n_rad_wvl
!!$              DO j=1,nxtrack_max
!!$                 DO k=1,4
!!$                    omi_fitspc(i,j,k,iloop) = fitspc_tmp(i,j,k)
!!$                 ENDDO
!!$              ENDDO
!!$           ENDDO

           ! ---------------------------------------------------------------
           ! Add fitted columns for possible removal from radiance reference
           ! ---------------------------------------------------------------
           IF ( yn_radiance_reference .AND. yn_remove_target ) THEN
              DO ipix = fpix, lpix
                 IF ( &
                      ( omi_fitconv_flag (ipix,iloop) > 0_i2       ) .AND. &
                      ( omi_column_amount(ipix,iloop) > r8_missval ) .AND. &
                      ( omi_column_amount(ipix,iloop) + &
                      2.0_r8*omi_column_uncert(ipix,iloop) >= 0.0_r8 )  ) THEN
                    targsum(1:n_fincol_idx,ipix) = &
                         targsum(1:n_fincol_idx,ipix) + target_var(1:n_fincol_idx,ipix)
                    targcnt(1:n_fincol_idx,ipix) = &
                         targcnt(1:n_fincol_idx,ipix) + 1.0_r8

                    target_col(ipix) = target_col(ipix) + omi_column_amount(ipix,iloop)
                 END IF
              END DO
           END IF

           ! --------------------------------------------
           ! Keeping the results of the fitting in memory
           ! --------------------------------------------
           mem_column_amount(fpix:lpix,omi_scanline_no)      = omi_column_amount(fpix:lpix,iloop)
           mem_column_uncertainty(fpix:lpix,omi_scanline_no) = omi_column_uncert(fpix:lpix,iloop)
           mem_rms(fpix:lpix,omi_scanline_no)                = omi_fit_rms(fpix:lpix,iloop)
           mem_latitude(fpix:lpix,omi_scanline_no)           = omi_latitude(fpix:lpix,iloop)
           mem_longitude(fpix:lpix,omi_scanline_no)          = omi_longitude(fpix:lpix,iloop)
           mem_sza(fpix:lpix,omi_scanline_no)                = omi_szenith(fpix:lpix,iloop)
           mem_vza(fpix:lpix,omi_scanline_no)                = omi_vzenith(fpix:lpix,iloop)
           mem_fit_flag(fpix:lpix,omi_scanline_no)           = omi_fitconv_flag(fpix:lpix,iloop)
           mem_xtrflg(fpix:lpix,omi_scanline_no)             = omi_xtrflg(fpix:lpix,iloop)
           mem_height(fpix:lpix,omi_scanline_no)             = REAL(omi_height(fpix:lpix,iloop), KIND = r4)

        END IF

        ! -----------------------
        ! Convert TAI to UTC time
        ! -----------------------
        CALL convert_tai_to_utc ( &
             nUTCdim, omi_time(iloop), omi_time_utc(1:nUTCdim,iloop) )

     END DO ScanLineBlock

     ! ----------------------------------------------------------------
     ! AMF calculation and update of fitting statistics only need to be
     ! done for the final round throught the common mode iteration loop
     ! ----------------------------------------------------------------
     IF ( .NOT. yn_commit ) THEN

        CALL he5_write_radfit_output (                            &
             pge_idx, iline, nx, nblock, fpix, lpix,              &
             all_fitted_columns (1:n_fitvar_rad,1:nx,0:nblock-1), &
             all_fitted_errors  (1:n_fitvar_rad,1:nx,0:nblock-1), &
             correlation_columns(1:n_fitvar_rad,1:nx,0:nblock-1), &
!!$             omi_fitspc,nt,
             locerrstat )
        errstat = MAX ( errstat, locerrstat )

     END IF

  END DO ScanLines

  ! -----------------------------------------
  ! Remove target gas from radiance reference
  ! -----------------------------------------
  IF ( yn_radiance_reference .AND. yn_remove_target ) THEN

     ! -----------------------------------------------
     ! Removing the target from the radiance reference
     ! -----------------------------------------------
     IF ( yn_remove_target ) THEN

        WHERE ( targcnt > 0.0_r8 )
           targsum = targsum / targcnt
        ELSEWHERE
           targsum = r8_missval
        END WHERE

        ! ----------------------------------------------------------------
        ! Removing the target gas from the radiance reference will alter
        ! OMI_RADREF_SPEC (1:NWVL,FPIX:LPIX). This is being passed to the
        ! subroutine via MODULE use rather than through the argument list.
        ! ----------------------------------------------------------------
        CALL remove_target_from_radiance (                              &
             nccd, fpix, lpix, n_fincol_idx, fincol_idx(1:2,1:n_fincol_idx),  &
             ntargpol, targsum(1:n_fincol_idx,fpix:lpix), target_fit(fpix:lpix) )
     END IF

  END IF

  RETURN

END SUBROUTINE omi_pge_swathline_loop_memory
