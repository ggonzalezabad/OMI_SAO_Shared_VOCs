MODULE OMSAO_destriping_module

  USE OMSAO_precision_module
  USE OMSAO_indices_module,    ONLY: pge_bro_idx, pge_o3_idx, xtrcor_didx
  USE OMSAO_parameters_module, ONLY: &
       r8_missval, normweight, downweight, elsunc_less_is_noise, elsunc_infloop_eval
  USE OMSAO_omidata_module,    ONLY: nxtrack_max, nlines_max
  USE OMSAO_variables_module,  ONLY: &
       radfit_latrange, num_fitfunc_calls, num_fitfunc_jacobi, yn_diagnostic_run
  USE OMSAO_median_module,     ONLY: median
  USE OMSAO_he5_module
  USE OMSAO_errstat_module

  IMPLICIT NONE


  ! ----------------------------------------------------------------------
  ! Variables local to this module but common to many of the routines here
  ! ----------------------------------------------------------------------
  REAL    (KIND=r8), PARAMETER :: max_corr_fac = 2.0_r8, max_corr_fac_inv = 0.5_r8

  ! ---------------------------------------------------------
  ! Common destriping arrays (passed to the fitting function)
  ! ---------------------------------------------------------
  REAL (KIND=r8), DIMENSION (nxtrack_max) :: &
       xtrack_striping_pos, xtrack_striping_col, xtrack_striping_wgt, xtrack_striping_pat, &
       xtrack_poly_pos, xtrack_poly_val, xtrack_poly_wgt

  ! -----------------------------------------
  ! Variables related to cross-track bias
  ! -----------------------------------------
  LOGICAL           :: yn_run_destriping
  LOGICAL           :: yn_remove_ctrbias
  INTEGER (KIND=i4) :: ctr_bias_pol

  ! -----------------------------------------
  ! Order of Baseline and Scaling Polynomials
  ! -----------------------------------------
  INTEGER (KIND=i4) :: ctr_pol_base, ctr_pol_scal, ctr_pol_patt
  INTEGER (KIND=i4) :: ctr_nloop, ctr_nblocks, ctr_fitfunc_calls

  ! --------------------------------------------------------------
  ! Latitude limits for the computation of the cross-track average
  ! --------------------------------------------------------------
  REAL (KIND=r4), DIMENSION (2) :: ctrdst_latrange

  ! --------------------------------------------------------------
  ! Absolute maximum column amount to include in averaging
  ! --------------------------------------------------------------
  REAL (KIND=r8) :: ctr_maxcol

CONTAINS

  SUBROUTINE xtrack_destriping (                           &
       pge_idx, ntimes, nxtrack, yn_process, xtrange,      &
       lat, saocol, saodco, saoamf, saofcf, saomqf, errstat )

    IMPLICIT NONE

    ! ---------------
    ! Input variables
    ! ---------------
    INTEGER (KIND=i4),                                 INTENT (IN) :: pge_idx, ntimes, nxtrack
    LOGICAL,           DIMENSION (0:ntimes-1),         INTENT (IN) :: yn_process
    INTEGER (KIND=i4), DIMENSION (0:ntimes-1,1:2),     INTENT (IN) :: xtrange
    REAL    (KIND=r4), DIMENSION (nxtrack,0:ntimes-1), INTENT (IN) :: lat
    REAL    (KIND=r8), DIMENSION (nxtrack,0:ntimes-1), INTENT (IN) :: saodco, saoamf
    INTEGER (KIND=i2), DIMENSION (nxtrack,0:ntimes-1), INTENT (IN) :: saofcf, saomqf

    ! ------------------
    ! Modified variables
    ! ------------------
    REAL    (KIND=r8), DIMENSION (nxtrack,0:ntimes-1), INTENT (INOUT) :: saocol
    INTEGER (KIND=i4),                                 INTENT (INOUT) :: errstat


    ! ---------------
    ! Local variables
    ! ---------------
    INTEGER (KIND=i4)                                 :: &
         i, it, k1, k2, iii, locerrstat, nl0, nl1, nbd2, midnum, ntmp
    REAL    (KIND=r8)                                 :: xtrack_norm, a_stripe
    REAL    (KIND=r8), DIMENSION (nxtrack)            :: xtrack_avg_norm, xtrack_pfit
    REAL    (KIND=r8), DIMENSION (nxtrack)            :: xtrack_cnt
    REAL    (KIND=r8), DIMENSION (0:ntimes-1)         :: xtrack_fit
    REAL    (KIND=r8), DIMENSION (nxtrack,0:ntimes-1) :: saodst, tmpcol, xtrack_cor
    INTEGER (KIND=i2), DIMENSION (nxtrack,0:ntimes-1) :: tmpmqf
    INTEGER (KIND=i4), DIMENSION (2)                  :: latitude_lines, avg_limits
    REAL    (KIND=r8)                                 :: tmp_norm
    REAL    (KIND=r4)                                 :: midlat
    REAL    (KIND=r8)                                 :: xtr_median_med
    REAL    (KIND=r8), DIMENSION (nxtrack)            :: xtr_median, xtr_weight
    LOGICAL,           DIMENSION (0:ntimes-1)         :: yn_dst_range

    locerrstat = pge_errstat_ok

    ! ---------------------------------------------------------------
    ! Initialize the destriped columns with the regular ones, just in
    ! case we miss processing some lines or pixels.
    ! ---------------------------------------------------------------
    saodst(1:nXtrack,0:nTimes-1) = saocol(1:nXtrack,0:nTimes-1)

    ! -------------------------------------------------------------------
    ! YN_RUN_DESTRIPING=F, negative order of polynomial, or negative
    ! number of iteratrions means that we can not/should not perform any
    ! destriping.
    ! -------------------------------------------------------------------
    IF ( ( .NOT. yn_run_destriping               ) .OR. &
         ( MAX( ctr_pol_base, ctr_pol_scal ) < 0 ) .OR. &
         ( ctr_nloop < 0                         ) ) THEN
       xtrack_cor(1:nxtrack,0:nTimes-1) = r8_missval
       xtrack_fit(0:nTimes-1)           = r8_missval
       GOTO 666
    END IF
   
    ! ------------------------------------------------------
    ! Determine any limiting range due to latitude selection
    ! ------------------------------------------------------
    yn_dst_range(0:ntimes-1) = .FALSE.
    CALL find_swathrange_by_latitude (                               &
         ntimes, nxtrack, ctrdst_latrange(1), ctrdst_latrange(2),    &
         lat(1:nxtrack,0:ntimes-1), xtrange(0:ntimes-1,1:2), midlat, &
         midnum, yn_dst_range(0:ntimes-1)                              )
   
    ! ----------------------------------------
    ! Set the range of swath lines to destripe
    ! ----------------------------------------
    nl0 = 0  ;  nl1 = nTimes-1

    ! ---------------------------------------------------------------
    ! Compute average cross-track pattern (difference of AVG and FIT)
    ! ---------------------------------------------------------------
    xtrack_striping_pos(1:nxtrack) = 0.1_r8 * &
         ( (/ (REAL(i,KIND=r8), i = 1, nxtrack) /) - REAL(nxtrack+1, KIND=r8)/2.0_r8 )
   

    ! ======
    ! FUDGE:
    ! ======
    IF ( ctr_nblocks >= 0 ) ctr_nblocks = -999
    k1 = -1 ; k2 = -2

    !k1 = avg_limits(1) ; k2 = avg_limits(2)

    ! ----------------------------------------
    ! Short-hand for half the number of blocks
    ! ----------------------------------------
    nbd2 = ctr_nblocks / 2

    ! ------------------------------------
    ! Set some variables to initial values
    ! ------------------------------------
    xtrack_cor = r8_missval ; xtrack_fit = r8_missval ; xtrack_striping_pat = 0.0_r8
   
    ! -----------------------------------
    ! Loop over scan lines for destriping
    ! -----------------------------------
    stripeloop: DO iii = 1, ctr_nloop
       ! --------------------------------------------------------------
       ! Copy the freshly destriped columns back on the original array.
       ! The initial "SAODST = SAOCOL" above assures that SAODST has
       ! proper values the first time around the loop.
       ! --------------------------------------------------------------
       saocol = saodst

       ! --------------------------------------------------------------
       ! Unless we are doing dynamic, block-based destriping, we don't
       ! have to compute the average for each and every scan line. Once
       ! before entering the line loop is sufficient.
       ! --------------------------------------------------------------
       IF ( ctr_nblocks < 0 ) THEN
       
          ! -------------------------------------------------------------
          ! For each cross-track position, compute the along-track median
          ! over the range of latitudes selected.
          ! -------------------------------------------------------------
          ! First condense the selected part for computation of median
          ! ----------------------------------------------------------
          CALL condense_columns_for_median (                                      &
               ntimes, nxtrack, yn_process(0:ntimes-1), yn_dst_range(0:ntimes-1), &
               saocol(1:nxtrack,0:ntimes-1), saomqf(1:nxtrack,0:ntimes-1),        &
               ntmp, tmpcol(1:nxtrack,0:ntimes-1), tmpmqf(1:nxtrack,0:ntimes-1)     )

          CALL xtrack_median_comp ( &
               ntmp, nxtrack, tmpcol(1:nxtrack,1:ntmp), tmpmqf(1:nxtrack,1:ntmp), &
               xtr_median(1:nxtrack), xtr_median_med, xtr_weight(1:nxtrack) )
          
          ! ------------------------------------------------------------------
          ! Compute the cross-track striping pattern: Median normalized to One
          ! ------------------------------------------------------------------
          IF ( (xtr_median_med /= 0.0_r8) .AND. (ANY(xtr_median /= 0.0_r8) ) ) THEN
             xtrack_avg_norm = xtr_median / xtr_median_med
          ELSE
             xtrack_avg_norm = 1.0_r8
          END IF
          tmp_norm = SUM(xtrack_avg_norm(1:nxtrack) * &
               xtr_weight(1:nxtrack)) /SUM(xtr_weight(1:nxtrack) )
          xtrack_striping_pat = xtrack_avg_norm / tmp_norm
                    
          CALL xtrack_pattern_polyfit (                                    &
               it, nxtrack, ctr_pol_patt,                                  &
               xtrack_avg_norm(1:nxtrack), xtrack_striping_pos(1:nxtrack), &
               xtr_weight(1:nxtrack), xtrack_striping_pat(1:nxtrack),      &
               xtrack_pfit(1:nxtrack)    )
          
       END IF

       DO it = nl0, nl1

          !DeadBranch: ! ------------------------------------------------------
          !DeadBranch: ! If CTR_NBLOCKS is 0 or greater, we are doing a dynamic 
          !DeadBranch: ! block-based destriping. Else we are using the fixed
          !DeadBranch: ! latitude limits set in the Fitting Control File. For
          !DeadBranch: ! the block-based destriping we need to compute the
          !DeadBranch: ! along-track average for each and every scan line.
          !DeadBranch: ! ------------------------------------------------------
          !DeadBranch: IF ( ctr_nblocks >= 0 ) THEN
          !DeadBranch:    IF      ( it <= nl0+(nbd2-1)                       ) THEN
          !DeadBranch:       k1 = nl0
          !DeadBranch:    ELSE IF ( it >  nl0+(nbd2-1) .AND. it <= nl1-nbd2  ) THEN
          !DeadBranch:       k1 = it-nbd2
          !DeadBranch:    ELSE IF ( it >  nl1-nbd2                           ) THEN
          !DeadBranch:       k1 = nl1-ctr_nblocks
          !DeadBranch:    END IF
          !DeadBranch:    k2 = k1 + ctr_nblocks
          !DeadBranch: 
          !DeadBranch: 
          !DeadBranch:    CALL xtrack_median_comp ( &
          !DeadBranch:         k2-k1+1, nxtrack, saocol(1:nxtrack,k1:k2), saomqf(1:nxtrack,k1:k2), &
          !DeadBranch:         xtr_median(1:nxtrack), xtr_median_med, xtr_weight(1:nxtrack) )
          !DeadBranch: 
          !DeadBranch:    ! ----------------------------------------
          !DeadBranch:    ! Normalize the cross-track stripe pattern
          !DeadBranch:    ! ----------------------------------------
          !DeadBranch:    xtrack_avg_norm = xtr_median / xtr_median_med
          !DeadBranch:    tmp_norm = SUM(xtrack_avg_norm(1:nxtrack) * &
          !DeadBranch:         xtr_weight(1:nxtrack)) /SUM(xtr_weight(1:nxtrack) )
          !DeadBranch:    xtrack_striping_pat = xtrack_avg_norm / tmp_norm
          !DeadBranch: 
          !DeadBranch:    CALL xtrack_pattern_polyfit (                                    &
          !DeadBranch:         it, nxtrack, ctr_pol_patt,                                  &
          !DeadBranch:         xtrack_avg_norm(1:nxtrack), xtrack_striping_pos(1:nxtrack), &
          !DeadBranch:         xtr_weight(1:nxtrack), xtrack_striping_pat(1:nxtrack),      &
          !DeadBranch:         xtrack_pfit(1:nxtrack)  )
          !DeadBranch: 
          !DeadBranch: ELSE
          !DeadBranch:    k1 = avg_limits(1) ; k2 = avg_limits(2)
          !DeadBranch: END IF

          ! -----------------------
          ! Copy the actual columns
          ! -----------------------
          xtrack_striping_col(1:nxtrack) = saocol(1:nxtrack,it)

          ! -----------------
          ! Check for weights
          ! -----------------
          WHERE ( saomqf (1:nxtrack,it) == 0_i2  )
             xtrack_striping_wgt(1:nxtrack) = normweight
             xtrack_cnt         (1:nxtrack) = 1.0_r8
          ELSEWHERE
             xtrack_striping_wgt(1:nxtrack) = downweight
             xtrack_cnt         (1:nxtrack) = 0.0_r8
          END WHERE
          xtrack_striping_wgt(1:1) = downweight
     
          xtrack_norm = 0.0_r8
          IF ( yn_process(it) .AND. SUM(xtrack_cnt(1:nxtrack)) > 0.0_r8 ) THEN
          
             ! ------------------------
             ! Compute cross-track norm
             ! ------------------------
             xtrack_norm = SUM ( xtrack_striping_col(1:nxtrack)*xtrack_cnt(1:nxtrack) ) / &
                  SUM( xtrack_cnt(1:nxtrack) )
          
             ! --------------------------------------------------------------------------
             ! Fit for contribution of the cross-track striping in the current swath line
             ! --------------------------------------------------------------------------
             num_fitfunc_calls = 0 ; num_fitfunc_jacobi = 0
             CALL xtrack_striping_fit (                                &
                  nxtrack, ctr_pol_base, ctr_pol_scal, xtrack_norm, &
                  a_stripe, xtrack_cor(1:nxtrack,it) )
          
             ! ----------------
             ! Apply correction
             ! ----------------
             saodst(1:nxtrack,it) = saocol(1:nxtrack,it) - xtrack_cor(1:nxtrack,it)
             xtrack_fit (it)      = a_stripe

          ELSE
          
             ! -----------------------------------------------------------------
             ! If we don't have anything to fit, assign MISSING VALUEs to the
             ! cross-track variables, and copy over the original columns to
             ! the destriped ones (cosmetic reasons only: plots better that way)
             ! -----------------------------------------------------------------
             saodst(1:nXtrack,it)      = saocol(1:nXtrack,it)
             xtrack_cor (1:nXtrack,it) = r8_missval
             xtrack_fit (it)           = r8_missval
             
          END IF
             
       END DO

    END DO stripeloop


    ! -----------------------------------------------------------
    ! After the stripe correction we correct for any cross-track bias.
    ! This is done either via a fitted cross-track polynomial or simply
    ! a median correction.
    ! -----------------------------------------------------------
    IF ( yn_remove_ctrbias ) THEN

       ! ----------------------------------------------------------
       ! First condense the selected part for computation of median
       ! ----------------------------------------------------------
       CALL condense_columns_for_median (                                      &
            ntimes, nxtrack, yn_process(0:ntimes-1), yn_dst_range(0:ntimes-1), &
            saodst(1:nxtrack,0:ntimes-1), saomqf(1:nxtrack,0:ntimes-1),        &
            ntmp, tmpcol(1:nxtrack,0:ntimes-1), tmpmqf(1:nxtrack,0:ntimes-1)     )
       CALL xtrack_median_comp ( &
            ntmp, nxtrack, tmpcol(1:nxtrack,1:ntmp), tmpmqf(1:nxtrack,1:ntmp), &
            xtr_median(1:nxtrack), xtr_median_med, xtr_weight(1:nxtrack) )
       IF ( xtr_median_med /= 0.0_r8 ) xtr_median = xtr_median / xtr_median_med
       
       IF ( ctr_bias_pol < 1 ) THEN

          WHERE ( xtr_median /= 0.0_r8 )
             xtr_median = 1.0_r8 / xtr_median
          ELSEWHERE
             xtr_median = 1.0_r8
          END WHERE
          DO iii = 0, ntimes-1
             IF ( .NOT. yn_process(iii) ) CYCLE
             xtrack_cor(1:nxtrack,iii) = xtr_median(1:nxtrack)
             xtrack_fit(          iii) = 0.0_r8
             WHERE ( saocol(1:nxtrack,iii) > r8_missval )
                saodst(1:nxtrack,iii) = saodst(1:nxtrack,iii) * xtr_median(1:nxtrack)
             ENDWHERE
          END DO
       ELSE

          ! ------------------------------------------------------------------
          ! Compute the cross-track striping pattern: Median normalized to One
          ! ------------------------------------------------------------------
          tmp_norm = SUM( xtr_median(1:nxtrack)*xtr_weight(1:nxtrack) ) / &
               SUM(xtr_weight(1:nxtrack) )

          xtrack_striping_pat = xtr_median / tmp_norm

          CALL xtrack_pattern_polyfit (                               &
               it, nxtrack, ctr_bias_pol,                             &
               xtr_median(1:nxtrack), xtrack_striping_pos(1:nxtrack), &
               xtr_weight(1:nxtrack), xtrack_striping_pat(1:nxtrack), &
               xtrack_pfit(1:nxtrack)    )          

          ! ------------------------------------------------------------
          ! Remove the low-order polynomial bias from the fitted columns
          ! ------------------------------------------------------------
          DO it = nl0, nl1
             IF ( .NOT. yn_process(it) ) CYCLE
             saocol(1:nxtrack,it) = saocol(1:nxtrack,it) - xtrack_pfit(1:nxtrack)
          END DO


          ! ----------------------------------------------------------
          ! First condense the selected part for computation of median
          ! ----------------------------------------------------------
          CALL condense_columns_for_median (                                      &
               ntimes, nxtrack, yn_process(0:ntimes-1), yn_dst_range(0:ntimes-1), &
               saocol(1:nxtrack,0:ntimes-1), saomqf(1:nxtrack,0:ntimes-1),        &
               ntmp, tmpcol(1:nxtrack,0:ntimes-1), tmpmqf(1:nxtrack,0:ntimes-1)     )
          CALL xtrack_median_comp ( &
               ntmp, nxtrack, tmpcol(1:nxtrack,1:ntmp), tmpmqf(1:nxtrack,1:ntmp), &
               xtr_median(1:nxtrack), xtr_median_med, xtr_weight(1:nxtrack) )
       
          ! ----------------------------------------
          ! Normalize the cross-track stripe pattern
          ! ----------------------------------------
          xtrack_avg_norm = xtr_median / xtr_median_med
          tmp_norm = SUM(xtrack_avg_norm(1:nxtrack) * &
               xtr_weight(1:nxtrack)) /SUM(xtr_weight(1:nxtrack) )
          xtrack_striping_pat = xtrack_avg_norm / tmp_norm

       END IF   ! on order of CTR Bial Polynomial

    END IF   ! on CTR Bias Removal

    ! -----------------------------------------------------------------
    ! Write the correction and the corrected columns to the output file
    ! -----------------------------------------------------------------
    CALL xtrack_destriping_writecol (                            &
         ntimes, nxtrack, saodst(1:nxtrack,0:ntimes-1),          &
         xtrack_cor(1:nxtrack,0:nTimes-1), xtrack_fit(0:nTimes-1)  )

666 CONTINUE

    RETURN
  END SUBROUTINE xtrack_destriping


  SUBROUTINE xtrack_median_comp ( &
       nt, nx, saocol, saomqf, xtrack_med, xtrack_med_med, xtrack_wgt )

    IMPLICIT NONE

    ! ---------------
    ! Input variables
    ! ---------------
    INTEGER (KIND=i4),                        INTENT (IN) :: nt, nx
    REAL    (KIND=r8), DIMENSION (nx,0:nt-1), INTENT (IN) :: saocol
    INTEGER (KIND=i2), DIMENSION (nx,0:nt-1), INTENT (IN) :: saomqf

    ! ----------------
    ! Output variables
    ! ----------------
    REAL    (KIND=r8),                 INTENT (OUT) :: xtrack_med_med
    REAL    (KIND=r8), DIMENSION (nx), INTENT (OUT) :: xtrack_med, xtrack_wgt

    ! ---------------
    ! Local variables
    ! ---------------
    INTEGER (KIND=i4) ::                     ix, it, ngood
    REAL    (KIND=r8), DIMENSION (0:nt-1) :: saocol_good

    xtrack_med = 0.0_r8 ; xtrack_med_med = 0.0_r8 ; xtrack_wgt = downweight

    ! -------------------------------------------
    ! Compute MEDIAN, but for "good" columns only
    ! -------------------------------------------
    DO ix = 1, nx
       saocol_good = 0.0_r8 ; ngood = 0

       DO it = 0, nt-1
          IF ( saomqf(ix,it) == 0 ) THEN
             ngood = ngood + 1
             saocol_good(ngood-1) = saocol(ix,it)
          END IF
       END DO

       IF ( ngood > 0 ) THEN
          xtrack_med(ix) = median( ngood, saocol_good(0:ngood-1) )
          xtrack_wgt(ix) = normweight
       ELSE
          xtrack_med(ix) = 0.0_r8
          xtrack_wgt(ix) = downweight
       END IF
    END DO

    saocol_good = 0.0_r8 ; ngood = 0
    DO ix = 1, nx
       IF ( xtrack_wgt(ix) == normweight ) THEN
          ngood = ngood + 1
          saocol_good(ngood-1) = xtrack_med(ix)
       END IF
    END DO
    IF ( ngood > 0 ) THEN
       xtrack_med_med = median ( ngood, saocol_good(0:ngood-1) )
    END IF

    RETURN
  END SUBROUTINE xtrack_median_comp



  SUBROUTINE xtrack_destriping_lat_limits (             &
       ntimes, fline, lline, nxtrack, lat, nblocks,     &
       latitude_lim, latitude_lines, avg_limits, errstat )

    IMPLICIT NONE

    ! ---------------
    ! Input variables
    ! ---------------
    INTEGER (KIND=i4),                INTENT (IN) :: ntimes, nxtrack, nblocks, fline, lline
    REAL    (KIND=r4), DIMENSION (2), INTENT (IN) :: latitude_lim
    REAL    (KIND=r4), DIMENSION (nxtrack,0:ntimes-1), INTENT (IN) :: lat

    ! -----------------
    ! Modified variable
    ! -----------------
    INTEGER (KIND=i4), INTENT (INOUT) :: errstat

    ! ----------------
    ! Output variables
    ! ----------------
    INTEGER (KIND=i4), DIMENSION (2), INTENT (OUT) :: latitude_lines, avg_limits


    ! ---------------
    ! Local variables
    ! ---------------
    INTEGER (KIND=i4) :: locerrstat, j, k
    INTEGER (KIND=i4) :: nxh, nth, lower_line, upper_line
    REAL    (KIND=r4) :: lat1


    ! -------------------------------------
    ! Initialize local and output variables
    ! -------------------------------------
    latitude_lines = -1 ; avg_limits = -1 ; locerrstat = pge_errstat_ok

    ! -----------------------------------------------------
    ! Check whether we actually have to do any processing.
    ! If first and last line don't differe by more than one
    ! line, we can easily set the return values without 
    ! checking in detail. Actually, this should be true for
    ! more than just 1 line difference, but for 2 or more
    ! the loops further down can handle all cases.
    ! -----------------------------------------------------
    IF ( fline >= lline-1 ) THEN
       latitude_lines(1) = fline
       latitude_lines(2) = lline
       avg_limits(1)     = fline
       avg_limits(2)     = lline
       RETURN
    END IF
       
    nxh = nxtrack/2 ; nth = (fline+lline)/2
    ! --------------------------------------------
    ! Find points of latitude wrap-around (if any)
    ! --------------------------------------------
    upper_line = -1 ; lat1 = lat(nxh,nth)
    LatUpper: DO j = nth, lline-1 !nth+1, ntimes-2
       IF ( ALL(lat(nxh,j:j+1) >= lat1) ) THEN
          upper_line = j+1 ; lat1 = lat(nxh,j+1)
       ELSE
          EXIT LatUpper
       END IF
    END DO LatUpper

    lower_line = -1 ; lat1 = lat(nxh,nth)
    LatLower: DO j = nth, fline+1, -1  !nth-1, 1, -1
       IF ( ALL(lat(nxh,j-1:j) <= lat1) ) THEN
          lower_line = j-1 ; lat1 = lat(nxh,j-1)
       ELSE
          EXIT LatLower
       END IF
    END DO LatLower

    GetLatLines1: DO j = lower_line, upper_line
       IF ( latitude_lines(2) == -1 .AND. ALL(lat(1:nxtrack,j) > latitude_lim(2)) ) THEN
          latitude_lines(2) = MAX(j-1,lower_line)
          EXIT GetLatLines1
       END IF
    END DO GetLatLines1

    GetLatLines2: DO j = upper_line, lower_line, -1
       IF ( latitude_lines(1) == -1 .AND. ALL(lat(1:nxtrack,j) < latitude_lim(1)) ) THEN
          latitude_lines(1) = MIN(j+1,upper_line)
          EXIT GetLatLines2
       END IF
    END DO GetLatLines2

    IF ( latitude_lines(1) == -1 ) latitude_lines(1) = lower_line
    IF ( latitude_lines(2) == -1 ) latitude_lines(2) = upper_line

    ! --------------------------------------------------------
    ! Find the range of swath lines that go into the averaging
    ! --------------------------------------------------------
    IF ( nblocks < 0 ) THEN
       GetAvgLines: DO j = nth+1, ntimes-2
          k = ntimes-1 - j
          IF ( ALL(lat(nxh,k-1:k) < MINVAL(ctrdst_latrange)) .AND. &
               avg_limits(1) == -1 ) avg_limits(1) = k + 1
          IF ( ALL(lat(nxh,j:j+1) > MAXVAL(ctrdst_latrange)) .AND. &
               avg_limits(2) == -1 ) avg_limits(2) = j - 1
          IF ( MINVAL(avg_limits) > -1 ) EXIT GetAvgLines
       END DO GetAvgLines
    END IF

    errstat = MAX ( errstat, locerrstat )

    RETURN
  END SUBROUTINE xtrack_destriping_lat_limits
    

  SUBROUTINE xtrack_pattern_polyfit ( &
       is_line, nxtrack, npol, xtrack_avg, xtrack_pos, xtrack_wgt, &
       xtrack_cor, xtrack_pfit )

    IMPLICIT NONE

    ! ---------------
    ! Input variables
    ! ---------------
    INTEGER (KIND=i4),                      INTENT (IN)    :: is_line, nxtrack, npol
    REAL    (KIND=r8), DIMENSION (nxtrack), INTENT (IN)    :: xtrack_avg, xtrack_pos, xtrack_wgt

    ! ---------------
    ! Output variable
    ! ---------------
    REAL    (KIND=r8), DIMENSION (nxtrack), INTENT(OUT) :: xtrack_cor, xtrack_pfit

    ! ----------------
    ! DPOLFt variables
    ! ----------------
    INTEGER (KIND=i4) :: ndeg, ierr
    REAL    (KIND=r8) :: eps
    REAL    (KIND=r8), DIMENSION (3*(nxtrack+npol+1)) :: a

    xtrack_cor (1:nxtrack) = 0.0_r8
    xtrack_pfit(1:nxtrack) = 0.0_r8

    eps =  0.0_r8  ! Fit the complete NPOL polynomial
    !eps = -1.0_r8  ! Choose the best-fitting order
    CALL dpolft (&
         nxtrack, xtrack_pos(1:nxtrack), xtrack_avg(1:nxtrack), xtrack_wgt(1:nxtrack), &
         npol, ndeg, eps, xtrack_pfit(1:nxtrack), ierr, a )

    ! ------------------------------------------
    ! Compute the cross-track correction factor;
    ! ------------------------------------------
    xtrack_cor(1:nxtrack) = xtrack_avg(1:nxtrack) - xtrack_pfit(1:nxtrack)
   
    RETURN
  END SUBROUTINE xtrack_pattern_polyfit


  SUBROUTINE xtrack_striping_fit ( nxtrack, npolb, npols, xnorm, a_stripe, xtrack_cor )

    USE OMSAO_parameters_module,     ONLY: elsunc_np, elsunc_nw
    USE OMSAO_elsunc_fitting_module, ONLY: elsunc

    IMPLICIT NONE

    ! ---------------
    ! Input variables
    ! ---------------
    INTEGER (KIND=i4), INTENT (IN) :: nxtrack, npolb, npols
    REAL    (KIND=r8), INTENT (IN) :: xnorm

    ! ---------------
    ! Output variable
    ! ---------------
    REAL (KIND=r8),                     INTENT (OUT) :: a_stripe
    REAL (KIND=r8), DIMENSION(nxtrack), INTENT (OUT) :: xtrack_cor

    ! ---------------
    ! Local variables
    ! ---------------
    INTEGER (KIND=i4) :: i, nfit, ipar
    REAL    (KIND=r8) :: chisq

    ! ----------------
    ! ELSUNC variables
    ! ----------------
    ! ------------------------------------------------------------------
    ! Note that we add "4" to NPOL. This is to prevent ELSUNC choking on
    ! not having a minimum of 4 parameters (NPOL=0), which seems to be a
    ! requirement somewhere inside the fitting code. For less than 4 
    ! parameters it is possible that arrays run out of bounds.
    ! ------------------------------------------------------------------
    INTEGER (KIND=i4)                                    :: exval, ctrl
    INTEGER (KIND=i4)                                    :: elbnd
    INTEGER (KIND=i4), DIMENSION (elsunc_np)             :: p
    REAL    (KIND=r8), DIMENSION (elsunc_nw)             :: w
    REAL    (KIND=r8), DIMENSION (npolb+npols+4)         :: blow, bupp
    REAL    (KIND=r8), DIMENSION (nxtrack)               :: f
    REAL    (KIND=r8), DIMENSION (nxtrack,npolb+npols+4) :: dfda
    REAL    (KIND=r8), DIMENSION (npolb+npols+4)         :: fitpar

    ! -------------------------------------------------------------
    ! Set the number of fitting parameters: Order of the baseline
    ! plus scaling polynomial, plus the XTrack Pattern.
    ! -------------------------------------------------------------
    nfit = (npolb+1) + (npols+1) + 1

    ! -------------------------------------------------------------
    ! Set ELSUNC fitting paramerters for unconstrained fit
    ! -------------------------------------------------------------
    ! ELBND: 0 = unconstrained
    !        1 = all variables have same lower bound
    !        else: lower and upper bounds must be supplied
    ! -------------------------------------------------------------
    elbnd = 2  ;  exval = 0 ; p = -1 ; p(1) = 0 ; p(3) = ctr_fitfunc_calls
    w = -1.0
    blow(1:nfit) = 0.0_r8  ;  bupp(1:nfit) = 0.0_r8

    ! ---------------------------------------------------------
    ! Initial guess for fitting parameters: For the polynomials
    ! we use alternating coefficients  of decreasing order of
    ! magnitude (this is just a test to check whether we get 
    ! anywhere). For the XTrack Pattern we use 0.1*xnorm.
    ! ---------------------------------------------------------    
    ! First the XTrack Pattern:
    ! -------------------------
    ipar = 1
    fitpar(ipar) =  0.10_r8*xnorm
    blow  (ipar) = -10.0_r8*ABS(xnorm) 
    bupp  (ipar) = +10.0_r8*ABS(xnorm) 
    ! -----------------------------------------------------------
    ! Now NPOLB baseline coefficients. The zero order is set to
    ! XNORM. We could place this inside the assignement loop,
    ! of course, but it is a bit more intuitive to do it outside.
    ! -----------------------------------------------------------
    ipar = ipar + 1
    fitpar(ipar) = xnorm
    blow  (ipar) = -10.0_r8*ABS(xnorm) 
    bupp  (ipar) = +10.0_r8*ABS(xnorm) 
    DO i = 1, npolb
       ipar = ipar + 1
       fitpar(ipar) = xnorm*(-0.1)**i
       blow  (ipar) = -ABS(xnorm)
       bupp  (ipar) = +ABS(xnorm)
    END DO
    ! ----------------------------------------------
    ! Now NPOLS scaling coefficients. The zero order
    ! coefficient is fixed at 1 in order to prevent
    ! strong correlations with the baseline offset.
    ! ----------------------------------------------
    ipar = ipar + 1
    fitpar(ipar) = 1.0_r8
    blow  (ipar) = 1.0_r8
    bupp  (ipar) = 1.0_r8
    DO i = 1, npols
       ipar = ipar + 1
       fitpar(ipar) = (-0.1)**(i-1)
       blow  (ipar) = -ABS(xnorm)
       bupp  (ipar) = +ABS(xnorm)
    END DO

    ! ---------------------------------------------------------------
    ! Here we check for the number of fitting parameters. IF < 4 then
    ! we need to change tack and add the required number. Since those
    ! are dummy parameters, we have to make sure that they stay Zero
    ! by specifying lower and upper bounds.
    ! ---------------------------------------------------------------
    IF ( nfit < 4 ) THEN

       ! -------------------------------------------------
       ! Set dummy fitting parameters and their range to 0
       ! -------------------------------------------------
       fitpar(nfit+1:4) = 0.0_r8 ; blow(nfit+1:4) = 0.0_r8 ; bupp(nfit+1:4) = 0.0_r8

       ! ------------------------------------------------
       ! Now we set the number of fitting parameters to 4
       ! ------------------------------------------------
       nfit = 4
    END IF

    CALL elsunc ( &
         fitpar(1:nfit), nfit, nxtrack, nxtrack, xtrack_striping_func, &
         elbnd, blow(1:nfit), bupp(1:nfit), p, w, exval, f(1:nxtrack), &
         dfda(1:nxtrack,1:nfit) )
    chisq = SUM  ( f(1:nxtrack)**2 ) ! This gives the same CHI**2 as the NR routines

    a_stripe = fitpar(1)

    ! ----------------------------------------------------------
    ! Calculate the cross-track stripe contribution. This is
    ! the stripe pattern multiplied by any scaling polynomial.
    ! We can use the fitting function to do this for us, but we
    ! have to make sure that all the parameters for the baseline
    ! polynomial are set to zero so that we just consider the
    ! stripe pattern and the scaling polynomial.
    ! ----------------------------------------------------------
    ctrl = 3
    fitpar(2:2+npolb) = 0.0_r8
    CALL xtrack_striping_func ( &
         fitpar(1:nfit), nfit, xtrack_cor(1:nxtrack), nxtrack, ctrl, &
         dfda(1:nxtrack,1:nfit), 0 )
    RETURN
  END SUBROUTINE xtrack_striping_fit


  SUBROUTINE xtrack_striping_func ( a, na, y, m, ctrl, dyda, mdy )

    IMPLICIT NONE
    
    ! ----------------
    ! Input parameters
    ! ----------------
    INTEGER (KIND=i4),                  INTENT (IN)  :: na, m, mdy
    REAL    (KIND=r8), DIMENSION (na),  INTENT (IN)  :: a

    ! -------------------
    ! Modified parameters
    ! -------------------
    INTEGER (KIND=i4), INTENT (INOUT) :: ctrl

    ! ----------------
    ! Output parameters
    ! ----------------
    REAL (KIND=r8), DIMENSION (m),    INTENT (OUT)  :: y
    REAL (KIND=r8), DIMENSION (m,na), INTENT (OUT)  :: dyda

    ! ----------------
    ! Local variables
    ! ----------------
    INTEGER (KIND=i4)                :: i, ipar
    REAL    (KIND=r8), DIMENSION (m) :: x, y0, xpow, scpol, blpol, xtr

    x   (1:m) = xtrack_striping_pos(1:m)
    xpow(1:m) = 1.0_r8
    
    ! -------------------------------------------------------------
    ! Compose the cross-track spectrum:
    !
    ! ---------    ----------
    ! Parameter    Represents
    ! ---------    ----------
    ! 1            Cross-Track Stripe Pattern
    ! 2  :i        Baseline Polynomial (of order i-2)
    ! i+1:k        Scaling  Polynmial  (of order k-i-1)
    ! -------------------------------------------------------------

    ! -------------------------
    ! First the XTrack Pattern:
    ! -------------------------
    ipar = 1
    xtr(1:m) = a(ipar)*xtrack_striping_pat(1:m)
    y0 (1:m) = xtr(1:m)

    ! --------------------------------
    ! Now add the baseline polynomial:
    ! --------------------------------
    ipar = ipar + 1
    blpol(1:m) = a(ipar)
    DO i = 1, ctr_pol_base
       ipar = ipar + 1
       blpol(1:m) = blpol(1:m) + a(ipar)*x(1:m)**i
    END DO
    y0(1:m) = y0(1:m) + blpol(1:m)

    ! ----------------------------------------------
    ! Now NPOLS scaling coefficients. The zero order
    ! coefficient is fixed at 1 in order to prevent
    ! strong correlations with the baseline offset.
    ! Hence we can se
    ! ----------------------------------------------
    ipar = ipar + 1
    scpol(1:m) = a(ipar)
    DO i = 1, ctr_pol_scal
       ipar = ipar + 1
       scpol(1:m) = scpol(1:m) + a(ipar)*x(1:m)**i
    END DO
    y0(1:m) = y0(1:m) * scpol(1:m)

    SELECT CASE ( ABS(ctrl) )
    CASE ( 1 )
       ! -------------------------------------------------
       ! Count the number of calls to the fitting function
       ! and terminate if we exceed the allowed maximum.
       ! -------------------------------------------------
       num_fitfunc_calls = num_fitfunc_calls + 1
       IF ( num_fitfunc_calls > ctr_fitfunc_calls ) THEN
          ctrl = INT(elsunc_infloop_eval, KIND=i4)
          RETURN
       END IF

       ! ------------------------------------------
       ! Return the residual between data and model
       ! ------------------------------------------
       y(1:m)  = ( y0(1:m) - xtrack_striping_col(1:m) ) * xtrack_striping_wgt(1:m)
    CASE ( 2 )
       ! -------------------------------------------------
       ! Count the number of calls to the fitting function
       ! with request for the Jacobian and terminate if we
       ! exceed the allowed maximum (just to be safe!).
       ! -------------------------------------------------
       num_fitfunc_jacobi = num_fitfunc_jacobi + 1
       IF ( num_fitfunc_jacobi > ctr_fitfunc_calls ) THEN
          ctrl = INT(elsunc_infloop_eval, KIND=i4)
          RETURN
       END IF

       ! -------------------------------------------
       ! Compute the Jacobian (we know the function)
       ! -------------------------------------------
       dyda(1:m,1:na) = 0.0_r8

       ! ----------------------------
       ! Cross-Track Stripe parameter
       ! ----------------------------
       ipar = 1
       dyda(1:m,ipar) = xtrack_striping_pat(1:m) * scpol(1:m)

       ! ------------------------------
       ! Baseline Polynomial parameters
       ! ------------------------------
       ipar = ipar + 1
       dyda(1:m,ipar) = a(ipar) * scpol(1:m)
       DO i = 1, ctr_pol_base
          ipar = ipar + 1
          dyda(1:m,ipar) = (a(ipar)*x(1:m)**i)*scpol(1:m)
       END DO

       ! -----------------------------
       ! Scaling Polynomial parameters
       ! -----------------------------
       blpol(1:m) = blpol(1:m) + xtr(1:m)
       ipar = ipar + 1
       dyda(1:m,ipar) = blpol(1:m) * a(ipar)
       DO i = 1, ctr_pol_scal
          ipar = ipar + 1
          dyda(1:m,ipar) = blpol(1:m) * (a(ipar)*x(1:m)**i)
       END DO

       ! --------------------------------------------------------------
       ! Any additional entries are taken care of by the initialization
       ! to 0 above. There could be two such entries since ELSUNC needs
       ! a minimum of 4 variables, and if both polynomials are of Zero
       ! order (and the scaling is fixed), this will lead to less than
       ! 4 fitting parameters.
       ! --------------------------------------------------------------

    CASE ( 3 )
       ! ---------------------------------------------------------
       ! This CASE is included to get the complete fitted spectrum
       ! ---------------------------------------------------------
       y(1:m)  = y0(1:m)
    CASE DEFAULT
       ! ----------------------------------------------------------
       ! This shouldn't happen, but just for the case we reach here
       ! ..... we do nothing
       ! ----------------------------------------------------------
       !WRITE (*, '(A,I3)') "Don't know how to handle CTRL = ", ctrl
    END SELECT

    RETURN
  END SUBROUTINE xtrack_striping_func

  SUBROUTINE xtrack_destriping_writecol (ntimes, nxtrack, saocol, xtrack_cor, xtrack_fit )

    IMPLICIT NONE

    ! ---------------
    ! Input variables
    ! ---------------
    INTEGER (KIND=i4),                                 INTENT (IN) :: ntimes, nxtrack
    REAL    (KIND=r8), DIMENSION (0:ntimes-1),         INTENT (IN) :: xtrack_fit
    REAL    (KIND=r8), DIMENSION (nxtrack,0:ntimes-1), INTENT (IN) :: saocol, xtrack_cor

    ! ---------------
    ! Local variables
    ! ---------------
    INTEGER (KIND=i4)            :: locerrstat, iline, ntimes_loop

    ! -------------------------------------------------------
    ! First write the corrected columns; loop over all lines
    ! in the file in blocks of NLINES_MAX
    ! -------------------------------------------------------
    DO iline = 0, ntimes-1, nlines_max

       ! --------------------------------------------------------
       ! Check if loop ends before n_times_loop max is exhausted.
       ! --------------------------------------------------------
       ntimes_loop = MIN ( nlines_max, ntimes-iline )

       ! --------------------------------------------
       ! Set Start, Stride, and Edge of writing block
       ! --------------------------------------------
       he5_start_2d  = (/       0,       iline /)
       he5_stride_2d = (/       1,           1 /)
       he5_edge_2d   = (/ nxtrack, ntimes_loop /)
    
       ! -------------
       ! Write columns
       ! -------------
       locerrstat = HE5_SWwrfld ( pge_swath_id, TRIM(ADJUSTL(dstrcol_field)),         &
            he5_start_2d, he5_stride_2d, he5_edge_2d, saocol(1:nxtrack,iline:iline+ntimes_loop-1) )

       ! ----------------------------------------------
       ! Write X-Track pattern ("diagnostic" runs only)
       ! ----------------------------------------------
       IF ( yn_diagnostic_run .AND. yn_output_diag(xtrcor_didx) )                 &
            locerrstat = HE5_SWwrfld ( pge_swath_id, TRIM(ADJUSTL(xtrcor_field)), &
            he5_start_2d, he5_stride_2d, he5_edge_2d, xtrack_cor(1:nxtrack,iline:iline+ntimes_loop-1) )
    END DO
    

    RETURN
  END SUBROUTINE xtrack_destriping_writecol


  SUBROUTINE condense_columns_for_median ( nt, nx, ynproc, ynrange, col, mqf, nc, ccol, cmqf)

    IMPLICIT NONE

    ! ---------------
    ! Input variables
    ! ---------------
    INTEGER (KIND=i4),                          INTENT (IN) :: nt, nx
    REAL    (KIND=r8), DIMENSION (1:nx,0:nt-1), INTENT (IN) :: col
    INTEGER (KIND=i2), DIMENSION (1:nx,0:nt-1), INTENT (IN) :: mqf
    LOGICAL,           DIMENSION (     0:nt-1), INTENT (IN) :: ynrange, ynproc


    ! ----------------
    ! Output variables
    ! ----------------
    INTEGER (KIND=i4),                          INTENT (OUT) :: nc
    REAL    (KIND=r8), DIMENSION (1:nx,0:nt-1), INTENT (OUT) :: ccol
    INTEGER (KIND=i2), DIMENSION (1:nx,0:nt-1), INTENT (OUT) :: cmqf

    ! --------------
    ! Local variable
    ! --------------
    INTEGER (KIND=i4) :: it

    ! ----------------------------------------------------------
    ! Initialize output variables
    ! ----------------------------------------------------------
    ccol = r8_missval  ;  cmqf = i2_missval  ;  nc = 0

    ! ------------------------------------------------------------
    ! Go through all swath lines and save the ones that are "live"
    ! ------------------------------------------------------------
    DO it = 0, nt-1
       IF ( ynproc(it) .AND. ynrange(it) ) THEN
          nc = nc + 1
          ccol(1:nx,nc) = col(1:nx,it)
          cmqf(1:nx,nc) = mqf(1:nx,it)
       END IF
    END DO

  END SUBROUTINE condense_columns_for_median

END MODULE OMSAO_destriping_module
