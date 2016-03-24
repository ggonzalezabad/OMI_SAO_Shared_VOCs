MODULE OMSAO_Reference_sector_module

  ! ------------------------------------------------------------------
  ! This module defines variables associated with the Reference Sector
  ! Correction and contains the subroutines needed to apply it to the
  ! HCHO retrieval.
  ! ------------------------------------------------------------------

  USE OMSAO_precision_module
  USE OMSAO_errstat_module
  USE OMSAO_parameters_module,   ONLY: maxchlen, r8_missval, r4_missval, &
       i4_missval, i2_missval, i1_missval
  USE OMSAO_radiance_ref_module, ONLY: l1b_radref_filename, yn_radiance_reference
  USE OMSAO_variables_module,    ONLY: l1b_rad_filename, common_latrange
  USE OMSAO_omidata_module,      ONLY: omi_radiance_swathname
  USE OMSAO_he5_module,          ONLY: granule_month

  IMPLICIT NONE

  ! ------------------------------------------------------------------
  ! maxngrid: Parameter to define maximum size of arrays; set to 1000
  ! grid_lat: Latitudes of GEOS-Chem background levels from the Pacif
  !            ic, Hawaii is not included
  ! background_correction: Background level-reference_sector_concentra
  !                        tion
  ! background_level: Total column median obtained from the Radiance
  !                   Reference granule (molecules/cm2)
  ! Reference_sector_concentration: Total column (molecules/cm2) obtai
  !                                 ned using the GEOS-Chem climatolog
  !                                 y by D. Millet
  ! ngridpoints: Number of points in grid_lat and Reference_Sector_con
  !              centration
  ! ------------------------------------------------------------------

  INTEGER (KIND=i2), PARAMETER               :: maxngrid = 1000
  REAL    (KIND=r8), DIMENSION (maxngrid)    :: grid_lat, background_level
  REAL    (KIND=r8), DIMENSION (maxngrid,12) :: Reference_sector_concentration
  REAL    (KIND=r8), DIMENSION (maxngrid)    :: Ref_column_month
  INTEGER (KIND=i2)                          :: ngridpoints

!!$  REAL    (KIND=r8), DIMENSION (:,:), ALLOCATABLE :: background_correction
  INTEGER (KIND=i2), DIMENSION (:,:), ALLOCATABLE :: refmqf

  CONTAINS
    
    SUBROUTINE Reference_Sector_correction (ntimes, nxtrack, xtrange,     &
               lat, saocol, saodco, saoamf, saomqf, saorms, saofcf,       &
               extr, pge_idx, n_max_rspec, &
               errstat)

      ! ---------------------------------------------------------------
      ! This subroutine is a wrapper for the Reference Background corre
      ! ction
      ! ---------------------------------------------------------------
      ! VARIABLE DECLARATION
      ! --------------------
      IMPLICIT NONE
      ! ---------------
      ! Input variables
      ! ---------------
      INTEGER (KIND=i4),                                   INTENT (IN)    :: ntimes, nxtrack
      INTEGER (KIND=i4), DIMENSION (0:ntimes-1,1:2),       INTENT (IN)    :: xtrange  
      REAL    (KIND=r4), DIMENSION (1:nxtrack,0:ntimes-1), INTENT (IN)    :: lat,extr
      REAL    (KIND=r8), DIMENSION (1:nxtrack,0:ntimes-1), INTENT (IN)    :: saoamf
      INTEGER (KIND=i2), DIMENSION (1:nxtrack,0:ntimes-1), INTENT (IN)    :: saofcf, saomqf
      INTEGER (KIND=i4),                                   INTENT (IN)    :: pge_idx, n_max_rspec
      REAL    (KIND=r8), DIMENSION (1:nxtrack,0:ntimes-1), INTENT (IN)    :: saocol, saodco, saorms

      ! ------------------
      ! Modified variables
      ! ------------------
      INTEGER (KIND=i4),                                   INTENT (INOUT) :: errstat

      ! ---------------
      ! Local variables
      ! ---------------
      INTEGER (KIND=i4)                                   :: nTimesRadRR, nXtrackRadRR, nWvlCCDrr, &
                                                             igrid, itrack, iline
      REAL    (KIND=r8), DIMENSION (1:nxtrack,0:ntimes-1) :: int_saocol, int_saodco
      INTEGER (KIND=i2), DIMENSION (1:nxtrack,0:ntimes-1) :: int_saomqf
      REAL    (KIND=r8), DIMENSION(1)                     :: background, out_lat, &
                                                             reference_sector, background_correction_old

      ! ------------------------
      ! Error handling variables
      ! ------------------------
      INTEGER (KIND=i4) :: version, locerrstat

      ! ------------------------------
      ! Name of this module/subroutine
      ! ------------------------------
      CHARACTER (LEN=27), PARAMETER :: modulename = 'Reference_sector_correction'
      
      locerrstat = pge_errstat_ok

      int_saocol   = saocol
      int_saodco   = saodco
      nTimesRadRR  = i4_missval
      nXtrackRadRR = i4_missval
      nWvlCCDrr    = i4_missval

      ! -------------------------------------------------
      ! Read the concentrations from the Reference Sector
      ! Total column molecules/cm2
      ! -------------------------------------------------
      CALL Read_reference_sector_concentration(errstat)
      
      ! ---------------------------------------------------
      ! Obtain dimensions of the Radiance Reference granule
      ! ---------------------------------------------------
      CALL omi_read_radiance_paras (                                  &
           l1b_radref_filename, nTimesRadRR, nXtrackRadRR, nWvlCCDrr, &
           omi_radiance_swathname, errstat                          )

!!$      ALLOCATE (background_correction(1:nXtrackRadRR,0:nTimesRadRR-1))
      ALLOCATE (refmqf(1:nXtrackRadRR,0:nTimesRadRR-1))
      ! --------------------------
      ! Initialize these variables
      ! --------------------------
!!$      background_correction           = 0.0_r8
      refmqf                          = i2_missval
      Ref_column_month(1:ngridpoints) = Reference_sector_concentration(1:ngridpoints,granule_month)


      ! -------------------------------------------------------------------------
      ! We need to retrieve the concentrations for the Radiance Reference Granule
      ! to be able to compute the latitude dependent correction due to the Refere
      ! ce sector hypothesis. This is a bit pedestrian, in the case of the Radian
      ! ce granule being the same that the Radiance Reference granule we could sk
      ! ip this step but since time is not an issue. Let us do it anyway.
      ! -------------------------------------------------------------------------
      CALL Reference_Sector_radref_retrieval_and_median(nTimesRadRR,  &
           nXTrackRadRR, nWvlCCDrr, pge_idx, n_max_rspec, errstat)

      DO itrack = 1, nxtrack
         DO iline = 0, ntimes-1
            int_saocol(itrack,iline) = int_saocol(itrack,iline) * saoamf(itrack,iline)
!!$            IF ( extr(itrack,iline) .NE. 0 .AND. saomqf(itrack,iline) .GE. 0 .AND. &
!!$                 saomqf(itrack,iline) .LE. 1 .AND.                                 &
!!$                 refmqf(itrack,iline) .GE. 0 .AND. refmqf(itrack,iline) .LE. 1) THEN
!!$               ! ----------------------------
!!$               ! Apply pixel based correction
!!$               ! ----------------------------
!!$               int_saocol(itrack,iline) = int_saocol(itrack,iline) - background_correction(itrack,iline)
!!$
!!$               ! --------------------------------
!!$               ! Apply smooth correction (median)
!!$               ! --------------------------------
!!$            ELSE
            IF (saomqf(itrack,iline) .EQ. 0 .AND. extr(itrack,iline) .EQ. 0 ) THEN               
               ! -------------------------------------------------------------
               ! Interpolate background_level, and reference_sector_correction
               ! ------------------------------------------------------------
               out_lat(1)                   = REAL(lat(itrack,iline), KIND=r8) 
               background_correction_old(1) = 0.0_r8
               background(1)                = 0.0_r8
               reference_sector(1)          = 0.0_r8
               ! ---------------------------------
               ! Enough room for the interpolation
               ! ---------------------------------
               IF (out_lat(1) .LT. grid_lat(2) .OR. out_lat(1) .GT. grid_lat(ngridpoints-1)) CYCLE
               CALL ezspline_1d_interpolation ( INT(ngridpoints, KIND=i4),  &
                    grid_lat(1:ngridpoints), Ref_column_month(1:ngridpoints),              &
                    1, out_lat(1), reference_sector(1), &
                    locerrstat )
               CALL ezspline_1d_interpolation ( INT(ngridpoints, KIND=i4),  &
                    grid_lat(1:ngridpoints), background_level(1:ngridpoints),              &
                    1, out_lat(1), background(1), &
                    locerrstat )

               background_correction_old(1) = background(1) - &
                    ( reference_sector(1) * saoamf(itrack,iline) )
               ! ---------------------------------------------------------
               ! And apply the result of the interpolation to the original
               ! saocol
               ! ---------------------------------------------------------
               int_saocol(itrack,iline) = int_saocol(itrack,iline) - &
                                          background_correction_old(1)
            ELSE
               ! Nothing to do here
            END IF
            ! -------------------
            ! Convert back to VCD
            ! -------------------
            int_saocol(itrack,iline) = int_saocol(itrack,iline) / saoamf(itrack,iline)
         END DO
      END DO
     
      ! ---------------------------------------------------
      ! Final step to output the results, the new reference
      ! sector corrected total columns
      ! ---------------------------------------------------
      CALL he5_write_reference_sector_corrected_column(pge_idx, &
           ntimes, nxtrack, int_saocol, int_saodco, locerrstat)

      ! ---------------------------------------------------
      ! Deallocate background_correction and mem_correction
      ! ---------------------------------------------------
!!$      IF ( ALLOCATED ( background_correction ) ) DEALLOCATE ( background_correction )
      IF ( ALLOCATED ( refmqf                ) ) DEALLOCATE ( refmqf                )

      errstat = MAX ( errstat, locerrstat )

    END SUBROUTINE Reference_Sector_correction

    SUBROUTINE Read_reference_sector_concentration(errstat)

      USE OMSAO_variables_module, ONLY: OMSAO_refseccor_filename
      USE OMSAO_indices_module,   ONLY: OMSAO_refseccor_lun
      
      IMPLICIT NONE
      ! ------------------
      ! Modified variables
      ! ------------------
      INTEGER (KIND=i4),                                  INTENT (INOUT) :: errstat
      ! ---------------------------------
      ! External OMI and Toolkit routines
      ! ---------------------------------
      INTEGER (KIND=i4), EXTERNAL :: &
           pgs_smf_teststatuslevel, pgs_io_gen_openf, pgs_io_gen_closef
      ! ---------------
      ! Local variables
      ! ---------------
      INTEGER (KIND=i4)           :: funit, igrid
      CHARACTER(LEN=1), PARAMETER :: hstr='#'
      LOGICAL                     :: file_header
      CHARACTER (LEN=maxchlen)    :: header_line
      ! ------------------------
      ! Error handling variables
      ! ------------------------
      INTEGER (KIND=i4) :: version, locerrstat, ios
      ! ------------------------------
      ! Name of this module/subroutine
      ! ------------------------------
      CHARACTER (LEN=35), PARAMETER :: modulename = 'Read_reference_sector_concentration'
      
      locerrstat = pge_errstat_ok

      ! --------------------
      ! Initialize variables
      ! --------------------
      ngridpoints = 0.0_r8
      grid_lat(1:maxngrid) =  0.0_r8
      Reference_sector_concentration(1:maxngrid,1:12) = 0.0_r8

      ! -----------------------------------------
      ! Open Reference Sector concentrations file
      ! -----------------------------------------
      version = 1
      locerrstat = PGS_IO_GEN_OPENF ( OMSAO_refseccor_lun, PGSd_IO_Gen_RSeqFrm, 0, funit, version )
      locerrstat = PGS_SMF_TESTSTATUSLEVEL(locerrstat)
      CALL error_check ( &
           locerrstat, pgs_smf_mask_lev_s, pge_errstat_error, OMSAO_E_OPEN_REFSECCOR_FILE, &
           modulename//f_sep//TRIM(ADJUSTL(OMSAO_refseccor_filename)), vb_lev_default, errstat )
      IF (  errstat /= pge_errstat_ok ) RETURN

      ! ------------
      ! Reading file
      ! -----------------------------------------------------------------------------------
      ! Skip header lines. The header is done when the first character of the line is not #
      ! -----------------------------------------------------------------------------------
      file_header = .TRUE.
      skip_header: DO WHILE (file_header .EQV. .TRUE.)
         READ (UNIT=funit, FMT='(A)', IOSTAT=ios) header_line
         IF ( ios /= 0 ) THEN
            CALL error_check ( &
                 ios, file_read_ok, pge_errstat_error, OMSAO_E_READ_REFSECCOR_FILE, &
                 modulename//f_sep//TRIM(ADJUSTL(OMSAO_refseccor_filename)), vb_lev_default, errstat )
            IF (  errstat /= pge_errstat_ok ) RETURN
         END IF
         IF (header_line(1:1) /= hstr) THEN
            file_header = .FALSE.
         ENDIF
      END DO skip_header

      ! -----------------------------------------------
      ! Read number of grid points. Variable via module
      ! -----------------------------------------------
      READ (UNIT=funit, FMT='(I5)', IOSTAT=ios) ngridpoints
      IF ( ios /= 0 ) THEN
         CALL error_check ( &
              ios, file_read_ok, pge_errstat_error, OMSAO_E_READ_REFSECCOR_FILE, &
              modulename//f_sep//TRIM(ADJUSTL(OMSAO_refseccor_filename)), vb_lev_default, errstat )
         IF (  errstat /= pge_errstat_ok ) RETURN
      END IF

      ! ---------------------------------------------------------
      ! Read reference sector concentrations. Variable via module
      ! ---------------------------------------------------------
      DO igrid = 1, ngridpoints
         READ (UNIT=funit, FMT='(F6.2,2x,12(1x,E14.7))', IOSTAT=ios) grid_lat(igrid), &
              Reference_sector_concentration(igrid, 1:12)
         IF ( ios /= 0 ) THEN
            CALL error_check ( &
                 ios, file_read_ok, pge_errstat_error, OMSAO_E_READ_REFSECCOR_FILE, &
                 modulename//f_sep//TRIM(ADJUSTL(OMSAO_refseccor_filename)), vb_lev_default, errstat )
            IF (  errstat /= pge_errstat_ok ) RETURN
         END IF
      END DO

      ! -----------------------------------------------
      ! Close monthly average file, report SUCCESS read
      ! -----------------------------------------------
      locerrstat = PGS_IO_GEN_CLOSEF ( funit )
      locerrstat = PGS_SMF_TESTSTATUSLEVEL(locerrstat)
      CALL error_check ( &
           locerrstat, pgs_smf_mask_lev_s, pge_errstat_warning, OMSAO_W_CLOSE_REFSECCOR_FILE, &
           modulename//f_sep//TRIM(ADJUSTL(OMSAO_refseccor_filename)), vb_lev_default, errstat )
      IF ( errstat >= pge_errstat_error ) RETURN
     
    END SUBROUTINE Read_reference_sector_concentration

    SUBROUTINE Reference_Sector_radref_retrieval_and_median(nTimesRadRR,   &
         nXTrackRadRR, nWvlCCDrr, pge_idx, n_max_rspec, errstat)

      USE OMSAO_wfamf_module,     ONLY: amf_calculation_bis
      USE OMSAO_variables_module, ONLY: OMSAO_refseccor_cld_filename, voc_amf_filenames
      USE OMSAO_indices_module,   ONLY: voc_omicld_idx


      IMPLICIT NONE

      ! ---------------
      ! Input variables
      ! ---------------
      INTEGER (KIND=i4), INTENT (IN) :: errstat, pge_idx, n_max_rspec
      INTEGER (KIND=i4), INTENT (IN) :: nTimesRadRR, nXtrackRadRR, nWvlCCDrr

      ! ---------------
      ! Local variables
      ! ---------------
      INTEGER (KIND=i4)                                           :: first_line, last_line
      INTEGER (KIND=i4), DIMENSION (0:nTimesRadRR-1,1:2)          :: omi_xtrpix_range_rr
      LOGICAL,           DIMENSION (0:nTimesRadRR-1)              :: yn_radfitref_range
      CHARACTER(LEN=maxchlen)                                     :: l1b_rad_save_filename 
      REAL    (KIND=r8), DIMENSION (nXtrackRadRR,0:nTimesRadRR-1) :: mem_column_amount, &
           mem_column_uncertainty, mem_amf, mem_rms
      REAL    (KIND=r4), DIMENSION (nXtrackRadRR,0:nTimesRadRR-1) :: mem_latitude, mem_longitude, &
           mem_sza, mem_vza, mem_height
      INTEGER (KIND=i2), DIMENSION (nXtrackRadRR,0:nTimesRadRR-1) :: mem_xtrflg, mem_fit_flag
      INTEGER (KIND=i2), DIMENSION (nXtrackRadRR,0:nTimesRadRR-1) :: mem_snow, mem_glint
      LOGICAL,           DIMENSION (0:nTimesRadRR-1)              :: yn_szoom_rs, yn_common_range
      INTEGER (KIND=i1), DIMENSION (0:nTimesRadRR-1)              :: binfac_rs
      REAl    (KIND=r8), DIMENSION (1)                            :: zerovec = 0.0_r8
      LOGICAL                                                     :: yn_write

      ! ------------------------
      ! Error handling variables
      ! ------------------------
      INTEGER (KIND=i4) :: version, locerrstat
      ! ------------------------------
      ! Name of this module/subroutine
      ! ------------------------------
      CHARACTER (LEN=64), PARAMETER :: modulename = &
           'Reference_Sector_radiance_reference_granule_retrieval'
      
      locerrstat = pge_errstat_ok

      ! ---------------------------------
      ! A bit of mess with the file names
      ! ---------------------------------
      l1b_rad_save_filename = l1b_rad_filename
      l1b_rad_filename      = l1b_radref_filename
      ! -----------------------
      ! Variable initialization
      ! -----------------------
      mem_column_amount      = r8_missval
      mem_column_uncertainty = r8_missval
      ! ---------------------------------------------------
      ! mem_correction needs to be initialized to 0.0 since
      ! we want no correction if any of the two reference
      ! sector retrievals does fail.
      ! ---------------------------------------------------
      mem_amf                = 1.0_r8
      mem_rms                = r8_missval
      mem_latitude           = r4_missval
      mem_longitude          = r4_missval
      mem_sza                = r4_missval
      mem_vza                = r4_missval
      mem_fit_flag           = i2_missval
      mem_xtrflg             = i2_missval
      mem_snow               = i2_missval
      mem_glint              = i2_missval
      mem_height             = r4_missval

      ! -----------------------------------------------------
      ! I want to perform the retrieval for the whole granule
      ! -----------------------------------------------------
      first_line = 0  ;  last_line = nTimesRadRR-1
      yn_radfitref_range = .TRUE.
      omi_xtrpix_range_rr(0:nTimesRadRR-1,1) = 1
      omi_xtrpix_range_rr(0:nTimesRadRR-1,2) = nXtrackRadRR

      ! ---------------------------------
      ! Read radiance reference latitudes
      ! ---------------------------------
      CALL read_latitude ( &
           TRIM(ADJUSTL(l1b_radref_filename)), TRIM(ADJUSTL(omi_radiance_swathname)), &
           nTimesRadRR, nXtrackRadRR, mem_latitude(1:nXtrackRadRR,0:nTimesRadRR-1) )
      ! --------------------------------------------------
      ! Compute the common mode for the Radiance Reference
      ! granule
      ! --------------------------------------------------
      yn_common_range(0:nTimesRadRR-1) = .FALSE.
      CALL find_swathline_range ( &
           TRIM(ADJUSTL(l1b_radref_filename)), TRIM(ADJUSTL(omi_radiance_swathname)),  &
           nTimesRadRR, nXtrackRadRR, mem_latitude(1:nXtrackRadRR,0:nTimesRadRR-1), &
           common_latrange(1:2), yn_common_range(0:nTimesRadRR-1), locerrstat        )

      ! ----------------------------------------------------------
      ! Interface to the loop over all swath lines for common mode
      ! ----------------------------------------------------------
      CALL omi_pge_swathline_loop (                                    &
           pge_idx, nTimesRadRR, nxtrackRadRR, nWvlCCDRR, n_max_rspec, &
           yn_common_range(0:nTimesRadRR-1),                           &
           omi_xtrpix_range_rr(0:nTimesRadRR-1,1:2),                   &
           yn_radiance_reference, .FALSE., -1,                         &
           .TRUE., locerrstat)

      ! ---------------------------------------------------
      ! Set the index value of the Common Mode spectrum and
      ! assign values to the fitting parameter arrays
      ! ---------------------------------------------------
      CALL compute_common_mode ( .FALSE., nXtrackRadRR, 1, zerovec, zerovec, .TRUE. )

      ! --------------------------------------
      ! Interface to loop over all swath lines
      ! --------------------------------------
      CALL omi_pge_swathline_loop_memory (                             &
           pge_idx, nTimesRadRR, nXtrackRadRR, nWvlCCDRR, n_max_rspec, &
           yn_radfitref_range(0:nTimesRadRR-1),                        &
           omi_xtrpix_range_rr(0:nTimesRadRR-1,1:2),                   &
           yn_radiance_reference, .FALSE., -1,                         &
           .TRUE., mem_column_amount(1:nXtrackRadRR,0:nTimesRadRR-1),  &
           mem_column_uncertainty(1:nXtrackRadRR,0:nTimesRadRR-1),     &
           mem_rms(1:nXtrackRadRR,0:nTimesRadRR-1),                    &
           mem_fit_flag(1:nXtrackRadRR,0:nTimesRadRR-1),               &
           mem_xtrflg(1:nXtrackRadRR,0:nTimesRadRR-1),                 &
           mem_latitude(1:nXtrackRadRR,0:nTimesRadRR-1),               &
           mem_longitude(1:nXtrackRadRR,0:nTimesRadRR-1),              &
           mem_sza(1:nXtrackRadRR,0:nTimesRadRR-1),                    &
           mem_vza(1:nXtrackRadRR,0:nTimesRadRR-1),                    &
           mem_height(1:nXtrackRadRR,0:nTimesRadRR-1),                 &
           locerrstat)

      ! ------------------------------------------------------
      ! mem_column_uncertainty, men_latitude, mem_fit_flag and
      ! mem_column_uncertainty are holding the radiance refere
      ! nce in memory (SCD).
      ! Now I need to work out VCD.
      ! ------------------------------------------------------
      ! Read snow and ice data
      ! ----------------------
      ! ----------------------------------
      ! Read L1b glint and snow/ice flags
      ! ----------------------------------
      CALL omi_read_glint_ice_flags ( &
           l1b_radref_filename, nXtrackRadRR, nTimesRadRR, mem_snow, &
           mem_glint, locerrstat )

      ! --------------------
      ! OMI zoom mode or not
      ! --------------------
      CALL omi_read_binning_factor ( &
           TRIM(ADJUSTL(l1b_radref_filename)), TRIM(ADJUSTL(omi_radiance_swathname)), &
           nTimesRadRR, binfac_rs(0:nTimesRadRR-1), yn_szoom_rs(0:nTimesRadRR-1),     &
           locerrstat )

      ! ------------------------------------------------------
      ! Compute AMF using internal variables to avoud conflict
      ! with main retrieval.
      ! ---------------------------------------------------------
      ! No output of the amf calculation for the reference sector
      ! Height equal to 0 everyehwere. Pacific Ocean
      ! ---------------------------------------------------------
      yn_write = .FALSE.
      ! --------------------------------------------------------
      ! To read clouds corresponding with the radiance reference
      ! granule
      ! --------------------------------------------------------
      voc_amf_filenames(voc_omicld_idx) = TRIM(ADJUSTL(OMSAO_refseccor_cld_filename))

      CALL amf_calculation_bis (                                            &
           pge_idx, nTimesRadRR, nXtrackRadRR, mem_latitude, mem_longitude, &
           mem_sza, mem_vza, mem_snow, mem_glint, omi_xtrpix_range_rr,      &
           yn_szoom_rs, mem_column_amount, mem_column_uncertainty, mem_amf, &
           mem_height, yn_write, locerrstat )

      ! --------------------------------------------------------
      ! Compute average fitting statistics and main quality flag
      ! --------------------------------------------------------
      CALL compute_fitting_statistics_nohe5 (                           &
           pge_idx, nTimesRadRR, nXtrackRadRR, omi_xtrpix_range_rr,     &
           mem_column_amount, mem_column_uncertainty, mem_rms,          &
           mem_fit_flag, refmqf, locerrstat )

      ! ------------------------------------------------
      ! Apply the background correction to Slant columns
      ! (new) for row anomaly pixels
      ! ------------------------------------------------
      mem_column_amount = mem_column_amount * mem_amf

!!$      CALL compute_background_correction_bis(mem_column_amount, mem_latitude, mem_amf, &
!!$                                         nXtrackRadRR, nTimesRadRR, refmqf,            &
!!$                                         mem_xtrflg, locerrstat)

      ! -----------------------------------------------------
      ! Once we have the radiance reference retrievals we can
      ! compute the background correction
      ! -----------------------------------------------------
      CALL compute_background_median(nTimesRadRR, nXtrackRadRR,     &
           mem_column_amount, mem_column_uncertainty, mem_amf,      &
           mem_latitude, mem_longitude, mem_xtrflg, refmqf, locerrstat)

      mem_column_amount = mem_column_amount / mem_amf

      ! -------------------------------
      ! No more mess with the filenames
      ! -------------------------------
      l1b_rad_filename = l1b_rad_save_filename
      
      ! -----------
      ! Error check
      ! -----------
      locerrstat = MAX ( locerrstat, errstat )
      IF ( locerrstat >= pge_errstat_error )  RETURN

    END SUBROUTINE Reference_Sector_radref_retrieval_and_median

!!$    SUBROUTINE compute_background_correction_bis(mem_column_amount, mem_latitude, mem_amf, &
!!$                                                 nXtrackRadRR, nTimesRadRR, refmqf, mem_xtrflg, locerrstat)
!!$
!!$      USE OMSAO_variables_module, ONLY: max_good_col
!!$
!!$      IMPLICIT NONE
!!$      
!!$      ! ---------------
!!$      ! Input variables
!!$      ! ---------------
!!$      INTEGER (KIND=i4),                                           INTENT(IN) :: &
!!$                         nTimesRadRR, nXtrackRadRR
!!$      REAL    (KIND=r8), DIMENSION (nXtrackRadRR,0:nTimesRadRR-1), INTENT(IN) :: &
!!$                         mem_column_amount, mem_amf
!!$      INTEGER (KIND=i2), DIMENSION (nXtrackRadRR,0:nTimesRadRR-1), INTENT(IN) :: refmqf, mem_xtrflg
!!$      REAL    (KIND=r4), DIMENSION (nXtrackRadRR,0:nTimesRadRR-1), INTENT(IN) :: mem_latitude
!!$      LOGICAL                                                                 :: yn_row
!!$
!!$      ! ------------------
!!$      ! Modified variables
!!$      ! ------------------
!!$      INTEGER (KIND=i4), INTENT(INOUT) :: locerrstat
!!$
!!$      ! ---------------
!!$      ! Local variables
!!$      ! ---------------
!!$      INTEGER (KIND=i4)               :: iline, itrack, igrid, npixels, ipixel
!!$      REAL    (KIND=r8), DIMENSION(1) :: Ref_column, latitude
!!$      REAL    (KIND=r8)               :: nTotal, Total, Average, nSted, Stddev
!!$
!!$      ! -------------------
!!$      ! Routine starts here
!!$      ! -------------------
!!$      locerrstat       = pge_errstat_ok
!!$      
!!$      ! -----------------------------------------------------------
!!$      ! Loop pixel by pixel to compute the background correction as
!!$      ! sao_colum - reference_colum only for row anomaly pixels
!!$      ! -----------------------------------------------------------
!!$      DO iline = 0,nTimesRadRR-1
!!$         DO itrack = 1,nXtrackRadRR
!!$
!!$            ! -----------------------------------------------------
!!$            ! If the main quality flag of the retrieval is not good
!!$            ! then cycle.
!!$            ! -----------------------------------------------------
!!$            IF (refmqf(itrack,iline) .EQ. 0) THEN
!!$               ! ---------------------------------------------------
!!$               ! For the rest of the pixels find out the latitude,
!!$               ! interpolate the reference column and substrack from
!!$               ! the retrieved colum
!!$               ! ---------------------------------------------------
!!$               latitude(1) = REAL(mem_latitude(itrack,iline), KIND=r8)
!!$               CALL ezspline_1d_interpolation ( INT(ngridpoints, KIND=i4),    &
!!$                    grid_lat(1:ngridpoints), Ref_column_month(1:ngridpoints), &
!!$                    1, latitude(1), Ref_column(1), locerrstat )
!!$               background_correction(itrack,iline) = mem_column_amount(itrack,iline) &
!!$                    - ( Ref_column(1) * mem_amf(itrack,iline) )
!!$            ENDIF
!!$         END DO
!!$      END DO      
!!$      
!!$    END SUBROUTINE compute_background_correction_bis    

    SUBROUTINE he5_write_reference_sector_corrected_column(pge_idx, &
         nt, nx, column, uncertainty, errstat)

      USE OMSAO_he5_module
      USE OMSAO_omidata_module,   ONLY: n_roff_dig
      USE OMSAO_indices_module,   ONLY: pge_hcho_idx, pge_gly_idx, pge_bro_idx

      IMPLICIT NONE

      ! ---------------
      ! Input variables
      ! ---------------
      INTEGER (KIND=i4),                          INTENT (IN) :: pge_idx, nt, nx
      REAL    (KIND=r8), DIMENSION (1:nx,0:nt-1), INTENT (IN) :: column, uncertainty

      ! ------------------
      ! Modified variables
      ! ------------------
      INTEGER (KIND=i4), INTENT (INOUT) :: errstat

      ! ------------------------------
      ! Name of this module/subroutine
      ! ------------------------------
      CHARACTER (LEN=43), PARAMETER :: modulename = &
           'he5_write_reference_sector_corrected_column'

      ! ---------------
      ! Local variables
      ! ---------------
      INTEGER (KIND=i4)                          :: locerrstat
      REAL    (KIND=r8), DIMENSION (1:nx,0:nt-1) :: colloc
  
      ! -----------------------------------------------
      ! Total column corrected by the reference sector.
      ! -----------------------------------------------
      ! Only implemented for HCHO.
      ! -----------------------------------------------
      
      locerrstat = pge_errstat_ok

      he5_start_2d  = (/ 0, 0 /) ;  he5_stride_2d = (/ 1, 1 /) ; he5_edge_2d = (/ nx, nt /)      

      ! ----------------------------------------------------------------------
      ! All PGEs: Output of columns and column uncertainties left as a comment
      ! for the future. So far only for HCHO
      ! ----------------------------------------------------------------------
      IF (pge_idx .EQ. pge_hcho_idx) THEN
         colloc = column
         CALL roundoff_2darr_r8 ( n_roff_dig, nx, nt, colloc(1:nx,0:nt-1) )
         locerrstat = HE5_SWWRFLD ( pge_swath_id, TRIM(ADJUSTL(rscol_field)), &
              he5_start_2d, he5_stride_2d, he5_edge_2d, colloc(1:nx,0:nt-1) )
         errstat = MAX ( errstat, locerrstat )
  
         !colloc = uncertainty
         !CALL roundoff_2darr_r8 ( n_roff_dig, nx, nt, colloc(1:nx,0:nt-1) )
         !locerrstat = HE5_SWWRFLD ( pge_swath_id, TRIM(ADJUSTL(rsdcol_field)), &
         !     he5_start_2d, he5_stride_2d, he5_edge_2d, colloc(1:nx,0:nt-1) )
         !errstat = MAX ( errstat, locerrstat )

      END IF

      ! ------------------
      ! Check error status
      ! ------------------
      CALL error_check ( locerrstat, HE5_STAT_OK, pge_errstat_error, OMSAO_E_HE5SWWRFLD, &
           modulename, vb_lev_default, errstat )


      RETURN
  
    END SUBROUTINE he5_write_reference_sector_corrected_column

    SUBROUTINE compute_background_median(nTimesRadRR, nXtrackRadRR, &
           mem_column_amount, mem_column_uncertainty, mem_amf,          &
           mem_latitude, mem_longitude, mem_xtrflg, refmqf, locerrstat)

      USE OMSAO_median_module, ONLY: median

      IMPLICIT NONE
      
      ! ---------------
      ! Input variables
      ! ---------------
      REAL    (KIND=r8), DIMENSION (nXtrackRadRR,0:nTimesRadRR-1), &
           INTENT(IN) :: mem_column_amount, mem_column_uncertainty,&
           mem_amf
      REAL    (KIND=r4), DIMENSION (nXtrackRadRR,0:nTimesRadRR-1), &
           INTENT(IN) :: mem_latitude, mem_longitude
      INTEGER (KIND=i2), DIMENSION (nXtrackRadRR,0:nTimesRadRR-1), &
           INTENT(IN) :: refmqf, mem_xtrflg
      INTEGER (KIND=i4),                                           &
           INTENT(IN) :: nTimesRadRR, nXtrackRadRR

      ! ------------------
      ! Modified variables
      ! ------------------
      INTEGER (KIND=i4), INTENT(INOUT) :: locerrstat
      
      ! ---------------
      ! Local variables
      ! ---------------
      INTEGER (KIND=i2)                           :: igrid, itrack, iline
      REAL    (KIND=r8), DIMENSION(ngridpoints+1) :: grid
      REAL    (KIND=r8), DIMENSION(1000)          :: into_median
      REAL    (KIND=r8), DIMENSION (nXtrackRadRR,0:nTimesRadRR-1) &
                         :: column_amounts
      INTEGER (KIND=i2), DIMENSION (nXtrackRadRR,0:nTimesRadRR-1) &
                         :: yn_median
      INTEGER (KIND=i4)                           :: npixels, ipixel, nearest
      INTEGER (KIND=i4), DIMENSION(ngridpoints)   :: nonzero

      ! -------------------
      ! Routine starts here
      ! -------------------

      locerrstat = pge_errstat_ok

      background_level = 0.0_r8

      ! ----------------------------------------------------------------
      ! Computing the median of all the pixels of the radiance reference
      ! granule within two consequtive latitudes of the background conce
      ! trations grid readed from the background_concentrations file.
      ! Only pixels between 180W and 140W are considered. Indeed because
      ! of the conditions impossed on the radiance reference granule, as
      ! close as possible to the 165W, they should be no pixels above
      ! 180W and below 140W.
      ! ----------------------------------------------------------------
      ! Generating grid
      ! ---------------
      DO igrid = 1, ngridpoints+1
         grid(igrid) = (180.0_r8 / (REAL(ngridpoints+1, KIND=r8)) * &
                       (REAL(igrid, KIND=r8))) - 90.0_r8
      END DO

      DO igrid = 1, ngridpoints

         column_amounts = 0.0_r8
         yn_median      = 0_i2
         into_median    = 0.0_r8
         npixels        = 0_i4

         ! -------------------------------------------------------------------
         ! Finging pixels in the granule between grid(igrid) and grid(igrid+1)
         ! Mid point is grid_lat(igrid).
         ! -------------------------------------------------------------------
         WHERE (mem_latitude .gt. grid(igrid) .AND. mem_latitude .lt. grid(igrid+1) &
                .AND. refmqf .EQ. 0 .AND. mem_xtrflg .EQ. 0 .AND. mem_longitude     &
                .GT. -180 .AND. mem_longitude .LT. -140)
            column_amounts = mem_column_amount
            yn_median      = 1
         END WHERE
         npixels = SUM(yn_median)
         ipixel  = 0_i4
         DO itrack = 1, nXtrackRadRR
            DO iline = 0, nTimesRadRR-1
               IF (yn_median(itrack,iline) .EQ. 1) THEN
                  ipixel = ipixel + 1
                  into_median(ipixel) = column_amounts(itrack,iline)
               END IF
            END DO
         END DO
         IF (npixels .GT. 2) THEN
            background_level(igrid) = median(npixels, into_median(1:ipixel))
         END IF

      END DO
     
    END SUBROUTINE COMPUTE_BACKGROUND_median

END MODULE OMSAO_Reference_sector_module
