MODULE OMSAO_pixelcorner_module

  ! =========================================================================== !
  !                                                                             !
  ! Compute geolocation of OMI pixel corners from latitude and longitude fields !
  !                                                                             !
  ! =========================================================================== !

  USE OMSAO_precision_module
  USE OMSAO_errstat_module
  USE L1B_Reader_class
  USE OMSAO_parameters_module, ONLY: r8_missval
  USE OMSAO_omidata_module,    ONLY: nlines_max, gzoom_spix, gzoom_epix
  USE OMSAO_he5_module

  IMPLICIT NONE

  ! ----------------
  ! Local Parameters
  ! ---------------------------------------------------------------------
  ! * Values for Pi (rad, deg) and Conversions between Degree and Radians
  ! ---------------------------------------------------------------------
  REAL (KIND=r8), PARAMETER, PRIVATE :: pi         = 3.14159265358979_r8  ! 2*ASIN(1.0_r8)
  REAL (KIND=r8), PARAMETER, PRIVATE :: pihalf     = 0.5_r8  * pi
  REAL (KIND=r8), PARAMETER, PRIVATE :: twopi      = 2.0_r8  * pi
  REAL (KIND=r8), PARAMETER, PRIVATE :: pi_deg     = 180.0_r8
  REAL (KIND=r8), PARAMETER, PRIVATE :: pihalf_deg =  90.0_r8
  REAL (KIND=r8), PARAMETER, PRIVATE :: twopi_deg  = 360.0_r8
  REAL (KIND=r8), PARAMETER, PRIVATE :: deg2rad    = pi       / 180.0_r8
  REAL (KIND=r8), PARAMETER, PRIVATE :: rad2deg    = 180.0_r8 / pi
  ! ---------------------------------------------------------------------
  ! * Precison for DEG <-> RAD conversion - anything less than EPS
  !   is effectively ZERO.
  ! ---------------------------------------------------------------------
  REAL (KIND=r8), PARAMETER, PRIVATE :: eps = 1.0E-10_r8
  ! ---------------------------------------------------------------------


CONTAINS


  SUBROUTINE compute_pixel_corners ( ntimes, nxtrack, lat, lon, yn_szoom, errstat )

    ! =======================================================
    ! Computes OMI pixel corner coordinates, start to finish:
    !
    !  * Reads L1b geolocation data
    !  * Computes the corners
    !  * Writes output to file
    ! =======================================================

    IMPLICIT NONE

    ! ---------------
    ! Input variables
    ! ---------------
    INTEGER (KIND=i4),                                    INTENT (IN) :: ntimes, nxtrack
    REAL    (KIND=r4), DIMENSION (1:nxtrack, 0:ntimes-1), INTENT (IN) :: lat, lon
    LOGICAL,           DIMENSION (           0:ntimes-1), INTENT (IN) :: yn_szoom

    ! ----------------
    ! Output variables
    ! ----------------
    INTEGER (KIND=i4), INTENT (OUT) :: errstat

    ! ---------------
    ! Local variables
    ! ---------------
    INTEGER (KIND=i4)                                  :: estat, iline, nchunk, nxtloc, spix, epix
    REAL (KIND=r8), DIMENSION (0:nxtrack,0:ntimes)     :: corner_lat, corner_lon
    REAL (KIND=r4), DIMENSION (:,:), ALLOCATABLE       :: cor_latlon

    errstat = pge_errstat_ok

    !! --------------------------------------------------------------------
    !! Read Latitudes and Longitudes from HE5 files (properly dimensioned!)
    !! --------------------------------------------------------------------
    !CALL pixcorner_latlon_read ( &
    !     ntimes, nxtrack, lon(1:nXtrack,0:nTimes-1), lat(1:nXtrack,0:nTimes-1), errstat )

    ! --------------------------------------------
    ! Initialize coner lat/lon with missing values
    ! --------------------------------------------
    corner_lat = r8_missval ; corner_lon = r8_missval
    
    ! ---------------------------
    ! Check for spatial zoom mode
    ! ----------------------------------------------------------
    ! YN_SZOOM is .TRUE. for rebinned "global" zoom mode only. 
    ! Both global and true zoom mode have the regular amount of
    ! cross-track positions (i.e., 60).
    ! ----------------------------------------------------------
    IF ( ANY(yn_szoom(0:nTimes-1)) ) THEN
       spix = gzoom_spix ; epix = gzoom_epix
    ELSE
       spix = 1          ; epix = nXtrack
    END IF
    nXtloc = epix - spix + 1

    ! ------------------------------
    ! Compute the corner coordinates
    ! ------------------------------
    corner_lat = r8_missval ;  corner_lon = r8_missval
    CALL sphgeo_comp_pixel_corners (                                                &
         nxtloc, ntimes, lon(spix:epix,0:ntimes-1), lat(spix:epix,0:ntimes-1),      &
         corner_lon(spix-1:epix,0:ntimes), corner_lat(spix-1:epix,0:ntimes), estat, &
         yn_omi_pixel_adjust_k=.FALSE. )

    ! --------------------------------------------------------------------
    ! Write corner coordinates to L2 output file. Do this in chunks of no
    ! more than 100 lines, since HE5 hiccups when slices become too large.
    ! NCHUNK carries the current chunk size, which is either NLINES_MAX or
    ! the remainder between the last multiple of 100 and NTIMES.
    ! --------------------------------------------------------------------
    DO iline = 0, ntimes, nlines_max

       ! --------------------------------------------------------------
       ! Make sure the chunk is not larger than the current array slice
       ! --------------------------------------------------------------
       nchunk = nlines_max
       IF ( (iline+nchunk) > ntimes  ) nchunk = ntimes - iline + 1

       ! ---------------------------------------------
       ! Set up Start, Stride, and Edge for HE5 output
       ! ---------------------------------------------
       he5_start_2d  = (/         0,  iline /)
       he5_stride_2d = (/         1,      1 /)
       he5_edge_2d   = (/ nxtrack+1, nchunk /)
    
       ALLOCATE ( cor_latlon(0:nxtrack,1:nchunk) )

       cor_latlon(0:nxtrack,1:nchunk) = REAL(corner_lat(0:nxtrack,iline:iline+nchunk-1), KIND=r4 )
       estat = HE5_SWWRFLD ( pge_swath_id, pxclat_field, he5_start_2d, he5_stride_2d, he5_edge_2d, &
            cor_latlon(0:nxtrack,1:nchunk) ) 
       cor_latlon(0:nxtrack,1:nchunk) = REAL(corner_lon(0:nxtrack,iline:iline+nchunk-1), KIND=r4 )
       estat = HE5_SWWRFLD ( pge_swath_id, pxclon_field, he5_start_2d, he5_stride_2d, he5_edge_2d, &
            cor_latlon(0:nxtrack,1:nchunk) ) 

       IF ( ALLOCATED (cor_latlon) ) DEALLOCATE (cor_latlon)

    END DO

    RETURN
  END SUBROUTINE compute_pixel_corners

  SUBROUTINE sphergeom_asa ( c, alp, bet, a, b )

    ! -------------------------------------------------------------------
    ! ASA: Angle-Side-Angle 
    ! 
    ! Computation of quantities in a spherical triangle: Given two angles
    ! ALP and BET, and the side C inbetween them, computes the sides A
    ! and B.
    ! -------------------------------------------------------------------
    IMPLICIT NONE

    ! ---------------
    ! Input variables
    ! ---------------
    REAL (KIND=r8),    INTENT (IN) :: c, alp, bet

    ! ----------------
    ! Output variables
    ! ----------------
    REAL (KIND=r8), INTENT (OUT) :: a, b

    ! ---------------
    ! Local variables
    ! ---------------
    REAL (KIND=r8) :: tmp1, tmp2, gam

    ! --------------------------------------------------
    ! Initialize output variables to keep compiler happy
    ! --------------------------------------------------
    a = 0.0_r8 ; b = 0.0_r8

    ! -------------------------------------
    ! Compute the angle BET opposite side B
    ! -------------------------------------
    tmp1 = ATAN(TAN(c/2.0_r8) * COS((alp-bet)/2.0_r8) / COS((alp+bet)/2.0_r8))
    tmp2 = ATAN(TAN(c/2.0_r8) * SIN((alp-bet)/2.0_r8) / SIN((alp+bet)/2.0_r8))
    
    a = tmp1 + tmp2
    b = tmp1 - tmp2

    gam = ACOS(-COS(alp)*COS(bet) + SIN(alp)*SIN(bet)*COS(c))
    a   = ASIN(SIN(c)*SIN(alp)/SIN(gam))
    b   = ASIN(SIN(c)*SIN(bet)/SIN(gam))

    RETURN
  END SUBROUTINE sphergeom_asa

  SUBROUTINE sphergeom_ssa ( a, b, alp, c, bet, gam )

    ! -------------------------------------------------------------------
    ! SSA: Side-Side-Angle 
    ! 
    ! Computation of quantities in a spherical triangle: Given two sides 
    ! A, B and an angle ALP opposite to one of the sides, computes the
    ! third side and the remaining two angles in the triangle.
    ! -------------------------------------------------------------------
    IMPLICIT NONE

    ! ---------------
    ! Input variables
    ! ---------------
    REAL (KIND=r8),    INTENT (IN) :: a, b, alp

    ! ----------------
    ! Output variables
    ! ----------------
    REAL (KIND=r8), INTENT (OUT) :: c, bet, gam

    ! ---------------
    ! Local variables
    ! ---------------
    REAL (KIND=r8) :: tmp1, tmp2

    ! --------------------------------------------------
    ! Initialize output variables to keep compiler happy
    ! --------------------------------------------------
    c = 0.0_r8 ; bet = 0.0_r8 ; gam = 0.0_r8

    ! -------------------------------------
    ! Compute the angle BET opposite side B
    ! -------------------------------------
    tmp1 = SIN(b) * SIN(alp) ; tmp2 = SIN(a)
    IF ( ABS(ABS(tmp1)-ABS(tmp2)) < eps ) THEN
       bet = pihalf
    ELSE IF ( ABS(tmp1) < eps ) THEN
       bet = 0.0_r8
    ELSE
       bet = ASIN(tmp1/tmp2)
    END IF

    ! -----------------------------------
    ! Now the length of C, the third side
    ! -----------------------------------
    tmp1 = TAN(b) * COS(alp) ; tmp2 = TAN(a) * COS(bet)
    c = ATAN(tmp1) + ATAN(tmp2)

    ! -------------------------------------------
    ! Finally the angle GAM between sides A and B
    ! -------------------------------------------
    ! Note: COT(x) = -TAN(x+pi/2)
    ! -----------------------------------------------------------
    ! COT(phi) = COS(b) * TAN(alp) ; COT(psi) = COS(a) * TAN(bet)
    ! gam = phi + psi
    ! -----------------------------------------------------------
    tmp1 = COS(b) * TAN(alp) ; tmp2 = COS(a) * TAN(bet)
    gam = ATAN(-tmp1) + ATAN(-tmp2) - pi
    IF ( ABS(ABS(gam)-pi) < eps ) gam = 0.0_r8

    RETURN
  END SUBROUTINE sphergeom_ssa

  SUBROUTINE sphergeom_sas ( a, b, gam, c, alp, bet )

    ! -------------------------------------------------------------------
    ! SAS: Side-Angle-Side 
    ! 
    ! Computation of quantities in a spherical triangle: Given two sides 
    ! A, B and the angle GAM inbetween them, computes the third side and
    ! one of the remaining two angles in the triangle.
    ! -------------------------------------------------------------------
    IMPLICIT NONE

    ! ---------------
    ! Input variables
    ! ---------------
    REAL (KIND=r8),    INTENT (IN) :: a, b, gam

    ! ----------------
    ! Output variables
    ! ----------------
    REAL (KIND=r8), INTENT (OUT) :: c, alp, bet

    ! ---------------
    ! Local variables
    ! ---------------
    REAL (KIND=r8) :: tmp1, tmp2

    ! --------------------------------------------------
    ! Initialize output variables to keep compiler happy
    ! --------------------------------------------------
    c = 0.0_r8 ; alp = 0.0_r8

    ! -----------------------------------------------------------
    ! Compute length of baseline c between the two points A and B
    ! -----------------------------------------------------------
    tmp1 = COS(a) * COS(b) + SIN(a) * SIN(b) * COS(gam)
    IF ( ABS(tmp1) < eps ) THEN
       c = pihalf
    ELSE
       c = ACOS(tmp1)
    END IF

    ! ----------------------------------------------------
    ! Now the angle ALP at point A (between sides a and c)
    ! ----------------------------------------------------
    tmp1 = COS(a) - COS(b) * COS(c) ; tmp2 = SIN(b) * SIN(c)
    IF ( tmp2 /= 0.0_r8 ) THEN
       IF ( ABS(tmp1/tmp2) < 1.0_r8 ) THEN
          alp = ACOS(tmp1/tmp2)
       ELSE
          alp = 0.0_r8
       END IF
    ELSE
       alp = 0.0_r8
    END IF

    ! ----------------------------------------------------
    ! Now the angle BET at point B (between sides b and c)
    ! ----------------------------------------------------
    tmp1 = COS(b) - COS(c) * COS(a) ; tmp2 = SIN(c) * SIN(a)
    IF ( tmp2 /= 0.0_r8 ) THEN
       IF ( ABS(tmp1/tmp2) < 1.0_r8 ) THEN
          bet = ACOS(tmp1/tmp2)
       ELSE
          bet = 0.0_r8
       END IF
    ELSE
       bet = 0.0_r8
    END IF

    RETURN
  END SUBROUTINE sphergeom_sas

  SUBROUTINE sphergeom_coordinates_at_point ( a0, b0, gam0, c_in, abs_or_frac, a, gam )

    ! -----------------------------------------------------------
    ! Finds the co-ordinates of the baseline extended from two
    ! lon/lat points on a sphere given the hypotenuse C_IN.
    ! -----------------------------------------------------------
    IMPLICIT NONE

    ! ---------------
    ! Input variables
    ! ---------------
    REAL (KIND=r8),    INTENT (IN) :: a0, b0, gam0, c_in
    CHARACTER (LEN=3), INTENT (IN) :: abs_or_frac

    ! ----------------
    ! Output variables
    ! ----------------
    REAL (KIND=r8), INTENT (OUT) :: a, gam

    ! ---------------
    ! Local variables
    ! ---------------
    REAL (KIND=r8) :: tmp1, tmp2, c0, c, alp0, gamsig

    ! --------------------------------------------------
    ! Initialize output variables to keep compiler happy
    ! --------------------------------------------------
    a = 0.0_r8 ; gam = 0.0_r8

    ! -------------------------------------------------
    ! Compute length of baseline between the two points
    ! -------------------------------------------------
    tmp1 = COS(a0) * COS(b0) + SIN(a0) * SIN(b0) * COS(gam0)
    IF ( ABS(tmp1) < eps ) THEN
       c0 = pihalf
    ELSE
       c0 = ACOS(tmp1)
    END IF

    tmp1 = COS(a0) - COS(b0) * COS(c0) ; tmp2 = SIN(b0) * SIN(c0)
    IF ( tmp2 /= 0.0_r8 ) THEN
       IF ( ABS(tmp1/tmp2) < 1.0_r8 ) THEN
          alp0 = ACOS(tmp1/tmp2)
       ELSE
          alp0 = 0.0_r8
       END IF
    ELSE
       alp0 = 0.0_r8
    END IF

    ! -------------------------------------------------------------
    ! Now that we have the basic parameters of the triangle, we can
    ! compute the lat/lon at the desired point.
    ! -------------------------------------------------------------
    SELECT CASE ( abs_or_frac )
    CASE ( 'abs' )
       c = c_in
    CASE ( 'frc' )
       c = c0 * c_in
    CASE DEFAULT
       ! ------------------------------------------------------------
       ! Temporary fix for unknown option: Complain and set c to c0/2
       ! ------------------------------------------------------------
       WRITE (*,'(A,A)') 'ERROR: Unknown option -- ', abs_or_frac
       c = c0 / 2.0_r8
    END SELECT

    ! -----------------------------
    ! The "latitude" of the mid-point
    ! -----------------------------
    tmp1 = COS(b0) * COS(c) + SIN(b0) * SIN(c) * COS(alp0)
    IF ( ABS(tmp1) < eps ) THEN
       a = pihalf
    ELSE
       a = ACOS(tmp1)
    END IF

    ! --------------------------------
    ! The "longitude" of the mid-point
    ! --------------------------------
    tmp1 = COS(c) - COS(a) * COS(b0)
    tmp2 = SIN(a) * SIN(b0)
    IF ( tmp2 /= 0.0_r8 ) THEN
       IF ( gam0 /= 0.0_r8 .AND. ABS(tmp1/tmp2) < 1.0_r8 ) THEN
          gam = ACOS(tmp1/tmp2)
       ELSE
          gam = 0.0_r8
       END IF
    ELSE
       gam = 0.0_r8
    END IF

    IF ( ((gam0 > 0.0_r8) .AND. (gam0 < pi)) .OR. (gam0 < -pi) ) THEN
       gamsig = 1.0_r8
    ELSE
       gamsig = -1.0_r8
    END IF

    ! --------------------------------------
    ! Return the angle with the correct sign
    ! --------------------------------------
    gam = ABS(gam) * gamsig

    RETURN
  END SUBROUTINE sphergeom_coordinates_at_point

  SUBROUTINE sphergeom_baseline_comp ( a0, b0, gam0, c0 )

    ! -------------------------------------------------------
    ! Finds the lengh of the baseline of a spherical triangle
    ! -------------------------------------------------------
    IMPLICIT NONE

    ! ---------------
    ! Input variables
    ! ---------------
    REAL (KIND=r8), INTENT (IN) :: a0, b0, gam0

    ! ---------------
    ! Output variable
    ! ---------------
    REAL (KIND=r8), INTENT (OUT) :: c0

    ! --------------
    ! Local variable
    ! --------------
    REAL (KIND=r8) :: tmp

    ! --------------------------------------------------
    ! Initialize output variables to keep compiler happy
    ! --------------------------------------------------
    c0 = 0.0_r8

    ! -------------------------------------------------
    ! Compute length of baseline between the two points
    ! -------------------------------------------------
    tmp = COS(a0) * COS(b0) + SIN(a0) * SIN(b0) * COS(gam0)
    IF ( ABS(tmp) < eps ) THEN
       c0 = pihalf
    ELSE
       c0 = ACOS(tmp)
    END IF

    RETURN
  END SUBROUTINE sphergeom_baseline_comp

  SUBROUTINE lonlat_to_pi ( lon, lat )

    IMPLICIT NONE

    ! ------------------
    ! Modified variables
    ! ------------------
    REAL (KIND=r8), INTENT (INOUT) :: lon, lat

    ! ------------------------------------
    ! Adjust longitude values to [-pi,+pi]
    ! ------------------------------------
    IF ( ABS(lon) > twopi ) lon = MOD(lon, twopi)
    IF ( lon >  pi ) lon = lon - twopi
    IF ( lon < -pi ) lon = lon + twopi
    ! ---------------------------------------
    ! Adjust latitude values to [-pi/2,+pi/2]
    ! ---------------------------------------
    IF ( ABS(lat) > pihalf ) lat = MOD(lat, pihalf)
    IF ( lat >  pihalf ) lat =   pi - lat
    IF ( lat < -pihalf ) lat = -(pi + lat)
    
    RETURN
  END SUBROUTINE lonlat_to_pi

  REAL (KIND=r8) FUNCTION angle_minus_twopi ( gamma0, pival ) RESULT ( gamma )

    IMPLICIT NONE

    ! ---------------
    ! Input variables
    ! ---------------
    REAL (KIND=r8), INTENT (IN) :: gamma0, pival

    IF ( ABS(gamma0) > pival ) THEN
       gamma = gamma0 - 2.0_r8 * pival !SIGN(2.0_r8*pival - gamma0, gamma0)
    ELSE
       gamma = gamma0
    END IF

    RETURN
  END FUNCTION angle_minus_twopi


  SUBROUTINE sphgeo_comp_pixel_corners ( &
       nxtrack, ntimes, lon, lat, clon, clat, estat, yn_omi_pixel_adjust_k )

    ! ------------------------------------------------------------------------
    ! Compute corner coordinates of ground pixels given only the pixel centers
    ! ------------------------------------------------------------------------

    IMPLICIT NONE

    ! ---------------
    ! Input variables
    ! ---------------
    INTEGER (KIND=i4),                                    INTENT (IN) :: nxtrack, ntimes
    REAL    (KIND=r4), DIMENSION (1:nxtrack, 0:ntimes-1), INTENT (IN) :: lon, lat
    LOGICAL, OPTIONAL,                                    INTENT (IN) :: yn_omi_pixel_adjust_k

    ! -------------------------
    ! Output/modified variables
    ! -------------------------------------------------------------------------
    ! NOTE: The corner coordinates are R8 rather that R4!
    !       (they have been initialized to MISSING earlier; hence the "INOUT")
    ! -------------------------------------------------------------------------
    INTEGER (KIND=i4),                                  INTENT (OUT)   :: estat
    REAL    (KIND=r8), DIMENSION (0:nxtrack, 0:ntimes), INTENT (INOUT) :: clon, clat

    ! ---------------
    ! Local variables
    ! ---------------
    INTEGER (KIND=i4)                                   :: i, j
    REAL    (KIND=r8)                                   :: a0, b0, c0, gam0, a, gam
    REAL    (KIND=r8), DIMENSION (1:nxtrack,0:ntimes-1) :: lonrad, latrad
    REAL    (KIND=r8), DIMENSION (0:nxtrack,0:ntimes)   :: tmplat, tmplon
    LOGICAL                                             :: yn_pixel_adjust_crosstrack

    estat = pge_errstat_ok

    ! --------------------------------------------------------------------------------
    ! Check whether we want to perform some pixel adjustment. This is an issue for
    ! OMI, where across-track pixels are of different size (larger towards the edges of
    ! the swath). Using the "center" method between adjacent across-track pixels will
    ! introduce a small error in the pixel boundaries, and the pixels at the sides of
    ! the swath will be underestimated. 
    !
    ! Note that this adjustment is still experimental until the distortions around the
    ! poles have been solved. Hence the default is not to make this adjustment.
    ! --------------------------------------------------------------------------------
    yn_pixel_adjust_crosstrack = .FALSE.
    IF ( PRESENT (yn_omi_pixel_adjust_k) ) yn_pixel_adjust_crosstrack = yn_omi_pixel_adjust_k

    ! ------------------------------------------------------------------
    ! Convert geolocation to radians; do everything in R8 rather than R4
    ! ------------------------------------------------------------------
    lonrad = REAL ( lon, KIND=r8 ) * deg2rad
    latrad = REAL ( lat, KIND=r8 ) * deg2rad

    ! -------------------------
    ! Initialize some variables
    ! -------------------------
    !tmplon = -999.9_r8 ; tmplat = -999.9_r8
    tmplon =     0.0_r8 ; tmplat =    0.0_r8

    a0   = 0.0_r8 ; b0 = 0.0_r8 ; c0  = 0.0_r8
    gam0 = 0.0_r8 ; a  = 0.0_r8 ; gam = 0.0_r8

    ! ---------------------------------------------------
    ! First interpolate between center points ALONG-TRACK
    ! ---------------------------------------------------
    ! Array points filled: [1:nxtrack, 0:ntimes]
    ! ---------------------------------------------------
    DO i = 1, nxtrack
       DO j = 0, ntimes -2
          a0   = pihalf - latrad(i,j  )
          b0   = pihalf - latrad(i,j+1) 
          gam0 = angle_minus_twopi ( lonrad(i,j) - lonrad(i,j+1), pi )

          CALL sphergeom_coordinates_at_point ( a0, b0, gam0, 0.5_r8, 'frc', a, gam )
          tmplat(i,j+1) = pihalf - a
          tmplon(i,j+1) = lonrad(i,j) - gam
          CALL lonlat_to_pi (tmplon(i,j+1), tmplat(i,j+1) )
       END DO

       j = 0
       a0 = pihalf - tmplat(i,j+1)
       b0 = pihalf - tmplat(i,j+2)
       gam0 = angle_minus_twopi ( tmplon(i,j+1) - tmplon(i,j+2), pi )
       CALL sphergeom_coordinates_at_point ( a0, b0, gam0, 2.0_r8, 'frc', a, gam )
       tmplat(i,j) = pihalf - a
       tmplon(i,j) = tmplon(i,j+2) + gam
       CALL lonlat_to_pi (tmplon(i,j), tmplat(i,j) )

       j = ntimes
       a0 = pihalf - tmplat(i,j-1)
       b0 = pihalf - tmplat(i,j-2)
       gam0 = angle_minus_twopi ( tmplon(i,j-1) - tmplon(i,j-2), pi )
       CALL sphergeom_coordinates_at_point ( a0, b0, gam0, 2.0_r8, 'frc', a, gam )
       tmplat(i,j) = pihalf - a
       tmplon(i,j) = tmplon(i,j-2) + gam
       CALL lonlat_to_pi (tmplon(i,j), tmplat(i,j) )
    END DO

    ! ----------------------------------------------------------
    ! Now interpolate ACROSS-TRACK whatever points are available
    ! ----------------------------------------------------------  
    DO j = 0, ntimes
       ! ----------------------------------
       ! First the center line of the swath
       ! ----------------------------------------------------------------------
       ! To be save, we compute ALL possible pixels using the center cross
       ! method. Only in favorable circumstances do we use a more sophisticated
       ! approach, and having all pixels precomputed saves us some IF THEN ELSE
       ! logic and several lines of code.
       ! ----------------------------------------------------------------------
       DO i = 1, nxtrack-1
          a0   = pihalf - tmplat(i+1,j)
          b0   = pihalf - tmplat(i,  j) 
          gam0 = angle_minus_twopi ( tmplon(i+1,j) - tmplon(i,j), pi )
          CALL sphergeom_coordinates_at_point ( a0, b0, gam0, 0.5_r8, 'frc', a, gam )
          clat(i,j) = pihalf - a
          clon(i,j) = tmplon(i+1,j) - gam
          CALL lonlat_to_pi (clon(i,j), clat(i,j) )

          IF ( i == 1 ) THEN
             a0 = pihalf - tmplat(i,j)
             b0 = pihalf -   clat(i,j)
             gam0 = angle_minus_twopi ( tmplon(i,j) - clon(i,j), pi )
             CALL sphergeom_coordinates_at_point ( a0, b0, gam0, 2.0_r8, 'frc', a, gam )
             clat(i-1,j) = pihalf - a
             clon(i-1,j) = clon(i,j) + gam
             CALL lonlat_to_pi (clon(i-1,j), clat(i-1,j) )
          END IF
          IF ( i == nxtrack-1 ) THEN
             a0 = pihalf - tmplat(i+1,j)
             b0 = pihalf -   clat(i,j)
             gam0 = angle_minus_twopi ( tmplon(i+1,j) - clon(i,j), pi )
             CALL sphergeom_coordinates_at_point ( a0, b0, gam0, 2.0_r8, 'frc', a, gam )
             clat(i+1,j) = pihalf - a
             clon(i+1,j) = clon(i,j) + gam
             CALL lonlat_to_pi (clon(i+1,j), clat(i+1,j) )
          END IF
       END DO
    END DO

    !!! ===============================================================================!!!
    !!! The following code is experimental and doesn't quite work yet. It is supposed  !!!
    !!! to adjust the corner coordinates according to the true distance between the    !!!
    !!! pixel centers, and it does so by working outwards from the center of the swath !!!
    !!! to the across-track edges. Known issues are:                                   !!!
    !!!   * At the along-track edges, we seem to be mis-aligning the pixels and are    !!!
    !!!     losing the proper last line of corners                                     !!!
    !!!   * At the poles, the distortion of pixels progresses outwards, resulting in   !!!
    !!!     heavy distortions past the poles (going outward). This can be ameliorated  !!!
    !!!     somewhat by excluding pole-most latitudes (within 3 deg of pole), but this !!!
    !!!     still leaves some undesirable distortions.                                 !!!
    !!! ===============================================================================!!!
    IF ( yn_pixel_adjust_crosstrack ) THEN
       DO j = 0, ntimes
          ! -------------------------------------------------------------------
          ! From the center of the swath to lower cross-track pixel numbers
          ! -------------------------------------------------------------------
          DO i = nxtrack/2-1, 1, -1
             a0 = pihalf -   clat(i+1,j)
             b0 = pihalf - tmplat(i+1,j) 
             gam0 = clon(i+1,j) - tmplon(i+1,j)
             gam0 = angle_minus_twopi ( gam0, pi )
             CALL sphergeom_baseline_comp ( a0, b0, gam0, c0 )
             a0 = pihalf - tmplat(i+1,j)
             b0 = pihalf - tmplat(i,  j) 
             gam0 = angle_minus_twopi ( tmplon(i+1,j) - tmplon(i,j), pi )
             CALL sphergeom_coordinates_at_point ( a0, b0, gam0, c0, 'abs', a, gam )
             clat(i,j) = pihalf - a
             clon(i,j) = tmplon(i+1,j) - gam
             CALL lonlat_to_pi (clon(i,j), clat(i,j) )
          END DO
          i = 0
          a0 = pihalf - tmplat(i+1,j)
          b0 = pihalf -   clat(i+1,j)
          gam0 = tmplon(i+1,j) - clon(i+1,j)
          CALL sphergeom_coordinates_at_point ( a0, b0, gam0, 2.0_r8, 'frc', a, gam )
          clat(i,j) = pihalf - a
          clon(i,j) = clon(i+1,j) + gam
          CALL lonlat_to_pi (clon(i,j), clat(i,j) )
       
       ! -----------------------------------------------------------------------
       ! From the center of the swath to higer cross-track pixel numbers
       ! -----------------------------------------------------------------------
          DO i = nxtrack/2+1, nxtrack-1
             a0 = pihalf - tmplat(i,  j)
             b0 = pihalf -   clat(i-1,j)
             gam0 = tmplon(i,j) - clon(i-1,j)
             gam0 = angle_minus_twopi ( gam0, pi )
             !CALL sphergeom_baseline_comp ( a0, b0, gam0, c0 )
             !IF ( j == ntimes .AND. i == nxtrack ) WRITE (*,'(4F10.6)') a0,b0,gam0,c0
             a0 = pihalf - tmplat(i,  j)
             b0 = pihalf - tmplat(i+1,j) 
             gam0 = angle_minus_twopi ( tmplon(i,j) - tmplon(i+1,j), pi )
             CALL sphergeom_coordinates_at_point ( a0, b0, gam0, c0, 'abs', a, gam )
             clat(i,j) = pihalf - a
             clon(i,j) = tmplon(i,j) - gam
             CALL lonlat_to_pi (clon(i,j), clat(i,j) )
          END DO
          i = nxtrack
          a0 = pihalf - tmplat(i,  j)
          b0 = pihalf - clat  (i-1,j)
          gam0 = tmplon(i,j) - clon(i-1,j)
          gam0 = angle_minus_twopi ( gam0, pi )
          CALL sphergeom_coordinates_at_point ( a0, b0, gam0, 2.0_r8, 'frc', a, gam )
          clat(i,j) = pihalf - a
          clon(i,j) = clon(i-1,j) + gam
          CALL lonlat_to_pi (clon(i,j), clat(i,j) )
       END DO
    END IF

    ! ------------------------------------
    ! Convert lons and lats back to Degree
    ! ------------------------------------
    clon = clon * rad2deg
    clat = clat * rad2deg

    RETURN
  END SUBROUTINE sphgeo_comp_pixel_corners

  !SUBROUTINE pixcorner_latlon_read ( ntimes, nxtrack, lons, lats, errstat )
  !
  !  IMPLICIT NONE
  !
  !  ! ---------------
  !  ! Input variables
  !  ! ---------------
  !  INTEGER (KIND=i4), INTENT (IN) :: ntimes, nxtrack
  !
  !  ! -----------------
  !  ! Modified variable
  !  ! -----------------
  !  INTEGER (KIND=i4), INTENT (INOUT) :: errstat
  !
  !  ! ----------------
  !  ! Output variables
  !  ! ----------------
  !  REAL    (KIND=r4), DIMENSION (nxtrack,0:ntimes-1), INTENT (OUT) :: lons, lats
  !
  !  ! ---------------
  !  ! Local variables
  !  ! ---------------
  !  INTEGER (KIND=i4)            :: locerrstat, iline, ntimes_loop
  !
  !  ! ---------------------------
  !  ! Initialize output variables
  !  ! ---------------------------
  !  lons = r4_missval ; lats = r4_missval
  !
  !  ! -------------------------------------------------------
  !  ! Loop over all lines in the file in blocks of NLINES_MAX
  !  ! -------------------------------------------------------
  !  ScanLines: DO iline = 0, ntimes-1, nlines_max
  !
  !     ! --------------------------------------------------------
  !     ! Check if loop ends before n_times_loop max is exhausted.
  !     ! --------------------------------------------------------
  !     ntimes_loop = MIN ( nlines_max, ntimes-iline )
  !
  !     ! ----------------------------------------------------
  !     ! Read current data block fitting output from HE5 file
  !     ! ----------------------------------------------------
  !     he5_start_2d  = (/       0,       iline /)
  !     he5_stride_2d = (/       1,           1 /)
  !     he5_edge_2d   = (/ nxtrack, ntimes_loop /)
  !  
  !     ! --------------------------------------
  !     ! Latitudes and longitudes
  !     ! --------------------------------------
  !     locerrstat = HE5_SWrdfld ( pge_swath_id, TRIM(ADJUSTL(lon_field)),         &
  !          he5_start_2d, he5_stride_2d, he5_edge_2d, lons(1:nxtrack,iline:iline+ntimes_loop-1) )
  !     locerrstat = HE5_SWrdfld ( pge_swath_id, TRIM(ADJUSTL(lat_field)),        &
  !          he5_start_2d, he5_stride_2d, he5_edge_2d, lats(1:nxtrack,iline:iline+ntimes_loop-1) )
  !
  !  END DO ScanLines
  !
  !  RETURN
  !END SUBROUTINE pixcorner_latlon_read

END MODULE OMSAO_pixelcorner_module
