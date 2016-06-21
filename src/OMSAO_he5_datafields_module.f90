MODULE OMSAO_he5_datafields_module

  USE OMSAO_precision_module
  USE OMSAO_indices_module,    ONLY: o3_t1_idx, o3_t3_idx
  USE OMSAO_parameters_module, ONLY: maxchlen
  USE OMSAO_he5_module

  IMPLICIT NONE

  !INCLUDE 'hdfeos5.inc'

  ! ******************************************************
  ! In this MODULE we define the HE5 Data Fields that will
  ! be written to the HE5 output swath.
  ! ******************************************************

  INTEGER (KIND=i4), PARAMETER, PRIVATE :: maxrank = 3

  ! -------------------------------------
  ! TYPE declaration for HE5 Swath fields
  ! -------------------------------------
  ! (Some fields not used)
  ! ----------------------
  TYPE, PUBLIC :: DataField_HE5
     REAL      (KIND=r8)                      :: FillValue, MissingValue, ScaleFactor, Offset, SpecTemp
     REAL      (KIND=r8), DIMENSION (2)       :: ValidRange
     CHARACTER (LEN=maxchlen)                 :: Name, Dimensions, Units, Title, UniqueFD
     INTEGER   (KIND=i4)                      :: LenUnits, LenTitle, LenUniqueFD
     INTEGER   (KIND=i4)                      :: Swath_ID
     INTEGER   (KIND=i4)                      :: HE5_DataType
     INTEGER   (KIND=i4)                      :: Rank
     INTEGER   (KIND=i4), DIMENSION (maxrank) :: Dims
  END TYPE DataField_HE5


  ! --------------------------
  ! Composite HE5 Swath Fields
  ! --------------------------
  ! * Geolocation (Scan Number, Time, Lat/Lon, Angles)
  TYPE (DataField_HE5), DIMENSION (n_gfields)             :: geo_he5fields
  ! * Diagnostic data fields
  TYPE (DataField_HE5), DIMENSION (n_diag_fields)         :: diagnostic_he5fields
  ! * Common data fields
  TYPE (DataField_HE5), DIMENSION (n_cdfields)            :: comdata_he5fields
  ! * Solar and radiance wavelength calibration
  TYPE (DataField_HE5), DIMENSION (n_solcal_fields)       :: sol_calfit_he5fields
  TYPE (DataField_HE5), DIMENSION (n_radcal_fields)       :: rad_calfit_he5fields
  TYPE (DataField_HE5), DIMENSION (n_radref_fields)       :: rad_reffit_he5fields
  ! * Special fields for OMHCHO and OMCHOCHO
  TYPE (DataField_HE5), DIMENSION (n_voc_fields)          :: voc_he5fields
  ! * Special fields for OMHCHO and OMCHOCHO
  TYPE (DataField_HE5), DIMENSION (n_cld_fields)          :: cld_he5fields
  ! * Special fields for wavelength-modified AMF fields
  TYPE (DataField_HE5), DIMENSION (n_wmamf_fields)        :: wmamf_he5fields
  ! * Special fields for OMSAO3
  TYPE (DataField_HE5), DIMENSION (o3_t1_idx:o3_t3_idx,2) :: o3_prefit_he5fields
  ! * Special fields for Reference Sector correction
  TYPE (DataField_He5), DIMENSION (n_rs_fields)           :: rs_he5fields
  ! * Special fields for Scattering Weights, Gas Profile
  !   Averaging Kernels and albedo
  TYPE (DataField_He5), DIMENSION (n_sw_fields)           :: sw_he5fields


  CONTAINS

    SUBROUTINE he5_initialize_datafields (  )
  
      USE OMSAO_parameters_module,   ONLY: r8_missval, r4_missval, i4_missval, i2_missval, str_missval
      USE OMSAO_radiance_ref_module, ONLY: yn_radiance_reference
      USe OMSAO_destriping_module,   ONLY: yn_run_destriping

      IMPLICIT NONE

      ! ---------------
      ! Input Variables
      ! ---------------

      ! ---------------
      ! Local Variables
      ! ---------------
      INTEGER (KIND=i4) :: i, j, numtype
      REAL    (KIND=r8) :: missval

      ! -----------------------------
      ! Initialize Geolocation Fields
      ! -----------------------------
      DO i = 1, n_gfields
         ! ----------------------------------------
         ! Find the numeric type of the data field.
         ! ----------------------------------------
         CALL he5_find_datatype ( TRIM(ADJUSTL(geo_field_specs(4,i))), numtype, missval )

         geo_he5fields(i)%ScaleFactor    = 1.0_r8
         geo_he5fields(i)%Offset         = 0.0_r8
         geo_he5fields(i)%Name           = TRIM(ADJUSTL(geo_field_names(1,i)))
         geo_he5fields(i)%Title          = TRIM(ADJUSTL(geo_field_names(2,i)))
         geo_he5fields(i)%Units          = TRIM(ADJUSTL(geo_field_specs(1,i)))
         geo_he5fields(i)%Dimensions     = TRIM(ADJUSTL(geo_field_specs(2,i)))
         geo_he5fields(i)%UniqueFD       = TRIM(ADJUSTL(geo_field_specs(3,i)))
         geo_he5fields(i)%FillValue      = missval
         geo_he5fields(i)%MissingValue   = missval
         geo_he5fields(i)%Swath_ID       = pge_swath_id
         geo_he5fields(i)%HE5_DataType   = numtype
         geo_he5fields(i)%Rank           = 1
         geo_he5fields(i)%Dims           = (/ -1, -1, -1 /)
         geo_he5fields(i)%LenTitle       = LEN_TRIM(ADJUSTL(geo_field_names(2,i)))
         geo_he5fields(i)%LenUnits       = LEN_TRIM(ADJUSTL(geo_field_specs(1,i)))
         geo_he5fields(i)%LenUniqueFD    = LEN_TRIM(ADJUSTL(geo_field_specs(3,i)))
         geo_he5fields(i)%ValidRange     = geo_valids(1:2,i)
      END DO

      ! ---------------------------------------------------------------------------
      ! Initialize TYPE fields fields for solar and radiance wavelength calibration
      ! as well as radiance reference parameters.
      ! ---------------------------------------------------------------------------
      DO i = 1, n_solcal_fields
         ! ----------------------------------------
         ! Find the numeric type of the data field.
         ! ----------------------------------------
         CALL he5_find_datatype ( TRIM(ADJUSTL(solcal_field_specs(4,i))), numtype, missval )

         sol_calfit_he5fields(i)%ScaleFactor    = 1.0_r8
         sol_calfit_he5fields(i)%Offset         = 0.0_r8
         sol_calfit_he5fields(i)%Name           = TRIM(ADJUSTL(solcal_field_names(1,i)))
         sol_calfit_he5fields(i)%Title          = TRIM(ADJUSTL(solcal_field_names(2,i)))
         sol_calfit_he5fields(i)%Units          = TRIM(ADJUSTL(solcal_field_specs(1,i)))
         sol_calfit_he5fields(i)%Dimensions     = TRIM(ADJUSTL(solcal_field_specs(2,i)))
         sol_calfit_he5fields(i)%UniqueFD       = TRIM(ADJUSTL(solcal_field_specs(3,i)))
         sol_calfit_he5fields(i)%FillValue      = missval
         sol_calfit_he5fields(i)%MissingValue   = missval
         sol_calfit_he5fields(i)%Swath_ID       = pge_swath_id
         sol_calfit_he5fields(i)%HE5_DataType   = numtype
         sol_calfit_he5fields(i)%Rank           = 1
         sol_calfit_he5fields(i)%Dims           = (/ -1, -1, -1 /)
         sol_calfit_he5fields(i)%LenTitle       = LEN_TRIM(ADJUSTL(solcal_field_names(2,i)))
         sol_calfit_he5fields(i)%LenUnits       = LEN_TRIM(ADJUSTL(solcal_field_specs(1,i)))
         sol_calfit_he5fields(i)%LenUniqueFD    = LEN_TRIM(ADJUSTL(solcal_field_specs(3,i)))
         sol_calfit_he5fields(i)%ValidRange     = solcal_valids(1:2,i)
      END DO
      DO i = 1, n_radcal_fields
         ! ----------------------------------------
         ! Find the numeric type of the data field.
         ! ----------------------------------------
         CALL he5_find_datatype ( TRIM(ADJUSTL(radcal_field_specs(4,i))), numtype, missval )

         rad_calfit_he5fields(i)%ScaleFactor    = 1.0_r8
         rad_calfit_he5fields(i)%Offset         = 0.0_r8
         rad_calfit_he5fields(i)%Name           = TRIM(ADJUSTL(radcal_field_names(1,i)))
         rad_calfit_he5fields(i)%Title          = TRIM(ADJUSTL(radcal_field_names(2,i)))
         rad_calfit_he5fields(i)%Units          = TRIM(ADJUSTL(radcal_field_specs(1,i)))
         rad_calfit_he5fields(i)%Dimensions     = TRIM(ADJUSTL(radcal_field_specs(2,i)))
         rad_calfit_he5fields(i)%UniqueFD       = TRIM(ADJUSTL(radcal_field_specs(3,i)))
         rad_calfit_he5fields(i)%FillValue      = missval
         rad_calfit_he5fields(i)%MissingValue   = missval
         rad_calfit_he5fields(i)%Swath_ID       = pge_swath_id
         rad_calfit_he5fields(i)%HE5_DataType   = numtype
         rad_calfit_he5fields(i)%Rank           = 1
         rad_calfit_he5fields(i)%Dims           = (/ -1, -1, -1 /)
         rad_calfit_he5fields(i)%LenTitle       = LEN_TRIM(ADJUSTL(radcal_field_names(2,i)))
         rad_calfit_he5fields(i)%LenUnits       = LEN_TRIM(ADJUSTL(radcal_field_specs(1,i)))
         rad_calfit_he5fields(i)%LenUniqueFD    = LEN_TRIM(ADJUSTL(radcal_field_specs(3,i)))
         rad_calfit_he5fields(i)%ValidRange     = radcal_valids(1:2,i)
      END DO
      IF (yn_radiance_reference) THEN
         DO i = 1, n_radref_fields
            ! ----------------------------------------
            ! Find the numeric type of the data field.
            ! ----------------------------------------
            CALL he5_find_datatype ( TRIM(ADJUSTL(radref_field_specs(4,i))), numtype, missval )
            
            rad_reffit_he5fields(i)%ScaleFactor    = 1.0_r8
            rad_reffit_he5fields(i)%Offset         = 0.0_r8
            rad_reffit_he5fields(i)%Name           = TRIM(ADJUSTL(radref_field_names(1,i)))
            rad_reffit_he5fields(i)%Title          = TRIM(ADJUSTL(radref_field_names(2,i)))
            rad_reffit_he5fields(i)%Units          = TRIM(ADJUSTL(radref_field_specs(1,i)))
            rad_reffit_he5fields(i)%Dimensions     = TRIM(ADJUSTL(radref_field_specs(2,i)))
            rad_reffit_he5fields(i)%UniqueFD       = TRIM(ADJUSTL(radref_field_specs(3,i)))
            rad_reffit_he5fields(i)%FillValue      = missval
            rad_reffit_he5fields(i)%MissingValue   = missval
            rad_reffit_he5fields(i)%Swath_ID       = pge_swath_id
            rad_reffit_he5fields(i)%HE5_DataType   = numtype
            rad_reffit_he5fields(i)%Rank           = 1
            rad_reffit_he5fields(i)%Dims           = (/ -1, -1, -1 /)
            rad_reffit_he5fields(i)%LenTitle       = LEN_TRIM(ADJUSTL(radref_field_names(2,i)))
            rad_reffit_he5fields(i)%LenUnits       = LEN_TRIM(ADJUSTL(radref_field_specs(1,i)))
            rad_reffit_he5fields(i)%LenUniqueFD    = LEN_TRIM(ADJUSTL(radref_field_specs(3,i)))
            rad_reffit_he5fields(i)%ValidRange     = radref_valids(1:2,i)
         END DO
      ENDIF

      ! -----------------------------
      ! Initialize Common Data Fields
      ! -----------------------------
      DO i = 1, n_cdfields
         ! ----------------------------------------------
         ! If we are not using destriping correction skip
         ! ColumnAmountDestriped field (number 12)
         ! ----------------------------------------------
         IF (i .EQ. 12 .AND. (.NOT. yn_run_destriping)) CYCLE

         ! ----------------------------------------
         ! Find the numeric type of the data field.
         ! ----------------------------------------
         CALL he5_find_datatype ( TRIM(ADJUSTL(comdata_field_specs(4,i))), numtype, missval )

         comdata_he5fields(i)%ScaleFactor    = 1.0_r8
         comdata_he5fields(i)%Offset         = 0.0_r8
         comdata_he5fields(i)%Name           = TRIM(ADJUSTL(comdata_field_names(1,i)))
         comdata_he5fields(i)%Title          = TRIM(ADJUSTL(comdata_field_names(2,i)))
         comdata_he5fields(i)%Units          = TRIM(ADJUSTL(comdata_field_specs(1,i)))
         comdata_he5fields(i)%Dimensions     = TRIM(ADJUSTL(comdata_field_specs(2,i)))
         comdata_he5fields(i)%UniqueFD       = TRIM(ADJUSTL(comdata_field_specs(3,i)))
         comdata_he5fields(i)%FillValue      = missval
         comdata_he5fields(i)%MissingValue   = missval
         comdata_he5fields(i)%Swath_ID       = pge_swath_id
         comdata_he5fields(i)%HE5_DataType   = numtype
         comdata_he5fields(i)%Rank           = 1
         comdata_he5fields(i)%Dims           = (/ -1, -1, -1 /)
         comdata_he5fields(i)%LenTitle       = LEN_TRIM(ADJUSTL(comdata_field_names(2,i)))
         comdata_he5fields(i)%LenUnits       = LEN_TRIM(ADJUSTL(comdata_field_specs(1,i)))
         comdata_he5fields(i)%LenUniqueFD    = LEN_TRIM(ADJUSTL(comdata_field_specs(3,i)))
         comdata_he5fields(i)%ValidRange     = comdata_valids(1:2,i)
      END DO

      ! ---------------------------------
      ! Initialize Diagnostic Data Fields
      ! ---------------------------------
      DO i = 1, n_diag_fields
         ! ----------------------------------------
         ! Find the numeric type of the data field.
         ! ----------------------------------------
         CALL he5_find_datatype ( TRIM(ADJUSTL(diagnostic_field_specs(4,i))), numtype, missval )

         diagnostic_he5fields(i)%ScaleFactor    = 1.0_r8
         diagnostic_he5fields(i)%Offset         = 0.0_r8
         diagnostic_he5fields(i)%Name           = TRIM(ADJUSTL(diagnostic_field_names(1,i)))
         diagnostic_he5fields(i)%Title          = TRIM(ADJUSTL(diagnostic_field_names(2,i)))
         diagnostic_he5fields(i)%Units          = TRIM(ADJUSTL(diagnostic_field_specs(1,i)))
         diagnostic_he5fields(i)%Dimensions     = TRIM(ADJUSTL(diagnostic_field_specs(2,i)))
         diagnostic_he5fields(i)%UniqueFD       = TRIM(ADJUSTL(diagnostic_field_specs(3,i)))
         diagnostic_he5fields(i)%FillValue      = missval
         diagnostic_he5fields(i)%MissingValue   = missval
         diagnostic_he5fields(i)%Swath_ID       = pge_swath_id
         diagnostic_he5fields(i)%HE5_DataType   = numtype
         diagnostic_he5fields(i)%Rank           = 1
         diagnostic_he5fields(i)%Dims           = (/ -1, -1, -1 /)
         diagnostic_he5fields(i)%LenTitle       = LEN_TRIM(ADJUSTL(diagnostic_field_names(2,i)))
         diagnostic_he5fields(i)%LenUnits       = LEN_TRIM(ADJUSTL(diagnostic_field_specs(1,i)))
         diagnostic_he5fields(i)%LenUniqueFD    = LEN_TRIM(ADJUSTL(diagnostic_field_specs(3,i)))
         diagnostic_he5fields(i)%ValidRange     = diagnostic_valids(1:2,i)
      END DO


      ! -------------------------------------------
      ! OMHCHO/OMCHOCHO special fields (AMF clouds)
      ! -------------------------------------------
      DO i = 1, n_voc_fields
         ! ----------------------------------------
         ! Find the numeric type of the data field.
         ! ----------------------------------------
         CALL he5_find_datatype ( TRIM(ADJUSTL(voc_field_specs(4,i))), numtype, missval )

         voc_he5fields(i)%ScaleFactor    = 1.0_r8
         voc_he5fields(i)%Offset         = 0.0_r8
         voc_he5fields(i)%Name           = TRIM(ADJUSTL(voc_field_names(1,i)))
         voc_he5fields(i)%Title          = TRIM(ADJUSTL(voc_field_names(2,i)))
         voc_he5fields(i)%Units          = TRIM(ADJUSTL(voc_field_specs(1,i)))
         voc_he5fields(i)%Dimensions     = TRIM(ADJUSTL(voc_field_specs(2,i)))
         voc_he5fields(i)%UniqueFD       = TRIM(ADJUSTL(voc_field_specs(3,i)))
         voc_he5fields(i)%FillValue      = missval
         voc_he5fields(i)%MissingValue   = missval
         voc_he5fields(i)%Swath_ID       = pge_swath_id
         voc_he5fields(i)%HE5_DataType   = numtype
         voc_he5fields(i)%Rank           = 1
         voc_he5fields(i)%Dims           = (/ -1, -1, -1 /)
         voc_he5fields(i)%LenTitle       = LEN_TRIM(ADJUSTL(voc_field_names(2,i)))
         voc_he5fields(i)%LenUnits       = LEN_TRIM(ADJUSTL(voc_field_specs(1,i)))
         voc_he5fields(i)%LenUniqueFD    = LEN_TRIM(ADJUSTL(voc_field_specs(3,i)))
         voc_he5fields(i)%ValidRange     = voc_valids(1:2,i)
      END DO


      ! ---------------------------------------------------------
      ! Wavelength-modified AMF fitting fields (SLANT quantities)
      ! ---------------------------------------------------------
      DO i = 1, n_wmamf_fields
         ! ----------------------------------------
         ! Find the numeric type of the data field.
         ! ----------------------------------------
         CALL he5_find_datatype ( TRIM(ADJUSTL(wmamf_field_specs(4,i))), numtype, missval )

         wmamf_he5fields(i)%ScaleFactor    = 1.0_r8
         wmamf_he5fields(i)%Offset         = 0.0_r8
         wmamf_he5fields(i)%Name           = TRIM(ADJUSTL(wmamf_field_names(1,i)))
         wmamf_he5fields(i)%Title          = TRIM(ADJUSTL(wmamf_field_names(2,i)))
         wmamf_he5fields(i)%Units          = TRIM(ADJUSTL(wmamf_field_specs(1,i)))
         wmamf_he5fields(i)%Dimensions     = TRIM(ADJUSTL(wmamf_field_specs(2,i)))
         wmamf_he5fields(i)%UniqueFD       = TRIM(ADJUSTL(wmamf_field_specs(3,i)))
         wmamf_he5fields(i)%FillValue      = missval
         wmamf_he5fields(i)%MissingValue   = missval
         wmamf_he5fields(i)%Swath_ID       = pge_swath_id
         wmamf_he5fields(i)%HE5_DataType   = numtype
         wmamf_he5fields(i)%Rank           = 1
         wmamf_he5fields(i)%Dims           = (/ -1, -1, -1 /)
         wmamf_he5fields(i)%LenTitle       = LEN_TRIM(ADJUSTL(wmamf_field_names(2,i)))
         wmamf_he5fields(i)%LenUnits       = LEN_TRIM(ADJUSTL(wmamf_field_specs(1,i)))
         wmamf_he5fields(i)%LenUniqueFD    = LEN_TRIM(ADJUSTL(wmamf_field_specs(3,i)))
         wmamf_he5fields(i)%ValidRange     = wmamf_valids(1:2,i)
      END DO

      ! -----------------------------------
      ! OMSAO3 has additional output fields
      ! -----------------------------------
      DO i = o3_t1_idx, o3_t3_idx
         DO j = 1, 2
            numtype = HE5T_NATIVE_DOUBLE
            missval = r8_missval
            !o3_prefit_he5fields(i,j)%SpecTemp       = missval ! this will be assigned elsewhere
            o3_prefit_he5fields(i,j)%ScaleFactor    = 1.0_r8
            o3_prefit_he5fields(i,j)%Offset         = 0.0_r8
            o3_prefit_he5fields(i,j)%Name           = TRIM(ADJUSTL(o3_prefit_fields       (i,j)))
            o3_prefit_he5fields(i,j)%Dimensions     = TRIM(ADJUSTL(o3_prefit_fields_dims  (i,j)))
            o3_prefit_he5fields(i,j)%Units          = TRIM(ADJUSTL(o3_prefit_fields_units (i,j)))
            o3_prefit_he5fields(i,j)%Title          = TRIM(ADJUSTL(o3_prefit_fields_titles(i,j)))
            o3_prefit_he5fields(i,j)%UniqueFD       = TRIM(ADJUSTL(o3_prefit_fields_ufd   (i,j)))
            o3_prefit_he5fields(i,j)%FillValue      = missval
            o3_prefit_he5fields(i,j)%MissingValue   = missval
            o3_prefit_he5fields(i,j)%Swath_ID       = pge_swath_id
            o3_prefit_he5fields(i,j)%HE5_DataType   = numtype
            o3_prefit_he5fields(i,j)%Rank           = 1
            o3_prefit_he5fields(i,j)%Dims           = (/ -1, -1, -1 /)
            o3_prefit_he5fields(i,j)%LenUnits       = LEN_TRIM(ADJUSTL(o3_prefit_fields_units (i,j)))
            o3_prefit_he5fields(i,j)%LenTitle       = LEN_TRIM(ADJUSTL(o3_prefit_fields_titles(i,j)))
            o3_prefit_he5fields(i,j)%LenUniqueFD    = LEN_TRIM(ADJUSTL(o3_prefit_fields_ufd   (i,j)))
            o3_prefit_he5fields(i,j)%ValidRange     = o3_prefit_valids(1:2,i,j)
         END DO
      END DO

      ! -----------------------
      ! Reference Sector fields
      ! -----------------------
      DO i = 1, n_rs_fields
         ! ----------------------------------------
         ! Find the numeric type of the data field.
         ! ----------------------------------------
         CALL he5_find_datatype ( TRIM(ADJUSTL(rs_field_specs(4,i))), numtype, missval )

         rs_he5fields(i)%ScaleFactor    = 1.0_r8
         rs_he5fields(i)%Offset         = 0.0_r8
         rs_he5fields(i)%Name           = TRIM(ADJUSTL(rs_field_names(1,i)))
         rs_he5fields(i)%Title          = TRIM(ADJUSTL(rs_field_names(2,i)))
         rs_he5fields(i)%Units          = TRIM(ADJUSTL(rs_field_specs(1,i)))
         rs_he5fields(i)%Dimensions     = TRIM(ADJUSTL(rs_field_specs(2,i)))
         rs_he5fields(i)%UniqueFD       = TRIM(ADJUSTL(rs_field_specs(3,i)))
         rs_he5fields(i)%FillValue      = missval
         rs_he5fields(i)%MissingValue   = missval
         rs_he5fields(i)%Swath_ID       = pge_swath_id
         rs_he5fields(i)%HE5_DataType   = numtype
         rs_he5fields(i)%Rank           = 1
         rs_he5fields(i)%Dims           = (/ -1, -1, -1 /)
         rs_he5fields(i)%LenTitle       = LEN_TRIM(ADJUSTL(rs_field_names(2,i)))
         rs_he5fields(i)%LenUnits       = LEN_TRIM(ADJUSTL(rs_field_specs(1,i)))
         rs_he5fields(i)%LenUniqueFD    = LEN_TRIM(ADJUSTL(rs_field_specs(3,i)))
         rs_he5fields(i)%ValidRange     = rs_valids(1:2,i)
      END DO

      ! ------------------------------------------------------------
      ! Scattering weights, Gas profile and Averaging Kernels fields
      ! ------------------------------------------------------------
      DO i = 1, n_sw_fields
         ! ----------------------------------------
         ! Find the numeric type of the data field.
         ! ----------------------------------------
         CALL he5_find_datatype ( TRIM(ADJUSTL(sw_field_specs(4,i))), numtype, missval )

         sw_he5fields(i)%ScaleFactor    = 1.0_r8
         sw_he5fields(i)%Offset         = 0.0_r8
         sw_he5fields(i)%Name           = TRIM(ADJUSTL(sw_field_names(1,i)))
         sw_he5fields(i)%Title          = TRIM(ADJUSTL(sw_field_names(2,i)))
         sw_he5fields(i)%Units          = TRIM(ADJUSTL(sw_field_specs(1,i)))
         sw_he5fields(i)%Dimensions     = TRIM(ADJUSTL(sw_field_specs(2,i)))
         sw_he5fields(i)%UniqueFD       = TRIM(ADJUSTL(sw_field_specs(3,i)))
         sw_he5fields(i)%FillValue      = missval
         sw_he5fields(i)%MissingValue   = missval
         sw_he5fields(i)%Swath_ID       = pge_swath_id
         sw_he5fields(i)%HE5_DataType   = numtype
         sw_he5fields(i)%Rank           = 1
         sw_he5fields(i)%Dims           = (/ -1, -1, -1 /)
         sw_he5fields(i)%LenTitle       = LEN_TRIM(ADJUSTL(sw_field_names(2,i)))
         sw_he5fields(i)%LenUnits       = LEN_TRIM(ADJUSTL(sw_field_specs(1,i)))
         sw_he5fields(i)%LenUniqueFD    = LEN_TRIM(ADJUSTL(sw_field_specs(3,i)))
         sw_he5fields(i)%ValidRange     = sw_valids(1:2,i)
      END DO

      RETURN
    END SUBROUTINE he5_initialize_datafields

    SUBROUTINE he5_find_datatype ( tyc, numtype, missval )

      USE OMSAO_parameters_module, ONLY: &
           r8_missval, r4_missval, i4_missval, i2_missval, i1_missval, str_missval
      IMPLICIT NONE

      CHARACTER (LEN=2), INTENT (IN) :: tyc

      INTEGER (KIND=i4), INTENT (OUT) :: numtype
      REAL    (KIND=r8), INTENT (OUT) :: missval
      
      numtype = 0  ;  missval = 0.0_r8

      SELECT CASE (tyc)
      CASE ("ch")
         numtype = HE5T_NATIVE_CHAR
         missval = REAL ( i1_missval, KIND=r8 )
      CASE ("i1")
         numtype = HE5T_NATIVE_INT8
         missval = REAL ( i1_missval, KIND=r8 )
      CASE ("i2")
         numtype = HE5T_NATIVE_INT16
         missval = REAL ( i2_missval, KIND=r8 )
      CASE ("i4")
         numtype = HE5T_NATIVE_INT
         missval = REAL ( i4_missval, KIND=r8 )
      CASE ("r4")
         numtype = HE5T_NATIVE_FLOAT
         missval = REAL ( r4_missval, KIND=r8 )
      CASE ("r8")
         numtype = HE5T_NATIVE_DOUBLE
         missval = r8_missval
      END SELECT

      RETURN
    END SUBROUTINE he5_find_datatype

END MODULE OMSAO_he5_datafields_module
