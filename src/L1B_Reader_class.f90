!! 
! MODULE L1B_reader_class contains the OMI L1B data structure and
! the functions for the creation, deletion and access of specific 
! information stored in the data structure.  A user can use the
! following functions:  L1Br_open, L1Br_close, L1Br_getDIMsize, 
! L1Br_getSWdims, L1Br_getGEOpix, L1Br_getGEOline, L1Br_getDATA,
! L1Br_getSIGpix, 1Br_getSIGline, L1Br_getSIGpixWL, L1Br_getSIGlineWL,
! and L1Br_getDATAFIELDline.
!!

MODULE L1B_Reader_class

   IMPLICIT NONE
   INTEGER (KIND = 4) :: OMI_SMF_setmsg !!this external function writes the PGS
   EXTERNAL              OMI_SMF_setmsg !!SMF error messages to the Log files
   PUBLIC  :: L1Br_open
   PUBLIC  :: L1Br_close
   PUBLIC  :: L1Br_getDIMsize
   PUBLIC  :: L1Br_getSWdims
   PUBLIC  :: L1Br_getGEOpix
   PUBLIC  :: L1Br_getGEOline
   PUBLIC  :: L1Br_getDATA
   PUBLIC  :: L1Br_getSIGpix
   PUBLIC  :: L1Br_getSIGline
   PUBLIC  :: L1Br_getSIGpixWL
   PUBLIC  :: L1Br_getSIGlineWL
   PUBLIC  :: L1Br_getDATAFIELDline
   PUBLIC  :: L1Br_getDATAline

   PRIVATE :: fill_l1b_blk
   PRIVATE :: init_l1b_blk
   PRIVATE :: calc_wl_pix
   PRIVATE :: calc_wl_line
   PRIVATE :: check_blk_line
   PRIVATE :: check_blk_pix
   INTEGER (KIND = 4), PARAMETER, PUBLIC :: MAX_NAME_LENGTH = 257
   INTEGER, PRIVATE :: ierr  !!error code returned from a function
   INTEGER, PARAMETER, PRIVATE :: zero = 0
   INTEGER, PARAMETER, PRIVATE :: one  = 1
   INTEGER, PARAMETER, PRIVATE :: two  = 2
   REAL (KIND = 4), PARAMETER, PRIVATE :: rfill  = -1.0*(2.0**100)
   REAL (KIND = 8), PARAMETER, PRIVATE :: r8fill  = -1.0*(2.0**100)
   INTEGER (KIND = 1), PARAMETER, PRIVATE :: int8fill  = -127
   INTEGER (KIND = 2), PARAMETER, PRIVATE :: int16fill  = -32767
   INTEGER (KIND = 4), PARAMETER, PRIVATE :: int32fill  = -2147483647

 !!
 ! L1B_block_type is designed to store a block of L1B swath data.
 ! A block contains all the geolocation, radiance and wavelength
 ! data fields for a number of "scan" lines or "exposure times"
 ! in a L1B swath. 
 !
 ! The fields in the L1B_block_type are defined as follows:
 !
 ! initialized: the logical variable to indicate whether the
 !              block has been initialized.  The block data structure
 !              has to be created (initialized) before it can be used.
 !              
 ! filename: the file name for the L1B file (including path)
 ! swathname: the swath name in the L1B file, either "Earth UV-1 Swath",
 !            "Earth UV-2 Swath", or "Earth VIS Swath"
 ! iLine: the starting line for the block.  It is a reference to the
 !        L1B file, where 0 is the first line, and (NumTimes - 1)
 !        is the last line.
 ! eLine: the ending line for the block
 ! nLine: the number of lines containing in the block with iLine
 !        as start and eLine as the end (inclusive)
 !
 ! NumTimes: the total number of lines contained in the L1B swath
 ! nXtrack: the number of cross-track pixels in the swath
 ! nWavel: the number of wavelengths in the swath
 ! nWavelCoef: the number of wavelength coefficients 
 ! NumTimesSmallPixel: the total number of small pixel columns contained in the L1B swath
 !
 ! The remaining elements in the L1B_block_type are directly associated with datasets in
 ! the L1B file.  They are listed here, associated with the L1B Dataset Name whose data
 ! they contain.
 !
 ! Element Name    L1B Dataset Name
 ! ============    ================
 ! GEOLOCATION FIELDS: (extractable with the getGEO functions)
 ! Time            Time
 ! SecInDay        SecondsInDay
 ! ScLat           SpacecraftLatitude
 ! ScLon           SpacecraftLongitude
 ! ScAlt           SpacecraftAltitude
 ! SolElevation    SolarElevation
 ! SolElevMin      SolarElevationMinimum
 ! SolElevMax      SolarElevationMaximum
 ! SolAzimuth      SolarAzimuth
 ! SolAziMin       SolarAzimuthMinimum
 ! SolAziMax       SolarAzimuthMaximum
 ! Lon             Longitude
 ! Lat             Latitude
 ! SolZenAng       SolarZenithAngle
 ! SolAziAng       SolarAzimuthAngle
 ! ViewZenAng      ViewingZenithAngle
 ! ViewAziAng      ViewingAzimuthAngle
 ! TerrainHeight   TerrainHeight
 ! GPQFlag         GroundPixelQualityFlags
 !
 ! RADIANCE RELATED FIELDS: (extractable with the getSIG functions)
 ! RadMantissa     RadianceMantissa or IrradianceMantissa or SignalMantissa
 ! RadPrecision    RadiancePrecision or IrradiancePrecision or SignalPrecision
 ! PQFlag          PixelQualityFlags
 ! RadExponent     RadianceExponent or IrradianceExponent or SignalExponent
 ! WavelenCoef     WavelengthCoefficient
 ! WavelenPrec     WavelengthCoefficientPrecision
 ! WavelenRefCol   WavelengthReferenceColumn
 ! SmPixRad        SmallPixelRadiance or SmallPixelIrradiance or SmallPixelSignal
 ! SmPixWavelen    SmallPixelWavelength
 ! NumSmPixCol     NumberSmallPixelColumns
 ! NumSmPix        NumberSmallPixelColumns for calculation
 ! SmPixCol        SmallPixelColumn
 !
 ! OTHER DATA FIELDS: (extractable with the getDATA function)
 ! Config          InstrumentConfigurationId
 ! MeasClass       MeasurementClass
 ! ExposureType    ExposureType
 ! ImBinFact       ImageBinningFactor
 ! GC1             GainCode1
 ! GC2             GainCode2
 ! GC3             GainCode3
 ! GC4             GainCode4
 ! DSGC            DSGainCode
 ! LSLABF          LowerStrayLightAreaBinningFactor
 ! USLABF          UpperStrayLightAreaBinningFactor
 ! LDABF           LowerDarkAreaBinningFactor
 ! UDABF           UpperDarkAreaBinningFactor
 ! MQFlag          MeasurementQualityFlags
 ! CalSet          CalibrationSettings
 ! GSC1            GainSwitchingColumn1
 ! GSC2            GainSwitchingColumn2
 ! GSC3            GainSwitchingColumn3
 ! SR1             SkipRows1
 ! SR2             SkipRows2
 ! SR3             SkipRows3
 ! SR4             SkipRows4
 ! BinImgRows      BinnedImageRows
 ! StopColumn      StopColumn
 ! MasterClkPer    MasterClockPeriod
 ! ExpTime         ExposureTime
 ! ReadTime        ReadoutTime
 ! DetcTemp        DetectorTemperature
 ! OptBenchTemp    OpticalBenchTemperature
 !
 ! The Geolocation fields can be accessed via the L1Br_getGEOpix and L1Br_getGEOline
 ! functions.  These functions simply return the requested values from the L1B file and swath
 ! for the input pixel (L1Br_getGEOpix) or line (L1Br_getGEOline).
 !
 ! The Radiance related fields can be accessed via the L1Br_getSIGpix, L1Br_getSIGline, 
 ! L1Br_getSIGpixWL, and L1Br_getSIGlineWL functions.  These functions will recreate the
 ! signal (Radiance, Irradiance, or Signal) and the wavelengths from the raw data.
 ! L1Br_getSIGpix and L1Br_getSIGline will also return the raw data for a specific
 ! pixel (L1Br_getSIGpix) or line (L1Br_getSIGline).  L1Br_getSIGpixWL and L1Br_getSIGlineWL
 ! will not return raw data, but will interpolate to give approximate values of the signal
 ! for specific input wavelengths, whether or not those wavelengths are exactly present
 ! in the L1B file (as long as the requested wavelengths fall within the overall spread
 ! of wavelengths in the pixel/line in question).
 !
 ! The Other Data fields can be accessed via the L1Br_getDATA function.  This function simply
 ! returns the requested values from the L1B file and swath for the input line (all values
 ! accessed by this function have only a single value per line in the L1B file).
 !
 ! The function L1Br_getDATAFIELDline can be used to access any datasets not covered
 ! by other functions and not contained in the L1B_block_type data block.
 !!

   TYPE, PUBLIC :: L1B_block_type
      PRIVATE
      CHARACTER (LEN = MAX_NAME_LENGTH) :: filename, swathname
      CHARACTER (LEN = MAX_NAME_LENGTH) :: fieldname
      INTEGER (KIND = 4) :: iLine
      INTEGER (KIND = 4) :: eLine
      INTEGER (KIND = 4) :: nLine
      INTEGER (KIND = 4) :: NumTimes
      INTEGER (KIND = 4) :: NumTimesSmallPixel
      INTEGER (KIND = 4) :: nXtrack
      INTEGER (KIND = 4) :: nWavel
      INTEGER (KIND = 4) :: nWavelCoef
      LOGICAL :: initialized = .FALSE. ! JED fix
      INTEGER (KIND = 4) :: SpaZoomIndex

      ! Data Type                                    Element Name    L1B Dataset Name
      ! ====================================         ============    ================
      ! Geolocation fields
      REAL (KIND = 8),    DIMENSION(:),   POINTER :: Time => NULL()         ! Time
      REAL (KIND = 4),    DIMENSION(:),   POINTER :: SecInDay => NULL()     ! SecondsInDay
      REAL (KIND = 4),    DIMENSION(:),   POINTER :: ScLat => NULL()        ! SpacecraftLatitude
      REAL (KIND = 4),    DIMENSION(:),   POINTER :: ScLon => NULL()        ! SpacecraftLongitude
      REAL (KIND = 4),    DIMENSION(:),   POINTER :: ScAlt => NULL()        ! SpacecraftAltitude
      REAL (KIND = 4),    DIMENSION(:),   POINTER :: SolElevation => NULL() ! SolarElevation
      REAL (KIND = 4),    DIMENSION(:),   POINTER :: SolAzimuth => NULL()   ! SolarAzimuth
      REAL (KIND = 4),    DIMENSION(:),   POINTER :: SolElevMin => NULL()   ! SolarElevationMinimum
      REAL (KIND = 4),    DIMENSION(:),   POINTER :: SolElevMax => NULL()   ! SolarElevationMaximum
      REAL (KIND = 4),    DIMENSION(:),   POINTER :: SolAziMin => NULL()    ! SolarAzimuthMinimum
      REAL (KIND = 4),    DIMENSION(:),   POINTER :: SolAziMax => NULL()    ! SolarAzimuthMaximum
      REAL (KIND = 4),    DIMENSION(:,:), POINTER :: Lon => NULL()          ! Longitude
      REAL (KIND = 4),    DIMENSION(:,:), POINTER :: Lat => NULL()          ! Latitude
      REAL (KIND = 4),    DIMENSION(:,:), POINTER :: SolZenAng => NULL()    ! SolarZenithAngle
      REAL (KIND = 4),    DIMENSION(:,:), POINTER :: SolAziAng => NULL()    ! SolarAzimuthAngle
      REAL (KIND = 4),    DIMENSION(:,:), POINTER :: ViewZenAng => NULL()   ! ViewingZenithAngle
      REAL (KIND = 4),    DIMENSION(:,:), POINTER :: ViewAziAng => NULL()   ! ViewingAzimuthAngle
      INTEGER (KIND = 2), DIMENSION(:,:), POINTER :: TerrainHeight => NULL()! TerrainHeight
      INTEGER (KIND = 2), DIMENSION(:,:), POINTER :: GPQFlag => NULL()      ! GroundPixelQualityFlags
      INTEGER (KIND = 1), DIMENSION(:,:), POINTER :: XTQFlag => NULL()      ! XTrackQualityFlags

      ! Data fields
      INTEGER (KIND = 2), DIMENSION(:,:,:), POINTER :: RadMantissa => NULL() ! RadianceMantissa or
                                                                             ! IrradianceMantissa or
                                                                             ! SignalMantissa
      INTEGER (KIND = 2), DIMENSION(:,:,:), POINTER :: RadPrecision => NULL()! RadiancePrecision or
                                                                             ! IrradiancePrecision or
                                                                             ! SignalPrecision
      INTEGER (KIND = 1), DIMENSION(:,:,:), POINTER :: RadExponent => NULL() ! RadianceExponent or
                                                                             ! IrradianceExponent or
                                                                             ! SignalExponent
      INTEGER (KIND = 2), DIMENSION(:,:,:), POINTER :: PQFlag => NULL()      ! PixelQualityFlags
      INTEGER (KIND = 2), DIMENSION(:,:),   POINTER :: PFlag => NULL()       ! PixelQualityFlags
      REAL (KIND = 4),    DIMENSION(:,:,:), POINTER :: WavelenCoef => NULL() ! WavelengthCoefficient
      REAL (KIND = 4),    DIMENSION(:,:,:), POINTER :: WavelenPrec => NULL() ! WavelengthCoefficientPrecision
      INTEGER (KIND = 2), DIMENSION(:), POINTER :: WavelenRefCol => NULL()   ! WavelengthReferenceColumn
      REAL (KIND = 4),  DIMENSION(:,:), POINTER :: SmPixRad => NULL()        ! SmallPixelRadiance or
                                                                             ! SmallPixelIrradiance or
                                                                             ! SmallPixelSignal
      REAL (KIND = 4),  DIMENSION(:,:), POINTER :: SmPixWavelen => NULL()! SmallPixelWavelength

      INTEGER (KIND = 1), DIMENSION(:), POINTER :: MeasClass => NULL()   ! MeasurementClass
      INTEGER (KIND = 1), DIMENSION(:), POINTER :: Config => NULL()      ! InstrumentConfigurationId
      INTEGER (KIND = 2), DIMENSION(:), POINTER :: MQFlag => NULL()      ! MeasurementQualityFlags
      INTEGER (KIND = 1), DIMENSION(:), POINTER :: NumSmPixCol => NULL() ! NumberSmallPixelColumns
      INTEGER (KIND = 1), DIMENSION(:), POINTER :: NumSmPix => NULL()    ! NumberSmallPixelColumns_cal
      INTEGER (KIND = 1), DIMENSION(:), POINTER :: ExposureType => NULL()! ExposureType
      REAL (KIND = 4),    DIMENSION(:), POINTER :: MasterClkPer => NULL()! MasterClockPeriod
      INTEGER (KIND = 2), DIMENSION(:), POINTER :: CalSet => NULL()      ! CalibrationSettings
      REAL (KIND = 4),    DIMENSION(:), POINTER :: ExpTime => NULL()     ! ExposureTime
      REAL (KIND = 4),    DIMENSION(:), POINTER :: ReadTime => NULL()    ! ReadoutTime
      INTEGER (KIND = 2), DIMENSION(:), POINTER :: SmPixCol => NULL()    ! SmallPixelColumn
      INTEGER (KIND = 2), DIMENSION(:), POINTER :: GSC1 => NULL()        ! GainSwitchingColumn1
      INTEGER (KIND = 2), DIMENSION(:), POINTER :: GSC2 => NULL()        ! GainSwitchingColumn2
      INTEGER (KIND = 2), DIMENSION(:), POINTER :: GSC3 => NULL()        ! GainSwitchingColumn3
      INTEGER (KIND = 1), DIMENSION(:), POINTER :: GC1 => NULL()         ! GainCode1
      INTEGER (KIND = 1), DIMENSION(:), POINTER :: GC2 => NULL()         ! GainCode2
      INTEGER (KIND = 1), DIMENSION(:), POINTER :: GC3 => NULL()         ! GainCode3
      INTEGER (KIND = 1), DIMENSION(:), POINTER :: GC4 => NULL()         ! GainCode4
      INTEGER (KIND = 1), DIMENSION(:), POINTER :: DSGC => NULL()        ! DSGainCode
      INTEGER (KIND = 1), DIMENSION(:), POINTER :: LSLABF => NULL()      ! LowerStrayLightAreaBinningFactor
      INTEGER (KIND = 1), DIMENSION(:), POINTER :: USLABF => NULL()      ! UpperStrayLightAreaBinningFactor
      INTEGER (KIND = 1), DIMENSION(:), POINTER :: LDABF => NULL()       ! LowerDarkAreaBinningFactor
      INTEGER (KIND = 1), DIMENSION(:), POINTER :: UDABF => NULL()       ! UpperDarkAreaBinningFactor
      INTEGER (KIND = 2), DIMENSION(:), POINTER :: SR1 => NULL()         ! SkipRows1
      INTEGER (KIND = 2), DIMENSION(:), POINTER :: SR2 => NULL()         ! SkipRows2
      INTEGER (KIND = 2), DIMENSION(:), POINTER :: SR3 => NULL()         ! SkipRows3
      INTEGER (KIND = 2), DIMENSION(:), POINTER :: SR4 => NULL()         ! SkipRows4
      REAL (KIND = 4),    DIMENSION(:), POINTER :: DetcTemp => NULL()    ! DetectorTemperature
      REAL (KIND = 4),    DIMENSION(:), POINTER :: OptBenchTemp => NULL()! OpticalBenchTemperature
      INTEGER (KIND = 1), DIMENSION(:), POINTER :: ImBinFact => NULL()   ! ImageBinningFactor
      INTEGER (KIND = 2), DIMENSION(:), POINTER :: BinImgRows => NULL()  ! BinnedImageRows
      INTEGER (KIND = 2), DIMENSION(:), POINTER :: StopColumn => NULL()  ! StopColumn

   END TYPE L1B_block_type

   CONTAINS

!! 1. L1Br_open
 !
 !    Functionality:
 !
 !       This function should be called first to initiate the
 !       the interface with the L1B swath.  This function allocated memory
 !       in the data block, based on metadata located within the L1B file
 !       and swath to be accessed.  It also sets parameter values within
 !       the data block with values read from the L1B file and swath.
 !
 !    Calling Arguments:
 !
 !    Inputs:
 !
 !       this: the block data structure
 !       fn  : the L1B file name
 !       swn : the swath name in the L1B file
 !       nL  : the number of lines the data block can store.  This is
 !             an optional input.  If it is not present, then a default
 !             value of 100 is used.
 !
 !    Outputs:
 !
 !       status: the return PGS_SMF status value
 !
 !    Change History:
 !
 !       Date            Author          Modifications
 !       ====            ======          =============
 !       January 2005    Jeremy Warner   Original Source (derived from original 
 !                                             L1B Reader code written by Kai Yang)
 !
!!
      FUNCTION L1Br_open(this, fn, swn, nL) RESULT (status)
        USE hdfeos4_parameters
        INCLUDE 'PGS_OMI_1900.f'  !defines the L1B PGS error codes
        TYPE (L1B_block_type), INTENT(INOUT) :: this
        CHARACTER (LEN = *), INTENT(IN) :: fn, swn
        INTEGER (KIND = 4), OPTIONAL, INTENT(IN) :: nL
        INTEGER (KIND = 4) :: swfid, swid, status
        INTEGER (KIND = 4) :: dum
        CHARACTER (LEN = 256) :: message
        INTEGER (KIND = 4) :: getl1bblk
        EXTERNAL              getl1bblk
        INTEGER (KIND = 4) :: rank
        INTEGER (KIND = 4) :: fldflg
        INTEGER (KIND = 4), DIMENSION(1:3) :: dims
        INTEGER (KIND = 4) ::  IposL1bIdf
        CHARACTER (LEN = 6) :: L1bIdf
        CHARACTER (LEN = 1) :: ZoomIdf
        
        ! open the L1B swath file
        swfid = swopen(fn, DFACC_READ)
        IF(swfid < zero) THEN
           status = OMI_E_FAILURE
           ierr = OMI_SMF_setmsg(OMI_E_FILE_OPEN, fn, "L1Br_open", zero)
           RETURN
        ENDIF 

        ! attach to the swath
        swid = swattach(swfid, swn)
        IF(swid < zero) THEN
           status = OMI_E_FAILURE
           ierr = OMI_SMF_setmsg(OMI_E_SWATH_ATTACH, swn, "L1Br_open", zero)
           RETURN
        ENDIF 

        ! intitialize the start and end line to invalid values
        this%iLine      = -1
        this%eLine      = -1

        ! retrieve the dimension info from the swath file.
        ! dimension names are obtained from file spec
        status = swrdattr(swid, "NumTimes", this%NumTimes)
        If(status < zero) THEN
           status = OMI_E_FAILURE
           ierr = OMI_SMF_setmsg(OMI_E_HDFEOS, "get NumTimes size failed", &
                                 "L1Br_open", zero)
           RETURN
        ENDIF 

        status = swrdattr(swid, "NumTimesSmallPixel", this%NumTimesSmallPixel)
        If(status < zero) THEN
           status = OMI_E_FAILURE
           ierr = OMI_SMF_setmsg(OMI_E_HDFEOS, "get NumTimesSmallPixel size failed", &
                                 "L1Br_open", zero)
           RETURN
        ENDIF 

        this%nXtrack    = swdiminfo(swid, "nXtrack")
        If(this%nXtrack < zero) THEN
           status = OMI_E_FAILURE
           ierr = OMI_SMF_setmsg(OMI_E_HDFEOS, "get xTrack size failed", &
                                  "L1Br_open", zero)
           RETURN
        ENDIF 

        this%nWavel     = swdiminfo(swid, "nWavel")
        If(this%nWavel < zero) THEN
           status = OMI_E_FAILURE
           ierr = OMI_SMF_setmsg(OMI_E_HDFEOS, "get nWavel size failed", &
                                  "L1Br_open", zero)
           RETURN
        ENDIF 

        this%nWavelCoef = swdiminfo(swid, "nWavelCoef")
        If(this%nWavelCoef < zero) THEN
           status = OMI_E_FAILURE
           ierr = OMI_SMF_setmsg(OMI_E_HDFEOS, "get nWavelCoef size failed", &
                                  "L1Br_open", zero)
           RETURN
        ENDIF 

        ! detach and close the L1B swath files.  No need for 
        ! error checking, for an error is unlikely to occur here.
        ierr = swdetach(swid)
        ierr = swclose(swfid)
          
        ! set the block size according to input, if present
        ! if not set it to 100
        IF(.not. PRESENT(nL)) THEN
           this%nLine = 100
           IF (this%nLine > this%NumTimes) THEN
              this%nLine = this%NumTimes
           ENDIF
        ELSE
           IF(nL <= 0) THEN
              this%nLine = 100
              IF (this%nLine > this%NumTimes) THEN
                 this%nLine = this%NumTimes
              ENDIF
           ELSE 
              this%nLine = nL
           ENDIF
        ENDIF
        
        ! the maximum number of lines cannot exceed the total
        ! number of lines in the L1B swath.
        IF(this%nLine > this%NumTimes) THEN
           status = OMI_E_FAILURE
           ierr = OMI_SMF_setmsg(OMI_E_INPUT, "input nL too large ", & 
                                 "L1Br_open", zero)
           RETURN
        ENDIF 
        
        this%filename = fn   
        this%swathname = swn

        L1bIdf='OML1BR'
        IposL1bIdf=index(fn,L1bIdf)
        ZoomIdf=fn(IposL1bIdf+7:IposL1bIdf+7)
        this%SpaZoomIndex = 0
        IF (ZoomIdf == 'Z') THEN
          this%SpaZoomIndex = 1
        ENDIF
        
        ! Allocate the memory for storage of geolocation and
        ! radiance fields.  First make sure pointers are not
        ! associated with any specific memory location.  If
        ! they are, deallocate the memory, then allocate the
        ! the proper amount of memory for each data field,
        ! error checking for each memory allocation to make
        ! sure memory is allocated successfully. 

        IF (ASSOCIATED(this%Time)) DEALLOCATE(this%Time)
        IF (ASSOCIATED(this%SecInDay)) DEALLOCATE(this%SecInDay)
        IF (ASSOCIATED(this%ScLat)) DEALLOCATE(this%ScLat)
        IF (ASSOCIATED(this%ScLon)) DEALLOCATE(this%ScLon)
        IF (ASSOCIATED(this%ScAlt)) DEALLOCATE(this%ScAlt)
        IF (ASSOCIATED(this%SolElevation)) DEALLOCATE(this%SolElevation)
        IF (ASSOCIATED(this%SolElevMin)) DEALLOCATE(this%SolElevMin)
        IF (ASSOCIATED(this%SolElevMax)) DEALLOCATE(this%SolElevMax)
        IF (ASSOCIATED(this%SolAzimuth)) DEALLOCATE(this%SolAzimuth)
        IF (ASSOCIATED(this%SolAziMin)) DEALLOCATE(this%SolAziMin)
        IF (ASSOCIATED(this%SolAziMax)) DEALLOCATE(this%SolAziMax)
        IF (ASSOCIATED(this%Lon)) DEALLOCATE(this%Lon)
        IF (ASSOCIATED(this%Lat)) DEALLOCATE(this%Lat)
        IF (ASSOCIATED(this%SolZenAng)) DEALLOCATE(this%SolZenAng)
        IF (ASSOCIATED(this%SolAziAng)) DEALLOCATE(this%SolAziAng)
        IF (ASSOCIATED(this%ViewZenAng)) DEALLOCATE(this%ViewZenAng)
        IF (ASSOCIATED(this%ViewAziAng)) DEALLOCATE(this%ViewAziAng)
        IF (ASSOCIATED(this%TerrainHeight)) DEALLOCATE(this%TerrainHeight)
        IF (ASSOCIATED(this%GPQFlag)) DEALLOCATE(this%GPQFlag)
        IF (ASSOCIATED(this%XTQFlag)) DEALLOCATE(this%XTQFlag)
        IF (ASSOCIATED(this%RadMantissa)) DEALLOCATE(this%RadMantissa)
        IF (ASSOCIATED(this%RadPrecision)) DEALLOCATE(this%RadPrecision)
        IF (ASSOCIATED(this%PQFlag)) DEALLOCATE(this%PQFlag)
        IF (ASSOCIATED(this%PFlag)) DEALLOCATE(this%PFlag)
        IF (ASSOCIATED(this%RadExponent)) DEALLOCATE(this%RadExponent)
        IF (ASSOCIATED(this%SmPixRad)) DEALLOCATE(this%SmPixRad)
        IF (ASSOCIATED(this%SmPixWavelen)) DEALLOCATE(this%SmPixWavelen)
        IF (ASSOCIATED(this%WavelenCoef)) DEALLOCATE(this%WavelenCoef)
        IF (ASSOCIATED(this%WavelenPrec)) DEALLOCATE(this%WavelenPrec)
        IF (ASSOCIATED(this%WavelenRefCol)) DEALLOCATE(this%WavelenRefCol)
        IF (ASSOCIATED(this%Config)) DEALLOCATE(this%Config)
        IF (ASSOCIATED(this%MeasClass)) DEALLOCATE(this%MeasClass)
        IF (ASSOCIATED(this%NumSmPixCol)) DEALLOCATE(this%NumSmPixCol)
        IF (ASSOCIATED(this%NumSmPix)) DEALLOCATE(this%NumSmPix)
        IF (ASSOCIATED(this%ExposureType)) DEALLOCATE(this%ExposureType)
        IF (ASSOCIATED(this%ImBinFact)) DEALLOCATE(this%ImBinFact)
        IF (ASSOCIATED(this%GC1)) DEALLOCATE(this%GC1)
        IF (ASSOCIATED(this%GC2)) DEALLOCATE(this%GC2)
        IF (ASSOCIATED(this%GC3)) DEALLOCATE(this%GC3)
        IF (ASSOCIATED(this%GC4)) DEALLOCATE(this%GC4)
        IF (ASSOCIATED(this%DSGC)) DEALLOCATE(this%DSGC)
        IF (ASSOCIATED(this%LSLABF)) DEALLOCATE(this%LSLABF)
        IF (ASSOCIATED(this%USLABF)) DEALLOCATE(this%USLABF)
        IF (ASSOCIATED(this%LDABF)) DEALLOCATE(this%LDABF)
        IF (ASSOCIATED(this%UDABF)) DEALLOCATE(this%UDABF)
        IF (ASSOCIATED(this%MQFlag)) DEALLOCATE(this%MQFlag)
        IF (ASSOCIATED(this%CalSet)) DEALLOCATE(this%CalSet)
        IF (ASSOCIATED(this%SmPixCol)) DEALLOCATE(this%SmPixCol)
        IF (ASSOCIATED(this%GSC1)) DEALLOCATE(this%GSC1)
        IF (ASSOCIATED(this%GSC2)) DEALLOCATE(this%GSC2)
        IF (ASSOCIATED(this%GSC3)) DEALLOCATE(this%GSC3)
        IF (ASSOCIATED(this%SR1)) DEALLOCATE(this%SR1)
        IF (ASSOCIATED(this%SR2)) DEALLOCATE(this%SR2)
        IF (ASSOCIATED(this%SR3)) DEALLOCATE(this%SR3)
        IF (ASSOCIATED(this%SR4)) DEALLOCATE(this%SR4)
        IF (ASSOCIATED(this%BinImgRows)) DEALLOCATE(this%BinImgRows)
        IF (ASSOCIATED(this%StopColumn)) DEALLOCATE(this%StopColumn)
        IF (ASSOCIATED(this%MasterClkPer)) DEALLOCATE(this%MasterClkPer)
        IF (ASSOCIATED(this%ExpTime)) DEALLOCATE(this%ExpTime)
        IF (ASSOCIATED(this%ReadTime)) DEALLOCATE(this%ReadTime)
        IF (ASSOCIATED(this%DetcTemp)) DEALLOCATE(this%DetcTemp)
        IF (ASSOCIATED(this%OptBenchTemp)) DEALLOCATE(this%OptBenchTemp)

        ! What follows is the tedious memory allocation

        ALLOCATE (this%Time(this%nLine), STAT = ierr)
           IF(ierr .NE. 0) THEN
              message = "Time"
              GOTO 999
           ENDIF
        ALLOCATE (this%SecInDay(this%nLine), STAT = ierr)
           IF(ierr .NE. 0) THEN
              message = "SecInDay"
              GOTO 999
           ENDIF
        ALLOCATE (this%ScLat(this%nLine), STAT = ierr)
           IF(ierr .NE. 0) THEN
              message = "ScLat"
              GOTO 999
           ENDIF
        ALLOCATE (this%ScLon(this%nLine), STAT = ierr)
           IF(ierr .NE. 0) THEN
              message = "ScLon"
              GOTO 999
           ENDIF
        ALLOCATE (this%ScAlt(this%nLine), STAT = ierr)
           IF(ierr .NE. 0) THEN
              message = "ScAlt"
              GOTO 999
           ENDIF
        ALLOCATE (this%SolElevation(this%nLine), STAT = ierr)
           IF(ierr .NE. 0) THEN
              message = "SolElevation"
              GOTO 999
           ENDIF
        ALLOCATE (this%SolAzimuth(this%nLine), STAT = ierr)
           IF(ierr .NE. 0) THEN
              message = "SolAzimuth"
              GOTO 999
           ENDIF
        ALLOCATE (this%SolElevMin(this%nLine), STAT = ierr)
           IF(ierr .NE. 0) THEN
              message = "SolElevMin"
              GOTO 999
           ENDIF
        ALLOCATE (this%SolElevMax(this%nLine), STAT = ierr)
           IF(ierr .NE. 0) THEN
              message = "SolElevMax"
              GOTO 999
           ENDIF
        ALLOCATE (this%SolAziMin(this%nLine), STAT = ierr)
           IF(ierr .NE. 0) THEN
              message = "SolAziMin"
              GOTO 999
           ENDIF
        ALLOCATE (this%SolAziMax(this%nLine), STAT = ierr)
           IF(ierr .NE. 0) THEN
              message = "SolAziMax"
              GOTO 999
           ENDIF
        ALLOCATE (this%Lat(this%nXtrack,this%nLine), STAT = ierr)
           IF(ierr .NE. 0) THEN
              message = "Lat"
              GOTO 999
           ENDIF
        ALLOCATE (this%Lon(this%nXtrack,this%nLine), STAT = ierr)
           IF(ierr .NE. 0) THEN
              message = "Lon"
              GOTO 999
           ENDIF
        ALLOCATE (this%SolZenAng(this%nXtrack,this%nLine), STAT = ierr)
           IF(ierr .NE. 0) THEN
              message = "SolZenAng"
              GOTO 999
           ENDIF
        ALLOCATE (this%SolAziAng(this%nXtrack,this%nLine), STAT = ierr)
           IF(ierr .NE. 0) THEN
              message = "SolAziAng"
              GOTO 999
           ENDIF
        ALLOCATE (this%ViewZenAng(this%nXtrack,this%nLine), STAT = ierr)
           IF(ierr .NE. 0) THEN
              message = "ViewZenAng"
              GOTO 999
           ENDIF
        ALLOCATE (this%ViewAziAng(this%nXtrack,this%nLine), STAT = ierr)
           IF(ierr .NE. 0) THEN
              message = "ViewAziAng"
              GOTO 999
           ENDIF
        ALLOCATE (this%TerrainHeight(this%nXtrack,this%nLine), STAT = ierr)
           IF(ierr .NE. 0) THEN
              message = "TerrainHeight"
              GOTO 999
           ENDIF
        ALLOCATE (this%GPQFlag(this%nXtrack,this%nLine), STAT = ierr)
           IF(ierr .NE. 0) THEN
              message = "GPQFlag"
              GOTO 999
           ENDIF
        ALLOCATE (this%XTQFlag(this%nXtrack,this%nLine), STAT = ierr)
           IF(ierr .NE. 0) THEN
             message = "XTQFlag"
             GOTO 999
           ENDIF
        ALLOCATE (this%RadMantissa(this%nWavel,this%nXtrack,this%nLine), STAT = ierr)
           IF(ierr .NE. 0) THEN
              message = "RadMantissa"
              GOTO 999
           ENDIF
        ALLOCATE (this%RadPrecision(this%nWavel,this%nXtrack,this%nLine), STAT = ierr)
           IF(ierr .NE. 0) THEN
              message = "RadPrecision"
              GOTO 999
           ENDIF
        ALLOCATE (this%RadExponent(this%nWavel,this%nXtrack,this%nLine), STAT = ierr)
           IF(ierr .NE. 0) THEN
              message = "RadExponent"
              GOTO 999
           ENDIF
        ALLOCATE (this%PQFlag(this%nWavel,this%nXtrack,this%nLine), STAT = ierr)
           IF(ierr .NE. 0) THEN
              message = "PQFlag"
              GOTO 999
           ENDIF

        ALLOCATE (this%PFlag(this%nXtrack,this%nLine), STAT = ierr)
           IF(ierr .NE. 0) THEN
              message = "PFlag"
              GOTO 999
           ENDIF

        ALLOCATE (this%WavelenCoef(this%nWavelCoef,this%nXtrack,this%nLine), &
              STAT = ierr)
           IF(ierr .NE. 0) THEN
              message = "WavelenCoef"
              GOTO 999
           ENDIF
        ALLOCATE (this%WavelenPrec(this%nWavelCoef,this%nXtrack,this%nLine), &
              STAT = ierr)
           IF(ierr .NE. 0) THEN
              message = "WavelenPrec"
              GOTO 999
           ENDIF
        ALLOCATE (this%WavelenRefCol(this%nLine), STAT = ierr)
           IF(ierr .NE. 0) THEN
              message = "WavelenRefCol"
              GOTO 999
           ENDIF
        ALLOCATE(this%SmPixRad(this%nXtrack,this%NumTimesSmallPixel), STAT = ierr)
           IF(ierr .NE. 0) THEN
              message = "SmPixRad"
              GOTO 999
           ENDIF
        ALLOCATE(this%SmPixWavelen(this%nXtrack,this%NumTimesSmallPixel), STAT = ierr)
           IF(ierr .NE. 0) THEN
              message = "SmPixWavelen"
              GOTO 999
           ENDIF
        ALLOCATE (this%MeasClass(this%nLine), STAT = ierr)
           IF(ierr .NE. 0) THEN
              message = "MeasClass"
              GOTO 999
           ENDIF
        ALLOCATE (this%Config(this%nLine), STAT = ierr)
           IF(ierr .NE. 0) THEN
              message = "Config"
              GOTO 999
           ENDIF
        ALLOCATE (this%MQFlag(this%nLine), STAT = ierr)
           IF(ierr .NE. 0) THEN
              message = "MQFlag"
              GOTO 999
           ENDIF
        ALLOCATE (this%NumSmPixCol(this%NumTimes), STAT = ierr) !tpk  ALLOCATE (this%NumSmPixCol(this%nLine), STAT = ierr)

           IF(ierr .NE. 0) THEN
              message = "NumSmPixCol"
              GOTO 999
           ENDIF
        ALLOCATE (this%NumSmPix(this%NumTimes), STAT = ierr)
           IF(ierr .NE. 0) THEN
              message = "NumSmPix"
              GOTO 999
           ENDIF

        ALLOCATE (this%ExposureType(this%nLine), STAT = ierr)
           IF(ierr .NE. 0) THEN
              message = "ExposureType"
              GOTO 999
           ENDIF
        ALLOCATE (this%MasterClkPer(this%nLine), STAT = ierr)
           IF(ierr .NE. 0) THEN
              message = "MasterClkPer"
              GOTO 999
           ENDIF
        ALLOCATE (this%CalSet(this%nLine), STAT = ierr)
           IF(ierr .NE. 0) THEN
              message = "CalSet"
              GOTO 999
           ENDIF
        ALLOCATE (this%ExpTime(this%nLine), STAT = ierr)
           IF(ierr .NE. 0) THEN
              message = "ExpTime"
              GOTO 999
           ENDIF
        ALLOCATE (this%ReadTime(this%nLine), STAT = ierr)
           IF(ierr .NE. 0) THEN
              message = "ReadTime"
              GOTO 999
           ENDIF
        ALLOCATE (this%SmPixCol(this%nLine), STAT = ierr)
           IF(ierr .NE. 0) THEN
              message = "SmPixCol"
              GOTO 999
           ENDIF
        ALLOCATE (this%GSC1(this%nLine), STAT = ierr)
           IF(ierr .NE. 0) THEN
              message = "GSC1"
              GOTO 999
           ENDIF
        ALLOCATE (this%GSC2(this%nLine), STAT = ierr)
           IF(ierr .NE. 0) THEN
              message = "GSC2"
              GOTO 999
           ENDIF
        ALLOCATE (this%GSC3(this%nLine), STAT = ierr)
           IF(ierr .NE. 0) THEN
              message = "GSC3"
              GOTO 999
           ENDIF
        ALLOCATE (this%GC1(this%nLine), STAT = ierr)
           IF(ierr .NE. 0) THEN
              message = "GC1"
              GOTO 999
           ENDIF
        ALLOCATE (this%GC2(this%nLine), STAT = ierr)
           IF(ierr .NE. 0) THEN
              message = "GC2"
              GOTO 999
           ENDIF
        ALLOCATE (this%GC3(this%nLine), STAT = ierr)
           IF(ierr .NE. 0) THEN
              message = "GC3"
              GOTO 999
           ENDIF
        ALLOCATE (this%GC4(this%nLine), STAT = ierr)
           IF(ierr .NE. 0) THEN
              message = "GC4"
              GOTO 999
           ENDIF
        ALLOCATE (this%DSGC(this%nLine), STAT = ierr)
           IF(ierr .NE. 0) THEN
              message = "DSGC"
              GOTO 999
           ENDIF
        ALLOCATE (this%LSLABF(this%nLine), STAT = ierr)
           IF(ierr .NE. 0) THEN
              message = "LSLABF"
              GOTO 999
           ENDIF
        ALLOCATE (this%USLABF(this%nLine), STAT = ierr)
           IF(ierr .NE. 0) THEN
              message = "USLABF"
              GOTO 999
           ENDIF
        ALLOCATE (this%LDABF(this%nLine), STAT = ierr)
           IF(ierr .NE. 0) THEN
              message = "LDABF"
              GOTO 999
           ENDIF
        ALLOCATE (this%UDABF(this%nLine), STAT = ierr)
           IF(ierr .NE. 0) THEN
              message = "UDABF"
              GOTO 999
           ENDIF
        ALLOCATE (this%SR1(this%nLine), STAT = ierr)
           IF(ierr .NE. 0) THEN
              message = "SR1"
              GOTO 999
           ENDIF
        ALLOCATE (this%SR2(this%nLine), STAT = ierr)
           IF(ierr .NE. 0) THEN
              message = "SR2"
              GOTO 999
           ENDIF
        ALLOCATE (this%SR3(this%nLine), STAT = ierr)
           IF(ierr .NE. 0) THEN
              message = "SR3"
              GOTO 999
           ENDIF
        ALLOCATE (this%SR4(this%nLine), STAT = ierr)
           IF(ierr .NE. 0) THEN
              message = "SR4"
              GOTO 999
           ENDIF
        ALLOCATE (this%DetcTemp(this%nLine), STAT = ierr)
           IF(ierr .NE. 0) THEN
              message = "DetcTemp"
              GOTO 999
           ENDIF
        ALLOCATE (this%OptBenchTemp(this%nLine), STAT = ierr)
           IF(ierr .NE. 0) THEN
              message = "OptBenchTemp"
              GOTO 999
           ENDIF
        ALLOCATE (this%ImBinFact(this%nLine), STAT = ierr)
           IF(ierr .NE. 0) THEN
              message = "ImBinFact"
              GOTO 999
           ENDIF
        ALLOCATE (this%BinImgRows(this%nLine), STAT = ierr)
           IF(ierr .NE. 0) THEN
              message = "BinImgRows"
              GOTO 999
           ENDIF
        ALLOCATE (this%StopColumn(this%nLine), STAT = ierr)
           IF(ierr .NE. 0) THEN
              message = "StopColumn"
              GOTO 999
           ENDIF

        status = getl1bblk( swid, &
                               "SmallPixelRadiance", DFNT_FLOAT32, &
                               fldflg, 0, 1, rank, dims, &
                                this%SmPixRad )
        IF( status .NE. OMI_S_SUCCESS ) THEN
           status = getl1bblk( swid, &
                                  "SmallPixelIrradiance", DFNT_FLOAT32, &
                                   fldflg, 0, 1, rank, dims, &
                                   this%SmPixRad )
           IF( status .NE. OMI_S_SUCCESS ) THEN
              status = getl1bblk( swid, &
                                     "SmallPixelSignal", DFNT_FLOAT32, &
                                      fldflg, 0, 1, rank, dims, &
                                      this%SmPixRad )
              IF( status .NE. OMI_S_SUCCESS ) THEN
                 ierr = OMI_SMF_setmsg( status, &
                                       "retrieve small pixel value name", &
                                       "SmPxr_open", one )
                 RETURN
              ELSE
                 this%fieldname = "SmallPixelSignal"
              ENDIF
           ELSE
              this%fieldname = "SmallPixelIrradiance"
           ENDIF
        ELSE
           this%fieldname = "SmallPixelRadiance"
        ENDIF

        this%initialized = .TRUE.
        status = OMI_S_SUCCESS
        RETURN      

999     CONTINUE
        status = OMI_E_FAILURE
        ierr = OMI_SMF_setmsg(OMI_E_MEM_ALLOC, message, "L1Br_open", zero)
        RETURN
        
      END FUNCTION L1Br_open 

!! 2. L1Br_getDIMsize
 !
 !    Functionality:
 !
 !       This function retrieves the dimension size for the named dimension.
 !       A user can use this function to find the size of each dimension in
 !       the L1B file and plan the amount of storage needed to in the calling 
 !       procedure.
 !
 !    Calling Arguments:
 !
 !    Inputs:
 !
 !       this:    the block data structure
 !       dimname: the dimension name 
 !
 !    Outputs:
 !
 !       size: the size for the input dimension, if error occurs, -1 is returned.
 !
 !    Change History:
 !
 !       Date            Author          Modifications
 !       ====            ======          =============
 !       January 2005    Jeremy Warner   Original Source (derived from original 
 !                                             L1B Reader code written by Kai Yang)
 !
!!
      FUNCTION L1Br_getDIMsize(this, dimname) RESULT(size) 
        INCLUDE 'PGS_OMI_1900.f'  !defines the L1B PGS error codes
        TYPE (L1B_block_type), INTENT(INOUT) :: this
        CHARACTER (LEN = *), INTENT(IN) :: dimname
        INTEGER (KIND = 4) :: size
        CHARACTER (LEN = MAX_NAME_LENGTH) :: msg

        IF(.NOT. this%initialized) THEN
           ierr = OMI_SMF_setmsg(OMI_E_INPUT, &
                                 "input block not initialized", &
                                 "L1Br_getDIMsize", zero)
           size = -1
           RETURN
        ENDIF

        SELECT CASE(dimname)
          CASE("NumTimes")
            size = this%NumTimes     
            RETURN
          CASE("nXtrack")
            size = this%nXtrack    
            RETURN
          CASE("nWavel")
            size = this%nWavel
            RETURN
          CASE("nWavelCoef")
            size = this%nWavelCoef 
            RETURN
          CASE("NumTimesSmallPixel")
            size = this%NumTimesSmallPixel     
            RETURN
          CASE DEFAULT
            msg = "unknown dimension name: " // dimname
            ierr = OMI_SMF_setmsg(OMI_E_INPUT, msg, "L1Br_getDIMsize", zero)
            size = -1
            RETURN      
        END SELECT 
      END FUNCTION L1Br_getDIMsize

!! 3. L1Br_getSWdims
 !
 !    Functionality:
 !
 !       This function retrieves the dimension sizes for the named dimension.
 !       A user can use this function to find the size of each dimension in
 !       the L1B file and plan the amount of storage needed to in the calling 
 !       procedure.  It is similar in function to L1Br_getDIMsize.
 !
 !    Calling Arguments:
 !
 !    Inputs:
 !
 !       this:    the block data structure
 !
 !    Outputs:
 !
 !       The following are keyword arguments.  Only those present in the argument 
 !       list will be set by the function.
 !
 !       Keyword                    Corresponding L1B Geolocation Data Field
 !       =======                    ========================================
 !       NumTimes_k                 total number of lines in the swath
 !       nXtrack_k                  number of cross-track pixels in the swath
 !       nWavel_k                   number of wavelengths for each line and pixel number
 !       nWavelCoef_k               number of wavelength coefficients
 !       NumTimesSmallPixel_k       total number of small pixel columns in the swath
 !
 !       status: the return PGS_SMF status value
 !
 !    Change History:
 !
 !       Date            Author          Modifications
 !       ====            ======          =============
 !       January 2005    Jeremy Warner   Original Source (derived from original 
 !                                             L1B Reader code written by Kai Yang)
 !
!!
      FUNCTION L1Br_getSWdims(this, NumTimes_k, nXtrack_k, NumTimesSmallPixel_k, &
                               nWavel_k, nWavelCoef_k,fieldname_k) RESULT(status)
        INCLUDE 'PGS_OMI_1900.f'  !defines the L1B PGS error codes
        TYPE (L1B_block_type), INTENT(INOUT) :: this
        INTEGER (KIND = 4), OPTIONAL, INTENT(OUT) :: NumTimes_k
        INTEGER (KIND = 4), OPTIONAL, INTENT(OUT) :: NumTimesSmallPixel_k
        INTEGER (KIND = 4), OPTIONAL, INTENT(OUT) :: nXtrack_k
        INTEGER (KIND = 4), OPTIONAL, INTENT(OUT) :: nWavel_k
        INTEGER (KIND = 4), OPTIONAL, INTENT(OUT) :: nWavelCoef_k
        CHARACTER (LEN = *), OPTIONAL, INTENT(OUT) :: fieldname_k
        INTEGER (KIND = 4) :: status

        status = OMI_S_SUCCESS
        IF(.NOT. this%initialized) THEN
           ierr = OMI_SMF_setmsg(OMI_E_INPUT, "input block not initialized", &
                                 "L1Br_getSWdims", zero)
           status = OMI_E_FAILURE
           RETURN
        ENDIF

        IF(PRESENT(NumTimes_k))           NumTimes_k = this%NumTimes
        IF(PRESENT(nXtrack_k))            nXtrack_k = this%nXtrack
        IF(PRESENT(nWavel_k))             nWavel_k = this%nWavel
        IF(PRESENT(nWavelCoef_k))         nWavelCoef_k = this%nWavelCoef
        IF(PRESENT(NumTimesSmallPixel_k)) NumTimesSmallPixel_k = this%NumTimesSmallPixel
        IF(PRESENT(fieldname_k))           fieldname_k=this%fieldname
        RETURN
      END FUNCTION L1Br_getSWdims

!! Private function: init_l1b_blk
 !
 !    Functionality:
 !
 !       This is a private function which is only used by other functions 
 !       in this MODULE to initialize the block data structure with fill values. 
 !
 !    Calling Arguments:
 !
 !    Inputs:
 !
 !       this:    the block data structure
 !
 !    Outputs:
 !
 !       status:  the return PGS_SMF status value; this function only returns OMI_S_SUCCESS
 !
 !    Change History:
 !
 !       Date            Author          Modifications
 !       ====            ======          =============
 !       January 2005    Jeremy Warner   Original Source
 !
!!
      FUNCTION init_l1b_blk(this) RESULT(status) 
        INCLUDE 'PGS_OMI_1900.f'  !defines the L1B PGS error codes
        TYPE (L1B_block_type), INTENT(INOUT) :: this
        INTEGER (KIND = 4) :: status

        status = OMI_S_SUCCESS

        this%Time(1:this%nLine) = r8fill
        this%SecInDay(1:this%nLine) = rfill
        this%ScLat(1:this%nLine) = rfill
        this%ScLon(1:this%nLine) = rfill
        this%ScAlt(1:this%nLine) = rfill
        this%SolElevation(1:this%nLine) = rfill
        this%SolElevMin(1:this%nLine) = rfill
        this%SolElevMax(1:this%nLine) = rfill
        this%SolAzimuth(1:this%nLine) = rfill
        this%SolAziMin(1:this%nLine) = rfill
        this%SolAziMax(1:this%nLine) = rfill
        this%Lon(1:this%nXtrack,1:this%nLine) = rfill
        this%Lat(1:this%nXtrack,1:this%nLine) = rfill
        this%SolZenAng(1:this%nXtrack,1:this%nLine) = rfill
        this%SolAziAng(1:this%nXtrack,1:this%nLine) = rfill
        this%ViewZenAng(1:this%nXtrack,1:this%nLine) = rfill
        this%ViewAziAng(1:this%nXtrack,1:this%nLine) = rfill
        this%TerrainHeight(1:this%nXtrack,1:this%nLine) = int16fill
        this%GPQFlag(1:this%nXtrack,1:this%nLine) = int16fill
        this%XTQFlag(1:this%nXtrack,1:this%nLine) = int8fill
        this%RadMantissa(1:this%nWavel,1:this%nXtrack,1:this%nLine) = int16fill
        this%RadPrecision(1:this%nWavel,1:this%nXtrack,1:this%nLine) = int16fill
        this%RadExponent(1:this%nWavel,1:this%nXtrack,1:this%nLine) = int8fill
        this%PQFlag(1:this%nWavel,1:this%nXtrack,1:this%nLine) = int16fill
        this%PFlag(1:this%nXtrack,1:this%nLine) = int16fill
        this%WavelenCoef(1:this%nWavelCoef,1:this%nXtrack,1:this%nLine) = rfill
        this%WavelenPrec(1:this%nWavelCoef,1:this%nXtrack,1:this%nLine) = rfill
        this%WavelenRefCol(1:this%nLine) = int16fill
        this%SmPixRad(1:this%nXtrack,1:this%NumTimesSmallPixel) = rfill
        this%SmPixWavelen(1:this%nXtrack,1:this%NumTimesSmallPixel) = rfill
        this%MeasClass(1:this%nLine) = int8fill
        this%Config(1:this%nLine) = int8fill
        this%MQFlag(1:this%nLine) = int16fill
        this%NumSmPixCol(1:this%nLine) = int8fill
        this%NumSmPix(1:this%NumTimes) = int8fill
        this%ExposureType(1:this%nLine) = int8fill
        this%MasterClkPer(1:this%nLine) = rfill
        this%CalSet(1:this%nLine) = int16fill
        this%ExpTime(1:this%nLine) = rfill
        this%ReadTime(1:this%nLine) = rfill
        this%SmPixCol(1:this%nLine) = int16fill
        this%GSC1(1:this%nLine) = int16fill
        this%GSC2(1:this%nLine) = int16fill
        this%GSC3(1:this%nLine) = int16fill
        this%GC1(1:this%nLine) = int8fill
        this%GC2(1:this%nLine) = int8fill
        this%GC3(1:this%nLine) = int8fill
        this%GC4(1:this%nLine) = int8fill
        this%DSGC(1:this%nLine) = int8fill
        this%LSLABF(1:this%nLine) = int8fill
        this%USLABF(1:this%nLine) = int8fill
        this%LDABF(1:this%nLine) = int8fill
        this%UDABF(1:this%nLine) = int8fill
        this%SR1(1:this%nLine) = int16fill
        this%SR2(1:this%nLine) = int16fill
        this%SR3(1:this%nLine) = int16fill
        this%SR4(1:this%nLine) = int16fill
        this%DetcTemp(1:this%nLine) = rfill
        this%OptBenchTemp(1:this%nLine) = rfill
        this%ImBinFact(1:this%nLine) = int8fill
        this%BinImgRows(1:this%nLine) = int16fill
        this%StopColumn(1:this%nLine) = int16fill

      END FUNCTION init_l1b_blk

!! Private function: fill_l1b_blk
 !
 !    Functionality:
 !
 !       This is a private function which is only used by other functions 
 !       in this MODULE to fill the block data structure with data from the L1B file.
 !
 !    Calling Arguments:
 !
 !    Inputs:
 !
 !       this:    the block data structure
 !       iLine:   the starting line for the data block
 !
 !    Outputs:
 !
 !       status:  the return PGS_SMF status value
 !
 !    Change History:
 !
 !       Date            Author          Modifications
 !       ====            ======          =============
 !       January 2005    Jeremy Warner   Original Source (derived from original 
 !                                             L1B Reader code written by Kai Yang)
 !
!!
      FUNCTION fill_l1b_blk(this, iLine) RESULT(status) 
        USE hdfeos4_parameters
        INCLUDE 'PGS_OMI_1900.f'  !defines the L1B PGS error codes
        TYPE (L1B_block_type), INTENT(INOUT) :: this
        INTEGER (KIND = 4), INTENT(IN) :: iLine
        INTEGER (KIND = 4) :: status
        INTEGER (KIND = 4) :: getl1bblk
        EXTERNAL              getl1bblk
        INTEGER (KIND = 4) :: rank
        INTEGER (KIND = 4) :: nL, swid, swfid, fldflg
        INTEGER (KIND = 4), DIMENSION(1:3) :: dims
        INTEGER (KIND = 4) :: i, j, k

        IF(.NOT. this%initialized) THEN
           ierr = OMI_SMF_setmsg(OMI_E_INPUT, &
                                 "input block not initialized", &
                                 "fill_l1b_blk", zero)
           status = OMI_E_FAILURE
           RETURN
        ENDIF

        this%iLine = iLine
        
        IF(iLine < 0 .OR. iLine >= this%NumTimes) THEN
           ierr = OMI_SMF_setmsg(OMI_E_INPUT, &
                                 "iLine out of range", "fill_l1b_blk", zero)
           status = OMI_E_FAILURE
           RETURN
        ENDIF

        ! If, from the starting line, the number of lines to be read into
        ! the block data structure goes past the final line in the L1B swath,
        ! then recalculate the number of lines to be read. 

        IF((iLine + this%nLine) > this%NumTimes) THEN          
           nL = this%NumTimes - iLine 
        ELSE
           nL = this%nLine
        ENDIF
        this%eLine = this%iLine + nL - 1

        ! Open the L1B swath file

        swfid = swopen(this%filename, DFACC_READ)
        IF(swfid < zero) THEN
           status = OMI_E_FAILURE
           ierr = OMI_SMF_setmsg(OMI_E_FILE_OPEN, this%filename, &
                 "fill_l1b_blk", zero)
           RETURN
        ENDIF 

        ! Attach to the swath

        swid = swattach(swfid, this%swathname)
        IF(swid < zero) THEN
           status = OMI_E_FAILURE
           ierr = OMI_SMF_setmsg(OMI_E_SWATH_ATTACH, this%swathname, &
                 "fill_l1b_blk", zero)
           ierr = swclose(swfid)
           RETURN
        ENDIF 
        
	! Initialize L1B block to fill values

	ierr = init_l1b_blk(this)

        ! Read all L1B datasets

        status = getl1bblk(swid, "Time", DFNT_FLOAT64, & 
                                fldflg, iLine, nL, rank, dims, this%Time)
        IF(status .NE. OMI_S_SUCCESS) goto 1001
 
        status = getl1bblk(swid, "SecondsInDay", DFNT_FLOAT32, & 
                                fldflg, iLine, nL, rank, dims, this%SecInDay)
        IF(status .NE. OMI_S_SUCCESS) goto 1001
 
        status = getl1bblk(swid, "SpacecraftLatitude", DFNT_FLOAT32, & 
                                fldflg, iLine, nL, rank, dims, this%ScLat)
        IF(status .NE. OMI_S_SUCCESS) goto 1001
 
        status = getl1bblk(swid, "SpacecraftLongitude", DFNT_FLOAT32, & 
                                fldflg, iLine, nL, rank, dims, this%ScLon)
        IF(status .NE. OMI_S_SUCCESS) goto 1001
 
        status = getl1bblk(swid, "SpacecraftAltitude", DFNT_FLOAT32, & 
                                fldflg, iLine, nL, rank, dims, this%ScAlt)
        IF(status .NE. OMI_S_SUCCESS) goto 1001
 
        status = getl1bblk(swid, "SolarElevation", DFNT_FLOAT32, & 
                                fldflg, iLine, nL, rank, dims, this%SolElevation)
        IF(status .NE. OMI_S_SUCCESS) goto 1001
 
        status = getl1bblk(swid, "SolarAzimuth", DFNT_FLOAT32, & 
                                fldflg, iLine, nL, rank, dims, this%SolAzimuth)
        IF(status .NE. OMI_S_SUCCESS) goto 1001
 
        status = getl1bblk(swid, "SolarElevationMinimum", DFNT_FLOAT32, & 
                                fldflg, iLine, nL, rank, dims, this%SolElevMin)
        IF(status .NE. OMI_S_SUCCESS) goto 1001
 
        status = getl1bblk(swid, "SolarElevationMaximum", DFNT_FLOAT32, & 
                                fldflg, iLine, nL, rank, dims, this%SolElevMax)
        IF(status .NE. OMI_S_SUCCESS) goto 1001
 
        status = getl1bblk(swid, "SolarAzimuthMinimum", DFNT_FLOAT32, & 
                                fldflg, iLine, nL, rank, dims, this%SolAziMin)
        IF(status .NE. OMI_S_SUCCESS) goto 1001
 
        status = getl1bblk(swid, "SolarAzimuthMaximum", DFNT_FLOAT32, & 
                                fldflg, iLine, nL, rank, dims, this%SolAziMax)
        IF(status .NE. OMI_S_SUCCESS) goto 1001
 
        status = getl1bblk(swid, "Latitude", DFNT_FLOAT32, & 
                                fldflg, iLine, nL, rank, dims, this%Lat)
        IF(status .NE. OMI_S_SUCCESS) goto 1001
 
        status = getl1bblk(swid, "Longitude", DFNT_FLOAT32, & 
                                fldflg, iLine, nL, rank, dims, this%Lon)
        IF(status .NE. OMI_S_SUCCESS) goto 1001

        status = getl1bblk(swid, "SolarZenithAngle", DFNT_FLOAT32, & 
                                fldflg, iLine, nL, rank, dims, this%SolZenAng)
        IF(status .NE. OMI_S_SUCCESS) goto 1001

        status = getl1bblk(swid, "SolarAzimuthAngle", DFNT_FLOAT32, & 
                                fldflg, iLine, nL, rank, dims, this%SolAziAng)
        IF(status .NE. OMI_S_SUCCESS) goto 1001

        status = getl1bblk(swid, "ViewingZenithAngle", DFNT_FLOAT32, & 
                                fldflg, iLine, nL, rank, dims, this%ViewZenAng)
        IF(status .NE. OMI_S_SUCCESS) goto 1001
        
        status = getl1bblk(swid, "ViewingAzimuthAngle", DFNT_FLOAT32, & 
                                fldflg, iLine, nL, rank, dims, this%ViewAziAng)
        IF(status .NE. OMI_S_SUCCESS) goto 1001

        status = getl1bblk(swid, "TerrainHeight", DFNT_INT16, & 
                                fldflg, iLine, nL, rank, dims, this%TerrainHeight)
        IF(status .NE. OMI_S_SUCCESS) goto 1001

        status = getl1bblk(swid, "GroundPixelQualityFlags", DFNT_UINT16, & 
                                fldflg, iLine, nL, rank, dims, this%GPQFlag)
        IF(status .NE. OMI_S_SUCCESS) goto 1001

        status = getl1bblk(swid, "XTrackQualityFlags", DFNT_UINT8, &
                                fldflg, iLine, nL, rank, dims, this%XTQFlag)
        IF(status .NE. OMI_S_SUCCESS) goto 1001

        status = getl1bblk(swid, "RadianceMantissa", DFNT_INT16, & 
                                fldflg, iLine, nL, rank, dims, this%RadMantissa)
        IF(status .NE. OMI_S_SUCCESS) goto 1001

        if (fldflg .ne. 0) then
           status = getl1bblk(swid, "IrradianceMantissa", DFNT_INT16, & 
                                fldflg, iLine, nL, rank, dims, this%RadMantissa)
           IF(status .NE. OMI_S_SUCCESS) goto 1001
           if (fldflg .ne. 0) then
              status = getl1bblk(swid, "SignalMantissa", DFNT_INT16, & 
                             fldflg, iLine, nL, rank, dims, this%RadMantissa)
              IF(status .NE. OMI_S_SUCCESS) goto 1001
           endif
        endif

        status = getl1bblk(swid, "RadiancePrecisionMantissa", DFNT_INT16, & 
                                fldflg, iLine, nL, rank, dims, this%RadPrecision)
        IF(status .NE. OMI_S_SUCCESS) goto 1001
        if (fldflg .ne. 0) then
           status = getl1bblk(swid, "IrradiancePrecisionMantissa", DFNT_INT16,&
                                fldflg, iLine, nL, rank, dims, this%RadPrecision)
           IF(status .NE. OMI_S_SUCCESS) goto 1001
           if (fldflg .ne. 0) then
              status = getl1bblk(swid,"SignalPrecisionMantissa",DFNT_INT16,& 
                             fldflg, iLine, nL, rank, dims, this%RadPrecision)
              IF(status .NE. OMI_S_SUCCESS) goto 1001
           endif
        endif

        status = getl1bblk(swid, "RadianceExponent", DFNT_INT8, & 
                                fldflg, iLine, nL, rank, dims, this%RadExponent)
        IF(status .NE. OMI_S_SUCCESS) goto 1001
        if (fldflg .ne. 0) then
           status = getl1bblk(swid, "IrradianceExponent", DFNT_INT8, & 
                                fldflg, iLine, nL, rank, dims, this%RadExponent)
           IF(status .NE. OMI_S_SUCCESS) goto 1001
           if (fldflg .ne. 0) then
              status = getl1bblk(swid, "SignalExponent", DFNT_INT8, & 
                             fldflg, iLine, nL, rank, dims, this%RadExponent)
              IF(status .NE. OMI_S_SUCCESS) goto 1001
           endif
        endif

        status = getl1bblk(swid, "PixelQualityFlags", DFNT_UINT16, & 
                                fldflg, iLine, nL, rank, dims, this%PQFlag)
        IF(status .NE. OMI_S_SUCCESS) goto 1001

        DO j = 1, this%nLine
           k = this%SmPixCol(j)
           IF (k .gt. 0) then
             DO i = 1, this%nXtrack
                this%PFlag(i,j) = this%PQFlag(k,i,j)
             ENDDO
           ENDIF
        ENDDO

        status = getl1bblk(swid, "WavelengthCoefficient", DFNT_FLOAT32, & 
                                fldflg, iLine, nL, rank, dims, this%WavelenCoef)
        IF(status .NE. OMI_S_SUCCESS) goto 1001

        status = getl1bblk(swid,"WavelengthCoefficientPrecision",DFNT_FLOAT32,&
                                fldflg, iLine, nL, rank, dims, this%WavelenPrec)
        IF(status .NE. OMI_S_SUCCESS) goto 1001

        status = getl1bblk(swid, "WavelengthReferenceColumn", DFNT_INT16, &
                                fldflg, iLine, nL, rank, dims, this%WavelenRefCol)
        IF(status .NE. OMI_S_SUCCESS) goto 1001

        status = getl1bblk(swid, "SmallPixelRadiance", DFNT_FLOAT32, & 
                                fldflg, 0, this%NumTimesSmallPixel, rank, dims, this%SmPixRad)
        IF(status .NE. OMI_S_SUCCESS) goto 1001
        if (fldflg .ne. 0) then
           status = getl1bblk(swid, "SmallPixelIrradiance", DFNT_FLOAT32, & 
                                fldflg, 0, this%NumTimesSmallPixel, rank, dims, this%SmPixRad)
           IF(status .NE. OMI_S_SUCCESS) goto 1001
           if (fldflg .ne. 0) then
              status = getl1bblk(swid, "SmallPixelSignal", DFNT_FLOAT32, & 
                             fldflg, 0, this%NumTimesSmallPixel, rank, dims, this%SmPixRad)
              IF(status .NE. OMI_S_SUCCESS) goto 1001
           endif
        endif

        status = getl1bblk(swid, "SmallPixelWavelength", DFNT_FLOAT32, fldflg, 0, & 
                                this%NumTimesSmallPixel, rank, dims, this%SmPixWavelen)
        IF(status .NE. OMI_S_SUCCESS) goto 1001

        status = getl1bblk(swid, "MeasurementClass", DFNT_UINT8, & 
                                fldflg, iLine, nL, rank, dims, this%MeasClass)
        IF(status .NE. OMI_S_SUCCESS) goto 1001

        status = getl1bblk(swid, "InstrumentConfigurationId", DFNT_UINT8, & 
                                fldflg, iLine, nL, rank, dims, this%Config)
        IF(status .NE. OMI_S_SUCCESS) goto 1001
 
        status = getl1bblk(swid, "MeasurementQualityFlags", DFNT_UINT16, & 
                                fldflg, iLine, nL, rank, dims, this%MQFlag)
        IF(status .NE. OMI_S_SUCCESS) goto 1001

        status = getl1bblk(swid, "NumberSmallPixelColumns", DFNT_INT8, & 
                                fldflg, iLine, nL, rank, dims, this%NumSmPixCol)
        IF(status .NE. OMI_S_SUCCESS) goto 1001

        status = getl1bblk(swid, "NumberSmallPixelColumns", DFNT_INT8, & 
                                fldflg, 0, this%NumTimes, rank, dims, this%NumSmPix)
        IF(status .NE. OMI_S_SUCCESS) goto 1001

        status = getl1bblk(swid, "ExposureType", DFNT_INT8, & 
                                fldflg, iLine, nL, rank, dims, this%ExposureType)
        IF(status .NE. OMI_S_SUCCESS) goto 1001

        status = getl1bblk(swid, "MasterClockPeriod", DFNT_FLOAT32, & 
                                fldflg, iLine, nL, rank, dims, this%MasterClkPer)
        IF(status .NE. OMI_S_SUCCESS) goto 1001

        status = getl1bblk(swid, "CalibrationSettings", DFNT_UINT16, & 
                                fldflg, iLine, nL, rank, dims, this%CalSet)
        IF(status .NE. OMI_S_SUCCESS) goto 1001

        status = getl1bblk(swid, "ExposureTime", DFNT_FLOAT32, & 
                                fldflg, iLine, nL, rank, dims, this%ExpTime)
        IF(status .NE. OMI_S_SUCCESS) goto 1001

        status = getl1bblk(swid, "ReadoutTime", DFNT_FLOAT32, & 
                                fldflg, iLine, nL, rank, dims, this%ReadTime)
        IF(status .NE. OMI_S_SUCCESS) goto 1001

        status = getl1bblk(swid, "SmallPixelColumn", DFNT_INT16, & 
                                fldflg, iLine, nL, rank, dims, this%SmPixCol)
        IF(status .NE. OMI_S_SUCCESS) goto 1001

        status = getl1bblk(swid, "GainSwitchingColumn1", DFNT_INT16, & 
                                fldflg, iLine, nL, rank, dims, this%GSC1)
        IF(status .NE. OMI_S_SUCCESS) goto 1001

        status = getl1bblk(swid, "GainSwitchingColumn2", DFNT_INT16, & 
                                fldflg, iLine, nL, rank, dims, this%GSC2)
        IF(status .NE. OMI_S_SUCCESS) goto 1001

        status = getl1bblk(swid, "GainSwitchingColumn3", DFNT_INT16, & 
                                fldflg, iLine, nL, rank, dims, this%GSC3)
        IF(status .NE. OMI_S_SUCCESS) goto 1001

        status = getl1bblk(swid, "GainCode1", DFNT_INT8, & 
                                fldflg, iLine, nL, rank, dims, this%GC1)
        IF(status .NE. OMI_S_SUCCESS) goto 1001

        status = getl1bblk(swid, "GainCode2", DFNT_INT8, & 
                                fldflg, iLine, nL, rank, dims, this%GC2)
        IF(status .NE. OMI_S_SUCCESS) goto 1001

        status = getl1bblk(swid, "GainCode3", DFNT_INT8, & 
                                fldflg, iLine, nL, rank, dims, this%GC3)
        IF(status .NE. OMI_S_SUCCESS) goto 1001

        status = getl1bblk(swid, "GainCode4", DFNT_INT8, & 
                                fldflg, iLine, nL, rank, dims, this%GC4)
        IF(status .NE. OMI_S_SUCCESS) goto 1001

        status = getl1bblk(swid, "DSGainCode", DFNT_INT8, & 
                                fldflg, iLine, nL, rank, dims, this%DSGC)
        IF(status .NE. OMI_S_SUCCESS) goto 1001

        status = getl1bblk(swid, "LowerStrayLightAreaBinningFactor", DFNT_INT8, & 
                                fldflg, iLine, nL, rank, dims, this%LSLABF)
        IF(status .NE. OMI_S_SUCCESS) goto 1001

        status = getl1bblk(swid, "UpperStrayLightAreaBinningFactor", DFNT_INT8, & 
                                fldflg, iLine, nL, rank, dims, this%USLABF)
        IF(status .NE. OMI_S_SUCCESS) goto 1001

        status = getl1bblk(swid, "LowerDarkAreaBinningFactor", DFNT_INT8, & 
                                fldflg, iLine, nL, rank, dims, this%LDABF)
        IF(status .NE. OMI_S_SUCCESS) goto 1001

        status = getl1bblk(swid, "UpperDarkAreaBinningFactor", DFNT_INT8, & 
                                fldflg, iLine, nL, rank, dims, this%UDABF)
        IF(status .NE. OMI_S_SUCCESS) goto 1001

        status = getl1bblk(swid, "SkipRows1", DFNT_INT16, & 
                                fldflg, iLine, nL, rank, dims, this%SR1)
        IF(status .NE. OMI_S_SUCCESS) goto 1001

        status = getl1bblk(swid, "SkipRows2", DFNT_INT16, & 
                                fldflg, iLine, nL, rank, dims, this%SR2)
        IF(status .NE. OMI_S_SUCCESS) goto 1001

        status = getl1bblk(swid, "SkipRows3", DFNT_INT16, & 
                                fldflg, iLine, nL, rank, dims, this%SR3)
        IF(status .NE. OMI_S_SUCCESS) goto 1001

        status = getl1bblk(swid, "SkipRows4", DFNT_INT16, & 
                                fldflg, iLine, nL, rank, dims, this%SR4)
        IF(status .NE. OMI_S_SUCCESS) goto 1001

        status = getl1bblk(swid, "DetectorTemperature", DFNT_FLOAT32, & 
                                fldflg, iLine, nL, rank, dims, this%DetcTemp)
        IF(status .NE. OMI_S_SUCCESS) goto 1001

        status = getl1bblk(swid, "OpticalBenchTemperature", DFNT_FLOAT32, & 
                                fldflg, iLine, nL, rank, dims, this%OptBenchTemp)
        IF(status .NE. OMI_S_SUCCESS) goto 1001

        status = getl1bblk(swid, "ImageBinningFactor", DFNT_INT8, & 
                                fldflg, iLine, nL, rank, dims, this%ImBinFact)
        IF(status .NE. OMI_S_SUCCESS) goto 1001

        status = getl1bblk(swid, "BinnedImageRows", DFNT_INT16, & 
                                fldflg, iLine, nL, rank, dims, this%BinImgRows)
        IF(status .NE. OMI_S_SUCCESS) goto 1001

        status = getl1bblk(swid, "StopColumn", DFNT_INT16, & 
                                fldflg, iLine, nL, rank, dims, this%StopColumn)
        IF(status .NE. OMI_S_SUCCESS) goto 1001

1001    continue
        ierr = swdetach(swid)
        ierr = swclose(swfid)
        IF(status .NE. OMI_S_SUCCESS) THEN
           ierr = OMI_SMF_setmsg(OMI_E_FAILURE, "Unrecoverable error reading L1B file", &
                                 "fill_l1b_blk", zero)
        ENDIF
        RETURN
      END FUNCTION fill_l1b_blk

!! Private function: check_blk_pix
 !
 ! Functionality:
 !
 !   This function checks to make sure the requested indecies are within
 !   range and returns 1-based indecies (from 0-based inputs).  Note that
 !   input indicies refer to the file, while output indicies refer to
 !   the structure in memory (only a part of the input L1B file is
 !   stored memory at a given time).
 !
 ! Calling Arguments:
 !
 ! Inputs:
 !
 !    this      L1B block
 !    iLine     requested line (0-indexed) of L1B file
 !    iPix      requested Xtrack position (0-indexed) of L1B file
 !
 ! Outputs:
 !
 !    j         Index into L1B structure, line number, (1-indexed)
 !    i         Index into L1B structure, Xtrack position, (1-indexed)
 !
 !    status    the return PGS_SMF status value
 !
 ! Change History:
 !
 !    Date            Author          Modifications
 !    ====            ======          =============
 !    January 2005    Jeremy Warner   Original Source
 !
!!
      FUNCTION check_blk_pix(this, iLine, iPix, j, i) RESULT(status)
        INCLUDE 'PGS_OMI_1900.f'  !defines the L1B PGS error codes
        TYPE (L1B_block_type), INTENT(INOUT) :: this
        INTEGER (KIND = 4), INTENT(IN) :: iLine, iPix
        INTEGER (KIND = 4), INTENT(OUT) :: i, j
        INTEGER (KIND = 4) :: status

        status = OMI_S_SUCCESS 
        IF(.NOT. this%initialized) THEN
           ierr = OMI_SMF_setmsg(OMI_E_INPUT, "input block not initialized", &
                                 "check_blk_pix", zero)
           status = OMI_E_FAILURE
           RETURN
        ENDIF

        IF(iLine < 0 .OR. iLine >= this%NumTimes) THEN
           ierr = OMI_SMF_setmsg(OMI_E_INPUT, "iLine out of range", &
                                 "check_blk_pix", zero)
           status = OMI_E_FAILURE
           RETURN
        ENDIF

        IF(iPix < 0 .OR. iPix >= this%nXtrack) THEN
           ierr = OMI_SMF_setmsg(OMI_E_INPUT, "iPix out of range", &
                                 "check_blk_pix", zero)
           status = OMI_E_FAILURE
           RETURN
        ENDIF

        IF(iLine < this%iLine .OR. iLine > this%eLine) THEN 
           status = fill_l1b_blk(this, iLine)
           IF(status .NE. OMI_S_SUCCESS) THEN
              ierr = OMI_SMF_setmsg(OMI_E_DATA_BLOCK, "retrieve data block", &
                                    "check_blk_pix", one)
              RETURN
           ENDIF
        ENDIF

        i = iPix + 1
        j = iLine - this%iLine + 1
        RETURN
      END FUNCTION check_blk_pix

!! Private function: check_blk_line
 !
 ! Functionality:
 !
 !   This function checks to make sure the requested index is within
 !   range and returns 1-based index (from 0-based input).  Note that
 !   input index refers to the file, while output index refers to
 !   the structure in memory (only a part of the input L1B file is
 !   stored memory at a given time).
 !
 ! Calling Arguments:
 !
 ! Inputs:
 !
 !    this      L1B block
 !    iLine     requested line (0-indexed) of the L1B file
 !
 ! Outputs:
 !
 !    j         Index into L1B structure, line number, (1-indexed)
 !
 !    status    the return PGS_SMF status value
 !
 ! Change History:
 !
 !    Date            Author          Modifications
 !    ====            ======          =============
 !    January 2005    Jeremy Warner   Original Source
 !
!!
      FUNCTION check_blk_line(this, iLine, j) RESULT(status)
        INCLUDE 'PGS_OMI_1900.f'  !defines the L1B PGS error codes
        TYPE (L1B_block_type), INTENT(INOUT) :: this
        INTEGER (KIND = 4), INTENT(IN) :: iLine
        INTEGER (KIND = 4), INTENT(OUT) :: j
        INTEGER (KIND = 4) :: ierr
        INTEGER (KIND = 4) :: status

        status = OMI_S_SUCCESS 
        IF(.NOT. this%initialized) THEN
           ierr = OMI_SMF_setmsg(OMI_E_INPUT, "input block not initialized", &
                                 "check_blk_line", zero)
           status = OMI_E_FAILURE
           RETURN
        ENDIF

        IF(iLine < 0 .OR. iLine >= this%NumTimes) THEN
           ierr = OMI_SMF_setmsg(OMI_E_INPUT, "iLine out of range", &
                                 "check_blk_line", zero)
           status = OMI_E_FAILURE
           RETURN
        ENDIF

        IF(iLine < this%iLine .OR. iLine > this%eLine) THEN 
           status = fill_l1b_blk(this, iLine)
           IF(status .NE. OMI_S_SUCCESS) THEN
              ierr = OMI_SMF_setmsg(OMI_E_DATA_BLOCK, "retrieve data block", &
                                    "check_blk_line", one)
              RETURN
           ENDIF
        ENDIF

        j = iLine - this%iLine + 1
        RETURN
      END FUNCTION check_blk_line

!! 4. L1Br_getGEOpix
 !
 !    Functionality:
 !
 !       This function gets one or more geolocation data field values from
 !       the data block.  This function returns only a single pixel value
 !       for each given keyword (i.e., a value for a single line and a single
 !       Xtrack position).
 !
 !    Calling Arguments:
 !
 !    Inputs:
 !
 !       this: the block data structure
 !       iLine: the line number in the L1B swath.  NOTE: this input is 0 based
 !              range from 0 to (NumTimes-1) inclusive.
 !       iPix: the pixel number in the L1B swath.  NOTE: this input is 0 based
 !             range from 0 to (nXtrack-1) inclusive.
 !
 !    Outputs:
 !
 !       The following are keyword arguments.  Only those present in the argument 
 !       list will be set by the function.
 !
 !       Keyword                    Corresponding L1B Geolocation Data Field
 !       =======                    ========================================
 !       Time_k                     Time
 !       SecondsInDay_k             SecondsInDay
 !       SpacecraftLatitude_k       SpacecraftLatitude
 !       SpacecraftLongitude_k      SpacecraftLongitude
 !       SpacecraftAltitude_k       SpacecraftAltitude
 !       SolarElevation_k           SolarElevation
 !       SolarAzimuth_k             SolarAzimuth
 !       SolarElevationMinimum_k    SolarElevationMinimum
 !       SolarElevationMaximum_k    SolarElevationMaximum
 !       SolarAzimuthMinimum_k      SolarAzimuthMinimum
 !       SolarAzimuthMaximum_k      SolarAzimuthMaximum
 !       Latitude_k                 Latitude
 !       Longitude_k                Longitude
 !       SolarZenithAngle_k         SolarZenithAngle
 !       SolarAzimuthAngle_k        SolarAzimuthAngle
 !       ViewingZenithAngle_k       ViewingZenithAngle
 !       ViewingAzimuthAngle_k      ViewingAzimuthAngle
 !       TerrainHeight_k            TerrainHeight
 !       GroundPixelQualityFlags_k  GroundPixelQualityFlags
 !       XTrackQualityFlags_k       XTrackQualityFlags
 !
 !       status: the return PGS_SMF status value
 !
 !    Change History:
 !
 !       Date            Author          Modifications
 !       ====            ======          =============
 !       January 2005    Jeremy Warner   Original Source ((derived from original 
 !                                             L1B Reader code written by Kai Yang)
 !
!!
      FUNCTION L1Br_getGEOpix(this, iLine, iPix, Time_k, SecondsInDay_k, &
           SpacecraftLatitude_k, SpacecraftLongitude_k, &
           SpacecraftAltitude_k, &
           SolarElevation_k, SolarAzimuth_k, &
           SolarElevationMinimum_k, SolarElevationMaximum_k, &
           SolarAzimuthMinimum_k, SolarAzimuthMaximum_k, &
           Latitude_k, Longitude_k, SolarZenithAngle_k, &
           SolarAzimuthAngle_k, ViewingZenithAngle_k, &
           ViewingAzimuthAngle_k, TerrainHeight_k, &
           GroundPixelQualityFlags_k, XTrackQualityFlags_k) RESULT (status)
      
        INCLUDE 'PGS_OMI_1900.f'  !defines the L1B PGS error codes
        TYPE (L1B_block_type), INTENT(INOUT) :: this
        INTEGER (KIND = 4), INTENT(IN) :: iLine, iPix
        REAL (KIND =8), OPTIONAL, INTENT(OUT) :: Time_k
        REAL (KIND =4), OPTIONAL, INTENT(OUT) :: Latitude_k, Longitude_k, &
             SpacecraftLatitude_k, SpacecraftLongitude_k, &
             SpacecraftAltitude_k, SolarElevation_k, &
             SolarAzimuth_k, SolarElevationMinimum_k, &
             SolarElevationMaximum_k, SolarAzimuthMinimum_k, &
             SolarAzimuthMaximum_k, SolarZenithAngle_k, &
             SolarAzimuthAngle_k, ViewingZenithAngle_k, &
             ViewingAzimuthAngle_k, SecondsInDay_k
        INTEGER (KIND = 2), OPTIONAL, INTENT(OUT) :: TerrainHeight_k
        INTEGER (KIND = 2), OPTIONAL, INTENT(OUT) :: GroundPixelQualityFlags_k
        INTEGER (KIND = 1), OPTIONAL, INTENT(OUT) :: XTrackQualityFlags_k
        INTEGER :: i,j
        INTEGER (KIND = 4) :: status
 
        status = check_blk_pix(this, iLine, iPix, j, i) 
        IF(status .NE. OMI_S_SUCCESS) THEN
           ierr = OMI_SMF_setmsg(OMI_E_GENERAL, "Failed retrieving geolocation data.", &
                                 "L1Br_getGEOpix", one)
           RETURN
        ENDIF

        IF(PRESENT(Time_k))                    Time_k = this%Time(j)
        IF(PRESENT(SecondsInDay_k))            SecondsInDay_k = this%SecInDay(j)
        IF(PRESENT(SpacecraftLatitude_k))      SpacecraftLatitude_k = this%ScLat(j)
        IF(PRESENT(SpacecraftLongitude_k))     SpacecraftLongitude_k = this%ScLon(j)
        IF(PRESENT(SpacecraftAltitude_k))      SpacecraftAltitude_k = this%ScAlt(j)
        IF(PRESENT(SolarElevation_k))          SolarElevation_k = this%SolElevation(j)
        IF(PRESENT(SolarAzimuth_k))            SolarAzimuth_k = this%SolAzimuth(j)
        IF(PRESENT(SolarElevationMinimum_k))   SolarElevationMinimum_k = this%SolElevMin(j)
        IF(PRESENT(SolarElevationMaximum_k))   SolarElevationMaximum_k = this%SolElevMax(j)
        IF(PRESENT(SolarAzimuthMinimum_k))     SolarAzimuthMinimum_k = this%SolAziMin(j)
        IF(PRESENT(SolarAzimuthMaximum_k))     SolarAzimuthMaximum_k = this%SolAziMax(j)
        IF(PRESENT(Latitude_k))                Latitude_k = this%Lat(i,j)
        IF(PRESENT(Longitude_k))               Longitude_k = this%Lon(i,j)
        IF(PRESENT(SolarZenithAngle_k))        SolarZenithAngle_k = this%SolZenAng(i,j)
        IF(PRESENT(SolarAzimuthAngle_k))       SolarAzimuthAngle_k = this%SolAziAng(i,j)
        IF(PRESENT(ViewingZenithAngle_k))      ViewingZenithAngle_k = this%ViewZenAng(i,j)
        IF(PRESENT(ViewingAzimuthAngle_k))     ViewingAzimuthAngle_k = this%ViewAziAng(i,j)
        IF(PRESENT(TerrainHeight_k))           TerrainHeight_k = this%TerrainHeight(i,j)
        IF(PRESENT(GroundPixelQualityFlags_k)) GroundPixelQualityFlags_k = this%GPQFlag(i,j)
        IF(PRESENT(XTrackQualityFlags_k))      XTrackQualityFlags_k = this%XTQFlag(i,j)

        RETURN
      END FUNCTION L1Br_getGEOpix
                              
!! 5. L1Br_getGEOline
 !
 !    Functionality:
 !
 !       This function gets one or more geolocation data field values from
 !       the data block.  This function returns a single line of values
 !       for each given keyword (i.e., for fields with only a single value
 !       per line of the L1B file (e.g., Time), a single value is returned, 
 !       but for fields with multiple values per line (e.g., Latitude), all values
 !       associated with that line are returned in an array; the caller has
 !       the responsibility of allocating sufficient space to hold the returned
 !       array).
 !
 !    Calling Arguments:
 !
 !    Inputs:
 !
 !       this: the block data structure
 !       iLine: the line number in the L1B swath.  NOTE: this input is 0 based
 !              range from 0 to (NumTimes-1) inclusive.
 !
 !    Outputs:
 !
 !       The following are keyword arguments.  Only those present in the argument 
 !       list will be set by the function.
 !
 !       Keyword                    Corresponding L1B Geolocation Data Field
 !       =======                    ========================================
 !       Time_k                     Time
 !       SecondsInDay_k             SecondsInDay
 !       SpacecraftLatitude_k       SpacecraftLatitude
 !       SpacecraftLongitude_k      SpacecraftLongitude
 !       SpacecraftAltitude_k       SpacecraftAltitude
 !       SolarElevation_k           SolarElevation
 !       SolarAzimuth_k             SolarAzimuth
 !       SolarElevationMinimum_k    SolarElevationMinimum
 !       SolarElevationMaximum_k    SolarElevationMaximum
 !       SolarAzimuthMinimum_k      SolarAzimuthMinimum
 !       SolarAzimuthMaximum_k      SolarAzimuthMaximum
 !       Latitude_k                 Latitude
 !       Longitude_k                Longitude
 !       SolarZenithAngle_k         SolarZenithAngle
 !       SolarAzimuthAngle_k        SolarAzimuthAngle
 !       ViewingZenithAngle_k       ViewingZenithAngle
 !       ViewingAzimuthAngle_k      ViewingAzimuthAngle
 !       TerrainHeight_k            TerrainHeight
 !       GroundPixelQualityFlags_k  GroundPixelQualityFlags
 !       XTrackQualityFlags_k       XTrackQualityFlags 
 !
 !       status: the return PGS_SMF status value
 !
 !    Change History:
 !
 !       Date            Author          Modifications
 !       ====            ======          =============
 !       January 2005    Jeremy Warner   Original Source ((derived from original 
 !                                             L1B Reader code written by Kai Yang)
 !
!!
      FUNCTION L1Br_getGEOline(this, iLine, Time_k, SecondsInDay_k, &
                SpacecraftLatitude_k, SpacecraftLongitude_k, &
			    SpacecraftAltitude_k, &
			    SolarElevation_k, SolarAzimuth_k, &
			    SolarElevationMinimum_k, SolarElevationMaximum_k, &
			    SolarAzimuthMinimum_k, SolarAzimuthMaximum_k, &
			    Latitude_k, Longitude_k, SolarZenithAngle_k, &
                SolarAzimuthAngle_k, ViewingZenithAngle_k, &
			    ViewingAzimuthAngle_k, TerrainHeight_k, &
                GroundPixelQualityFlags_k, XTrackQualityFlags_k ) RESULT (status)
                !tpk TrueZoom_k ) RESULT (status)   ! tpk addition

        USE Szoom_Parameter_Module
        INCLUDE 'PGS_OMI_1900.f'  !defines the L1B PGS error codes
        TYPE (L1B_block_type), INTENT(INOUT) :: this
        INTEGER, PARAMETER :: i1 = 1
        INTEGER (KIND=i1), PARAMETER :: global_mode = 8_i1, szoom_mode = 4_i1
        INTEGER (KIND = 4), INTENT(IN) :: iLine
        !tpk INTEGER (KIND = 1), INTENT(IN), OPTIONAL :: TrueZoom_k    !tpk
        REAL (KIND = 8), OPTIONAL, INTENT(OUT) :: Time_k
        REAL (KIND = 4), OPTIONAL, INTENT(OUT) :: &
             SpacecraftLatitude_k, SpacecraftLongitude_k, &
             SpacecraftAltitude_k, SolarElevation_k, &
             SolarAzimuth_k, SolarElevationMinimum_k, &
             SolarElevationMaximum_k, SolarAzimuthMinimum_k, &
             SolarAzimuthMaximum_k, SecondsInDay_k
        REAL (KIND = 4), OPTIONAL, DIMENSION(:), INTENT(OUT) :: &
             Latitude_k, Longitude_k, SolarZenithAngle_k, &
             SolarAzimuthAngle_k, ViewingZenithAngle_k, &
             ViewingAzimuthAngle_k
        REAL (KIND = 4), DIMENSION(this%nXtrack) :: &
             tmp_Latitude_k, tmp_Longitude_k, tmp_SolarZenithAngle_k, &
             tmp_SolarAzimuthAngle_k, tmp_ViewingZenithAngle_k, &
             tmp_ViewingAzimuthAngle_k
        INTEGER (KIND = 2), OPTIONAL, DIMENSION(:), INTENT(OUT) :: TerrainHeight_k
        INTEGER (KIND = 2), OPTIONAL, DIMENSION(:), INTENT(OUT) :: GroundPixelQualityFlags_k
        INTEGER (KIND = 1), OPTIONAL, DIMENSION(:), INTENT(OUT) :: XTrackQualityFlags_k
        INTEGER (KIND = 2),  DIMENSION(this%nXtrack) :: tmp_TerrainHeight_k
        INTEGER (KIND = 2),  DIMENSION(this%nXtrack) :: tmp_GroundPixelQualityFlags_k
        INTEGER (KIND = 1),  DIMENSION(this%nXtrack) :: tmp_XTrackQualityFlags_k

        INTEGER :: i
        INTEGER :: itmp
        INTEGER (KIND = 4) :: status

        !tpk LOGICAL :: yn_truezoom   !tpk


        status = check_blk_line(this, iLine, i) 
        IF(status .NE. OMI_S_SUCCESS) THEN
           ierr = OMI_SMF_setmsg(OMI_E_GENERAL, "Failed retrieving geolocation data.", &
                                 "L1Br_getGEOline", one)
           RETURN
        ENDIF

        ! These first ones have only 1 value per line

        IF(PRESENT(Time_k))                    Time_k = this%Time(i)
        IF(PRESENT(SecondsInDay_k))            SecondsInDay_k = this%SecInDay(i)
        IF(PRESENT(SpacecraftLatitude_k))      SpacecraftLatitude_k = this%ScLat(i)
        IF(PRESENT(SpacecraftLongitude_k))     SpacecraftLongitude_k = this%ScLon(i)
        IF(PRESENT(SpacecraftAltitude_k))      SpacecraftAltitude_k = this%ScAlt(i)
        IF(PRESENT(SolarElevation_k))          SolarElevation_k = this%SolElevation(i)
        IF(PRESENT(SolarAzimuth_k))            SolarAzimuth_k = this%SolAzimuth(i)
        IF(PRESENT(SolarElevationMinimum_k))   SolarElevationMinimum_k = this%SolElevMin(i)
        IF(PRESENT(SolarElevationMaximum_k))   SolarElevationMaximum_k = this%SolElevMax(i)
        IF(PRESENT(SolarAzimuthMinimum_k))     SolarAzimuthMinimum_k = this%SolAziMin(i)
        IF(PRESENT(SolarAzimuthMaximum_k))     SolarAzimuthMaximum_k = this%SolAziMax(i)


        !tpk !tpk -----------------------------------------------------------------
        !tpk !tpk Check for TrueZoom flag (60 cross-track positions instead of 30)
        !tpk !tpk -----------------------------------------------------------------
        !tpk yn_truezoom = .FALSE.
        !tpk IF ( PRESENT (TrueZoom_k) ) THEN
        !tpk    IF (TrueZoom_k > 0_i1 ) yn_truezoom = .TRUE.
        !tpk END IF

        ! These next ones have multiple values per line

        IF(PRESENT(Latitude_k)) THEN
           IF(SIZE(Latitude_k) < this%nXtrack) THEN
              ierr = OMI_SMF_setmsg(OMI_E_INPUT, &
                                    "input latitude array too small", &
                                    "L1Br_getGEOline", zero)
              status = OMI_E_FAILURE
              RETURN
           ENDIF
           Latitude_k = this%Lat(1:this%nXtrack, i)
             IF (this%ImBinFact(i) == szoom_mode) THEN
                IF (this%SpaZoomIndex == 0) THEN
             !tpk IF (this%ImBinFact(i) == szoom_mode .AND. (.NOT. yn_truezoom) ) THEN  !tpk
                  DO itmp = 1, this%nXtrack
                     tmp_Latitude_k(itmp) = Latitude_k(itmp)
                  ENDDO
                  Latitude_k(gfpix:zspix) = rfill 
                  Latitude_k(zspix+1:zspix+zlpix) =  tmp_Latitude_k(zfpix:zlpix)
                  Latitude_k(zlpix+zspix+1:glpix) = rfill
                ENDIF
             ENDIF
        ENDIF
     
        IF(PRESENT(Longitude_k)) THEN
           IF(SIZE(Longitude_k) < this%nXtrack) THEN
              ierr = OMI_SMF_setmsg(OMI_E_INPUT, &
                                    "input longitude array too small", &
                                    "L1Br_getGEOline", zero)
              status = OMI_E_FAILURE
              RETURN
           ENDIF
           Longitude_k = this%Lon(1:this%nXtrack, i)

             IF (this%ImBinFact(i) == szoom_mode) THEN
                IF (this%SpaZoomIndex == 0) THEN
             !tpk IF (this%ImBinFact(i) == szoom_mode .AND. (.NOT. yn_truezoom) ) THEN !tpk
                  DO itmp = 1, this%nXtrack
                     tmp_Longitude_k(itmp) = Longitude_k(itmp)
                  ENDDO
                  Longitude_k(gfpix:zspix) =  rfill
                  Longitude_k(zspix+1:zspix+zlpix) =  tmp_Longitude_k(zfpix:zlpix)
                  Longitude_k(zlpix+zspix+1:glpix) =  rfill
                ENDIF
             ENDIF
        ENDIF
            
        IF(PRESENT(SolarZenithAngle_k)) THEN
           IF(SIZE(SolarZenithAngle_k) < this%nXtrack) THEN
              ierr = OMI_SMF_setmsg(OMI_E_INPUT, &
                                    "input solarZenith array too small", &
                                    "L1Br_getGEOline", zero)
              status = OMI_E_FAILURE
              RETURN
           ENDIF
           SolarZenithAngle_k = this%SolZenAng(1:this%nXtrack, i)

             IF (this%ImBinFact(i) == szoom_mode) THEN
                IF (this%SpaZoomIndex == 0) THEN
             !tpk IF (this%ImBinFact(i) == szoom_mode .AND. (.NOT. yn_truezoom) ) THEN  !tpk
                  DO itmp = 1, this%nXtrack
                     tmp_SolarZenithAngle_k(itmp) = SolarZenithAngle_k(itmp)
                  ENDDO
                  SolarZenithAngle_k(gfpix:zspix) = rfill
                  SolarZenithAngle_k(zspix+1:zspix+zlpix) =  tmp_SolarZenithAngle_k(zfpix:zlpix)
                  SolarZenithAngle_k(zlpix+zspix+1:glpix) =  rfill
                ENDIF
             ENDIF
        ENDIF

        IF(PRESENT(SolarAzimuthAngle_k)) THEN
           IF(SIZE(SolarAzimuthAngle_k) <  this%nXtrack) THEN
              ierr = OMI_SMF_setmsg(OMI_E_INPUT, &
                                    "input solarAzimuth array too small", &
                                    "L1Br_getGEOline", zero)
              status = OMI_E_FAILURE
              RETURN
           ENDIF
           SolarAzimuthAngle_k= this%SolAziAng(1:this%nXtrack, i)

             IF (this%ImBinFact(i) == szoom_mode) THEN
                IF (this%SpaZoomIndex == 0) THEN
             !tpk IF (this%ImBinFact(i) == szoom_mode .AND. (.NOT. yn_truezoom) ) THEN  !tpk
                   DO itmp = 1, this%nXtrack
                      tmp_SolarAzimuthAngle_k(itmp) = SolarAzimuthAngle_k(itmp)
                   ENDDO
                   SolarAzimuthAngle_k(gfpix:zspix) = rfill
                   SolarAzimuthAngle_k(zspix+1:zspix+zlpix) =  tmp_SolarAzimuthAngle_k(zfpix:zlpix)
                   SolarAzimuthAngle_k(zlpix+zspix+1:glpix) = rfill 
                ENDIF
             ENDIF

        ENDIF

        IF(PRESENT(ViewingZenithAngle_k)) THEN
           IF(SIZE(ViewingZenithAngle_k) < this%nXtrack) THEN
              ierr = OMI_SMF_setmsg(OMI_E_INPUT, &
                                    "input viewZenith array too small", &
                                    "L1Br_getGEOline", zero)
              status = OMI_E_FAILURE
              RETURN
           ENDIF
           ViewingZenithAngle_k  = this%ViewZenAng(1:this%nXtrack, i)

             IF (this%ImBinFact(i) == szoom_mode) THEN
                IF (this%SpaZoomIndex == 0) THEN
             !tpkIF (this%ImBinFact(i) == szoom_mode .AND. (.NOT. yn_truezoom) ) THEN  !tpk
                   DO itmp = 1, this%nXtrack
                      tmp_ViewingZenithAngle_k(itmp) = ViewingZenithAngle_k(itmp)
                   ENDDO
                   ViewingZenithAngle_k(gfpix:zspix) = rfill 
                   ViewingZenithAngle_k(zspix+1:zspix+zlpix) =  tmp_ViewingZenithAngle_k(zfpix:zlpix)
                   ViewingZenithAngle_k(zlpix+zspix+1:glpix) = rfill
                ENDIF
             ENDIF
        ENDIF

        IF(PRESENT(ViewingAzimuthAngle_k)) THEN
           IF(SIZE(ViewingAzimuthAngle_k) < this%nXtrack) THEN
              ierr = OMI_SMF_setmsg(OMI_E_INPUT, &
                                    "input viewAzimuth array too small", &
                                    "L1Br_getGEOline", zero)
              status = OMI_E_FAILURE
              RETURN
           ENDIF
           ViewingAzimuthAngle_k = this%ViewAziAng(1:this%nXtrack, i)
           IF (this%ImBinFact(i) == szoom_mode) THEN
              IF (this%SpaZoomIndex == 0) THEN
           !tpk IF (this%ImBinFact(i) == szoom_mode .AND. (.NOT. yn_truezoom) ) THEN !tpk 
                 DO itmp = 1, this%nXtrack
                    tmp_ViewingAzimuthAngle_k(itmp) = ViewingAzimuthAngle_k(itmp)
                 ENDDO
                 ViewingAzimuthAngle_k(gfpix:zspix) = rfill
                 ViewingAzimuthAngle_k(zspix+1:zspix+zlpix) =  tmp_ViewingAzimuthAngle_k(zfpix:zlpix)
                 ViewingAzimuthAngle_k(zlpix+zspix+1:glpix) = rfill
              ENDIF
           ENDIF
        ENDIF

        IF(PRESENT(TerrainHeight_k)) THEN
           IF(SIZE(TerrainHeight_k) < this%nXtrack) THEN
              ierr = OMI_SMF_setmsg(OMI_E_INPUT, &
                                    "input height array too small", &
                                    "L1Br_getGEOline", zero)
              status = OMI_E_FAILURE
              RETURN
           ENDIF
           TerrainHeight_k  = this%TerrainHeight(1:this%nXtrack, i)
           IF (this%ImBinFact(i) == szoom_mode) THEN
              IF (this%SpaZoomIndex == 0) THEN
          !tpk IF (this%ImBinFact(i) == szoom_mode .AND. (.NOT. yn_truezoom) ) THEN  !tpk
                 DO itmp = 1, this%nXtrack
                    tmp_TerrainHeight_k(itmp) = TerrainHeight_k(itmp)
                 ENDDO
                 TerrainHeight_k(gfpix:zspix) = int16fill
                 TerrainHeight_k(zspix+1:zspix+zlpix) =  tmp_TerrainHeight_k(zfpix:zlpix)
                 TerrainHeight_k(zlpix+zspix+1:glpix) = int16fill
              ENDIF
             ENDIF
        ENDIF

        IF(PRESENT(GroundPixelQualityFlags_k)) THEN
           IF(SIZE(GroundPixelQualityFlags_k) < this%nXtrack ) THEN
              ierr = OMI_SMF_setmsg(OMI_E_INPUT, &
                                    "input geolocation flag array too small", &
                                    "L1Br_getGEOline", zero)
              status = OMI_E_FAILURE
              RETURN
           ENDIF
           GroundPixelQualityFlags_k = this%GPQFlag(1:this%nXtrack, i)
           IF (this%ImBinFact(i) == szoom_mode) THEN
              IF (this%SpaZoomIndex == 0) THEN
           !tpk IF (this%ImBinFact(i) == szoom_mode .AND. (.NOT. yn_truezoom) ) THEN  !tpk
                 DO itmp = 1, this%nXtrack
                    tmp_GroundPixelQualityFlags_k(itmp) = GroundPixelQualityFlags_k(itmp)
                 ENDDO
                 GroundPixelQualityFlags_k(gfpix:zspix) = int16fill 
                 GroundPixelQualityFlags_k(zspix+1:zspix+zlpix) =  tmp_GroundPixelQualityFlags_k(zfpix:zlpix)
                 GroundPixelQualityFlags_k(zlpix+zspix+1:glpix) = int16fill
              ENDIF
           ENDIF
        ENDIF
  
        IF(PRESENT(XTrackQualityFlags_k)) THEN
           IF(SIZE(XTrackQualityFlags_k) < this%nXtrack ) THEN
              ierr = OMI_SMF_setmsg(OMI_E_INPUT, &
                   "input geolocation flag array too small", &
                   "L1Br_getGEOline", zero)
              status = OMI_E_FAILURE
              RETURN
           ENDIF
           XTrackQualityFlags_k = this%XTQFlag(1:this%nXtrack, i)
           IF (this%ImBinFact(i) == szoom_mode) THEN
              IF (this%SpaZoomIndex == 0) THEN
                 DO itmp = 1, this%nXtrack
                    tmp_XTrackQualityFlags_k(itmp) = XTrackQualityFlags_k(itmp)
                 ENDDO
                 XTrackQualityFlags_k(gfpix:zspix) = int8fill
                 XTrackQualityFlags_k(zspix+1:zspix+zlpix) =  tmp_XTrackQualityFlags_k(zfpix:zlpix)
                 XTrackQualityFlags_k(zlpix+zspix+1:glpix) = int8fill
              ENDIF
           ENDIF
        ENDIF
      
        RETURN
      END FUNCTION L1Br_getGEOline

!! 6. L1Br_getDATA
 !
 !    Functionality:
 !
 !       This function gets one or more data field values from
 !       the data block.  This function returns a single line of values
 !       for each given keyword.  All of the data fields accessed by this 
 !       function have a single value per line, so only a single value is
 !       ever passed back for any given keyword.
 !
 !    Calling Arguments:
 !
 !    Inputs:
 !
 !       this: the block data structure
 !       iLine: the line number in the L1B swath.  NOTE: this input is 0 based
 !              range from 0 to (NumTimes-1) inclusive.
 !
 !    Outputs:
 !
 !       The following are keyword arguments.  Only those present in the argument 
 !       list will be set by the function.
 !
 !       Keyword                             Corresponding L1B Data Field
 !       =======                             ============================
 !       MeasurementClass_k                  MeasurementClass
 !       InstrumentConfigurationId_k         InstrumentConfigurationId
 !       MeasurementQualityFlags_k           MeasurementQualityFlags
 !       ExposureType_k                      ExposureType
 !       MasterClockPeriod_k                 MasterClockPeriod
 !       CalibrationSettings_k               CalibrationSettings
 !       ExposureTime_k                      ExposureTime
 !       ReadoutTime_k                       ReadoutTime
 !       GainSwitchingColumn1_k              GainSwitchingColumn1
 !       GainSwitchingColumn2_k              GainSwitchingColumn2
 !       GainSwitchingColumn3_k              GainSwitchingColumn3
 !       GainCode1_k                         GainCode1
 !       GainCode2_k                         GainCode2
 !       GainCode3_k                         GainCode3
 !       GainCode4_k                         GainCode4
 !       DSGainCode_k                        DSGainCode
 !       LowerStrayLightAreaBinningFactor_k  LowerStrayLightAreaBinningFactor
 !       UpperStrayLightAreaBinningFactor_k  UpperStrayLightAreaBinningFactor
 !       LowerDarkAreaBinningFactor_k        LowerDarkAreaBinningFactor
 !       UpperDarkAreaBinningFactor_k        UpperDarkAreaBinningFactor
 !       SkipRows1_k                         SkipRows1
 !       SkipRows2_k                         SkipRows2
 !       SkipRows3_k                         SkipRows3
 !       SkipRows4_k                         SkipRows4
 !       DetectorTemperature_k               DetectorTemperature
 !       OpticalBenchTemperature_k           OpticalBenchTemperature
 !       ImageBinningFactor_k                ImageBinningFactor
 !       BinnedImageRows_k                   BinnedImageRows
 !       StopColumn_k                        StopColumn
 !
 !       status: the return PGS_SMF status value
 !
 !    Change History:
 !
 !       Date            Author          Modifications
 !       ====            ======          =============
 !       January 2005    Jeremy Warner   Original Source
 !
!!
      FUNCTION L1Br_getDATA(this, iLine, MeasurementClass_k, InstrumentConfigurationId_k, &
           MeasurementQualityFlags_k, ExposureType_k, MasterClockPeriod_k, &
           CalibrationSettings_k, ExposureTime_k, ReadoutTime_k, &
           GainSwitchingColumn1_k, GainSwitchingColumn2_k, &
           GainSwitchingColumn3_k, GainCode1_k, GainCode2_k, GainCode3_k, &
           GainCode4_k, DSGainCode_k, LSLABF_k, USLABF_k, LDABF_k, &
           UDABF_k, SkipRows1_k, SkipRows2_k, &
           SkipRows3_k, SkipRows4_k, DetectorTemperature_k, &
           OpticalBenchTemperature_k, ImageBinningFactor_k, &
           BinnedImageRows_k, StopColumn_k) RESULT (status)

        INCLUDE 'PGS_OMI_1900.f'  !defines the L1B PGS error codes
        TYPE (L1B_block_type), INTENT(INOUT) :: this
        INTEGER (KIND = 4), INTENT(IN) :: iLine
        INTEGER (KIND = 1), OPTIONAL, INTENT(OUT) :: MeasurementClass_k, &
			      InstrumentConfigurationId_k, ExposureType_k, GainCode1_k, &
			      GainCode2_k, GainCode3_k, GainCode4_k, DSGainCode_k, &
			      LSLABF_k, USLABF_k, LDABF_k, UDABF_k, ImageBinningFactor_k
        INTEGER (KIND = 2), OPTIONAL, INTENT(OUT) :: MeasurementQualityFlags_k, &
                              CalibrationSettings_k, GainSwitchingColumn1_k, &
			      GainSwitchingColumn2_k, GainSwitchingColumn3_k, &
                              SkipRows1_k, SkipRows2_k, SkipRows3_k, SkipRows4_k, &
                              BinnedImageRows_k, StopColumn_k
        REAL (KIND =4), OPTIONAL, INTENT(OUT) :: MasterClockPeriod_k, ExposureTime_k, &
	                      ReadoutTime_k, DetectorTemperature_k, OpticalBenchTemperature_k
        INTEGER :: i
        INTEGER (KIND = 4) :: status

        status = check_blk_line(this, iLine, i) 
        IF(status .NE. OMI_S_SUCCESS) THEN
           ierr = OMI_SMF_setmsg(OMI_E_GENERAL, "Failed retrieving L1B data.", &
                                 "L1Br_getDATA", one)
           RETURN
        ENDIF

        IF(PRESENT(MeasurementClass_k)) MeasurementClass_k = this%MeasClass(i)
        IF(PRESENT(InstrumentConfigurationId_k)) InstrumentConfigurationId_k = this%Config(i)
        IF(PRESENT(MeasurementQualityFlags_k)) MeasurementQualityFlags_k = this%MQFlag(i)
        IF(PRESENT(ExposureType_k)) ExposureType_k = this%ExposureType(i)
        IF(PRESENT(MasterClockPeriod_k)) MasterClockPeriod_k = this%MasterClkPer(i)
        IF(PRESENT(CalibrationSettings_k)) CalibrationSettings_k = this%CalSet(i)
        IF(PRESENT(ExposureTime_k)) ExposureTime_k = this%ExpTime(i)
        IF(PRESENT(ReadoutTime_k)) ReadoutTime_k = this%ReadTime(i)
        IF(PRESENT(GainSwitchingColumn1_k)) GainSwitchingColumn1_k = this%GSC1(i)
        IF(PRESENT(GainSwitchingColumn2_k)) GainSwitchingColumn2_k = this%GSC2(i)
        IF(PRESENT(GainSwitchingColumn3_k)) GainSwitchingColumn3_k = this%GSC3(i)
        IF(PRESENT(GainCode1_k)) GainCode1_k = this%GC1(i)
        IF(PRESENT(GainCode2_k)) GainCode2_k = this%GC2(i)
        IF(PRESENT(GainCode3_k)) GainCode3_k = this%GC3(i)
        IF(PRESENT(GainCode4_k)) GainCode4_k = this%GC4(i)
        IF(PRESENT(DSGainCode_k)) DSGainCode_k = this%DSGC(i)
        IF(PRESENT(LSLABF_k)) LSLABF_k = this%LSLABF(i)
        IF(PRESENT(USLABF_k)) USLABF_k = this%USLABF(i)
        IF(PRESENT(LDABF_k)) LDABF_k = this%LDABF(i)
        IF(PRESENT(UDABF_k)) UDABF_k = this%UDABF(i)
        IF(PRESENT(SkipRows1_k)) SkipRows1_k = this%SR1(i)
        IF(PRESENT(SkipRows2_k)) SkipRows2_k = this%SR2(i)
        IF(PRESENT(SkipRows3_k)) SkipRows3_k = this%SR3(i)
        IF(PRESENT(SkipRows4_k)) SkipRows4_k = this%SR4(i)
        IF(PRESENT(DetectorTemperature_k)) DetectorTemperature_k = this%DetcTemp(i)
        IF(PRESENT(OpticalBenchTemperature_k)) OpticalBenchTemperature_k = this%OptBenchTemp(i)
        IF(PRESENT(ImageBinningFactor_k)) ImageBinningFactor_k = this%ImBinFact(i)
        IF(PRESENT(BinnedImageRows_k)) BinnedImageRows_k = this%BinImgRows(i)
        IF(PRESENT(StopColumn_k)) StopColumn_k = this%StopColumn(i)

        RETURN
      END FUNCTION L1Br_getDATA

!! Private function: calc_wl_pix
 !
 ! Functionality:
 !
 !   This function calculates the wavelengths values for a specific pixel
 !   and returns those values, plus the limits of the specified wavelength
 !   range.
 !
 ! Calling Arguments:
 !
 ! Inputs:
 !
 !    j         Index into L1B structure, line number
 !    i         Index into L1B structure, Xtrack position
 !    this      L1B block
 !    minwl     minimum wavelength requested
 !    maxwl     maximum wavelength requested
 !
 ! Outputs:
 !
 !    wl_local  wavelength values calculated from L1B coefficients
 !    il        index of lower bound of wavelength range (in wl_local)
 !    ih        index of upper bound of wavelength range (in wl_local)
 !    Nwl_l     number of wavelengths in range
 !
 !    status    the return PGS_SMF status value
 !
 ! Change History:
 !
 !    Date            Author          Modifications
 !    ====            ======          =============
 !    January 2005    Jeremy Warner   Original Source
 !
!!
      FUNCTION calc_wl_pix(j, i, this, minwl, maxwl, wl_local, il, ih, &
	                   Nwl_l) RESULT (status)
      INCLUDE 'PGS_OMI_1900.f'  !defines the L1B PGS error codes
      TYPE (L1B_block_type), INTENT(INOUT) :: this
      INTEGER (KIND = 4), INTENT(IN) :: i, j
      REAL (KIND = 4), INTENT(INOUT) :: minwl, maxwl
      INTEGER (KIND = 4), INTENT(OUT) :: il, ih, Nwl_l
      REAL (KIND = 4), DIMENSION(:), INTENT(OUT) :: wl_local

      INTEGER (KIND = 4) :: status
      INTEGER (KIND = 4) :: fflag, k, q, i_foo(1)
      INTEGER (KIND = 4) :: fac

      ! First calculate the wavelength values based on the wavelength
      ! coefficients from the L1B file.  Set fflag if we find any
      ! fill values.

      status = OMI_S_SUCCESS
      fflag = 0
      DO k = 1, this%nWavel
         wl_local(k) = 1.0
         DO q = 1, this%nWavelCoef
            IF(this%WavelenCoef(q,i,j) .eq. rfill) wl_local(k) = rfill
         ENDDO
         fac = k-1-this%WavelenRefCol(j)
         IF (wl_local(k) > 0.0) THEN
            wl_local(k) = fac*this%WavelenCoef(this%nWavelCoef,i,j)
            DO q = this%nWavelCoef -1, 2, -1
               wl_local(k) = fac*(wl_local(k) + this%WavelenCoef(q,i,j))
            ENDDO
            wl_local(k) = wl_local(k) + this%WavelenCoef(1,i,j)
         ELSE
            fflag = 1
         ENDIF
      ENDDO

      ! Now find the limits of the wavelength range.

      IF (fflag .EQ. 0) THEN  ! We found no fill values, so do it the fast way.

        ! Check to see if minwl and maxwl are consistent (i.e., minwl < maxwl)

        IF(minwl .gt. 0.0 .AND. maxwl .gt. 0.0) THEN
           IF((minwl - maxwl) > 0.0) THEN
              ierr = OMI_SMF_setmsg(OMI_E_INPUT, "input Wlmin_k greater than Wlmax_k", &
                                    "calc_wl_pix", zero)
              status = OMI_E_FAILURE
              RETURN
           ENDIF
        ENDIF

        ! Check to see if minwl is compatable with the wavelength range
        ! If minwl was not given, set the lower bound to the first element.

        IF(minwl .gt. 0.0) THEN
           IF(minwl > MAXVAL(wl_local))THEN
              ierr = OMI_SMF_setmsg(OMI_E_INPUT, "input Wlmin_k out of bound", &
                                    "calc_wl_pix", zero)
              status = OMI_E_FAILURE
              RETURN
           ENDIF
 
	   ! Set the lower bound based on wl_local and minwl

           IF(minwl < MINVAL(wl_local)) THEN
              il = 1
           ELSE
              i_foo = MAXLOC(wl_local, MASK = wl_local-minwl <= 0.0)
              il    = i_foo(1)
           ENDIF            
        ELSE
           il = 1
        ENDIF

        ! Check to see if maxwl is compatable with the wavelength range
        ! If maxwl was not given, set the upper bound to the last element.

        IF(maxwl .gt. 0.0) THEN
           IF(maxwl < MINVAL(wl_local)) THEN
              ierr = OMI_SMF_setmsg(OMI_E_INPUT, "input Wlmax_k out of bound", &
                                    "calc_wl_pix", zero)
              status = OMI_E_FAILURE
              RETURN
           ENDIF
 
	   ! Set the upper bound based on wl_local and maxwl

           IF(maxwl > MAXVAL(wl_local)) THEN
              ih = this%nWavel
           ELSE
              i_foo = MINLOC(wl_local, MASK = wl_local-maxwl >= 0.0)
              ih    = i_foo(1)
           ENDIF       
        ELSE
           ih = this%nWavel
        ENDIF

      ELSE  ! => we had a fill value

        ! Initialize il and ih

        ih = 0
        il = 0

        ! Loop through all wavelengths, and defind bounds based on minwl and maxwl

	if (minwl < 0.0) minwl = 0.0
	if (maxwl < 0.0) maxwl = 1000.0
        DO k = 1, this%nWavel
           IF (wl_local(k) >= minwl .AND. wl_local(k) <= maxwl) THEN
              IF (il .EQ. 0) il = k
              IF (ih .EQ. 0) ih = k
              IF (ih .LT. k) ih = k
           ENDIF
        ENDDO
        IF (il > 1) il = il - 1
        IF (ih < this%nWavel .AND. ih .NE. 0) ih = ih + 1
      ENDIF

      ! Set the number of wavelengths in the span to Nwl_l

      IF (il .EQ. 0 .AND. ih .EQ. 0) THEN
         Nwl_l = 0
      ELSE
         Nwl_l = ih-il+1
      ENDIF

      RETURN
      END FUNCTION calc_wl_pix

!! Private function: calc_wl_line
 !
 ! Functionality:
 !
 !   This function calculates the wavelengths values for a specific line
 !   and returns those values, plus the limits of the specified wavelength
 !   range.
 !
 ! Calling Arguments:
 !
 ! Inputs:
 !
 !    i         Index into L1B structure, line number
 !    this      L1B block
 !    minwl     minimum wavelength requested
 !    maxwl     maximum wavelength requested
 !
 ! Outputs:
 !
 !    wl_local  wavelength values calculated from L1B coefficients
 !    il        index of lower bound of wavelength range (in wl_local)
 !    ih        index of upper bound of wavelength range (in wl_local)
 !    Nwl_l     number of wavelengths in range
 !
 !    status    the return PGS_SMF status value
 !
 ! Change History:
 !
 !    Date            Author          Modifications
 !    ====            ======          =============
 !    January 2005    Jeremy Warner   Original Source
 !
!!
      FUNCTION calc_wl_line(j, this, minwl, maxwl, wl_local, il, ih, &
	                   Nwl_l) RESULT (status)
      INCLUDE 'PGS_OMI_1900.f'  !defines the L1B PGS error codes
      TYPE (L1B_block_type), INTENT(INOUT) :: this
      INTEGER (KIND = 4), INTENT(IN) :: j
      REAL (KIND = 4), INTENT(INOUT) :: minwl, maxwl
      INTEGER (KIND = 4), INTENT(OUT) :: il, ih, Nwl_l
      REAL (KIND = 4), DIMENSION(:,:), INTENT(OUT) :: wl_local

      INTEGER (KIND = 4) :: status
      INTEGER (KIND = 4) :: fflag, i, k, q, i_foo(1)
      INTEGER (KIND = 4) :: fac


      ! First calculate the wavelength values based on the wavelength
      ! coefficients from the L1B file.  Set fflag if we find any
      ! fill values.

      status = OMI_S_SUCCESS
      fflag = 0
      DO i = 1, this%nXtrack
        DO k = 1, this%nWavel
           wl_local(k,i) = 1.0
           DO q = 1, this%nWavelCoef
              IF(this%WavelenCoef(q,i,j) .EQ. rfill) wl_local(k,i) = rfill
           ENDDO
         fac = k-1-this%WavelenRefCol(j)
         IF (wl_local(k,i) > 0.0) THEN
            wl_local(k,i) = fac*this%WavelenCoef(this%nWavelCoef,i,j)
            DO q = this%nWavelCoef -1, 2, -1
               wl_local(k,i) = fac*(wl_local(k,i) + this%WavelenCoef(q,i,j))
            ENDDO
            wl_local(k,i) = wl_local(k,i) + this%WavelenCoef(1,i,j)
           ELSE
              fflag = 1
           ENDIF
        ENDDO
      ENDDO

      ! Now find the limits of the wavelength range.

      IF (fflag .EQ. 0) THEN  ! We found no fill values, so do it the fast way.

        ! Check to see if minwl and maxwl are consistent (i.e., minwl < maxwl)

        IF(minwl .gt. 0.0 .AND. maxwl .gt. 0.0) THEN
           IF((minwl-maxwl) > 0.0) THEN
              ierr = OMI_SMF_setmsg(OMI_E_INPUT, &
                                    "input Wlmin_k greater than Wlmax_k", &
                                    "calc_wl_line", zero)
              status = OMI_E_FAILURE
              RETURN
           ENDIF
        ENDIF

        ! Check to see if minwl is compatable with the wavelength range
        ! If minwl was not given, set the lower bound to the first element.

        IF(minwl .gt. 0.0) THEN
           IF(minwl > MAXVAL(wl_local))THEN
              ierr = OMI_SMF_setmsg(OMI_E_INPUT, &
                                    "input Wlmin_k out of bound", &
                                    "calc_wl_line", zero)
              status = OMI_E_FAILURE
              RETURN
           ENDIF
 
	   ! Set the lower bound based on wl_local and minwl

           IF(minwl < MINVAL(wl_local)) THEN
              il = 1
           ELSE
              il = 0
              DO i = 1, this%nXtrack
                 DO k = 1, this%nWavel
                    if (wl_local(k,i) > minwl) THEN
                       if (il .eq. 0 .or. il .gt. k) il = k
                       EXIT
                    ENDIF
                 ENDDO
              ENDDO
              IF (il > 1) THEN 
                 il = il - 1
              ELSE 
                 il = 1
              ENDIF
           ENDIF            
        ELSE
           il = 1
        ENDIF

        ! Check to see if maxwl is compatable with the wavelength range
        ! If maxwl was not given, set the upper bound to the last element.

        IF(maxwl .gt. 0.0) THEN
           IF(maxwl < MINVAL(wl_local)) THEN
              ierr = OMI_SMF_setmsg(OMI_E_INPUT, &
                                    "input Wlmax_k out of bound", &
                                    "calc_wl_line", zero)
              status = OMI_E_FAILURE
              RETURN
           ENDIF

	   ! Set the upper bound based on wl_local and maxwl

           IF(maxwl > MAXVAL(wl_local)) THEN
              ih = this%nWavel
           ELSE
              ih = this%nWavel+1
              DO i = 1, this%nXtrack
                 DO k = 1, this%nWavel
                    if (wl_local(k,i) > maxwl) THEN
                       if (ih .eq. this%nWavel+1 .or. ih .lt. k) ih = k
                       EXIT
                    ENDIF
                 ENDDO
              ENDDO
              IF (ih < this%nWavel) THEN 
                 ih = ih + 1
              ELSE 
                 ih = this%nWavel
              ENDIF
           ENDIF       
        ELSE
           ih = this%nWavel
        ENDIF

      ELSE  ! => we had a fill value

        ! Initialize il and ih

        ih = 0
        il = 0

        ! Loop through all wavelengths, and defind bounds based on minwl and maxwl

	if (minwl < 0.0) minwl = 0.0
	if (maxwl < 0.0) maxwl = 1000.0
        DO k = 1, this%nWavel
           DO i = 1, this%nXtrack
              IF (wl_local(k,i) >= minwl .AND. wl_local(k,i) <= maxwl) THEN
                 IF (il .EQ. 0) il = k
                 IF (ih .EQ. 0) ih = k
                 IF (ih .LT. k) ih = k
              ENDIF
           ENDDO
        ENDDO
        IF (il .EQ. 0 .AND. ih .EQ. 0) THEN
        ELSE
           IF (il > 1) il = il - 1
           IF (ih < this%nWavel) ih = ih + 1
        ENDIF
      ENDIF

      ! Set the number of wavelengths in the span to Nwl_l

      IF (il .EQ. 0 .AND. ih .EQ. 0) THEN
         Nwl_l = 0
      ELSE
         Nwl_l = ih-il+1
      ENDIF
      RETURN
      END FUNCTION calc_wl_line

!! 7. L1Br_getSIGpix
 !
 !    Functionality:
 !
 !       This function gets one or more radiance data field values from
 !       the data block.  This function returns values for each given output
 !       keyword for the each wavelength within the input wavelength range.
 !       If no wavelength range is input, values are returned for all
 !       wavelengths.
 !
 !    Calling Arguments:
 !
 !    Inputs:
 !
 !       this: the block data structure
 !       iLine: the line number in the L1B swath.  NOTE: this input is 0 based
 !              range from 0 to (NumTimes-1) inclusive.
 !       iPix: the pixel number in the L1B swath.  NOTE: this input is 0 based
 !             range from 0 to (nXtrack-1) inclusive.
 !       Wlmin_k: input wavelength value for the band beginning (in nm) (Optional)
 !       Wlmax_k: input wavelength value for the band end (in nm) (Optional)
 !
 !       Note that Wlmin_k < Wlmax_k.
 !       If Wlmin_k is omitted, all wavelengths < Wlmax_k are returned.
 !       If Wlmax_k is omitted, all wavelengths > Wlmin_k are returned.
 !       If Wlmin_k and Wlmax_k are omitted, all wavelengths are returned.
 !
 !    Outputs:
 !
 !       The following are keyword arguments.  Only those present in the argument 
 !       list will be set by the function.  The first set of keywords is for
 !       combined data, whereas the second set is for raw data; the two sets are
 !       not mutually exclusive, i.e., raw and combined data can be extracted
 !       in the same function call.
 !
 !       Keyword                    Corresponding L1B Data
 !       =======                    ======================
 !       Signal_k                   Radiance or Irradiance or Signal
 !                                     Mantissa*10^(Exponent)
 !       SignalPrecision_k          RadiancePrecision or IrradiancePrecision or
 !                                     SignalPrecision
 !                                     Precision*10^(Exponent)
 !       PixelQualityFlags_k        PixelQualityFlags
 !       Wavelength_k               Wavelengths
 !                                     WavelengthCoefficient aggregated into Wavelengths
 !       Nwl_k                      None - just the number of wavelengths returned
 !
 !       If the user wants the raw L1B data, rather than combined data, that can be
 !       retrieved using the following keywords.
 !
 !       Keyword                     Corresponding L1B Data Field
 !       =======                     ============================
 !       Mantissa_k                  RadianceMantissa or IrradianceMantissa or SignalMantissa
 !       Precision_k                 RadiancePrecisionMantissa or IrradiancePrecisionMantissa
 !                                      or SignalPrecisionMantissa
 !       Exponent_k                  RadianceExponent or IrradianceExponent or SignalExponent
 !       WavelengthCoefficient_k     WavelengthCoefficient
 !       WavelengthCoefficientPrec_k WavelengthCoefficientPrecision
 !       WavelengthReferenceColumn_k WavelengthReferenceColumn
 !       SmallPixelSignal_k          SmallPixelRadiance or SmallPixelIrradiance or
 !                                      SmallPixelSignal
 !       SmallPixelWavelength_k      SmallPixelWavelength
 !       NumberSmallPixelColumns_k   NumberSmallPixelColumns
 !       SmallPixelColumn_k          SmallPixelColumn
 !
 !       status: the return PGS_SMF status value
 !
 !    Change History:
 !
 !       Date            Author          Modifications
 !       ====            ======          =============
 !       January 2005    Jeremy Warner   Original Source (derived from original 
 !                                             L1B Reader code written by Kai Yang)
 !
!!
      FUNCTION L1Br_getSIGpix(this, iLine, iPix, Wlmin_k, Wlmax_k, &
                              Signal_k, SignalPrecision_k, &
                              PixelQualityFlags_k, Wavelength_k, &
                              Nwl_k, Mantissa_k, Precision_k, Exponent_k, &
                              WavelengthCoefficient_k, WavelengthCoefficientPrec_k, &
                              WavelengthReferenceColumn_k, SmallPixelSignal_k, &
                              SmallPixelWavelength_k, NumberSmallPixelColumns_k, &
	                      SmallPixelColumn_k) RESULT (status)
        INCLUDE 'PGS_OMI_1900.f'  !defines the L1B PGS error codes

        TYPE (L1B_block_type), INTENT(INOUT) :: this
        INTEGER (KIND = 4), INTENT(IN) :: iLine, iPix
        REAL (KIND = 4), OPTIONAL, INTENT(IN) :: Wlmin_k, Wlmax_k    ! wavelength range
        INTEGER (KIND = 4), OPTIONAL, INTENT(OUT) :: Nwl_k
        INTEGER (KIND = 2), OPTIONAL, INTENT(OUT) :: WavelengthReferenceColumn_k
        INTEGER (KIND = 2), OPTIONAL, INTENT(OUT) :: SmallPixelColumn_k
        INTEGER (KIND = 1), OPTIONAL, INTENT(OUT) :: NumberSmallPixelColumns_k
        REAL (KIND = 4), OPTIONAL, DIMENSION(:), INTENT(OUT) :: Signal_k, &
                                             SignalPrecision_k, Wavelength_k, &
					     WavelengthCoefficient_k, &
					     WavelengthCoefficientPrec_k, &
					     SmallPixelSignal_k, SmallPixelWavelength_k 
        INTEGER (KIND = 2), OPTIONAL, DIMENSION(:), INTENT(OUT) :: &
                                             PixelQualityFlags_k, Mantissa_k, Precision_k
        INTEGER (KIND = 1), OPTIONAL, DIMENSION(:), INTENT(OUT) :: Exponent_k

        REAL (KIND = 4), DIMENSION(1:this%nWavel) :: wl_local
        INTEGER (KIND = 4) :: Nwl_l
        INTEGER :: i, j, k, q, l
        INTEGER (KIND = 4) :: il, ih, i_foo(1)
        INTEGER (KIND = 4) :: status
        REAL (KIND = 4) :: minWl, maxWl

      ! Check requested indecies for consistancy

      status = check_blk_pix(this, iLine, iPix, j, i) 
      IF(status .NE. OMI_S_SUCCESS) THEN
         ierr = OMI_SMF_setmsg(OMI_E_GENERAL, "Failed retrieving L1B signal data.", &
                                 "L1Br_getSIGpix", one)
         RETURN
      ENDIF

      ! Set min and max wavelengths

      IF(PRESENT(Wlmin_k)) THEN
         minWl = Wlmin_k
      ELSE
         minWl = -1.0
      ENDIF

      IF(PRESENT(Wlmax_k)) THEN
         maxWl = Wlmax_k
      ELSE
         maxWl = -1.0
      ENDIF

      ! Calculate wavelengths from L1B data

      status = calc_wl_pix(j, i, this, minWl, maxWl, wl_local, il, ih, Nwl_l)
      IF(status .NE. OMI_S_SUCCESS) THEN
         ierr = OMI_SMF_setmsg(OMI_E_GENERAL, "Failed calculating L1B wavelengths.", &
                                 "L1Br_getSIGpix", one)
         RETURN
      ENDIF

      ! Retrieve requested data; first do data for possible multiple wavelengths

      IF (Nwl_l > 0) THEN
        IF(PRESENT(Signal_k)) THEN
           IF(SIZE(Signal_k) < Nwl_l) THEN
              ierr = OMI_SMF_setmsg(OMI_E_INPUT, "input Signal_k array too small", &
                                    "L1Br_getSIGpix", zero)
              status = OMI_E_FAILURE
              RETURN
           ENDIF
           Signal_k(1:Nwl_l) = this%RadMantissa(il:ih,i,j) &
                               * 10.0**this%RadExponent(il:ih,i,j)
           DO q = il, ih
              IF(this%RadMantissa(q,i,j) .eq. -32767) Signal_k(q-il+1) = rfill
           ENDDO
           
        ENDIF

        IF(PRESENT(SignalPrecision_k)) THEN
           IF(SIZE(SignalPrecision_k) < Nwl_l) THEN
              ierr = OMI_SMF_setmsg(OMI_E_INPUT, "input SignalPrecision_k array too small",&
                                     "L1Br_getSIGpix", zero)
              status = OMI_E_FAILURE
              RETURN
           ENDIF
           SignalPrecision_k(1:Nwl_l) = this%RadPrecision(il:ih,i,j) &
                                        * 10.0**this%RadExponent(il:ih,i,j)
           DO q = il, ih
              IF(this%RadPrecision(q,i,j).eq.-32767) SignalPrecision_k(q-il+1) =rfill
           ENDDO
        ENDIF

        IF(PRESENT(PixelQualityFlags_k)) THEN
           IF(SIZE(PixelQualityFlags_k) < Nwl_l) THEN
              ierr = OMI_SMF_setmsg(OMI_E_INPUT, "input PixelQualityFlags_k array too small",&
                                    "L1Br_getSIGpix", zero)
              status = OMI_E_FAILURE
              RETURN
           ENDIF
           PixelQualityFlags_k(1:Nwl_l) = this%PQFlag(il:ih,i,j)
        ENDIF

        IF(PRESENT(Wavelength_k)) THEN
           IF(SIZE(Wavelength_k) < Nwl_l) THEN 
              ierr = OMI_SMF_setmsg(OMI_E_INPUT, "input Wavelength_k array too small",&
                                    "L1Br_getSIGpix", zero)
              status = OMI_E_FAILURE
              RETURN
           ENDIF
           Wavelength_k(1:Nwl_l) = wl_local(il:ih)
        ENDIF

	! Here is the raw data extraction

        IF(PRESENT(Mantissa_k)) THEN
           IF(SIZE(Mantissa_k) < Nwl_l) THEN 
              ierr = OMI_SMF_setmsg(OMI_E_INPUT, "input Mantissa_k array too small",&
                                    "L1Br_getSIGpix", zero)
              status = OMI_E_FAILURE
              RETURN
           ENDIF
           Mantissa_k(1:Nwl_l) = this%RadMantissa(il:ih,i,j)
        ENDIF

        IF(PRESENT(Precision_k)) THEN
           IF(SIZE(Precision_k) < Nwl_l) THEN 
              ierr = OMI_SMF_setmsg(OMI_E_INPUT, "input Precision_k array too small",&
                                    "L1Br_getSIGpix", zero)
              status = OMI_E_FAILURE
              RETURN
           ENDIF
           Precision_k(1:Nwl_l) = this%RadPrecision(il:ih,i,j)
        ENDIF

        IF(PRESENT(Exponent_k)) THEN
           IF(SIZE(Exponent_k) < Nwl_l) THEN 
              ierr = OMI_SMF_setmsg(OMI_E_INPUT, "input Exponent_k array too small",&
                                    "L1Br_getSIGpix", zero)
              status = OMI_E_FAILURE
              RETURN
           ENDIF
           Exponent_k(1:Nwl_l) = this%RadExponent(il:ih,i,j)
        ENDIF

      ENDIF  ! Nwl_l > 0

      IF(PRESENT(Nwl_k)) Nwl_k = Nwl_l

      ! The rest of the outputs do not depend on the value of Nwl_l

      IF(PRESENT(WavelengthCoefficient_k)) THEN
         IF(SIZE(WavelengthCoefficient_k) < this%nWavelCoef) THEN
            ierr = OMI_SMF_setmsg(OMI_E_INPUT, &
		                    "input WavelengthCoefficient_k array too small", &
				    "L1Br_getSIGpix", zero)
            status = OMI_E_FAILURE
            RETURN
         ENDIF
         WavelengthCoefficient_k(1:this%nWavelCoef) = this%WavelenCoef(1:this%nWavelCoef,i,j)
      ENDIF

      IF(PRESENT(WavelengthCoefficientPrec_k)) THEN
         IF(SIZE(WavelengthCoefficientPrec_k) < this%nWavelCoef) THEN
            ierr = OMI_SMF_setmsg(OMI_E_INPUT, &
		                    "input WavelengthCoefficientPrec_k array too small",&
                                    "L1Br_getSIGpix", zero)
            status = OMI_E_FAILURE
            RETURN
         ENDIF
         WavelengthCoefficientPrec_k(1:this%nWavelCoef) = &
                                        this%WavelenPrec(1:this%nWavelCoef,i,j)
      ENDIF

      IF(PRESENT(WavelengthReferenceColumn_k)) &
                                   WavelengthReferenceColumn_k = this%WavelenRefCol(j)

      ! Here is the Small Pixel data

      IF (this%NumTimesSmallPixel > 0) THEN
         IF(PRESENT(NumberSmallPixelColumns_k)) NumberSmallPixelColumns_k = this%NumSmPixCol(j)
         IF(PRESENT(SmallPixelColumn_k)) SmallPixelColumn_k = this%SmPixCol(j)

         IF(PRESENT(SmallPixelSignal_k) .OR. PRESENT(SmallPixelWavelength_k)) THEN

            ! Calculate index into small pixel arrays

            l = 1
            do k = 1, iLine
               l = l + this%NumSmPixCol(k)
            enddo
            k = l + this%NumSmPixCol(iLine+1)

            IF(PRESENT(SmallPixelSignal_k)) THEN
               IF(SIZE(SmallPixelSignal_k) < this%NumSmPixCol(iLine+1)) THEN
                  ierr = OMI_SMF_setmsg(OMI_E_INPUT, &
                                   "input SmallPixelSignal_k array too small", &
                                   "L1Br_getSIGpix", zero)
                  status = OMI_E_FAILURE
                  RETURN
               ENDIF
               SmallPixelSignal_k(1:this%NumSmPixCol(iLine+1)) = this%SmPixRad(i,l:k)
            ENDIF

            IF(PRESENT(SmallPixelWavelength_k)) THEN
               IF(SIZE(SmallPixelWavelength_k) < this%NumSmPixCol(iLine+1)) THEN 
                  ierr = OMI_SMF_setmsg(OMI_E_INPUT, &
                                   "input SmallPixelWavelength_k array too small",&
                                   "L1Br_getSIGpix", zero)
                  status = OMI_E_FAILURE
                  RETURN
               ENDIF
               SmallPixelWavelength_k(1:this%NumSmPixCol(iLine+1)) = this%SmPixWavelen(i,l:k)
            ENDIF
         ENDIF  ! PRESENT(SmallPixelSignal_k) .OR. PRESENT(SmallPixelWavelength_k)
      ENDIF  ! this%NumTimesSmallPixel > 0

      RETURN
      END FUNCTION L1Br_getSIGpix

!! 8. L1Br_getSIGline
 !
 !    Functionality:
 !
 !       This function gets one or more radiance data field values from
 !       the data block.  This function returns values for each given output
 !       keyword for the each wavelength within the input wavelength range.
 !       If no wavelength range is input, values are returned for all
 !       wavelengths.
 !
 !    Calling Arguments:
 !
 !    Inputs:
 !
 !       this: the block data structure
 !       iLine: the line number in the L1B swath.  NOTE: this input is 0 based
 !              range from 0 to (NumTimes-1) inclusive.
 !       iPix: the pixel number in the L1B swath.  NOTE: this input is 0 based
 !             range from 0 to (nXtrack-1) inclusive.
 !       Wlmin_k: input wavelength value for the band beginning (in nm) (Optional)
 !       Wlmax_k: input wavelength value for the band end (in nm) (Optional)
 !
 !       Note that Wlmin_k < Wlmax_k.
 !       If Wlmin_k is omitted, all wavelengths < Wlmax_k are returned.
 !       If Wlmax_k is omitted, all wavelengths > Wlmin_k are returned.
 !       If Wlmin_k and Wlmax_k are omitted, all wavelengths are returned.
 !
 !    Outputs:
 !
 !       The following are keyword arguments.  Only those present in the argument 
 !       list will be set by the function.  The first set of keywords is for
 !       combined data, whereas the second set is for raw data; the two sets are
 !       not mutually exclusive, i.e., raw and combined data can be extracted
 !       in the same function call.
 !
 !       Keyword                    Corresponding L1B Data
 !       =======                    ======================
 !       Signal_k                   Radiance or Irradiance or Signal
 !                                     Mantissa*10^(Exponent)
 !       SignalPrecision_k          RadiancePrecision or IrradiancePrecision or
 !                                     SignalPrecision
 !                                     Precision*10^(Exponent)
 !       PixelQualityFlags_k        PixelQualityFlags
 !       Wavelength_k               Wavelengths
 !                                     WavelengthCoefficient aggregated into Wavelengths
 !       Nwl_k                      None - just the number of wavelengths returned
 !
 !       If the user wants the raw L1B data, rather than combined data, that can be
 !       retrieved using the following keywords.
 !
 !       Keyword                     Corresponding L1B Data Field
 !       =======                     ============================
 !       Mantissa_k                  RadianceMantissa or IrradianceMantissa or SignalMantissa
 !       Precision_k                 RadiancePrecisionMantissa or IrradiancePrecisionMantissa
 !                                      or SignalPrecisionMantissa
 !       Exponent_k                  RadianceExponent or IrradianceExponent or SignalExponent
 !       WavelengthCoefficient_k     WavelengthCoefficient
 !       WavelengthCoefficientPrec_k WavelengthCoefficientPrecision
 !       WavelengthReferenceColumn_k WavelengthReferenceColumn
 !       SmallPixelSignal_k          SmallPixelRadiance or SmallPixelIrradiance or
 !                                      SmallPixelSignal
 !       SmallPixelWavelength_k      SmallPixelWavelength
 !
 !       status: the return PGS_SMF status value
 !
 !    Change History:
 !
 !       Date            Author          Modifications
 !       ====            ======          =============
 !       January 2005    Jeremy Warner   Original Source (derived from original 
 !                                             L1B Reader code written by Kai Yang)
 !
!!
      FUNCTION L1Br_getSIGline(this, iLine, Wlmin_k, Wlmax_k, Signal_k, &
                             SignalPrecision_k, PixelQualityFlags_k, Wavelength_k, &
                             Nwl_k, Mantissa_k, Precision_k, Exponent_k, &
			     WavelengthCoefficient_k, WavelengthCoefficientPrec_k, &
			     WavelengthReferenceColumn_k, SmallPixelSignal_k, &
			     SmallPixelWavelength_k, SmallPixelColumn_k, &
			     NumberSmallPixelColumns_k) RESULT (status)
!tpk       TrueZoom_k ) RESULT (status)     !tpk
        USE Szoom_Parameter_Module
        INCLUDE 'PGS_OMI_1900.f'  !defines the L1B PGS error codes
        TYPE (L1B_block_type), INTENT(INOUT) :: this
        INTEGER, PARAMETER :: i1 = 1
        INTEGER (KIND=i1), PARAMETER :: global_mode = 8_i1, szoom_mode = 4_i1
        INTEGER (KIND = 4), INTENT(IN) :: iLine 
        !tpk INTEGER (KIND = 1), INTENT(IN), OPTIONAL :: TrueZoom_k    !tpk 
        REAL (KIND = 4), OPTIONAL, INTENT(IN) :: Wlmin_k, Wlmax_k    ! wavelength range
        INTEGER (KIND = 4), OPTIONAL, INTENT(OUT) :: Nwl_k
        REAL (KIND = 4), OPTIONAL, DIMENSION(:,:), INTENT(OUT) :: Signal_k, &
                                             SignalPrecision_k, Wavelength_k, &
					     WavelengthCoefficient_k, &
					     WavelengthCoefficientPrec_k, &
					     SmallPixelSignal_k, SmallPixelWavelength_k 
        INTEGER (KIND = 2), OPTIONAL, DIMENSION(:,:), INTENT(OUT) :: &
                                             PixelQualityFlags_k, Mantissa_k, Precision_k
        INTEGER (KIND = 2), OPTIONAL, INTENT(OUT) :: WavelengthReferenceColumn_k
        INTEGER (KIND = 2), OPTIONAL, INTENT(OUT) :: SmallPixelColumn_k
        INTEGER (KIND = 1), OPTIONAL, INTENT(OUT) :: NumberSmallPixelColumns_k
        INTEGER (KIND = 1), OPTIONAL, DIMENSION(:,:), INTENT(OUT) :: Exponent_k


        REAL (KIND = 4),  DIMENSION(1:this%nWavel,1:this%nXtrack) :: tmp_Signal_k, &
                                             tmp_SignalPrecision_k, tmp_Wavelength_k
        REAL (KIND = 4),  DIMENSION(1:this%nWavelCoef, 1:this%nXtrack) :: &
					     tmp_WavelengthCoefficient_k, &
					     tmp_WavelengthCoefficientPrec_k
        REAL (KIND = 4),  DIMENSION(1:this%nXtrack, 1:this%NumSmPixCol(iLine+1)) :: &
					     tmp_SmallPixelSignal_k, tmp_SmallPixelWavelength_k 
        INTEGER (KIND = 2),  DIMENSION(1:this%nWavel,1:this%nXtrack) :: &
                                             tmp_PixelQualityFlags_k, tmp_Mantissa_k, tmp_Precision_k
        INTEGER (KIND = 1), DIMENSION(1:this%nWavel,1:this%nXtrack) :: tmp_Exponent_k

        REAL (KIND = 4), DIMENSION(1:this%nWavel,1:this%nXtrack) :: wl_local
        INTEGER :: i, j, k, q, l, iw
        INTEGER (KIND = 4) :: Nwl_l
        INTEGER (KIND = 4) :: il, ih
        INTEGER (KIND = 4) :: status
        REAL (KIND = 4) :: minWl, maxWl

        !tpk LOGICAL :: yn_truezoom   !tpk
        REAL (KIND = 8),  DIMENSION(1:this%nWavel,1:this%nXtrack) :: &
             tmp_Signal_k8, tmp_SignalPrecision_k8

      status = check_blk_line(this, iLine, j) 
      IF(status .NE. OMI_S_SUCCESS) THEN
         ierr = OMI_SMF_setmsg(OMI_E_GENERAL, "Failed retrieving L1B signal data.", &
                                 "L1Br_getSIGline", one)
         RETURN
      ENDIF

      !tpk !tpk -----------------------------------------------------------------
      !tpk !tpk Check for TrueZoom flag (60 cross-track positions instead of 30)
      !tpk !tpk -----------------------------------------------------------------
      !tpk yn_truezoom = .FALSE.
      !tpk IF ( PRESENT (TrueZoom_k) ) THEN
      !tpk    IF (TrueZoom_k > 0_i1 ) yn_truezoom = .TRUE.
      !tpk END IF


      ! Set min and max wavelengths

      IF(PRESENT(Wlmin_k)) THEN
         minWl = Wlmin_k
      ELSE
         minWl = -1.0
      ENDIF

      IF(PRESENT(Wlmax_k)) THEN
         maxWl = Wlmax_k
      ELSE
         maxWl = -1.0
      ENDIF

      ! Calculate wavelengths from L1B data

      status = calc_wl_line(j, this, minWl, maxWl, wl_local, il, ih, Nwl_l)
      IF(status .NE. OMI_S_SUCCESS) THEN
         ierr = OMI_SMF_setmsg(OMI_E_GENERAL, "Failed calculating L1B wavelengths.", &
                                 "L1Br_getSIGline", one)
         RETURN
      ENDIF

      IF (Nwl_l > 0) THEN
        IF(PRESENT(Signal_k)) THEN
           IF(SIZE(Signal_k, 1) < Nwl_l  .OR. &
               SIZE(Signal_k, 2) < this%nXtrack) THEN
              ierr = OMI_SMF_setmsg(OMI_E_INPUT, &
                                    "input Signal_k array too small", &
                                    "L1Br_getSIGline", zero)
              status = OMI_E_FAILURE
              RETURN
           ENDIF

           tmp_Signal_k8(1:Nwl_l, 1:this%nXtrack) =  &
                REAL(this%RadMantissa(il:ih,1:this%nXtrack,j),KIND=8) * &
                10.0_8**this%RadExponent(il:ih,1:this%nXtrack,j)

           Signal_k(1:Nwl_l, 1:this%nXtrack) =  &
                REAL(tmp_Signal_k8(1:Nwl_l, 1:this%nXtrack), KIND=4)
           !tpkSignal_k(1:Nwl_l, 1:this%nXtrack) =  &
           !tpk                  this%RadMantissa(il:ih,1:this%nXtrack,j) &
           !tpk                * 10.0**this%RadExponent(il:ih,1:this%nXtrack,j)
           DO q = il, ih
              DO i = 1, this%nXtrack
                 IF(this%RadMantissa(q,i,j) .eq. -32767) Signal_k(q-il+1,i) = rfill
              ENDDO
           ENDDO

           IF (this%ImBinFact(j) == szoom_mode) THEN
              IF (this%SpaZoomIndex == 0) THEN
           !tpk IF (this%ImBinFact(j) == szoom_mode .AND. (.NOT. yn_truezoom) ) THEN  !tpk
                 DO iw = 1, Nwl_l
                    DO i = 1, this%nXtrack
                       tmp_Signal_k(iw,i) = Signal_k(iw,i)
                    ENDDO
                 ENDDO
                 Signal_k(1:Nwl_l,gfpix:zspix) = rfill
                 Signal_k(1:Nwl_l,zspix+1:zspix+zlpix) =  tmp_Signal_k(1:Nwl_l,zfpix:zlpix)
                 Signal_k(1:Nwl_l,zlpix+zspix+1:glpix) =  rfill
              ENDIF
           ENDIF

        ENDIF

        IF(PRESENT(SignalPrecision_k)) THEN
           IF(SIZE(SignalPrecision_k, 1) < Nwl_l .OR. &
               SIZE(SignalPrecision_k, 2) < this%nXtrack) THEN
              ierr = OMI_SMF_setmsg(OMI_E_INPUT, &
                                   "input SignalPrecision_k array too small",&
                                    "L1Br_getSIGline", zero)
              status = OMI_E_FAILURE
              RETURN
           ENDIF

           tmp_SignalPrecision_k8(1:Nwl_l, 1:this%nXtrack) =              &
                REAL(this%RadPrecision(il:ih,1:this%nXtrack,j),KIND=8) *  &
                10.0_8**this%RadExponent(il:ih,1:this%nXtrack,j)
           SignalPrecision_k(1:Nwl_l, 1:this%nXtrack) =                   &
                REAL ( tmp_SignalPrecision_k8(1:Nwl_l, 1:this%nXtrack), KIND=4 )

          !tpkSignalPrecision_k(1:Nwl_l, 1:this%nXtrack) = &
          !tpk     this%RadPrecision(il:ih,1:this%nXtrack,j) &
          !tpk     * 10.0**this%RadExponent(il:ih,1:this%nXtrack,j)

           DO q = il, ih
            DO i = 1, this%nXtrack
              IF(this%RadPrecision(q,i,j) .eq. -32767) SignalPrecision_k(q-il+1,i) = rfill
            ENDDO
           ENDDO

           IF (this%ImBinFact(j) == szoom_mode) THEN
              IF (this%SpaZoomIndex == 0) THEN
           !tpk IF (this%ImBinFact(j) == szoom_mode .AND. (.NOT. yn_truezoom) ) THEN  !tpk
                 DO iw = 1, Nwl_l
                    DO i = 1, this%nXtrack
                       tmp_SignalPrecision_k(iw,i) = SignalPrecision_k(iw,i)
                    ENDDO
                 ENDDO
                 SignalPrecision_k(1:Nwl_l,gfpix:zspix) =  rfill
                 SignalPrecision_k(1:Nwl_l,zspix+1:zspix+zlpix) =  tmp_SignalPrecision_k(1:Nwl_l,zfpix:zlpix)
                 SignalPrecision_k(1:Nwl_l,zlpix+zspix+1:glpix) = rfill 
              ENDIF
           ENDIF
        ENDIF

        IF(PRESENT(PixelQualityFlags_k)) THEN
           IF(SIZE(PixelQualityFlags_k, 1) < Nwl_l .OR. &
               SIZE(PixelQualityFlags_k, 2) < this%nXtrack) THEN
              ierr = OMI_SMF_setmsg(OMI_E_INPUT, &
                                    "input PixelQualityFlags_k array too small",&
                                    "L1Br_getSIGline", zero)
              status = OMI_E_FAILURE
              RETURN
           ENDIF
           PixelQualityFlags_k(1:Nwl_l, 1:this%nXtrack) = this%PQFlag(il:ih,1:this%nXtrack,j)
           IF (this%ImBinFact(j) == szoom_mode) THEN
              IF (this%SpaZoomIndex == 0) THEN
           !tpk IF (this%ImBinFact(j) == szoom_mode .AND. (.NOT. yn_truezoom) ) THEN  !tpk
                 DO iw = 1, Nwl_l
                    DO i = 1, this%nXtrack
                       tmp_PixelQualityFlags_k(iw,i) = PixelQualityFlags_k(iw,i)
                    ENDDO
                 ENDDO
                 PixelQualityFlags_k(1:Nwl_l,gfpix:zspix) =  int16fill
                                   
                 PixelQualityFlags_k(1:Nwl_l,zspix+1:zspix+zlpix) =  &
                      tmp_PixelQualityFlags_k(1:Nwl_l,zfpix:zlpix)
                 PixelQualityFlags_k(1:Nwl_l,zlpix+zspix+1:glpix) =  int16fill
              ENDIF
           ENDIF
        ENDIF

        IF(PRESENT(Wavelength_k)) THEN
           IF(SIZE(Wavelength_k, 1) < Nwl_l .OR. &
               SIZE(Wavelength_k, 2) < this%nXtrack) THEN 
              ierr = OMI_SMF_setmsg(OMI_E_INPUT, &
                                    "input Wavelength_k array too small",&
                                    "L1Br_getSIGline", zero)
              status = OMI_E_FAILURE
              RETURN
           ENDIF
           Wavelength_k(1:Nwl_l, 1:this%nXtrack) = wl_local(il:ih,1:this%nXtrack)

           IF (this%ImBinFact(j) == szoom_mode) THEN
                IF (this%SpaZoomIndex == 0) THEN
           !tpk IF (this%ImBinFact(j) == szoom_mode .AND. (.NOT. yn_truezoom) ) THEN  !tpk
                   DO iw = 1, Nwl_l
                      DO i = 1, this%nXtrack
                         tmp_Wavelength_k(iw,i) = Wavelength_k(iw,i)
                      ENDDO
                   ENDDO
                   Wavelength_k(1:Nwl_l,gfpix:zspix) = rfill 
                   Wavelength_k(1:Nwl_l,zspix+1:zspix+zlpix) =  tmp_Wavelength_k(1:Nwl_l,zfpix:zlpix)
                   Wavelength_k(1:Nwl_l,zlpix+zspix+1:glpix) = rfill 
                ENDIF
             ENDIF
        ENDIF

        ! Here is the raw data

        IF(PRESENT(Mantissa_k)) THEN
           IF(SIZE(Mantissa_k, 1) < Nwl_l  .OR. &
              SIZE(Mantissa_k, 2) < this%nXtrack) THEN
              ierr = OMI_SMF_setmsg(OMI_E_INPUT, "input Mantissa_k array too small", &
                                    "L1Br_getSIGline", zero)
              status = OMI_E_FAILURE
              RETURN
           ENDIF
           Mantissa_k(1:Nwl_l, 1:this%nXtrack) = this%RadMantissa(il:ih,1:this%nXtrack,j)
           IF (this%ImBinFact(j) == szoom_mode) THEN
              IF (this%SpaZoomIndex == 0) THEN
           !tpk IF (this%ImBinFact(j) == szoom_mode .AND. (.NOT. yn_truezoom) ) THEN  !tpk
                 DO iw = 1, Nwl_l
                    DO i = 1, this%nXtrack
                       tmp_Mantissa_k(iw,i) = Mantissa_k(iw,i)
                    ENDDO
                 ENDDO
                 Mantissa_k(1:Nwl_l,gfpix:zspix) = int16fill
                 Mantissa_k(1:Nwl_l,zspix+1:zspix+zlpix) =  tmp_Mantissa_k(1:Nwl_l,zfpix:zlpix)
                 Mantissa_k(1:Nwl_l,zlpix+zspix+1:glpix) = int16fill
              ENDIF
           ENDIF
        ENDIF

        IF(PRESENT(Precision_k)) THEN
           IF(SIZE(Precision_k, 1) < Nwl_l  .OR. &
              SIZE(Precision_k, 2) < this%nXtrack) THEN
              ierr = OMI_SMF_setmsg(OMI_E_INPUT, "input Precision_k array too small", &
                                    "L1Br_getSIGline", zero)
              status = OMI_E_FAILURE
              RETURN
           ENDIF
           Precision_k(1:Nwl_l, 1:this%nXtrack) = this%RadPrecision(il:ih,1:this%nXtrack,j)
           IF (this%ImBinFact(j) == szoom_mode) THEN
              IF (this%SpaZoomIndex == 0) THEN
           !tpk IF (this%ImBinFact(j) == szoom_mode .AND. (.NOT. yn_truezoom) ) THEN  !tpk
                 DO iw = 1, Nwl_l
                    DO i = 1, this%nXtrack
                       tmp_Precision_k(iw,i) = Precision_k(iw,i)
                    ENDDO
                 ENDDO
                 Precision_k(1:Nwl_l,gfpix:zspix) = int16fill
                 Precision_k(1:Nwl_l,zspix+1:zspix+zlpix) =  tmp_Precision_k(1:Nwl_l,zfpix:zlpix)
                 Precision_k(1:Nwl_l,zlpix+zspix+1:glpix) = int16fill
              ENDIF
           ENDIF
        ENDIF

        IF(PRESENT(Exponent_k)) THEN
           IF(SIZE(Exponent_k, 1) < Nwl_l  .OR. &
                SIZE(Exponent_k, 2) < this%nXtrack) THEN
              ierr = OMI_SMF_setmsg(OMI_E_INPUT, "input Exponent_k array too small", &
                   "L1Br_getSIGline", zero)
              status = OMI_E_FAILURE
              RETURN
           ENDIF
           Exponent_k(1:Nwl_l, 1:this%nXtrack) = this%RadExponent(il:ih,1:this%nXtrack,j)
           IF (this%ImBinFact(j) == szoom_mode) THEN
              IF (this%SpaZoomIndex == 0) THEN
           !tpk IF (this%ImBinFact(j) == szoom_mode .AND. (.NOT. yn_truezoom) ) THEN  !tpk
                 DO iw = 1, Nwl_l
                    DO i = 1, this%nXtrack
                       tmp_Exponent_k(iw,i) = Exponent_k(iw,i)
                    ENDDO
                 ENDDO
                 Exponent_k(1:Nwl_l,gfpix:zspix) = int8fill
                 Exponent_k(1:Nwl_l,zspix+1:zspix+zlpix) =  tmp_Exponent_k(1:Nwl_l,zfpix:zlpix)
                 Exponent_k(1:Nwl_l,zlpix+zspix+1:glpix) = int8fill
              ENDIF
           ENDIF
        ENDIF

     ENDIF  ! Nwl_l > 0

     IF(PRESENT(Nwl_k)) Nwl_k = Nwl_l

     ! The rest of the outputs do not depend on the value of Nwl_l

     IF(PRESENT(WavelengthCoefficient_k)) THEN
        IF(SIZE(WavelengthCoefficient_k, 1) < this%nWavelCoef  .OR. &
             SIZE(WavelengthCoefficient_k, 2) < this%nXtrack) THEN
           ierr = OMI_SMF_setmsg(OMI_E_INPUT, &
                "input WavelengthCoefficient_k array too small", &
                "L1Br_getSIGline", zero)
           status = OMI_E_FAILURE
           RETURN
        ENDIF
        WavelengthCoefficient_k(1:this%nWavelCoef, 1:this%nXtrack) = &
             this%WavelenCoef(1:this%nWavelCoef,1:this%nXtrack,j)


        IF (this%ImBinFact(j) == szoom_mode) THEN
           IF (this%SpaZoomIndex == 0) THEN
        !tpk IF (this%ImBinFact(j) == szoom_mode .AND. (.NOT. yn_truezoom) ) THEN  !tpk
              DO iw = 1, this%nWavelCoef
                 DO i = 1, this%nXtrack
                    tmp_WavelengthCoefficient_k(iw,i) = WavelengthCoefficient_k(iw,i)
                 ENDDO
              ENDDO
              WavelengthCoefficient_k(1:this%nWavelCoef,gfpix:zspix) = rfill   
              
              WavelengthCoefficient_k(1:this%nWavelCoef,zspix+1:zspix+zlpix) =  &   
                   tmp_WavelengthCoefficient_k(1:this%nWavelCoef,zfpix:zlpix)
              WavelengthCoefficient_k(1:this%nWavelCoef,zlpix+zspix+1:glpix) = rfill   
           ENDIF
        ENDIF
     ENDIF

     IF(PRESENT(WavelengthCoefficientPrec_k)) THEN
        IF(SIZE(WavelengthCoefficientPrec_k, 1) < this%nWavelCoef  .OR. &
             SIZE(WavelengthCoefficientPrec_k, 2) < this%nXtrack) THEN
           ierr = OMI_SMF_setmsg(OMI_E_INPUT, &
                "input WavelengthCoefficientPrec_k array too small", &
                "L1Br_getSIGline", zero)
           status = OMI_E_FAILURE
           RETURN
        ENDIF
        WavelengthCoefficientPrec_k(1:this%nWavelCoef, 1:this%nXtrack) = &
             this%WavelenPrec(1:this%nWavelCoef,1:this%nXtrack,j)


        IF (this%ImBinFact(j) == szoom_mode) THEN
           IF (this%SpaZoomIndex == 0) THEN
              !IF (this%ImBinFact(j) == szoom_mode .AND. (.NOT. yn_truezoom) ) THEN  !tpk
              DO iw = 1, this%nWavelCoef
                 DO i = 1, this%nXtrack
                    tmp_WavelengthCoefficientPrec_k(iw,i) = WavelengthCoefficientPrec_k(iw,i)
                 ENDDO
              ENDDO
              WavelengthCoefficientPrec_k(1:this%nWavelCoef,gfpix:zspix) = rfill     
              WavelengthCoefficientPrec_k(1:this%nWavelCoef,zspix+1:zspix+zlpix) =  &   
                   tmp_WavelengthCoefficientPrec_k(1:this%nWavelCoef,zfpix:zlpix)
              WavelengthCoefficientPrec_k(1:this%nWavelCoef,zlpix+zspix+1:glpix) = &
                   rfill  
           ENDIF
        ENDIF
     ENDIF

      IF(PRESENT(WavelengthReferenceColumn_k)) &
                                   WavelengthReferenceColumn_k = this%WavelenRefCol(j)

      ! Here is the Small Pixel data

      IF (this%NumTimesSmallPixel > 0) THEN
         IF(PRESENT(NumberSmallPixelColumns_k)) NumberSmallPixelColumns_k = this%NumSmPixCol(j)
         IF(PRESENT(SmallPixelColumn_k)) SmallPixelColumn_k = this%SmPixCol(j)

         IF(PRESENT(SmallPixelSignal_k) .OR. PRESENT(SmallPixelWavelength_k)) THEN

            ! Calculate index into small pixel arrays

            l = 1
            do k = 1, iLine
               l = l + this%NumSmPixCol(k)
            enddo
            k = l + this%NumSmPixCol(iLine+1)

            IF(PRESENT(SmallPixelSignal_k)) THEN
              IF(SIZE(SmallPixelSignal_k, 1) < this%nXtrack  .OR. &
                 SIZE(SmallPixelSignal_k, 2) < this%NumSmPixCol(iLine+1)) THEN
                 ierr = OMI_SMF_setmsg(OMI_E_INPUT, &
		                    "input SmallPixelSignal_k array too small", &
                                    "L1Br_getSIGline", zero)
                 status = OMI_E_FAILURE
                 RETURN
              ENDIF
              SmallPixelSignal_k(1:this%nXtrack, 1:this%NumSmPixCol(iLine+1)) = &
                             this%SmPixRad(1:this%nXtrack,l:k)

              IF (this%ImBinFact(j) == szoom_mode) THEN
                 IF (this%SpaZoomIndex == 0) THEN
                    !tpk IF (this%ImBinFact(j) == szoom_mode .AND. (.NOT. yn_truezoom) ) THEN  !tpk
                    DO i = 1, this%nXtrack
                       DO iw = 1, this%NumSmPixCol(iLine+1)
                          tmp_SmallPixelSignal_k(i,iw) = SmallPixelSignal_k(i,iw)
                       ENDDO
                    ENDDO
                    SmallPixelSignal_k(gfpix:zspix,1:this%NumSmPixCol(iLine+1)) =  & 
                         rfill
                    SmallPixelSignal_k(zspix+1:zspix+zlpix,1:this%NumSmPixCol(iLine+1)) =  & 
                         tmp_SmallPixelSignal_k(zfpix:zlpix,1:this%NumSmPixCol(iLine+1))
                    SmallPixelSignal_k(zlpix+zspix+1:glpix,1:this%NumSmPixCol(iLine+1)) =  &
                         rfill
                 ENDIF
              ENDIF
           ENDIF

            IF(PRESENT(SmallPixelWavelength_k)) THEN
              IF(SIZE(SmallPixelWavelength_k, 1) < this%nXtrack  .OR. &
                 SIZE(SmallPixelWavelength_k, 2) < this%NumSmPixCol(iLine+1)) THEN
                 ierr = OMI_SMF_setmsg(OMI_E_INPUT, &
		                    "input SmallPixelWavelength_k array too small", &
                                    "L1Br_getSIGline", zero)
                 status = OMI_E_FAILURE
                 RETURN
              ENDIF
              SmallPixelWavelength_k(1:this%nXtrack, 1:this%NumSmPixCol(iLine+1)) = &
                             this%SmPixWavelen(1:this%nXtrack,l:k)
              IF (this%ImBinFact(j) == szoom_mode) THEN
                 IF (this%SpaZoomIndex == 0) THEN
                    !tpk IF (this%ImBinFact(j) == szoom_mode .AND. (.NOT. yn_truezoom) ) THEN  !tpk
                    DO i = 1, this%nXtrack
                       DO iw = 1, this%NumSmPixCol(iLine+1)
                          tmp_SmallPixelWavelength_k(i,iw) = SmallPixelWavelength_k(i,iw)
                       ENDDO
                    ENDDO
                    SmallPixelWavelength_k(gfpix:zspix,1:this%NumSmPixCol(iLine+1)) =  &
                         rfill
                    SmallPixelWavelength_k(zspix+1:zspix+zlpix,1:this%NumSmPixCol(iLine+1)) = & 
                         tmp_SmallPixelWavelength_k(zfpix:zlpix,1:this%NumSmPixCol(iLine+1))
                    SmallPixelWavelength_k(zlpix+zspix+1:glpix,1:this%NumSmPixCol(iLine+1)) = & 
                         rfill
                 ENDIF
              ENDIF
           ENDIF

        ENDIF  ! PRESENT(SmallPixelSignal_k) .OR. PRESENT(SmallPixelWavelength_k)
     ENDIF  ! this%NumTimesSmallPixel > 0
     RETURN
   END FUNCTION L1Br_getSIGline

!! 9. L1Br_getSIGpixWL
 !
 !    Functionality:
 !
 !       This function gets one or more radiance data field values from
 !       the data block.  This function returns values for each given output
 !       keyword for the each wavelength specified in the input wavelength list.
 !       This function only returns data from a single pixel in the L1B file.
 !       If the exact wavelength is not found, the data is interpolated via a
 !       simple linear interpolation to find the desired approximate value, which
 !       is then returned.
 !
 !    Calling Arguments:
 !
 !    Inputs:
 !
 !       this: the block data structure
 !       iLine: the line number in the L1B swath.  NOTE: this input is 0 based
 !              range from 0 to (NumTimes-1) inclusive.
 !       iPix: the pixel number in the L1B swath.  NOTE: this input is 0 based
 !              range from 0 to (nXtrack-1) inclusive.
 !       Wavelength_k:  list of input wavelengths
 !       Nwl_k is optional.  If present, then it indicates the number of valid
 !              wavelengths in Wavelength_k.  If not, then all the values in the
 !              Wavelength_k are assumed to be valid.
 !
 !    Outputs:
 !
 !       The following are keyword arguments.  Only those present in the argument 
 !       list will be set by the function.
 !
 !       Keyword                    Corresponding L1B Data
 !       =======                    ======================
 !       Signal_k                   Radiance or Irradiance or Signal
 !                                     Mantissa*10^(Exponent)
 !       SignalPrecision_k          RadiancePrecision or IrradiancePrecision or
 !                                     SignalPrecision
 !                                     Precision*10^(Exponent)
 !       PixelQualityFlags_k        PixelQualityFlags
 !
 !       status: the return PGS_SMF status value
 !
 !    Change History:
 !
 !       Date            Author          Modifications
 !       ====            ======          =============
 !       January 2005    Jeremy Warner   Original Source (derived from original 
 !                                             L1B Reader code written by Kai Yang)
 !
!!
      FUNCTION L1Br_getSIGpixWL(this, iLine, iPix, Wavelength_k, Nwl_k, Signal_k, &
                               SignalPrecision_k, PixelQualityFlags_k) RESULT (status)
        INCLUDE 'PGS_OMI_1900.f'  !defines the L1B PGS error codes
        TYPE (L1B_block_type), INTENT(INOUT) :: this
        INTEGER (KIND = 4), INTENT(IN) :: iLine, iPix
        REAL (KIND = 4), DIMENSION(:), INTENT(IN) :: Wavelength_k
        INTEGER (KIND = 4), OPTIONAL, INTENT(IN) :: Nwl_k
        REAL (KIND = 4), OPTIONAL, DIMENSION(:), INTENT(OUT) :: Signal_k, SignalPrecision_k
        INTEGER (KIND=2), OPTIONAL, DIMENSION(:), INTENT(OUT) :: PixelQualityFlags_k
        REAL (KIND = 4), DIMENSION(1:this%nWavel) :: wl_local
        INTEGER (KIND = 4) :: Nwl_l
        INTEGER :: i, j, k, q
        INTEGER :: fac
        INTEGER :: iw, last_good
        INTEGER (KIND = 4) :: il, ih, i_foo(1)
        REAL (KIND = 4) :: rad_il, rad_ih, frac, loc_max, loc_min
        INTEGER (KIND = 4) :: status
 
        status = check_blk_pix(this, iLine, iPix, j, i) 
        IF(status .NE. OMI_S_SUCCESS) THEN
           ierr = OMI_SMF_setmsg(OMI_E_GENERAL, "Failed retrieving L1B signal data.", &
                                 "L1Br_getSIGpixWL", one)
           RETURN
        ENDIF

        loc_max = -1.0
        loc_min = -1.0
        DO k = 1, this%nWavel
           wl_local(k) = 1.0
           DO q = 1, this%nWavelCoef
              IF(this%WavelenCoef(q,i,j) .EQ. rfill) wl_local(k) = rfill
           ENDDO
           fac = k-1-this%WavelenRefCol(j)
           IF (wl_local(k) .gt. 0) THEN

              wl_local(k) = fac*this%WavelenCoef(this%nWavelCoef,i,j)
              DO q = this%nWavelCoef -1, 2, -1
                 wl_local(k) = fac*(wl_local(k) + this%WavelenCoef(q,i,j))
              ENDDO
              wl_local(k) = wl_local(k) + this%WavelenCoef(1,i,j)
              if (loc_min .lt. 0 .or. loc_min .gt. wl_local(k)) loc_min = wl_local(k)
              if (loc_max .lt. 0 .or. loc_max .lt. wl_local(k)) loc_max = wl_local(k)
           ENDIF
        ENDDO

        IF(PRESENT(Nwl_k)) THEN
           Nwl_l = Nwl_k
           IF(SIZE(Wavelength_k) < Nwl_l) THEN  
              ierr = OMI_SMF_setmsg(OMI_E_INPUT, &
                                    "input Wavelength_k array too small",&
                                    "L1Br_getSIGpixWL", zero)
              status = OMI_E_FAILURE
              RETURN
           ENDIF
        ELSE
           Nwl_l = SIZE(Wavelength_k)
        ENDIF

        IF(PRESENT(Signal_k)) THEN
           IF(SIZE(Signal_k) < Nwl_l) THEN
               ierr = OMI_SMF_setmsg(OMI_E_INPUT, &
                                     "input Signal_k array too small", &
                                     "L1Br_getSIGpixWL", zero)
              status = OMI_E_FAILURE
              RETURN
           ENDIF
        ENDIF

        IF(PRESENT(SignalPrecision_k)) THEN
           IF(SIZE(SignalPrecision_k) < Nwl_l) THEN
              ierr = OMI_SMF_setmsg(OMI_E_INPUT, &
                                  "input SignalPrecision_k array too small",&
                                    "L1Br_getSIGpixWL", zero)
              status = OMI_E_FAILURE
              RETURN
           ENDIF
        ENDIF

        IF(PRESENT(PixelQualityFlags_k)) THEN
           IF(SIZE(PixelQualityFlags_k) < Nwl_l) THEN
              ierr = OMI_SMF_setmsg(OMI_E_INPUT, &
                                    "input PixelQualityFlags_k array too small",&
                                    "L1Br_getSIGpixWL", zero)
              status = OMI_E_FAILURE
              RETURN
           ENDIF
        ENDIF

        DO iw = 1, Nwl_l
           IF(Wavelength_k(iw) > loc_max  .OR. &
               Wavelength_k(iw) < loc_min) THEN
              IF(PRESENT(Signal_k)) Signal_k(iw) =  0
              IF(PRESENT(SignalPrecision_k)) &
                                       SignalPrecision_k(iw) = 0
              IF(PRESENT(PixelQualityFlags_k)) PixelQualityFlags_k(iw) = 0
              CYCLE
           ENDIF

           il = 0
           ih = 0
           last_good = 0
           DO k = 1, this%nWavel
              if (wl_local(k) .gt. Wavelength_k(iw) .and. il .eq. 0) then
                 il = last_good
                 ih = k
                 exit
              else if (wl_local(k) .eq. Wavelength_k(iw)) then
                 il = k
                 ih = k
                 exit
              else if (wl_local(k) .ne. rfill) then
                 last_good = k
              end if
           ENDDO

           ! Interpolate

           if (il .eq. 0 .or. ih .eq. 0) then
              IF(PRESENT(Signal_k)) Signal_k(iw) = rfill
              IF(PRESENT(SignalPrecision_k)) SignalPrecision_k(iw) = rfill
              IF(PRESENT(PixelQualityFlags_k)) PixelQualityFlags_k(iw) = &
                 int16fill
           else IF(il .EQ. ih) THEN
              IF(PRESENT(Signal_k)) Signal_k(iw) = &
                 this%RadMantissa(il,i,j) * 10.0**this%RadExponent(il,i,j)
              IF(PRESENT(SignalPrecision_k)) SignalPrecision_k(iw) = &
                 this%RadPrecision(il,i,j) * 10.0**this%RadExponent(il,i,j)
              IF(PRESENT(PixelQualityFlags_k)) PixelQualityFlags_k(iw) = &
                 this%PQFlag(il,i,j)
           ELSE  
              IF(PRESENT(Signal_k)) THEN
                 frac = (Wavelength_k(iw) - wl_local(il))/ &
                        (wl_local(ih) - wl_local(il))
                 rad_il = this%RadMantissa(il,i,j) * 10.0**this%RadExponent(il,i,j) 
                 rad_ih = this%RadMantissa(ih,i,j) * 10.0**this%RadExponent(ih,i,j) 
                 Signal_k(iw) = rad_il*(1.0+frac*(rad_ih/rad_il-1.0))
                 if (rad_il .eq. 0.0) Signal_k(iw) = rfill
              ENDIF

              IF(PRESENT(SignalPrecision_k)) THEN
                 frac = (Wavelength_k(iw) - wl_local(il))/ &
                        (wl_local(ih) - wl_local(il))
                 rad_il = this%RadPrecision(il,i,j) * 10.0**this%RadExponent(il,i,j) 
                 rad_ih = this%RadPrecision(ih,i,j) * 10.0**this%RadExponent(ih,i,j) 
                 SignalPrecision_k(iw) = rad_il*(1.0+frac*(rad_ih/rad_il-1.0))
              ENDIF

              IF(PRESENT(PixelQualityFlags_k)) THEN
                 PixelQualityFlags_k(iw) = this%PQFlag(il,i,j)
              ENDIF
           ENDIF
        ENDDO
        RETURN
      END FUNCTION L1Br_getSIGpixWL

!! 10. L1Br_getSIGlineWL
 !
 !    Functionality:
 !
 !       This function gets one or more radiance data field values from
 !       the data block.  This function returns values for each given output
 !       keyword for the each wavelength specified in the input wavelength list.
 !       This function returns data from an entire line of the L1B file.
 !       If the exact wavelength is not found, the data is interpolated via a
 !       simple linear interpolation to find the desired approximate value, which
 !       is then returned.
 !
 !    Calling Arguments:
 !
 !    Inputs:
 !
 !       this: the block data structure
 !       iLine: the line number in the L1B swath.  NOTE: this input is 0 based
 !              range from 0 to (NumTimes-1) inclusive.
 !       iPix: the pixel number in the L1B swath.  NOTE: this input is 0 based
 !              range from 0 to (nXtrack-1) inclusive.
 !       Wavelength_k:  list of input wavelengths
 !       Nwl_k is optional.  If present, then it indicates the number of valid
 !              wavelengths in Wavelength_k.  If not, then all the values in the
 !              Wavelength_k are assumed to be valid.
 !
 !    Outputs:
 !
 !       The following are keyword arguments.  Only those present in the argument 
 !       list will be set by the function.
 !
 !       Keyword                    Corresponding L1B Data
 !       =======                    ======================
 !       Signal_k                   Radiance or Irradiance or Signal
 !                                     Mantissa*10^(Exponent)
 !       SignalPrecision_k          RadiancePrecision or IrradiancePrecision or
 !                                     SignalPrecision
 !                                     Precision*10^(Exponent)
 !       PixelQualityFlags_k        PixelQualityFlags
 !
 !       status: the return PGS_SMF status value
 !
 !    Change History:
 !
 !       Date            Author          Modifications
 !       ====            ======          =============
 !       January 2005    Jeremy Warner   Original Source (derived from original 
 !                                             L1B Reader code written by Kai Yang)
 !
!!
      FUNCTION L1Br_getSIGlineWL(this, iLine, Wavelength_k, Nwl_k, Signal_k, &
                                SignalPrecision_k, PixelQualityFlags_k) RESULT (status)
       !tpk                                TrueZoom_k ) RESULT (status)        !tpk
        USE Szoom_Parameter_Module
        INCLUDE 'PGS_OMI_1900.f'  !defines the L1B PGS error codes
        TYPE (L1B_block_type), INTENT(INOUT) :: this
        INTEGER, PARAMETER :: i1 = 1
        INTEGER (KIND=i1), PARAMETER :: global_mode = 8_i1, szoom_mode = 4_i1
        INTEGER (KIND = 4), INTENT(IN) :: iLine
        REAL (KIND = 4), DIMENSION(:), INTENT(IN) :: Wavelength_k
        !tpk INTEGER (KIND = 1), INTENT(IN), OPTIONAL :: TrueZoom_k    !tpk
        INTEGER (KIND = 4), OPTIONAL, INTENT(IN) :: Nwl_k
        REAL (KIND = 4), OPTIONAL, DIMENSION(:,:), INTENT(OUT) :: Signal_k, SignalPrecision_k
        INTEGER (KIND=2), OPTIONAL, DIMENSION(:,:), INTENT(OUT) :: PixelQualityFlags_k
        REAL (KIND = 4), DIMENSION(1:this%nWavel,1:this%nXtrack) :: wl_local
        INTEGER (KIND = 4) :: Nwl_l
        INTEGER :: i, j, k, q
        INTEGER :: fac
        INTEGER :: iw
        INTEGER (KIND = 4) :: il, ih, last_good
        REAL (KIND = 4) :: rad_il, rad_ih, frac
        REAL (KIND = 4), DIMENSION(1:this%nXtrack) :: loc_min, loc_max
        INTEGER (KIND = 4) :: status

        REAL (KIND = 4), DIMENSION(1:this%nWavel,1:this%nXtrack)  :: tmp_Signal_k
        REAL (KIND = 4), DIMENSION(1:this%nWavel,1:this%nXtrack)  :: tmp_SignalPrecision_k
        INTEGER (KIND=2), DIMENSION(1:this%nWavel,1:this%nXtrack) :: tmp_PixelQualityFlags_k

        status = check_blk_line(this, iLine, j) 
        IF(status .NE. OMI_S_SUCCESS) THEN
           ierr = OMI_SMF_setmsg(OMI_E_GENERAL, "Failed retrieving L1B signal data.", &
                                 "L1Br_getSIGlineWL", one)
           RETURN
        ENDIF


        !tpk !tpk -----------------------------------------------------------------
        !tpk !tpk Check for TrueZoom flag (60 cross-track positions instead of 30)
        !tpk !tpk -----------------------------------------------------------------
        !tpk yn_truezoom = .FALSE.
        !tpk IF ( PRESENT (TrueZoom_k) ) THEN
        !tpk    IF (TrueZoom_k > 0_i1 ) yn_truezoom = .TRUE.
        !tpk END IF


        loc_min(1:this%nXtrack) = -1.0
        loc_max (1:this%nXtrack)= -1.0
        DO i = 1, this%nXtrack
           DO k = 1, this%nWavel
              wl_local(k,i) = 1.0
              DO q = 1, this%nWavelCoef
                 IF(this%WavelenCoef(q,i,j) >  1.0E30) wl_local(k,i) = rfill
                 IF(this%WavelenCoef(q,i,j) < -1.0E28) wl_local(k,i) = rfill
              ENDDO
              IF (wl_local(k,i) .gt. 0) THEN
                 
                 fac = k-1-this%WavelenRefCol(j)
                 wl_local(k,i) = fac*this%WavelenCoef(this%nWavelCoef,i,j)
                 DO q = this%nWavelCoef -1, 2, -1
                    wl_local(k,i) = fac*(wl_local(k,i) + this%WavelenCoef(q,i,j))
                 ENDDO
                 wl_local(k,i) = wl_local(k,i) + this%WavelenCoef(1,i,j)
                 if (loc_min(i) .lt. 0 .or. loc_min(i) .gt. wl_local(k,i)) &
                      loc_min(i) = wl_local(k,i)
                 if (loc_max(i) .lt. 0 .or. loc_max(i) .lt. wl_local(k,i)) &
                      loc_max(i) = wl_local(k,i)
              ENDIF
           ENDDO
        ENDDO

        IF(PRESENT(Nwl_k)) THEN
           Nwl_l = Nwl_k
           IF(SIZE(Wavelength_k) < Nwl_l) THEN
              ierr = OMI_SMF_setmsg(OMI_E_INPUT, "input Nwl_l too large",&
                                    "L1Br_getSIGlineWL", zero)
              status = OMI_E_FAILURE
              RETURN
           ENDIF
        ELSE
           Nwl_l = SIZE(Wavelength_k)
        ENDIF

        IF(PRESENT(Signal_k)) THEN
           IF(SIZE(Signal_k, 1) < Nwl_l .OR. &
               SIZE(Signal_k, 2) < this%nXtrack) THEN
               ierr = OMI_SMF_setmsg(OMI_E_INPUT, &
                                     "input Signal_k array too small", &
                                     "L1Br_getSIGlineWL", zero)
              status = OMI_E_FAILURE
              RETURN
           ENDIF
        ENDIF

        IF(PRESENT(SignalPrecision_k)) THEN
           IF(SIZE(SignalPrecision_k, 1) < Nwl_l .OR. &
               SIZE(SignalPrecision_k, 2) <  this%nXtrack) THEN
              ierr = OMI_SMF_setmsg(OMI_E_INPUT, &
                                  "input SignalPrecision_k array too small",&
                                    "L1Br_getSIGlineWL", zero)
              status = OMI_E_FAILURE
              RETURN
           ENDIF
        ENDIF

        IF(PRESENT(PixelQualityFlags_k)) THEN
           IF(SIZE(PixelQualityFlags_k, 1) < Nwl_l .OR. &
               SIZE(PixelQualityFlags_k, 2) < this%nXtrack) THEN
              ierr = OMI_SMF_setmsg(OMI_E_INPUT, &
                                    "input PixelQualityFlags_k array too small",&
                                    "L1Br_getSIGlineWL", zero)
              status = OMI_E_FAILURE
              RETURN
           ENDIF
        ENDIF

        DO iw = 1, Nwl_l
        DO i = 1, this%nXtrack
           IF(Wavelength_k(iw) > loc_max(i)  .OR. &
               Wavelength_k(iw) < loc_min(i)) THEN
              IF(PRESENT(Signal_k)) Signal_k(iw,i) =  0
              IF(PRESENT(SignalPrecision_k)) &
                                       SignalPrecision_k(iw,i) = 0
              IF(PRESENT(PixelQualityFlags_k)) PixelQualityFlags_k(iw,i) = 0
           ELSE
              il = 0
              ih = 0
              last_good = 0
              DO k = 1, this%nWavel
                 if (wl_local(k,i) .gt. Wavelength_k(iw) .and. il .eq. 0) then
                    il = last_good
                    ih = k
                    exit
                 else if (wl_local(k,i) .eq. Wavelength_k(iw)) then
                    il = k
                    ih = k
                    exit
                 else if (wl_local(k,i) .ne. rfill) then
                    last_good = k
                 end if
              ENDDO

              ! Interpolate

              IF (il .eq. 0 .or. ih .eq. 0) THEN
                 IF(PRESENT(Signal_k)) Signal_k(iw,i) = rfill
                 IF(PRESENT(SignalPrecision_k)) SignalPrecision_k(iw,i) =rfill
                 IF(PRESENT(PixelQualityFlags_k)) PixelQualityFlags_k(iw,i) = &
                    int16fill
              ELSE IF(il .EQ. ih) THEN
                 IF(PRESENT(Signal_k)) Signal_k(iw,i) = &
                    this%RadMantissa(il,i,j) * 10.0**this%RadExponent(il,i,j)
                 IF(PRESENT(SignalPrecision_k)) SignalPrecision_k(iw,i) =&
                    this%RadPrecision(il,i,j) * 10.0**this%RadExponent(il,i,j)
                 IF(PRESENT(PixelQualityFlags_k)) PixelQualityFlags_k(iw,i) = &
                    this%PQFlag(il,i,j)
              ELSE  
                 IF(PRESENT(Signal_k)) THEN
                    frac = (Wavelength_k(iw) - wl_local(il,i))/ &
                           (wl_local(ih,i)   - wl_local(il,i))
                    rad_il = this%RadMantissa(il,i,j) * 10.0**this%RadExponent(il,i,j)
                    rad_ih = this%RadMantissa(ih,i,j) * 10.0**this%RadExponent(ih,i,j)
                    Signal_k(iw,i) = rad_il*(1.0+frac*(rad_ih/rad_il-1.0))
                 ENDIF

                 IF(PRESENT(SignalPrecision_k)) THEN
                    rad_il = this%RadPrecision(il,i,j) * 10.0**this%RadExponent(il,i,j)
                    rad_ih = this%RadPrecision(ih,i,j) * 10.0**this%RadExponent(ih,i,j)
                    SignalPrecision_k(iw,i) = MAX(rad_il, rad_ih)
                 ENDIF

                 IF(PRESENT(PixelQualityFlags_k)) THEN
                    PixelQualityFlags_k(iw,i) = this%PQFlag(il,i,j)
                 ENDIF
              ENDIF
           ENDIF
        ENDDO
        ENDDO

        IF (this%ImBinFact(j) == szoom_mode) THEN
           IF (this%SpaZoomIndex == 0) THEN
              !tpk IF (this%ImBinFact(j) == szoom_mode .AND. (.NOT. yn_truezoom) ) THEN  !tpk
              DO iw = 1, Nwl_l
                 DO i = 1, this%nXtrack
                    tmp_Signal_k(iw,i) = Signal_k(iw,i)
                    tmp_SignalPrecision_k(iw,i) =  SignalPrecision_k(iw,i)
                    tmp_PixelQualityFlags_k(iw,i) = PixelQualityFlags_k(iw,i)
                 ENDDO
              ENDDO

           DO iw = 1, Nwl_l
               Signal_k(iw,gfpix:zspix) = rfill
               Signal_k(iw,zspix+1:zspix+zlpix) =  tmp_Signal_k(iw,zfpix:zlpix)
               Signal_k(iw,zlpix+zspix+1:glpix) =  rfill
               SignalPrecision_k(iw,gfpix:zspix) = rfill
               SignalPrecision_k(iw,zspix+1:zspix+zlpix) = tmp_SignalPrecision_k(iw,zfpix:zlpix)
               SignalPrecision_k(iw,zlpix+zspix+1:glpix) = rfill
               PixelQualityFlags_k(iw,gfpix:zspix) = int16fill
               PixelQualityFlags_k(iw,zspix+1:zspix+zlpix) = tmp_PixelQualityFlags_k(iw,zfpix:zlpix)
               PixelQualityFlags_k(iw,zlpix+zspix+1:glpix) = int16fill
           ENDDO
        ENDIF
     ENDIF
     RETURN
   END FUNCTION L1Br_getSIGlineWL

!! 11: L1Br_getDATAFIELDline
 !
 !    Functionality:
 !
 !       This function will return the data from a single line of any specified dataset
 !       in the L1B file.  This function exists in case any modifications are made to the
 !       L1B file, or in case a PGE requires the use of some of the more obscure datasets
 !       which are not accommodated by the rest of the L1B Reader code.  In particular,
 !       this includes a few of the calibration datasets.  This function is much less
 !       efficient than other L1B Reader functions, so it is advised that it is only
 !       used for datasets not otherwise included in the L1B Reader functions.
 !
 !    Calling Arguments:
 !
 !    Inputs:
 !
 !       filename      the L1B filename
 !       swathname     the L1B swathname
 !       fieldname     the L1B dataset name
 !       Line          the L1B line number to be retrieved
 !
 !    Outputs:
 !
 !       buf           data space to hold L1B values
 !       status        the return PGS_SMF status value
 !
 !    Change History:
 !
 !       Date            Author          Modifications
 !       ====            ======          =============
 !       January 2005    Jeremy Warner   Original Source
 !
!!
      FUNCTION L1Br_getDATAFIELDline(filename, swathname, fieldname, Line, numtype, buf) &
                RESULT(status) 
        USE hdfeos4_parameters
        INCLUDE 'PGS_OMI_1900.f'  !defines the L1B PGS error codes
        INTEGER (KIND = 4), INTENT(IN) :: Line, numtype
        CHARACTER (LEN = *), INTENT(IN) :: filename, swathname, fieldname
        REAL (KIND = 4), INTENT(OUT) :: buf
        CHARACTER (LEN = 256) :: message
        INTEGER (KIND = 4) :: status
        INTEGER (KIND = 4) :: getl1bblk
        EXTERNAL              getl1bblk
        INTEGER (KIND = 4) :: rank
        INTEGER (KIND = 4) :: nL, swid, swfid, fldflg, rLine
        INTEGER (KIND = 4), DIMENSION(1:3) :: dims

        nL = 1
        rLine = Line + 1

        ! Open the L1B swath file

        swfid = swopen(filename, DFACC_READ)
        IF(swfid < zero) THEN
           status = OMI_E_FAILURE
           ierr = OMI_SMF_setmsg(OMI_E_FILE_OPEN, filename, &
                 "L1Br_getDATAFIELDline", zero)
           RETURN
        ENDIF 

        ! Attach to the swath

        swid = swattach(swfid, swathname)
        IF(swid < zero) THEN
           status = OMI_E_FAILURE
           ierr = OMI_SMF_setmsg(OMI_E_SWATH_ATTACH, swathname, &
                 "L1Br_getDATAFIELDline", zero)
           ierr = swclose(swfid)
           RETURN
        ENDIF 
        
        ! Read requested L1B datasets

        status = getl1bblk(swid, fieldname, numtype, fldflg, Line, nL, rank, dims, buf)
        IF(status .NE. OMI_S_SUCCESS) THEN
           ierr = OMI_SMF_setmsg(OMI_E_FAILURE, "Unrecoverable error reading L1B file", &
                                 "L1Br_getDATAFIELDline", zero)
        ENDIF
        IF(fldflg .NE. 0) THEN
           status = OMI_E_FAILURE
	   message = "Unable to read "//fieldname
           ierr = OMI_SMF_setmsg(OMI_E_FAILURE, message, "L1Br_getDATAFIELDline", zero)
        ENDIF

        ! Close, detach, and go home

        ierr = swdetach(swid)
        ierr = swclose(swfid)

        RETURN
      END FUNCTION L1Br_getDATAFIELDline


!! 7. L1Br_getDATAline
 !    This function gets one or more small pixel data field values from
 !    the data block.
 !    this: the block data structure
 !    iLine: the line number in the SmPx swath.  NOTE: this input is 0 based
 !           range from 0 to (nTimes-1) inclusive.
 !
 !    Data_K, Quality_k, and Wavelength_k
 !    are keyword arguments.  Only those present in the argument 
 !    list will be set by the function.
 !
 !    status: the return PGS_SMF status value
!!
      FUNCTION L1Br_getDATAline( this, iLine, nPix, Data_k, &
           Wavelength_k, Quality_k ) & !tpk,       & TrueZoom_k  )   &  !tpk
           RESULT (status)
        USE Szoom_Parameter_Module
        INCLUDE 'PGS_OMI_1900.f'  !!this external file defines the SmPx PGS error codes
        TYPE (L1B_block_type), INTENT( INOUT ) :: this
        INTEGER, PARAMETER :: i1 = 1
        INTEGER (KIND=i1), PARAMETER :: global_mode = 8_i1, szoom_mode = 4_i1
        INTEGER (KIND = 4), INTENT( IN ) :: iLine 
        !tpk INTEGER (KIND = 1), INTENT(IN), OPTIONAL :: TrueZoom_k    !tpk
        REAL (KIND = 4), OPTIONAL, DIMENSION(:,:), INTENT( OUT ) :: & 
                                         Data_k, Wavelength_k
        REAL (KIND = 4), DIMENSION(1:this%nXtrack,1:this%NumTimesSmallPixel) :: & 
                                         tmp_Data_k, tmp_Wavelength_k
        INTEGER (KIND = 2), OPTIONAL, DIMENSION(:), INTENT( OUT ) :: Quality_k
        INTEGER (KIND = 2), DIMENSION(1:this%nXtrack) :: tmp_Quality_k
        INTEGER (KIND = 2), INTENT( OUT ) :: nPix
        INTEGER (KIND = 4) :: i, j, k
        INTEGER :: itmp, ip
        INTEGER (KIND = 4) :: status

        !tpk LOGICAL :: yn_truezoom   !tpk

        status = OMI_S_SUCCESS
        IF( .NOT. this%initialized ) THEN
           ierr = OMI_SMF_setmsg( OMI_E_INPUT, &
                                 "input block not initialized", &
                                 "L1Br_getDATAline", zero )
           status = OMI_E_FAILURE
           RETURN
        ENDIF


        !tpk !tpk -----------------------------------------------------------------
        !tpk !tpk Check for TrueZoom flag (60 cross-track positions instead of 30)
        !tpk !tpk -----------------------------------------------------------------
        !tpk yn_truezoom = .FALSE.
        !tpk IF ( PRESENT (TrueZoom_k) ) THEN
        !tpk    IF (TrueZoom_k > 0_i1 ) yn_truezoom = .TRUE.
        !tpk END IF


        IF( iLine < 0 .OR. iLine >= this%NumTimes) THEN
           ierr = OMI_SMF_setmsg( OMI_E_INPUT, "iLine out of range", &
                                 "L1Br_getDATAline", zero )
           status = OMI_E_FAILURE
           RETURN
        ENDIF

        IF( iLine < this%iLine .OR. iLine > this%eLine ) THEN
           status = fill_l1b_blk( this, iLine )
           IF( status .NE. OMI_S_SUCCESS ) THEN
              ierr = OMI_SMF_setmsg( OMI_E_DATA_BLOCK, "retrieve data block", &
                                    "L1Br_getDATAline", zero )
              RETURN
           ENDIF
        ENDIF

        status = check_blk_line(this, iLine, ip)
        IF(status .NE. OMI_S_SUCCESS) THEN
           ierr = OMI_SMF_setmsg(OMI_E_GENERAL, "Failed retrieving getDATAline data.", &
                                 "L1Br_getDATAline", one)
           RETURN
        ENDIF

        nPix = this%NumSmPix(iLine+1)


        ! Calculate index into small pixel arrays

        i = 1
        do k = 1, iLine
           i = i + this%NumSmPix(k)
        enddo
        j = i + nPix

        IF( PRESENT( Quality_k ) ) THEN
           IF( SIZE( Quality_k ) < this%nXtrack ) THEN 
              ierr = OMI_SMF_setmsg( OMI_E_INPUT, &
                                    "input Quality_k array too small", &
                                    "L1Br_getDATAline", zero )
              status = OMI_E_FAILURE
              RETURN
           ENDIF
           Quality_k(1:this%nXtrack) = &
                             this%PFlag(1:this%nXtrack, iLine-this%iLine+1)
           IF (this%ImBinFact(ip) == szoom_mode) THEN
              IF (this%SpaZoomIndex == 0) THEN
                 !tpk IF (this%ImBinFact(ip) == szoom_mode .AND. (.NOT. yn_truezoom) ) THEN  !tpk
                 DO itmp = 1, this%nXtrack
                    tmp_Quality_k(itmp) = Quality_k(itmp)
                 ENDDO
                 Quality_k(gfpix:zspix) = int16fill
                 Quality_k(zspix+1:zspix+zlpix) =  tmp_Quality_k(zfpix:zlpix)
                 Quality_k(zlpix+zspix+1:glpix) = int16fill
              ENDIF
           ENDIF           
        ENDIF


        IF( PRESENT( Data_k ) ) THEN
           IF( SIZE( Data_k, 2 ) < nPix  .OR. &
               SIZE( Data_k, 1 ) < this%nXtrack ) THEN
              ierr = OMI_SMF_setmsg( OMI_E_INPUT, &
                                    "input Data_k array too small", &
                                    "L1Br_getDATAline", zero )
              status = OMI_E_FAILURE
              RETURN
           ENDIF
           Data_k(1:this%nXtrack, 1:nPix) =  &
                             this%SmPixRad(1:this%nXtrack,i:j)
           IF (this%ImBinFact(ip) == szoom_mode) THEN
              IF (this%SpaZoomIndex == 0) THEN
                 !tpk IF (this%ImBinFact(ip) == szoom_mode .AND. (.NOT. yn_truezoom) ) THEN  !tpk
                 DO itmp = 1, this%nXtrack
                    tmp_Data_k(itmp, 1:nPix) = Data_k(itmp, 1:nPix)
                 ENDDO
                 Data_k(gfpix:zspix, 1:nPix) = rfill
                 Data_k(zspix+1:zspix+zlpix, 1:nPix) =  tmp_Data_k(zfpix:zlpix, 1:nPix)
                 Data_k(zlpix+zspix+1:glpix, 1:nPix) = rfill
              ENDIF
           ENDIF
        ENDIF

        IF( PRESENT( Wavelength_k ) ) THEN
           IF( SIZE( Wavelength_k, 2 ) < nPix .OR. &
               SIZE( Wavelength_k, 1 ) < this%nXtrack ) THEN 
              ierr = OMI_SMF_setmsg( OMI_E_INPUT, &
                                    "input Wavelength_k array too small",&
                                    "L1Br_getDATAline", zero )
              status = OMI_E_FAILURE
              RETURN
           ENDIF
           Wavelength_k(1:this%nXtrack, 1:nPix) = &
                                         this%SmPixWavelen(1:this%nXtrack,i:j)
           IF (this%ImBinFact(ip) == szoom_mode) THEN
              IF (this%SpaZoomIndex == 0) THEN
                 !tpk IF (this%ImBinFact(ip) == szoom_mode .AND. (.NOT. yn_truezoom) ) THEN  !tpk
                 DO itmp = 1, this%nXtrack
                    tmp_Wavelength_k(itmp, 1:nPix) = Wavelength_k(itmp, 1:nPix)
                 ENDDO
                 Wavelength_k(gfpix:zspix, 1:nPix) = rfill
                 Wavelength_k(zspix+1:zspix+zlpix, 1:nPix) =  tmp_Wavelength_k(zfpix:zlpix, 1:nPix)
                 Wavelength_k(zlpix+zspix+1:glpix, 1:nPix) = rfill
              ENDIF
           ENDIF
        ENDIF

        RETURN
      END FUNCTION L1Br_getDATAline

!! 12. L1Br_close
 !
 !    Functionality:
 !
 !       This function should be called when the data block is no longer
 !       needed.  It deallocates all the allocated memory, and sets
 !       all the parameters to invalid values.
 !
 !    Calling Arguments:
 !
 !    Inputs:
 !
 !       this: the block data structure
 !
 !    Outputs:
 !
 !       status: the return PGS_SMF status value
 !
 !    Change History:
 !
 !       Date            Author          Modifications
 !       ====            ======          =============
 !       January 2005    Jeremy Warner   Original Source (derived from original 
 !                                             L1B Reader code written by Kai Yang)
 !
!!
      FUNCTION L1Br_close(this) RESULT(status)
        INCLUDE 'PGS_OMI_1900.f'  !defines the L1B PGS error codes
        TYPE (L1B_block_type), INTENT(INOUT) :: this
        INTEGER (KIND = 4) :: status
        status = OMI_S_SUCCESS

        IF(.NOT. this%initialized) THEN
           ierr = OMI_SMF_setmsg(OMI_E_INPUT, "input block not initialized", &
                                 "L1Br_close", zero)
           status = OMI_E_FAILURE
           RETURN
        ENDIF

        ! DEALLOCATE all the memory
        IF (ASSOCIATED(this%Time)) DEALLOCATE(this%Time)
        IF (ASSOCIATED(this%SecInDay)) DEALLOCATE(this%SecInDay)
        IF (ASSOCIATED(this%ScLat)) DEALLOCATE(this%ScLat)
        IF (ASSOCIATED(this%ScLon)) DEALLOCATE(this%ScLon)
        IF (ASSOCIATED(this%ScAlt)) DEALLOCATE(this%ScAlt)
        IF (ASSOCIATED(this%SolElevation)) DEALLOCATE(this%SolElevation)
        IF (ASSOCIATED(this%SolElevMin)) DEALLOCATE(this%SolElevMin)
        IF (ASSOCIATED(this%SolElevMax)) DEALLOCATE(this%SolElevMax)
        IF (ASSOCIATED(this%SolAzimuth)) DEALLOCATE(this%SolAzimuth)
        IF (ASSOCIATED(this%SolAziMin)) DEALLOCATE(this%SolAziMin)
        IF (ASSOCIATED(this%SolAziMax)) DEALLOCATE(this%SolAziMax)
        IF (ASSOCIATED(this%Lon)) DEALLOCATE(this%Lon)
        IF (ASSOCIATED(this%Lat)) DEALLOCATE(this%Lat)
        IF (ASSOCIATED(this%SolZenAng)) DEALLOCATE(this%SolZenAng)
        IF (ASSOCIATED(this%SolAziAng)) DEALLOCATE(this%SolAziAng)
        IF (ASSOCIATED(this%ViewZenAng)) DEALLOCATE(this%ViewZenAng)
        IF (ASSOCIATED(this%ViewAziAng)) DEALLOCATE(this%ViewAziAng)
        IF (ASSOCIATED(this%TerrainHeight)) DEALLOCATE(this%TerrainHeight)
        IF (ASSOCIATED(this%GPQFlag)) DEALLOCATE(this%GPQFlag)
        IF (ASSOCIATED(this%XTQFlag)) DEALLOCATE(this%XTQFlag)
        IF (ASSOCIATED(this%RadMantissa)) DEALLOCATE(this%RadMantissa)
        IF (ASSOCIATED(this%RadPrecision)) DEALLOCATE(this%RadPrecision)
        IF (ASSOCIATED(this%PQFlag)) DEALLOCATE(this%PQFlag)
        IF (ASSOCIATED(this%RadExponent)) DEALLOCATE(this%RadExponent)
        IF (ASSOCIATED(this%SmPixRad)) DEALLOCATE(this%SmPixRad)
        IF (ASSOCIATED(this%SmPixWavelen)) DEALLOCATE(this%SmPixWavelen)
        IF (ASSOCIATED(this%WavelenCoef)) DEALLOCATE(this%WavelenCoef)
        IF (ASSOCIATED(this%WavelenPrec)) DEALLOCATE(this%WavelenPrec)
        IF (ASSOCIATED(this%WavelenRefCol)) DEALLOCATE(this%WavelenRefCol)
        IF (ASSOCIATED(this%Config)) DEALLOCATE(this%Config)
        IF (ASSOCIATED(this%MeasClass)) DEALLOCATE(this%MeasClass)
        IF (ASSOCIATED(this%NumSmPixCol)) DEALLOCATE(this%NumSmPixCol)
        IF (ASSOCIATED(this%ExposureType)) DEALLOCATE(this%ExposureType)
        IF (ASSOCIATED(this%ImBinFact)) DEALLOCATE(this%ImBinFact)
        IF (ASSOCIATED(this%GC1)) DEALLOCATE(this%GC1)
        IF (ASSOCIATED(this%GC2)) DEALLOCATE(this%GC2)
        IF (ASSOCIATED(this%GC3)) DEALLOCATE(this%GC3)
        IF (ASSOCIATED(this%GC4)) DEALLOCATE(this%GC4)
        IF (ASSOCIATED(this%DSGC)) DEALLOCATE(this%DSGC)
        IF (ASSOCIATED(this%LSLABF)) DEALLOCATE(this%LSLABF)
        IF (ASSOCIATED(this%USLABF)) DEALLOCATE(this%USLABF)
        IF (ASSOCIATED(this%LDABF)) DEALLOCATE(this%LDABF)
        IF (ASSOCIATED(this%UDABF)) DEALLOCATE(this%UDABF)
        IF (ASSOCIATED(this%MQFlag)) DEALLOCATE(this%MQFlag)
        IF (ASSOCIATED(this%CalSet)) DEALLOCATE(this%CalSet)
        IF (ASSOCIATED(this%SmPixCol)) DEALLOCATE(this%SmPixCol)
        IF (ASSOCIATED(this%GSC1)) DEALLOCATE(this%GSC1)
        IF (ASSOCIATED(this%GSC2)) DEALLOCATE(this%GSC2)
        IF (ASSOCIATED(this%GSC3)) DEALLOCATE(this%GSC3)
        IF (ASSOCIATED(this%SR1)) DEALLOCATE(this%SR1)
        IF (ASSOCIATED(this%SR2)) DEALLOCATE(this%SR2)
        IF (ASSOCIATED(this%SR3)) DEALLOCATE(this%SR3)
        IF (ASSOCIATED(this%SR4)) DEALLOCATE(this%SR4)
        IF (ASSOCIATED(this%BinImgRows)) DEALLOCATE(this%BinImgRows)
        IF (ASSOCIATED(this%StopColumn)) DEALLOCATE(this%StopColumn)
        IF (ASSOCIATED(this%MasterClkPer)) DEALLOCATE(this%MasterClkPer)
        IF (ASSOCIATED(this%ExpTime)) DEALLOCATE(this%ExpTime)
        IF (ASSOCIATED(this%ReadTime)) DEALLOCATE(this%ReadTime)
        IF (ASSOCIATED(this%DetcTemp)) DEALLOCATE(this%DetcTemp)
        IF (ASSOCIATED(this%OptBenchTemp)) DEALLOCATE(this%OptBenchTemp)

        this%iLine              = -1
        this%nLine              = -1
        this%eLine              = -1
        this%NumTimes           = -1
        this%NumTimesSmallPixel = -1
        this%nXtrack            = -1
        this%nWavel             = -1     
        this%nWavelCoef         = -1     
        this%filename           = ""
        this%swathname          = ""
        this%initialized        = .FALSE.
        RETURN
      END FUNCTION L1Br_close

END MODULE L1B_Reader_class
