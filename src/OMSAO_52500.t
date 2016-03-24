# Status Message Text for Toolkit --- OMI SAO PGE
#
# Error messages for input module
#
# Author:                    T.P. Kurosu
# Co-Author:                 G. Gonzalez Abad
# Date of last modification: 29 September 2012
# Version:                   1.0
# History: 
#
%INSTR = OMI
%LABEL = OMSAO
%SEED  = 52500

OMSAO_W_GETLUN    WARNING...failed to read LUN from PCF
OMSAO_E_GETLUN    ERROR...failed to read LUN from PCF
OMSAO_F_GETLUN    FATAL ERROR...failed to read LUN from PCF. PGE aborting with exit code 1

OMSAO_F_GET_MOLINDEX      FATAL ERROR...failed to convert molecule ID to index. PGE aborting with exit code 1
OMSAO_S_GET_MOLINDEX      SUCCESS...identified molecule ID for current PGE run

OMSAO_F_OPEN_FITCTRL_FILE  FATAL ERROR...failed to open fitting control file. PGE aborting with exit code 1
OMSAO_F_READ_FITCTRL_FILE  FATAL ERROR...failed during READ from fitting control file. PGE aborting with exit code 1
OMSAO_W_CLOSE_FITCTRL_FILE WARNING...failed to close fitting control file
OMSAO_S_READ_FITCTRL_FILE  SUCCESS...completed READ from fitting control file

OMSAO_F_GET_MOLFITNAME     FATAL ERROR...failed to convert mol-name to fit-index. PGE aborting with exit code 1

OMSAO_E_OPEN_REFSPEC_FILE  ERROR...failed to open reference spectrum file
OMSAO_E_READ_REFSPEC_FILE  ERROR...failed to READ from reference spectrum file
OMSAO_W_CLOSE_REFSPEC_FILE WARNING...failed to close reference spectrum file
OMSAO_S_READ_REFSPEC_FILE  SUCCESS...completed READ of all reference spectra
OMSAO_E_REFSPEC_MAXPTS     ERROR...reference spectrum exceeds maximum number of spectral points

OMSAO_W_READ_AMFTABLE_NAME  WARNING...failed to read name of AMF table file from PCF
OMSAO_W_OPEN_AMFTABLE_FILE  WARNING...failed to open AMF table file
OMSAO_W_READ_AMFTABLE_FILE  WARNING...failed during READ from AMF table file
OMSAO_W_CLOSE_AMFTABLE_FILE WARNING...failed to close AMF table file
OMSAO_S_READ_AMFTABLE_FILE  SUCCESS...completed READ of AMF table file
OMSAO_W_AMFTABLE_MEMALLOC   WARNING...failed to allocate memory for AMF table

OMSAO_E_OPEN_RANCIL_FILE   ERROR...failed to open ancillary file for reading
OMSAO_S_OPEN_RANCIL_FILE   SUCCESS...opened ancillary file for reading
OMSAO_E_OPEN_WANCIL_FILE   ERROR...failed to open ancillary file for writing
OMSAO_S_OPEN_WANCIL_FILE   SUCCESS...openened ancillary file for writing
OMSAO_E_OPEN_UANCIL_FILE   ERROR...unknown action in OPEN of ancillary file
OMSAO_E_UACTION_ANCIL_FILE ERROR...unknown request to ancillary file
OMSAO_W_CLOSE_ANCIL_FILE   WARNING...failed to close ancillary file
OMSAO_S_CLOSE_ANCIL_FILE   SUCCESS...closed ancillary file

OMSAO_F_OPEN_L1B_FILE      FATAL ERROR...failed to open L1b file. PGE aborting with exit code 1
OMSAO_E_READ_L1B_FILE      ERROR...failed to read from L1b file
OMSAO_W_READ_L1B_FILE      WARNING...failed to read from L1b file
OMSAO_W_CLOS_L1B_FILE      WARNING...failed to close OMI L1B file

OMSAO_F_MIXMODE            FATAL ERROR...unable to handle mixed GLOBAL and SPATIAL ZOOM modes. PGE aborting with exit code 1
OMSAO_F_SPAZOOM_PIX        FATAL ERROR...cross-track pixel out of bounds for spatial zoom mode. PGE aborting with exit code 1

OMSAO_F_MEM_ALLOC          FATAL ERROR...memory allocation failure. PGE aborting with exit code 1
OMSAO_E_MEM_DALLOC         ERROR...memory de-allocation failure
OMSAO_W_MEMALLOC           WARNING...memory allocation failure

OMSAO_A_INTERPOL           ALERT...interpolation failed. 
OMSAO_E_INTERPOL           ERROR...interpolation failed
OMSAO_W_INTERPOL           WARNING...interpolation failed
OMSAO_E_INTERPOL_REFSPEC   ERROR...interpolation failed for reference spectum
OMSAO_W_INTERPOL_RANGE     WARNING...interpolation found incomplete wavelength range

OMSAO_E_HE5SWOPEN          ERROR...failed to open HE5 output file
OMSAO_F_HE5SWOPEN          ERROR...failed to open HE5 output file. PGE aborting with exit code 1
OMSAO_F_HE5SWCREATE        FATAL ERROR...failed to create HE5 swath. PGE aborting with exit code 1
OMSAO_F_HE5SWDEFDIM        FATAL ERROR...failed to define HE5 swath dimensions. PGE aborting with exit code 1
OMSAO_E_HE5SWDEFFLD        ERROR...failed to define HE5 swath geo/data fields
OMSAO_E_HE5SWATTACH        ERROR...failed to attach to existing HE5 swath
OMSAO_E_HE5SWLOCATE        ERROR...failed to locate HE5 swath
OMSAO_E_HE5SWWRFLD         ERROR...failed to write to HE5 data field
OMSAO_E_HE5SWRDFLD         ERROR...failed to read from HE5 data field
OMSAO_W_HE5SWCLOSE         WARNING...failed to close HE5 output file
OMSAO_E_HE5SWWRATTR        ERROR...failed to write HE5 swath attribute
OMSAO_W_HE5SWRLATTR        WARNING...failed to write HE5 local attribute
OMSAO_W_HE5EHWRGLATT       WARNING...failed to write HE5 global attribute

OMSAO_F_METINIT            FATAL ERROR...failed to initialize MCF. PGE aborting with exit code 1
OMSAO_E_GETATTR            ERROR...failed to get MetaData Attribute
OMSAO_E_METINIT            ERROR...failed to initialize MCF
OMSAO_E_GETREF             ERROR...failed to get L2 Input Pointer
OMSAO_E_MDMISMATCH         ERROR...mismatch for Metadata value field
OMSAO_E_LGIDCOMP           ERROR...could not compose LocalGranuleID
OMSAO_E_MDL2INV            ERROR...failed to set L2 inventory Metadata
OMSAO_E_MDL2ARC            ERROR...failed to set L2 archived Metadata
OMSAO_E_MDINIT             ERROR...failed to initiate L2 Metadata
OMSAO_W_AMDWRT             WARNING...problems encountered during write of ArchivedMetaData
OMSAO_W_CMDWRT             WARNING...problems encountered during write of CoreMetaData
OMSAO_W_MDMISMATCH         WARNING...mismatch for Metadata value field
OMSAO_W_MDL2INV            WARNING...failed to set L2 inventory Metadata
OMSAO_W_MDL2ARC            WARNING...failed to set L2 archived Metadata
OMSAO_W_GETATTR            WARNING...failed to get MetaData Attribute

OMSAO_W_SKIPPIX            WARNING...unable to process cross-track pixel
OMSAO_W_NOPIXEL            WARNING...no valid cross-track positions to process

OMSAO_E_PREFITCOL          ERROR...failed to access ancillary L2 infomation
OMSAO_E_PREFITDIM          ERROR...dimension mismatch for prefitted columns
OMSAO_W_PREFITDIM          WARNING...dimension mismatch for prefitted columns

OMSAO_W_TAI93              WARNING...failed to compute TAI time at 0z of granule

OMSAO_S_PROGRESS           PGE PROGRESS Message...
OMSAO_S_ENDOFRUN           END_OF_RUN_S...PGE execution completed normally, exit code 0
OMSAO_W_ENDOFRUN           END_OF_RUN_W...PGE execution completed with warnings
OMSAO_E_ENDOFRUN           END_OF_RUN_E...PGE execution completed with errors
OMSAO_F_ENDOFRUN           END_OF_RUN_F...PGE execution terminated abnormally, exit code 1
OMSAO_U_ENDOFRUN           END_OF_RUN_U...PGE execution completed with unknown exit code

OMSAO_A_SUBROUTINE         ALERT...non-zero exit status returned from subroutine
OMSAO_W_SUBROUTINE         WARNING... warning-status returned from subroutine

OMSAO_F_XTRMISRAD          FATAL ERROR...Cross-track position different between L1b radiance granules

OMSAO_F_SOLCOM_VS_SOLAVE   FATAL ERROR...Check control file, both solar composite and solar average set true

OMSAO_E_OPEN_SOLMONAVE_FILE  ERROR...failed to open solar monthly average spectrum file
OMSAO_E_READ_SOLMONAVE_FILE  ERROR...failed to READ from solar monthly average spectrum file
OMSAO_W_CLOSE_SOLMONAVE_FILE WARNING...failed to close solar monthly average spectrum file
OMSAO_S_READ_SOLMONAVE_FILE  SUCCESS...completed READ of solar monthly average spectrum file

OMSAO_E_OPEN_REFSECCOR_FILE  ERROR...failed to open reference sector corrrection file
OMSAO_E_READ_REFSECCOR_FILE  ERROR...failed to READ from reference sector corrrection file
OMSAO_W_CLOSE_REFSECCOR_FILE WARNING...failed to close reference sector corrrection file
OMSAO_S_READ_REFSECCOR_FILE  SUCCESS...completed READ of reference sector corrrection file

OMSAO_E_HE5GDOPEN    ERROR...failed to open GRID file
OMSAO_W_HE5GDCLOSE   WARNING...failed to close GRID file
OMSAO_E_HE5GDATTACH  ERROR...failed to attach to GRID
