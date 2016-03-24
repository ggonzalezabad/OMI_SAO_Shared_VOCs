#include "PGS_OMI_1900.h"
#include "l1b_reader.h"
#include "hdf.h"
#include "cfortHdf.h"
#include "HdfEosDef.h"
#include "omi_smf.h"

#define FUNCTION_NAME "OMI_Get_L1B_Data_Block()"
#define MAX_LEN_DIM_LIST 256
#define MAX_RANK         4

PGSt_SMF_status
OMI_Get_L1B_Data_Block( int   SWid_he4,       /* input  */
                        char *fieldname,      /* input  */
                        int   numberType,     /* input  */
                        int  *fldflag,        /* output  */
                        int   Line_start,     /* input  */
                        int   nLines,         /* input  */
                        int  *rank,           /* output */
                        int  *dims,           /* output */
                        void *data_buffer )   /* output */
{  intn   status_he4= FAIL; 
   int32  numberType_local;
   int32  rank_local;
   int32  dims_local[MAX_RANK];
   int32  start[MAX_RANK] = {0,0,0,0};
   int32  edge[MAX_RANK];
   int    di;
   static PGSt_integer s_code = 9;
   static PGSt_integer e_code = 0;
   char   dimlist[MAX_LEN_DIM_LIST] = {'\0'};
   char   msg[PGS_SMF_MAX_MSG_SIZE];  /* PGS_SMF_MAX_MSG_SIZE is
                                         defined in PGS_SMF.h */

   *fldflag = 0;
   if( Line_start < 0 || nLines < 0 )
   {  sprintf( msg, "Line_start = %d, nLines = %d.", Line_start, nLines );
      OMI_SMF_setmsg( OMI_E_VALID_RANGE, msg, FUNCTION_NAME, e_code );
      return OMI_E_FAILURE;
   }
   
   status_he4 = SWfieldinfo( SWid_he4, fieldname, 
                            &rank_local, dims_local,
                            &numberType_local, dimlist );
   if( status_he4 == FAIL )
   {  sprintf( msg, 
              "%s: retrieve field info for %s failed.",
              "SWfieldinfo()", fieldname);
      OMI_SMF_setmsg( OMI_E_HDFEOS, msg, FUNCTION_NAME, s_code );
      *fldflag = 1;
      return OMI_S_SUCCESS;
   }

   if( numberType_local != (int32) numberType )
   {  sprintf( msg, "pre-set number type different from that in L1B file." ); 
      OMI_SMF_setmsg( OMI_E_INPUT, msg, FUNCTION_NAME, e_code );
      return OMI_E_FAILURE;
   }

   *rank = rank_local;

   for( di = 0; di < rank_local; di++ )
   {  dims[di] = (int) dims_local[di];
      edge[di] = (int) dims_local[di];
   }

   if( Line_start >= dims[0] || (Line_start + nLines) > dims[0] )
   {  sprintf( msg, "Line_start = %d, Line_start+ nLines = %d, dims[0] = %d.", 
               Line_start, Line_start + nLines, dims[0] );
      OMI_SMF_setmsg( OMI_E_VALID_RANGE, msg, FUNCTION_NAME, e_code );
      return OMI_E_FAILURE;
   }

   if( nLines > 0 )
   {  start[0] = Line_start;
      edge[0]  = nLines;
      status_he4 = SWreadfield( SWid_he4, fieldname, 
                                start, NULL, edge,
                                (VOIDP) data_buffer );
      if( status_he4 == FAIL )
      {  sprintf( msg, 
              "%s: read field data for %s failed", "SWreadfield()", fieldname);
         OMI_SMF_setmsg( OMI_E_HDFEOS, msg, FUNCTION_NAME, s_code );
         return OMI_S_SUCCESS;
   }  }

   return OMI_S_SUCCESS;
}

/* HDF types used in FORTRAN bindings */

#if defined(DEC_ALPHA) || defined(IRIX) || defined(UNICOS)

#define INT32  INT
#define INT32V INTV
#define PINT32 PINT

#else

#define INT32  LONG
#define INT32V LONGV
#define PINT32 PLONG

#endif

/* FORTRAN bindings */

FCALLSCFUN9(INT, OMI_Get_L1B_Data_Block, GETL1BBLK, getl1bblk, INT, STRING,
            INT, PINT, INT, INT, PINT, INTV, PVOID )
