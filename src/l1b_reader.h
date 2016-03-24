#include <PGS_SMF.h>

PGSt_SMF_status
OMI_Get_L1B_Data_Block( int   SWid_he4,       /* input  */
                        char *fieldname,      /* input  */
                        int   numberType,     /* input  */
                        int  *fldflag,        /* output  */
                        int   Line_start,     /* input  */
                        int   nLines,         /* input  */
                        int  *rank,           /* output */
                        int  *dims,           /* output */
                        void *data_buffer );  /* output */
