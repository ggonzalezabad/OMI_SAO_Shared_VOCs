#ifndef omismf_h
#define omismf_h

/* SDP Toolkit Status Message Facility header file */
#include "PGS_SMF.h"
#include "PGS_PC.h"

#define SEVERITYCODE_LUN         200100
#define SCODE_THRESHOLD_DEFAULT  3
#define SCODE_THRESHOLD_MAX      256

/* Prototypes */
PGSt_integer
OMI_SMF_setmsg( PGSt_SMF_code  mnemonicstring,
                char          *infostring,
                char          *functionstring,
                PGSt_integer   severity_code );

#endif 
