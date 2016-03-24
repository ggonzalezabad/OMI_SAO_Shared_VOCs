#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <cfortran.h>
#include "omi_smf.h"  /* include PGS_PC.h PGS_SMF.h */

/* severity_code threshold is initialized to -1, so it will be  
   set when this function is first called.
   0 = minimal reporting  
   1 = reporting appropriate for I&T debug 
   2 and higher - reporting for debug by STM
*/
static PGSt_integer Scode_threshold = -1; 

PGSt_integer
OMI_SMF_setmsg( PGSt_SMF_code  mnemonicstring, 
                char          *infostring, 
                char          *functionstring,
                PGSt_integer   severity_code )
/*****************************************************************************
!C

!Description:

  SMF module,
  this module writes errors and messages through the SDP toolkit to
  LogStatus and causes program termination upon encountering a fatal
  error (as determined by the error mnemonic).  As input it acccepts an
  error mnemonic, string containing user information, a string
  containing the function name and module where the error/message
  occurs, and an integer denotes the importance of the error/message.

!Input Parameters:
  PGSt_SMF_code mnemonicstring      SMF message mnemonic
  char *infostring                  String containing info to be written to 
                                    LogStatus
  char *functionstring              String containing function/module name
                                    to be written to LogStatus
  PGSt_integer severity_code        The higher the value of severity_code, 
                                    the lower the importance of the
                                    error/message.
!Output Parameters:
  None
  
!Return
  0  successful return
 
!Revision History:
  Revision 0.1  11/26/2001  Kai Yang/SSAI
  Adopted smf.c from MODIS code.

!Team-unique Header:
  This software was developed by the OMI Science Team Support
  Group for the National Aeronautics and Space Administration, Goddard
  Space Flight Center, under NASA Task 916-003-1

!References and Credits
  Written by 
  Kai Yang 
  Science Systems and Applications, Inc.
  email: kyang@ltpmail.gsfc.nasa.gov
  
!Design Notes
  Scode_threshold is read from the PCF file though logical number
  SEVERITYCODE_LUN defined in omi_smf.h 

  The determination of time is through native C function(s).

  Upon determination of a fatal error, exit code 1 is used.

  0 return is supplied here to complete Fortran call, no need
  for checking of this return.

  Designed to run with seed file format of:
    OMI_X_MNEMONIC STRING %sStatic Message.  %s

  Example of function call:

  C:
  OMI_SMF_setmsg(OMI_X_MNEMONIC_STRING, "user message string", 
                                        "function, module.c string",
                                         severity_code );

  FORTRAN:
  ierr = OMI_SMF_setmsg(OMI_X_MNEMONIC_STRING, "user message string", 
                                               "function, module.c string",
                                                severity_code );

  -- severity_code is an integer reflecting the importance of the message.
  -- The mnemonic string must be identical to one contained in an included 
       seed file.
  -- The user message string can contain anything that needs to be written
       to the log(s).
  -- The function, module.c string should contain "the name of the function,
       the name of the module".

!END
*****************************************************************************/
{
  struct tm *local    = NULL;/* Time functions */
  time_t t;

  PGSt_SMF_status ret = PGS_S_SUCCESS;  /* Holds string of message from 
                                           PGS_SMF_GetMsgByCode */
  char message[PGS_SMF_MAX_MSGBUF_SIZE]={'\0'};    
                                        /* holds the message string associated 
                                           with the error code returned by 
                                           GetMsg */
  char buf[PGS_SMF_MAX_MSGBUF_SIZE] = {'\0'};  
                                        /* Holds a concatenated dynamic 
                                           message string */
  char pval[PGSd_PC_VALUE_LENGTH_MAX];  /* used as buffer to get params 
                                           from pcf */

  /* Scode_threshold is not set, i.e., Scode_threshold = -1,
     need to read the Scode_threshold value from PCF */
  if( Scode_threshold == -1 )
  {  if( PGS_PC_GetConfigData( SEVERITYCODE_LUN, pval ) != PGS_S_SUCCESS )
     { Scode_threshold = SCODE_THRESHOLD_DEFAULT;
       sprintf( buf, "Scode_threshold is set to default %d", 
                      SCODE_THRESHOLD_DEFAULT);
       PGS_SMF_SetDynamicMsg( PGSPC_W_NO_CONFIG_VALUE, buf, "OMI_SMF_setmsg" );
     }
     else
       sscanf( pval,"%d",&Scode_threshold );
     
     if( Scode_threshold < 0 || Scode_threshold > SCODE_THRESHOLD_MAX )
       Scode_threshold = SCODE_THRESHOLD_DEFAULT;
  }

  if( severity_code >= Scode_threshold )  /* no message required, success */
     return 0;
   
  t     = time(NULL);                   /* Calculate timestamp */
  local = (struct tm *) localtime( &t );

  /* Get message from seed file based on error mnemonic */
  ret = PGS_SMF_GetMsgByCode( mnemonicstring, message );
  if( ret != PGS_S_SUCCESS ) 
     exit(1);

  /* Concatenate timestamp, and user massage with buffer message */
  sprintf( buf, "%s %s. %s", asctime(local), message, infostring );
  
  /* Write message to LogStatus */
  ret = PGS_SMF_SetDynamicMsg( mnemonicstring, buf, functionstring );
  if( ret != PGS_S_SUCCESS ) 
     exit(1);

  /* Test for fatal error */
  if( PGS_SMF_TestFatalLevel( mnemonicstring ) == PGS_TRUE )
     exit(1);  
  else 
     return 0;
}

/* FORTRAN bindings */

FCALLSCFUN4( INT, OMI_SMF_setmsg, OMI_SMF_SETMSG, omi_smf_setmsg, INT, STRING, STRING, INT )
