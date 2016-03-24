! ================================================================================
!
! Collection of subroutines that relate to PGE identification and processing aids.
!
! ================================================================================

SUBROUTINE get_pge_ident ( tmpchar, pge_name, pge_idx, errstat )

  ! ====================================================
  ! Find name and index of current PGE from input string
  ! ====================================================

  USE OMSAO_indices_module, ONLY: n_sao_pge, sao_pge_names, sao_pge_min_idx, sao_pge_max_idx
  USE OMSAO_errstat_module

  IMPLICIT NONE

  ! --------------
  ! Input variable
  ! --------------
  CHARACTER (LEN=*), INTENT (IN) :: tmpchar

  ! ----------------
  ! Output variables
  ! ----------------
  CHARACTER (LEN=*),   INTENT (OUT)   :: pge_name
  INTEGER   (KIND=i4), INTENT (OUT)   :: pge_idx
  INTEGER   (KIND=i4), INTENT (INOUT) :: errstat

  ! --------------
  ! Local variable
  ! --------------
  INTEGER (KIND=i4) :: i, locerrstat

  ! ---------------------------
  ! Initialize output variables
  ! ---------------------------
  pge_idx = -1  ;  pge_name(1:1) = '?' ; locerrstat = pge_errstat_ok

  ! -------------------------------------------------
  ! Find name and index by looping over all SAO PGEs.
  ! Not very elegant but simple and effective.
  ! -------------------------------------------------
  getpge: DO i = sao_pge_min_idx, sao_pge_max_idx
     IF ( TRIM(ADJUSTL(tmpchar)) == TRIM(ADJUSTL(sao_pge_names(i))) ) THEN
        pge_name = sao_pge_names(i)  ;  pge_idx = i
        EXIT getpge
     END IF
  END DO getpge
  IF ( pge_idx == -1 .OR.  pge_name(1:1) == '?' ) locerrstat = pge_errstat_error

  errstat = MAX ( errstat, locerrstat )
  RETURN
END SUBROUTINE get_pge_ident


SUBROUTINE pge_error_message ( errstat, ok_msg, warn_msg, err_msg )

  ! ====================================================================
  ! Check PGE error status and report appropriate message. STOP on error
  ! ====================================================================

  USE OMSAO_errstat_module

  IMPLICIT NONE

  ! ---------------
  ! Input variables
  ! ---------------
  INTEGER   (KIND=i4), INTENT (IN) :: errstat
  CHARACTER (LEN=*),   INTENT (IN) :: ok_msg, warn_msg, err_msg

  SELECT CASE ( errstat )
  CASE ( pge_errstat_ok )
     WRITE (*, '(A)') TRIM(ADJUSTL(ok_msg))
  CASE ( pge_errstat_warning )
     WRITE (*, '(A,A)') 'WARNING: ', TRIM(ADJUSTL(warn_msg))
  CASE ( pge_errstat_error )
     WRITE (*, '(A,A)') 'ERROR: ', TRIM(ADJUSTL(err_msg))
     STOP 1
  END SELECT

  RETURN
END SUBROUTINE pge_error_message

FUNCTION int2string ( int_num, ni ) RESULT ( int_str )

  ! ===============================================================
  ! Converts INTEGER number INT_NUM to STRING INT_STR of length NI
  ! or the number of digits in INT_NUM, whatever is larger.
  ! ===============================================================

  USE OMSAO_precision_module, ONLY: i4
  IMPLICIT NONE

  ! ---------------
  ! Input variables
  ! ---------------
  INTEGER (KIND=i4),  INTENT (IN)  :: int_num, ni

  ! ---------------
  ! Result variable
  ! ---------------
  CHARACTER (LEN=*) :: int_str

  ! ------------------------------
  ! Local variables and parameters
  ! ------------------------------
  ! * Arrays containing indices for number strings in ASCII table
  INTEGER (KIND=i4),                   PARAMETER :: n = 10
  INTEGER (KIND=i4), DIMENSION(0:n-1)            :: aidx
  CHARACTER (LEN=1), DIMENSION(0:n-1), PARAMETER :: astr = (/ "0","1","2","3","4","5","6","7","8","9" /)
  ! * Temporary and loop variables
  INTEGER (KIND=i4)                              :: i, k, nd, tmpint, ld

  ! ----------------------------------------------------
  ! Compute the index entries of ASTR in the ASCII table
  ! ----------------------------------------------------
  aidx = IACHAR(astr)

  ! ----------------------------------------------------------
  ! Find the number of digits in INT_NUM. This is equal to the 
  ! truncated integer of LOG10(INT_NUM) plus 1.
  ! ----------------------------------------------------------
  SELECT CASE ( int_num )
  CASE ( 0:9 ) 
     nd = 1
  CASE ( 10: )
     nd = INT ( LOG10(REAL(int_num)) ) + 1
  CASE DEFAULT
     int_str = '?' ; RETURN
  END SELECT

  ! -------------------------------------------------------------
  ! We may want to create a string that is longer than the number
  ! of digits in INT_NUM. This will create leading "0"s.
  ! -------------------------------------------------------------
  nd = MAX ( nd, ni )

  ! ----------------------------------------------
  ! Convert the integer to string "digit by digit"
  ! ----------------------------------------------
  int_str = "" ; tmpint = int_num
  DO i = 1, nd
     ld = MOD ( tmpint, 10 )        ! Current last digit
     tmpint = ( tmpint - ld ) / 10  ! Remove current last digit from INT_STR
     k = nd - i + 1                 ! Position of current digit in INT_STR
     int_str(k:k) = ACHAR(aidx(ld)) ! Convert INTEGER digit to CHAR
  END DO

  RETURN
END FUNCTION int2string


SUBROUTINE utc_julian_date_and_time ( year, month, day, julday, hour, minute, second )

  ! -------------------------------------------------------------
  ! Converts current date and time to UTC time and "Julian" date,
  ! i.e., day of the year.
  ! -------------------------------------------------------------

  USE OMSAO_precision_module, ONLY: i4, r4
  IMPLICIT NONE

  ! ---------------
  ! RESULT
  ! ---------------
  INTEGER (KIND=i4), INTENT (OUT) :: year, month, day, julday, hour, minute, second

  !INTEGER (KIND=i4), DIMENSION(8), INTENT (OUT) :: utc_juldate
  
  ! ---------------
  ! Local variables
  ! ---------------
  ! * Arguments for DATE_AND_TIME, and position indices for Year, Month, Day, Zone, Hour, Minute, Seconds
  CHARACTER (LEN= 8)                 :: date
  CHARACTER (LEN=10)                 :: time
  CHARACTER (LEN= 5)                 :: zone
  INTEGER   (KIND=i4), DIMENSION (8) :: date_vector
  INTEGER   (KIND=i4), PARAMETER     :: y_idx=1, m_idx=2, d_idx=3, z_idx=4, hh_idx=5, mm_idx=6, ss_idx=7
  ! * Some MAX values for Minutes and Hours (not that we would expect those to change :-)
  INTEGER   (KIND=i4), PARAMETER     :: max_mi = 60, max_hr = 24, max_dy = 31, max_mo = 12
  ! * Other local variables
  INTEGER   (KIND=i4)                :: del_mm, del_hh, max_julday, yyy, mmm, ddd

  ! ------------------
  ! External functions
  ! ------------------
  INTEGER (KIND=i4), EXTERNAL :: day_of_year

  ! ----------------------------------------------------------
  ! First copy the input variable to the output variable. That
  ! saves work on premature return.
  ! ----------------------------------------------------------
  CALL DATE_AND_TIME ( date, time, zone, date_vector )

  ! -----------------------------------------------
  ! Initialize all output variables (except JULDAY)
  ! -----------------------------------------------
  year   = date_vector(y_idx )
  month  = date_vector(m_idx )
  day    = date_vector(d_idx )
  hour   = date_vector(hh_idx)
  minute = date_vector(mm_idx)
  second = date_vector(ss_idx)

  ! ----------------------------------------------------------------------
  ! Next, find the day of the year. At this point, we also compute the
  ! maximum number of days per year. This will come in handy further down.
  ! ----------------------------------------------------------------------
  julday     = day_of_year ( year, month,  day    )
  max_julday = day_of_year ( year, max_mo, max_dy )

  ! -----------------------------------------------
  ! Check if there is any difference to UTC at all.
  ! -----------------------------------------------
  IF ( ABS(date_vector(z_idx)) == 0 ) RETURN

  ! ------------------------------------------------------------
  ! Apply difference in UTC time to local time vector. This is
  ! given in minutes, thus apply to minutes index. 
  !
  ! REMEMBER: This is the DIFFERENCE to UTC, so we have to
  !           SUBTRACT whatever difference from the current time.
  ! -------------------------------------------------------------
  del_mm = MOD ( date_vector(z_idx), max_mi )
  minute = minute - del_mm

  ! ----------------------------------------
  ! Compute hour difference, return if none.
  ! ----------------------------------------
  del_hh =   ( date_vector(z_idx) - del_mm ) / max_mi
  IF ( del_hh == 0 ) RETURN

  ! --------------------------------------------------------
  ! Apply hour difference and check for/apply day difference
  ! --------------------------------------------------------
  hour = hour - del_hh

  ! -------------------------------------------------------------
  ! If we reach here, we may have had a day difference. We could
  ! make our lives miserable by checking for all the months of 
  ! the year. But what we ultimately want is the day of the year
  ! (a proxy for Julian day), and that makes life a lot easier.
  ! -------------------------------------------------------------
  SELECT CASE ( hour )
  CASE ( 0:23 )
     RETURN
  CASE ( :-1 )
     hour = max_hr + hour
     day  = day - 1
  CASE ( 24: )
     hour = hour - max_hr
     day  = day  + 1
  END SELECT
  
  ! --------------------------------------------------------------
  ! Finally, if we ever reach here, we've had a year difference.
  ! Apply correction, compute the new day of the year, and return.
  ! --------------------------------------------------------------
  yyy = year ; mmm = month ; ddd = day                     ! save old values
  CALL year_month_day ( yyy, mmm, ddd, year, month, day )  ! overwrite with new ones
  julday = day_of_year ( year, month, day )                ! recompute day-of-year

  RETURN
END SUBROUTINE utc_julian_date_and_time

FUNCTION day_of_year ( year, month, day ) RESULT ( jday )

  USE OMSAO_precision_module, ONLY: i4
  IMPLICIT NONE

  ! ---------------
  ! Input variables
  ! ---------------
  INTEGER (KIND=i4), INTENT (IN) :: year, month, day

  ! ---------------
  ! Result variable
  ! ---------------
  INTEGER (KIND=i4) :: jday

  ! ------------------------------
  ! Local variables and parameters
  ! ------------------------------
  INTEGER (KIND=i4), PARAMETER :: n_month = 12
  INTEGER (KIND=i4), DIMENSION (n_month), PARAMETER :: &
       days_per_month = (/ 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 /)

  ! -------------------------------------------------------
  ! First find day of the year in a regular (non leap) year
  ! -------------------------------------------------------
  SELECT CASE ( month )
  CASE ( 1 ) 
     jday = day
  CASE (2:12)
     jday = SUM ( days_per_month(1:month-1) ) + day
  CASE DEFAULT
     jday = -9999
  END SELECT

  ! -------------------------
  ! Now apply leap year rules
  ! ------------------------------------
  ! * Divisible by 4:   leap year
  ! * Divisible by 100: not a leap year
  ! * Divisible by 400: leap year
  ! ------------------------------------
  IF ( MOD(year,4) == 0 .AND. ( MOD(year,100) /= 0 .OR. MOD(year,400) == 0 ) ) jday = jday+1

  RETURN
END FUNCTION day_of_year

SUBROUTINE year_month_day ( year, month, day, newyear, newmonth, newday )

  USE OMSAO_precision_module, ONLY: i4
  IMPLICIT NONE

  ! ----------------------
  ! Input/Output variables
  ! ----------------------
  INTEGER (KIND=i4), INTENT (IN)    :: year, month, day
  INTEGER (KIND=i4), INTENT (INOUT) :: newyear, newmonth, newday

  ! ------------------------------
  ! Local variables and parameters
  ! ------------------------------
  INTEGER (KIND=i4), PARAMETER :: n_month = 12
  INTEGER (KIND=i4), DIMENSION (n_month), PARAMETER :: &
       days_per_month    = (/ 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 /), &
       ly_days_per_month = (/ 31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 /)

  INTEGER (KIND=i4), DIMENSION (n_month) :: dpm

  ! -------------------------
  ! Check for leap year
  ! ------------------------------------
  ! * Divisible by 4:   leap year
  ! * Divisible by 100: not a leap year
  ! * Divisible by 400: leap year
  ! ------------------------------------
  IF ( MOD(year,4) == 0 .AND. ( MOD(year,100) /= 0 .OR. MOD(year,400) == 0 ) ) THEN
     dpm = days_per_month
  ELSE
     dpm = ly_days_per_month
  END IF

  newyear = year  ;  newmonth = month  ;  newday = day

  ! -------------------------------------------------------
  ! If DAY is within bounds there's nothing to do
  ! (MONTH hasn't been adjusted yet in the calling program)
  ! -------------------------------------------------------
  IF ( day > 0 .AND. day <= dpm(month) ) RETURN

  ! ---------------------------------------------------------------------------
  ! If we reach here, DAY must be either too large or <= 0. The crucial
  ! variable now is MONTH, since it might lead to an additional change in YEAR.
  ! ---------------------------------------------------------------------------
  SELECT CASE ( month )
  CASE ( 2:11 )   ! The easiest case: change only MONTH and DAY
     IF ( day <= 0 ) THEN
        newmonth = month - 1
        newday   = dpm(newmonth) + day  ! DAY is negative
     ELSE
        newmonth = month + 1
        newday   = day - dpm(month)
     END IF
  CASE ( 1 )      ! January requires a YEAR change if DAY <= 0
     IF ( day <= 0 ) THEN
        newmonth = 12
        newday   = dpm(newmonth) + day  ! DAY is negative
        newyear  = year - 1
     ELSE
        newmonth = month + 1
        newday   = day - dpm(month)
     END IF
  CASE ( 12 )      ! December requires a YEAR change if DAY > 31
     IF ( day <= 0 ) THEN
        newmonth = month - 1
        newday   = dpm(newmonth) + day  ! DAY is negative
     ELSE
        newmonth = 1
        newday   = day - dpm(month)
        newyear  = year + 1
     END IF
  END SELECT

  RETURN
END SUBROUTINE year_month_day


CHARACTER(LEN=*) FUNCTION upper_case ( mixstring ) RESULT ( upcase )

  ! =============================================
  ! Function to convert strings to all upper case
  ! =============================================

  USE OMSAO_precision_module, ONLY: i4
  IMPLICIT NONE

  ! ------------------------------------------------------------
  CHARACTER (LEN=*), INTENT (IN) :: mixstring   ! Input string
  ! ------------------------------------------------------------

  ! -------------------------------------------------------------
  !CHARACTER (LEN=*)              :: upcase      ! Result string
  ! -------------------------------------------------------------

  ! -------------------------------------------------------------
  INTEGER (KIND=i4) :: i                         ! Local variables
  ! -------------------------------------------------------------

  DO i = 1, LEN(mixstring)
     SELECT CASE ( ICHAR(mixstring(i:i)) )
     CASE ( 97:122 )
        upcase(i:i) = ACHAR(ICHAR(mixstring(i:i))-32)
     CASE DEFAULT
        upcase(i:i) = mixstring(i:i)
     END SELECT
  END DO

  RETURN
END FUNCTION upper_case


SUBROUTINE remove_quotes ( quotestring ) 

  ! ==========================================================
  ! Subroutine to remove any quotation marks (") from a string
  ! ==========================================================

  USE OMSAO_precision_module, ONLY: i4
  IMPLICIT NONE

  ! ---------------------
  ! Input/modified string
  ! ---------------------
  CHARACTER (LEN=*), INTENT (INOUT) :: quotestring

  ! ---------------
  ! Local variables
  ! ---------------
  INTEGER (KIND=i4) :: i, j                 
  CHARACTER (LEN=LEN(quotestring)) :: tmpstring

  IF ( INDEX(quotestring, '"') == 0 .OR. LEN(quotestring) == 0 ) RETURN

  tmpstring = ""; j = 1
  DO i = 1, LEN(quotestring)
     SELECT CASE ( quotestring(i:i) )
     CASE ( '"' )
        ! skip; nothing else to do
     CASE DEFAULT
        tmpstring(j:j) = quotestring(i:i)
        j = j+1
     END SELECT
  END DO

  quotestring = tmpstring

  RETURN
END SUBROUTINE remove_quotes

SUBROUTINE check_for_endofinput ( iostring, yn_eoi )

  USE OMSAO_indices_module, ONLY : eoi_str
  IMPLICIT NONE

  ! ==============
  ! Input variable
  ! ==============
  CHARACTER (LEN=*), INTENT (IN) :: iostring

  ! ===============
  ! Output variable
  ! ===============
  LOGICAL, INTENT (OUT) :: yn_eoi

  yn_eoi = .FALSE.
  IF ( TRIM(ADJUSTL(iostring)) == eoi_str ) yn_eoi = .TRUE.

  RETURN
END SUBROUTINE check_for_endofinput

SUBROUTINE gome_check_read_status ( ios, file_read_stat )

  USE OMSAO_precision_module, ONLY: i4
  USE OMSAO_errstat_module,   ONLY: file_read_ok, file_read_failed, file_read_eof
  IMPLICIT NONE

  INTEGER (KIND=i4), INTENT (IN)  :: ios
  INTEGER (KIND=i4), INTENT (OUT) :: file_read_stat

  SELECT CASE ( ios )
  CASE ( :-1 )  ;  file_read_stat = file_read_eof
  CASE (   0 )  ;  file_read_stat = file_read_ok
  CASE DEFAULT  ;  file_read_stat = file_read_failed
  END SELECT

  RETURN
END SUBROUTINE gome_check_read_status

SUBROUTINE skip_to_filemark ( funit, lm_string, lastline, file_read_stat )

  USE OMSAO_precision_module, ONLY: i4
  USE OMSAO_variables_module, ONLY: maxchlen
  IMPLICIT NONE

  ! ===============
  ! Input variables
  ! ===============
  INTEGER   (KIND=i4), INTENT (IN) :: funit
  CHARACTER (LEN=*),   INTENT (IN) :: lm_string

  ! ================
  ! Output variables
  ! ================
  INTEGER   (KIND=i4), INTENT (OUT) :: file_read_stat
  CHARACTER (LEN=*),   INTENT (OUT) :: lastline

  ! ===============
  ! Local variables
  ! ===============
  INTEGER   (KIND=i4)      :: lmlen, ios
  CHARACTER (LEN=maxchlen) :: tmpline

  ! -------------------------------------------
  ! Determine the length of the string landmark
  ! -------------------------------------------
  lmlen = LEN(TRIM(ADJUSTL(lm_string)))

  ! -------------------------------------------------------
  ! Read lines in the file until we either find the string,
  ! reach the end of the file, or reading fails otherwise.
  ! ----------------------------------------------------
  ios = 0
  getlm: DO WHILE ( ios == 0 )
     READ (UNIT=funit, FMT='(A)', IOSTAT=ios) tmpline
     tmpline = TRIM(ADJUSTL(tmpline))
     IF ( ios /= 0 .OR. tmpline(1:lmlen) == lm_string ) EXIT getlm
  END DO getlm
  
  ! ---------------------------------------------------
  ! Return the last line read for the case that we need 
  ! to extract further information from it
  ! ---------------------------------------------------
  lastline = TRIM(ADJUSTL(tmpline))

  CALL gome_check_read_status ( ios, file_read_stat )

  RETURN
END SUBROUTINE skip_to_filemark


SUBROUTINE get_substring ( string, sstart, substring, nsubstring, eostring )

  USE OMSAO_precision_module, ONLY: i4
  IMPLICIT NONE

  ! ==========================================================================
  ! Given string STRING, extract first space- or comma-delimited substring
  ! beginning on or after position SSTART, and return its value, SUBSTRING and
  ! its length, NSUBSTRING; update STRING and SSTART to remove this substring.
  ! 
  ! F90 version of the original GET_TOKEN by J. Lavanigno
  ! ==========================================================================

  ! NSUBSTRING represents the token's length without trailing blanks. It's
  ! set to 0 if no substring was found in STRING: this can be used to determine
  ! when to stop looking.
  !
  ! This routine makes no attempt to detect or handle the case in which
  ! the token is bigger than the token buffer.  It's assumed that the
  ! caller will declare line and token to be the same size so that this
  ! won't ever happen.

  ! ==================
  ! Modified arguments.
  ! ==================
  CHARACTER (LEN = *), INTENT (INOUT) :: string
  INTEGER   (KIND=i4), INTENT (INOUT) :: sstart

  ! =================
  ! Output arguments.
  ! =================
  CHARACTER (LEN =LEN(string)), INTENT (OUT) :: substring
  INTEGER   (KIND=i4),          INTENT (OUT) :: nsubstring
  LOGICAL,                      INTENT (OUT) :: eostring

  ! ================
  ! Local arguments.
  ! ================
  CHARACTER (LEN=1)   :: char
  INTEGER   (KIND=i4) :: lstart, lend, nline
  

  ! ------------------
  ! Initialize outputs
  ! ------------------
  substring  = ' '  ;  nsubstring = 0  ;  eostring = .FALSE.

  nline = LEN(TRIM(ADJUSTL(string)) )


  ! ------------------------------------------------------
  ! We are working on character variables, i.e., positions 
  ! are always larger than 0
  ! ------------------------------------------------------
  IF ( sstart <= 0 ) sstart = 1

  ! ----------------------------------------------
  ! Check first whether we have to any work at all
  ! ----------------------------------------------
  IF ( sstart >= nline ) THEN
     eostring = .TRUE.  ;  RETURN
  END IF

  ! --------------------------------------------------------
  ! Find first character in line that's not a blank or comma
  ! --------------------------------------------------------
  findchar: DO lstart = sstart, nline
     char = string(lstart:lstart)
     SELECT CASE ( char )
     CASE ( ' ' )
        IF ( lstart == nline ) THEN
           sstart = nline + 1; RETURN  ! there are no further substrings in STRING
        END IF
     CASE ( ',' )
        IF ( lstart == nline ) THEN
           sstart = nline + 1; RETURN  ! there are no further substrings in STRING
        END IF
     CASE DEFAULT
        EXIT findchar
     END SELECT
  END DO findchar


  ! --------------------------------------------------
  ! Start of the next substring is at position LSTART.
  ! Now find separator that ends the substring.
  ! --------------------------------------------------
  findsep: DO lend = lstart + 1, nline
     char = string (lend:lend)
     If ( char == ' ' .OR. char == ',' ) EXIT findsep
  END DO findsep
  IF ( (lend == nline) .AND. (char /= ' ' .AND. char /= ',') ) lend = lend + 1
  lend = lend - 1

  ! --------------------
  ! Output the substring
  ! --------------------
  nsubstring = lend - lstart + 1

  substring(1:nsubstring) = string (lstart:lend)

  ! -----------------------------------------------------------------
  ! The next substring, if any, starts at least two characters beyond
  ! the end of thelast (we have to skip over the comma or space that 
  ! marks the substring's end).
  ! -----------------------------------------------------------------
  sstart = lend + 2

  ! ---------------------------------------------
  ! Final check, whether we have to any more work
  ! ---------------------------------------------
  IF ( sstart >= nline ) eostring = .TRUE.

  RETURN
END SUBROUTINE get_substring


SUBROUTINE string2index ( table, ntable, string, stridx )

  ! =====================================================
  ! Looks up STRING in character table TABLE of dimension
  ! NTABLE, and returns position STRIDX. Defaults to
  ! STRIDX = -1 if STRING is not found in TABLE.
  ! =====================================================

  USE OMSAO_precision_module, ONLY: i4
  IMPLICIT NONE

  ! ===============
  ! Input variables
  ! ===============
  INTEGER   (KIND=i4),                     INTENT (IN) :: ntable
  CHARACTER (LEN=*),   DIMENSION (ntable), INTENT (IN) :: table
  CHARACTER (LEN=*),                       INTENT (IN) :: string

  ! ================
  ! Output variables
  ! ================
  INTEGER (KIND=i4), INTENT (OUT) :: stridx

  ! ===============
  ! Local variables
  ! ===============
  INTEGER (KIND=i4) :: i

  stridx = -1

  getidx: DO i = 1,  ntable
     IF ( TRIM(ADJUSTL(string)) == TRIM(ADJUSTL(table(i))) ) THEN
        stridx = i
        EXIT getidx
     END IF
  END DO getidx

  RETURN
END SUBROUTINE string2index


REAL (KIND=KIND(1.0D0)) FUNCTION signdp ( x )

  USE OMSAO_precision_module, ONLY: r8
  IMPLICIT NONE

  REAL (KIND=r8), INTENT (IN) :: x

  signdp = 0.0_r8
  IF ( x < 0.0_r8 ) THEN
     signdp = -1.0_r8
  ELSE
     signdp = +1.0_r8
  END IF

  RETURN
END FUNCTION signdp

SUBROUTINE interpolation ( &
     calledfrom,           &
     n1, x1, y1, n2, x2, y2, filltype, fillval, yn_full_range, errstat )

  USE OMSAO_precision_module
  USE OMSAO_errstat_module, ONLY: pge_errstat_ok, pge_errstat_warning, pge_errstat_error
  IMPLICIT NONE

  ! ---------------
  ! Input variables
  ! ---------------
  CHARACTER (LEN=*),                 INTENT (IN) :: calledfrom, filltype
  INTEGER (KIND=i4),                 INTENT (IN) :: n1, n2
  REAL    (KIND=r8),                 INTENT (IN) :: fillval
  REAL    (KIND=r8), DIMENSION (n1), INTENT (IN) :: x1, y1
  REAL    (KIND=r8), DIMENSION (n2), INTENT (IN) :: x2

  ! ----------------
  ! Output variables
  ! ----------------
  INTEGER (KIND=i4),                 INTENT (INOUT) :: errstat
  REAL    (KIND=r8), DIMENSION (n2), INTENT (OUT)   :: y2
  LOGICAL,                           INTENT (OUT)   :: yn_full_range

  ! --------------
  ! Local variable
  ! --------------
  INTEGER (KIND=i4)                 :: locerrstat, imin, imax, nloc
  REAL    (KIND=r8), DIMENSION (n2) :: xtmp, ytmp

  locerrstat = pge_errstat_ok

  ! -------------------------------
  ! Initialize interpolation output
  ! -------------------------------
  y2(1:n2) = 0.0_r8

  ! --------------------------------------------------------------------------
  ! Find indices in radiance wavelength spectrum that cover reference spectrum
  ! --------------------------------------------------------------------------
  imin = -1 ; imax = -1
  !imin = MAXVAL ( MINLOC ( x2(1:n2), MASK=(x2(1:n2) >= x1( 1)) ) )
  !imax = MAXVAL ( MAXLOC ( x2(1:n2), MASK=(x2(1:n2) <= x1(n1)) ) )
  CALL array_locate_r8 ( n2, x2(1:n2), x1( 1), 'GE', imin )
  CALL array_locate_r8 ( n2, x2(1:n2), x1(n1), 'LE', imax )

  ! -------------------------------------------------------------------------------
  ! Check whether we have the whole wavelength range. We don't set the error status
  ! variable for this case, because it is too insignificant to set non-Zero exit at
  ! PGE termination. Instead, we check in the calling subroutine for YN_FULL_RANGE
  ! and report this at higher verbosity thresholds.
  ! ------------------------------------------------------------------------------
  yn_full_range = .TRUE.
  IF ( imin /= 1 .OR. imax /= n2 ) yn_full_range = .FALSE.

  ! --------------------------------------------------------------------------
  ! Now that we know the first and last index to cover with the interpolation,
  ! we can treat all cases alike. Only we need make sure that the number of
  ! interpolation points is consistent. And if we don't have *any* points,
  ! then we must return without calling the interpolation routine.
  ! --------------------------------------------------------------------------
  nloc = imax - imin + 1

  ! ------------------------------------------------------------------
  ! Anything less than 4 data points will make the interpolation fail.
  ! So before we choose to interpolate or not, we set the interpolated
  ! array to Zero. That is easier than a IF_THEN_ELSEing based on the
  ! IMIN and IMAX indices.
  ! ------------------------------------------------------------------
  SELECT CASE ( nloc )
  CASE ( :3 ) ! Less than 4 data points available for interpolation
     errstat = pge_errstat_error
     RETURN
  CASE DEFAULT
     xtmp(1:nloc) = x2(imin:imax)
     CALL ezspline_1d_interpolation (                &
          n1,   x1  (1:n1),   y1  (1:n1),            &
          nloc, xtmp(1:nloc), ytmp(1:nloc), locerrstat )
     y2(imin:imax) = ytmp(1:nloc)
     CALL fill_nonoverlap ( n2, y2(1:n2), imin, imax, filltype, fillval )
     ! -----------------------------------------------------------------
     ! If we have non-zero exit status, something must have gone wrong
     ! in the interpolation. Set PGE_ERROR_STATUS to ERROR in this case.
     ! -----------------------------------------------------------------
     errstat = MAX ( errstat, locerrstat )
  END SELECT

  RETURN
END SUBROUTINE interpolation

SUBROUTINE array_sort_r8 (n, x, y )

  USE OMSAO_precision_module, ONLY: r8
  IMPLICIT NONE

  ! --------------
  ! Input variable
  ! --------------
  INTEGER, INTENT (IN) :: n

  ! ------------------
  ! Modified variables
  ! ------------------
  REAL (KIND=r8), DIMENSION (n), INTENT (INOUT) :: x, y

  ! ---------------
  ! Local variables
  ! ---------------
  REAL (KIND=r8) :: xfirst, xtemp, yfirst, ytemp
  INTEGER        :: index, i, j

  ! -----------------------------------
  ! Loop to sort N-1 entries into order
  ! -----------------------------------
  DO i = 1, n-1

     ! ----------------------------------------------------------------
     ! Initialize earliest so far to be in the first place in this pass
     ! ----------------------------------------------------------------
     xfirst = x(i) ; yfirst = y(i)
     index = i

     ! ----------------------------------------------
     ! If X(I+1) > X(I) we can skip to the next index
     ! ----------------------------------------------
     IF ( x(i+1) > x(i) ) CYCLE

     ! --------------------------------------------------
     ! Search remaining (unsorted) items for earliest one
     ! --------------------------------------------------
     DO j = i+1, n
        IF ( x(j) < xfirst ) THEN
           xfirst = x(j) ; yfirst = y(j)
           index = j
        END IF
     END DO

     IF ( index /= i ) THEN
        xtemp    = x(i)     ; ytemp    = y(i)
        x(i)     = x(index) ; y(i)     = y(index)
        x(index) = xtemp    ; y(index) = ytemp
     END IF

  END DO

  RETURN
END SUBROUTINE array_sort_r8

SUBROUTINE find_overlap ( n1, arr1, n2, arr2, i1, i2 )

  ! ----------------------------------------------------------------
  !
  ! Given arrays ARR1 of dimension N1 and ARR2 of dimension N2, find
  ! array positions I1 and I2 in ARR2 such that ARR2(I1:i2) fully
  ! overlaps with ARR1(1:N1):
  !
  !        ARR1(1) <= ARR2(I1) <= ARR2(I2) <= ARR1(N1)
  !
  ! ----------------------------------------------------------------

  USE OMSAO_precision_module, ONLY: i4, r8
  IMPLICIT NONE

  ! ---------------
  ! Input variables
  ! ---------------
  INTEGER (KIND=i4), INTENT (IN) :: n1, n2
  REAL    (KIND=r8), DIMENSION (n1), INTENT (IN) :: arr1
  REAL    (KIND=r8), DIMENSION (n2), INTENT (IN) :: arr2

  ! ----------------
  ! Output variables
  ! ----------------
  INTEGER (KIND=i4), INTENT (OUT) :: i1, i2

  ! -----------------------------------------------------------------
  ! ARRAY_LOCATE is the preferred way of determining the indices, but
  ! this branch of the code has not yet been tested. tpk 09 Feb 2007
  ! -----------------------------------------------------------------
  i1 = -1 ; i2 = -1
  i1 = MAXVAL ( MINLOC ( arr2(1:n2), MASK=(arr2(1:n2) >= arr1(1 )) ) )
  i2 = MAXVAL ( MAXLOC ( arr2(1:n2), MASK=(arr2(1:n2) <= arr1(n1)) ) )
  !CALL array_locate_r8 ( n2, arr2(1:n2), arr1( 1), 'GE', i1 )
  !CALL array_locate_r8 ( n2, arr2(1:n2), arr1(n1), 'LE', i2 )

  RETURN
END SUBROUTINE find_overlap

SUBROUTINE fill_nonoverlap ( n, arr, i1, i2, filltype, fillval )

  ! ----------------------------------------------------------------
  !
  ! Given array ARR of dimension N and array positions I1 and I2 in
  ! ARR as returned from FIND_OVERLAP, fill in any non-overlapping
  ! segments either by constant extrapolation of endpoints or by
  ! using the value of FILLVAL.
  !
  ! ----------------------------------------------------------------

  USE OMSAO_precision_module, ONLY: i4, r8

  ! ---------------
  ! Input variables
  ! ---------------
  CHARACTER (LEN=9),                INTENT (IN)    :: filltype
  INTEGER (KIND=i4),                INTENT (IN)    :: n, i1, i2
  REAL    (KIND=r8),                INTENT (IN)    :: fillval

  ! ------------------
  ! Modified variables
  ! ------------------
  REAL    (KIND=r8), DIMENSION (n), INTENT (INOUT) :: arr

  ! ---------------
  ! Local variables
  ! ---------------
  REAL (KIND=r8) :: low_end, hig_end

  low_end = 0.0_r8 ; hig_end = 0.0_r8

  SELECT CASE ( filltype )
  CASE ( 'endpoints' )
     low_end = arr(i1) ; hig_end = arr(i2)
  CASE ( 'fillvalue' )
     low_end = fillval ; hig_end = fillval
  END SELECT

  IF ( i1 > 1 ) arr (   1:i1-1) = low_end
  IF ( i2 < n ) arr (i2+1:n   ) = hig_end

  RETURN
END SUBROUTINE fill_nonoverlap

FUNCTION roundoff_r8 ( ndecim, r8value ) RESULT ( r8rounded )

  USE OMSAO_precision_module, ONLY: i4, r8
  IMPLICIT NONE

  ! ---------------------------------------------------------------------
  ! Explanation of FUNCTION arguments:
  !
  !    ndecim ........... number of significant digits to keep
  !    r8value .......... DOUBLE PRECISION number to be truncated
  !    r8rounded ........ truncated value
  ! ---------------------------------------------------------------------

  ! ---------------
  ! Input variables
  ! ---------------
  INTEGER (KIND=i4), INTENT (IN)    :: ndecim
  REAL    (KIND=r8), INTENT (INOUT) :: r8value

  ! --------------
  ! Returned value
  ! --------------
  REAL    (KIND=r8) :: r8rounded

  ! ---------------
  ! Local variables
  ! ---------------
  INTEGER (KIND=i4) :: pow
  REAL    (KIND=r8) :: tmpval


  ! ------------------------------------------
  ! Nothing to be done if we have a Zero value
  ! ------------------------------------------
  IF (  r8value == 0.0_r8 ) THEN
     r8rounded = 0.0_r8 ; RETURN
  END IF

  ! ---------------------------------
  ! Power of 10 of the original value
  ! ---------------------------------
  pow = INT ( LOG10 (ABS(r8value)), KIND=i4 )

  ! ------------------------------------------------------
  ! Remove original power of 10 and shift by NDECIM powers
  ! ------------------------------------------------------
  tmpval = r8value * 10.0_r8**(ndecim-pow)

  ! ---------------------------------------------
  ! Find the nearest INTEGER and undo power-shift
  ! ---------------------------------------------
  r8rounded = ANINT ( tmpval )  * 10.0_r8**(pow-ndecim)

  RETURN
END FUNCTION roundoff_r8

FUNCTION roundoff_r4 ( ndecim, r4value ) RESULT ( r4rounded )

  USE OMSAO_precision_module, ONLY: i4, r4
  IMPLICIT NONE

  ! ---------------
  ! Input variables
  ! ---------------
  INTEGER (KIND=i4), INTENT (IN)    :: ndecim
  REAL    (KIND=r4), INTENT (INOUT) :: r4value


  INTEGER (KIND=i4) :: pow
  REAL    (KIND=r4) :: r4rounded, tmpval

  ! ------------------------------------------
  ! Nothing to be done if we have a Zero value
  ! ------------------------------------------
  IF (  r4value == 0.0_r4 ) THEN
     r4rounded = 0.0_r4 ; RETURN
  END IF

  ! ---------------------------------
  ! Power of 10 of the original value
  ! ---------------------------------
  pow = INT ( LOG10 (ABS(r4value)), KIND=i4 )

  ! ------------------------------------------------------
  ! Remove original power of 10 and shift by NDECIM powers
  ! ------------------------------------------------------
  tmpval = r4value * 10.0_r4**(ndecim-pow)

  ! ---------------------------------------------
  ! Find the nearest INTEGER and undo power-shift
  ! ---------------------------------------------
  r4rounded = ANINT ( tmpval ) * 10.0_r4**(pow-ndecim)

  RETURN
END FUNCTION roundoff_r4


SUBROUTINE roundoff_1darr_r8 ( ndecim, ndim, r8value )

  USE OMSAO_precision_module, ONLY: i4, r8
  IMPLICIT NONE

  ! ---------------
  ! Input variables
  ! ---------------
  INTEGER (KIND=i4),                   INTENT (IN)    :: ndecim, ndim
  REAL    (KIND=r8), DIMENSION (ndim), INTENT (INOUT) :: r8value


  INTEGER (KIND=i4) :: i
  REAL    (KIND=r8) :: roundoff_r8

  DO i = 1, ndim
     r8value(i) = ROUNDOFF_R8 ( ndecim, r8value(i) )
  END DO

  RETURN
END SUBROUTINE roundoff_1darr_r8


SUBROUTINE roundoff_1darr_r4 ( ndecim, ndim, r4value )

  USE OMSAO_precision_module, ONLY: i4, r4
  IMPLICIT NONE

  ! ---------------
  ! Input variables
  ! ---------------
  INTEGER (KIND=i4),                   INTENT (IN)    :: ndecim, ndim
  REAL    (KIND=r4), DIMENSION (ndim), INTENT (INOUT) :: r4value


  INTEGER (KIND=i4) :: i
  REAL    (KIND=r4) :: roundoff_r4

  DO i = 1, ndim
     r4value(i) = ROUNDOFF_R4 ( ndecim, r4value(i) )
  END DO

  RETURN
END SUBROUTINE roundoff_1darr_r4


SUBROUTINE roundoff_2darr_r8 ( ndecim, n1, n2, r8value )

  USE OMSAO_precision_module, ONLY: i4, r8
  IMPLICIT NONE

  ! ---------------
  ! Input variables
  ! ---------------
  INTEGER (KIND=i4),                    INTENT (IN)    :: ndecim, n1, n2
  REAL    (KIND=r8), DIMENSION (n1,n2), INTENT (INOUT) :: r8value


  INTEGER (KIND=i4) :: i, j
  REAL    (KIND=r8) :: roundoff_r8

  DO i = 1, n1
     DO j = 1, n2
        r8value(i,j) = ROUNDOFF_R8 ( ndecim, r8value(i,j) )
     END DO
  END DO

  RETURN
END SUBROUTINE roundoff_2darr_r8

SUBROUTINE roundoff_2darr_r4 ( ndecim, n1, n2, r4value )

  USE OMSAO_precision_module, ONLY: i4, r4
  IMPLICIT NONE

  ! ---------------
  ! Input variables
  ! ---------------
  INTEGER (KIND=i4),                    INTENT (IN)    :: ndecim, n1, n2
  REAL    (KIND=r4), DIMENSION (n1,n2), INTENT (INOUT) :: r4value


  INTEGER (KIND=i4) :: i, j
  REAL    (KIND=r4) :: roundoff_r4

  DO i = 1, n1
     DO j = 1, n2
        r4value(i,j) = ROUNDOFF_R4 ( ndecim, r4value(i,j) )
     END DO
  END DO

  RETURN
END SUBROUTINE roundoff_2darr_r4

SUBROUTINE array_locate_r4 ( n, x, x0, psel, ipos )

  ! -------------------------------------------------------------
  ! Given an ordered array X of size N and a point X0, return the
  ! position IPOS of X0 in X based on the selection PSEL.
  ! -------------------------------------------------------------

  USE OMSAO_precision_module, ONLY: i4, r4
  IMPLICIT NONE

  ! ---------------
  ! Input variables
  ! ---------------
  INTEGER (KIND=i4),                INTENT (IN) :: n
  REAL    (KIND=r4), DIMENSION (n), INTENT (IN) :: x
  REAL    (KIND=r4),                INTENT (IN) :: x0
  CHARACTER (LEN=2),                INTENT (IN) :: psel

  ! ----------------
  ! Output variables
  ! ----------------
  INTEGER (KIND=i4), INTENT (OUT) :: ipos

  ! ---------------
  ! Local variables
  ! ---------------
  INTEGER (KIND=i4)  :: ilow, iupp, imid

  ! --------------------
  ! Initialize variables
  ! --------------------
  ipos = -1
  ilow = 0
  iupp = n+1
  imid = 0

  ! -------------------------------------
  ! The actual bisection. Short and sweet
  ! -------------------------------------
  DO WHILE ( iupp-ilow > 1 )
     imid = ( iupp + ilow ) / 2
     IF ( x(imid) >= x0 ) THEN
        iupp = imid
     ELSE
        ilow = imid
     END IF
  END DO

  ! -------------------------------------------------------
  ! The final value of IPOS depends on what we are actually
  ! looking for. Depending on our selection, we may have to
  ! return IPOS=-1 for cases that don't fit the selection.
  ! -------------------------------------------------------
  SELECT CASE ( psel )
  CASE ( 'LT' )
     ! X(IPOS)  < X0
     ipos = ilow
     IF ( x(ipos) == x0 ) ipos = ipos - 1
  CASE ( 'LE' )
     ! X(IPOS) <= X0
     ipos = ilow
     IF ( x(iupp) == x0 ) ipos = iupp
  CASE ( 'GE' )
     ! X(IPOS) >= X0
     ipos = iupp
  CASE ( 'GT' )
     ! X(IPOS)  > X0
     ipos = iupp
     IF ( x(ipos) == x0 ) ipos = ipos + 1
  CASE DEFAULT
     ! Return IPOS=IUP
     ipos = iupp
  END SELECT

  IF ( ipos <= 0 .OR. ipos > n ) ipos = -1

  RETURN
END SUBROUTINE array_locate_r4

SUBROUTINE array_locate_r8 ( n, x, x0, psel, ipos )

  ! -------------------------------------------------------------
  ! Given an ordered array X of size N and a point X0, return the
  ! position IPOS of X0 in X based on the selection PSEL.
  ! -------------------------------------------------------------

  USE OMSAO_precision_module, ONLY: i4, r8
  IMPLICIT NONE

  ! ---------------
  ! Input variables
  ! ---------------
  INTEGER (KIND=i4),                INTENT (IN) :: n
  REAL    (KIND=r8), DIMENSION (n), INTENT (IN) :: x
  REAL    (KIND=r8),                INTENT (IN) :: x0
  CHARACTER (LEN=2),                INTENT (IN) :: psel

  ! ----------------
  ! Output variables
  ! ----------------
  INTEGER (KIND=i4), INTENT (OUT) :: ipos

  ! ---------------
  ! Local variables
  ! ---------------
  INTEGER (KIND=i4)  :: ilow, iupp, imid

  ! --------------------
  ! Initialize variables
  ! --------------------
  ipos = -1
  ilow = 0
  iupp = n+1
  imid = 0

  ! -------------------------------------
  ! The actual bisection. Short and sweet
  ! -------------------------------------
  DO WHILE ( iupp-ilow > 1 )
     imid = ( iupp + ilow ) / 2
     IF ( x(imid) >= x0 ) THEN
        iupp = imid
     ELSE
        ilow = imid
     END IF
  END DO
  
  ! -------------------------------------------------------------
  ! Force the Upper and Lower indices into the range of the array
  ! -------------------------------------------------------------
  ilow = MAX ( ilow, 1 )
  iupp = MIN ( iupp, n )

  ! -------------------------------------------------------
  ! The final value of IPOS depends on what we are actually
  ! looking for. Depending on our selection, we may have to
  ! return IPOS=-1 for cases that don't fit the selection.
  ! -------------------------------------------------------
  SELECT CASE ( psel )
  CASE ( 'LT' )
     ! X(IPOS)  < X0
     ipos = ilow
     IF ( x(ipos) == x0 ) ipos = ipos - 1
  CASE ( 'LE' )
     ! X(IPOS) <= X0
     ipos = ilow
     IF ( x(iupp) == x0 ) ipos = iupp
  CASE ( 'GE' )
     ! X(IPOS) >= X0
     ipos = iupp
  CASE ( 'GT' )
     ! X(IPOS)  > X0
     ipos = iupp
     IF ( x(ipos) == x0 ) ipos = ipos + 1
  CASE DEFAULT
     ! Return IPOS=IUP
     ipos = iupp
  END SELECT

  IF ( ipos <= 0 .OR. ipos > n ) ipos = -1

  RETURN
END SUBROUTINE array_locate_r8

SUBROUTINE roundoff_3darr_r8 ( ndecim, n1, n2, n3, r8value )

  USE OMSAO_precision_module, ONLY: i4, r8
  IMPLICIT NONE

  ! ---------------
  ! Input variables
  ! ---------------
  INTEGER (KIND=i4),                       INTENT (IN)    :: ndecim, n1, n2, n3
  REAL    (KIND=r8), DIMENSION (n1,n2,n3), INTENT (INOUT) :: r8value


  INTEGER (KIND=i4) :: i, j, k
  REAL    (KIND=r8) :: roundoff_r8

  DO i = 1, n1
     DO j = 1, n2
        DO k = 1, n3
           r8value(i,j,k) = ROUNDOFF_R8 ( ndecim, r8value(i,j,k) )
        END DO
     END DO
  END DO

  RETURN
END SUBROUTINE roundoff_3darr_r8

