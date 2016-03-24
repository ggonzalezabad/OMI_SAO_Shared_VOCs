MODULE OMSAO_casestring_module

  ! ---------------------------------------------------------
  ! STRING conversion from upper to lower case and vice versa
  ! (based on the FDAS_STRING module by Richard Maine)
  ! ---------------------------------------------------------

  IMPLICIT NONE

  INTEGER, PARAMETER, PRIVATE :: mascii = 127
  INTEGER,            PRIVATE :: i_do
  INTEGER, PARAMETER, PRIVATE :: up_map_ascii(0:mascii) = (/ &
       (i_do,    i_do =   0,  96), &  ! Anything up to "Z"
       (i_do-32, i_do =  97, 122), &  ! "a" to "z"
       (i_do,    i_do = 123, 127) /)  ! Anything after "z"
  INTEGER, PARAMETER, PRIVATE :: lo_map_ascii(0:mascii) = (/ &      
       (i_do,    i_do =  0,  64),  &  ! Anything before "A"
       (i_do+32, i_do = 65,  96),  &  ! "A" to "Z"
       (i_do,    i_do = 97, 122),  &  ! "a" to "z"
       (i_do, i_do=123,127)       /)  ! Anything after "z"

CONTAINS

  FUNCTION upper_case (string) RESULT (upcase)

    ! ---------------------------------
    ! lowercase to UPPERCASE conversion
    ! ---------------------------------

    ! --------------
    ! Input variable
    ! --------------
    CHARACTER (LEN=*), INTENT(IN) :: string  !-- An arbitrary string.

    ! ---------------
    ! Result variable
    ! ---------------
    CHARACTER (LEN=LEN(string)) :: upcase  !-- The converted string.
    
    ! --------------
    ! Local variable
    ! --------------
    INTEGER :: i

    ! --------------------------------------------------
    ! Peform the conversion iff STRING(i:i) is within
    ! the accept ASCII code. Else copy character as is.
    ! --------------------------------------------------
    DO i = 1 , LEN(string)
       SELECT CASE ( IACHAR(string(i:i)) )
       CASE ( 0:mascii )
          upcase(i:i) = ACHAR(up_map_ascii(IACHAR(string(i:i))))
       CASE DEFAULT
          upcase(i:i) = string(i:i)
       END SELECT
    END DO
    RETURN

  END FUNCTION upper_case

  FUNCTION lower_case (string) RESULT (locase)

    ! ---------------------------------
    ! UPPERCASE to lowercase conversion
    ! ---------------------------------

    ! --------------
    ! Input variable
    ! --------------
    CHARACTER (LEN=*), INTENT(IN) :: string  !-- An arbitrary string.

    ! ---------------
    ! Result variable
    ! ---------------
    CHARACTER (LEN=LEN(string)) :: locase  !-- The converted string.
    
    ! --------------
    ! Local variable
    ! --------------
    INTEGER :: i

    ! --------------------------------------------------
    ! Peform the conversion iff STRING(i:i) is within
    ! the accept ASCII code. Else copy character as is.
    ! --------------------------------------------------
    DO i = 1 , LEN(string)
       SELECT CASE ( IACHAR(string(i:i)) )
       CASE ( 0:mascii )
          locase(i:i) = ACHAR(lo_map_ascii(IACHAR(string(i:i))))
       CASE DEFAULT
          locase(i:i) = string(i:i)
       END SELECT
    END DO
    RETURN

  END FUNCTION lower_case


END MODULE OMSAO_casestring_module
