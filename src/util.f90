MODULE util
  use decimal
IMPLICIT NONE
CONTAINS
  SUBROUTINE util_get_unit (iunit)
    !
    !*******************************************************************************
    !
    !! GET_UNIT returns a free FORTRAN unit number.
    !
    !
    !  Discussion:
    !
    !    A "free" FORTRAN unit number is an integer between 1 and 99 which
    !    is not currently associated with an I/O device.  A free FORTRAN unit
    !    number is needed in order to open a file with the OPEN command.
    !
    !  Modified:
    !
    !    02 March 1999
    !
    !  Author:
    !
    !    John Burkardt
    !
    !  Parameters:
    !
    !    Output, integer IUNIT.
    !
    !    If IUNIT = 0, then no free FORTRAN unit could be found, although
    !    all 99 units were checked (except for units 5 and 6).
    !
    !    Otherwise, IUNIT is an integer between 1 and 99, representing a
    !    free FORTRAN unit.  Note that GET_UNIT assumes that units 5 and 6
    !    are special, and will never return those values.
    !
    !
    INTEGER i, ios, iunit
    LOGICAL lopen
    !  
    iunit = 0
    !
    DO i = 1, 99    
       IF ( i /= 5 .AND. i /= 6 ) THEN        
          INQUIRE ( unit = i, opened = lopen, iostat = ios )        
          IF ( ios == 0 ) THEN
             IF ( .NOT. lopen ) THEN
                iunit = i
                RETURN
             END IF
          END IF
       END IF
    END DO
    !
  END SUBROUTINE util_get_unit

END MODULE util