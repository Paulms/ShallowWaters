MODULE decimal
  IMPLICIT NONE
  !
  ! Definicion de la precision
  !
  INTEGER, PARAMETER :: dp=KIND(1.0d0)   
  !
  !
CONTAINS
  SUBROUTINE print_Rmatrix(mat)
    !========================================================
    ! Rutina para imprimir la matriz en formato completo
    !========================================================
    REAL(kind=dp)   :: mat(:,:)
    INTEGER         :: cols, rows, i, j
    rows = SIZE(mat, 1)
    cols = SIZE(mat, 2)
    do i=1,rows
      write (*,"("//trim(str(cols))//"(F10.4))") ( mat(i,j), j=1,cols )
    end do
  END SUBROUTINE

    CHARACTER(len=20) FUNCTION str(k)
  !   "Convert an integer to string."
    integer, intent(in) :: k
    write (str, *) k
    str = adjustl(str)
  END FUNCTION str
END MODULE decimal
