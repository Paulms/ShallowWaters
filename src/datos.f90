MODULE datos
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Modulo para lectura de datos del problema            !!
  !!                                                      !!
  !! METODOS       :                                      !!
  !! leer_archivo                                         !!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  USE decimal
  USE util
  IMPLICIT NONE
  PUBLIC  :: leer_archivo
  PRIVATE 
  CONTAINS
  SUBROUTINE leer_archivo (file_input_name, cellsize, cellnumber, nt, dt, name, dims, &
  order, ejemplo, errors)
    !======================================================
    ! Leemos datos de un problema lineal almacenados
    ! en un solo archivo file_input_name
    !====================================================== 
    CHARACTER(32), INTENT(in)       ::  file_input_name 
    REAL(kind = dp)                 :: cellsize       !Tamaño de celda [m]
    INTEGER                         :: cellnumber     !numero de celdas
    INTEGER                         :: nt             !# pasos temporales
    REAL(kind =dp)                  :: dt             ! Tamaño del paso temporal
    CHARACTER(LEN=32)               :: name
    INTEGER                         :: ejemplo        ! 1 agua desde la esquina, 2 gota de agua
    INTEGER                         :: dims           ! dimensiones del problema
    INTEGER                         :: order          ! orden del esquema
    INTEGER                         :: errors         ! estimar errores
    INTEGER                         :: iunit          ! id para archivo
    CALL util_get_unit(iunit)
    OPEN(iunit,file=file_input_name, status='old', action='read')
      READ(iunit,*) cellsize
      READ(iunit,*) cellnumber
      READ(iunit,*) nt
      READ(iunit,*) dt
      READ(iunit,*) name
      READ(iunit,*) dims
      READ(iunit,*) order
      READ(iunit,*) ejemplo
      READ(iunit,*) errors
    CLOSE(iunit)
   END SUBROUTINE leer_archivo
END MODULE datos