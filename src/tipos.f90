MODULE tipos
  USE decimal
  IMPLICIT NONE
  REAL(kind = dp), PARAMETER      :: grav = 9.8_dp
  TYPE SWSolution
    REAL(kind = dp), ALLOCATABLE    :: hh(:,:)        ! profundidad
    REAL(kind = dp), ALLOCATABLE    :: uu(:,:,:)      ! velocidad en x, (y)
    INTEGER                         :: dims           ! dimensiones
  END TYPE
END MODULE tipos