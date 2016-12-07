MODULE tipos
  USE decimal
  IMPLICIT NONE
  REAL(kind = dp), PARAMETER      :: grav = 9.8_dp
  TYPE SWSolution
    REAL(kind = dp), ALLOCATABLE    :: hh(:,:)        !profundidad
    REAL(kind = dp), ALLOCATABLE    :: uu(:,:)        !velocidad en x
    REAL(kind = dp), ALLOCATABLE    :: vv(:,:)        !velocidad en y
  END TYPE
END MODULE tipo