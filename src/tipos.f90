MODULE tipos
  USE decimal
  IMPLICIT NONE
  REAL(kind = dp), PARAMETER      :: grav = 9.8_dp
  TYPE SWSolution
    REAL(kind = dp), ALLOCATABLE    :: hh(:,:)        ! profundidad
    REAL(kind = dp), ALLOCATABLE    :: uu(:,:,:)      ! velocidad en x, (y)
    INTEGER                         :: dims           ! dimensiones
  END TYPE
  TYPE SWBed
    REAL(kind = dp), ALLOCATABLE   :: elev(:,:)        ! elevacion centros
    REAL(kind = dp), ALLOCATABLE    :: hc(:,:)        ! elevacion centros
    REAL(kind = dp), ALLOCATABLE    :: dx(:,:)      ! pendiente en x
    REAL(kind = dp), ALLOCATABLE    :: dy(:,:)      ! pendiente en y
  END TYPE
END MODULE tipos