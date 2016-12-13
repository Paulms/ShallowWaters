MODULE tipos
  USE decimal
  IMPLICIT NONE
  REAL(kind = dp), PARAMETER      :: grav = 9.8_dp
  TYPE SWSolution
    REAL(kind = dp), ALLOCATABLE    :: hh(:,:)        ! profundidad
    REAL(kind = dp), ALLOCATABLE    :: eta(:,:)        ! profundidad
    REAL(kind = dp), ALLOCATABLE    :: deta(:,:,:)        ! profundidad
    REAL(kind = dp), ALLOCATABLE    :: uu(:,:,:)      ! velocidad en x, (y)
    REAL(kind = dp), ALLOCATABLE    :: du(:,:,:,:)      ! cell averaged with limiter
    INTEGER                         :: dims           ! dimensiones
  END TYPE
  TYPE SWBed
    REAL(kind = dp), ALLOCATABLE   :: elev(:,:)        ! elevacion centros
    REAL(kind = dp), ALLOCATABLE    :: hc(:,:)        ! elevacion centros
    REAL(kind = dp), ALLOCATABLE    :: dz(:,:,:)      ! pendiente en x (y)
  END TYPE
END MODULE tipos